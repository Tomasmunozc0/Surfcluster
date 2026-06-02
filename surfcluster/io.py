import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional
import pickle
import string

# --- Probe mapping -----------------------------------------------------------
# Source: pyMDMix solvent config files (github.com/CBDD/pyMDmix)
# Key: (residue_name, atom_name) as they appear in the raw MDMix hotspot PDB
# Value: pharmacophore type used by the clustering algorithm

PROBE_MAP = {
    'ETA':  {'CT': 'HYD', 'OH': 'POL'},
    'MOH':  {'CT': 'HYD', 'OH': 'POL'},
    'ISO':  {'CT': 'HYD', 'OH': 'POL'},
    'ISO5': {'CT': 'HYD', 'OH': 'POL'},
    'MAM':  {'CT': 'HYD', 'N':  'DON', 'O': 'ACC'},
    'PYR':  {'C':  'HYD', 'N':  'ACC'},
    'ANT':  {'C':  'HYD', 'N':  'ACC'},
    'WAT':  {'O':  'POL'},
    # DG0 is the MDMix energy grid label used as residue name in older outputs
    # atom name carries the actual atom type in this case
    'DG0':  {'CT': 'HYD', 'OH': 'POL', 'C': 'HYD', 'N': 'DON', 'O': 'ACC'},
}

# Element-based VdW radii (Angstrom) — used when mol2 is not provided
ELEMENT_VDW = {
    'C': 1.7, 'N': 1.64, 'O': 1.46, 'S': 1.78,
    'P': 1.87, 'F': 1.56, 'CL': 1.74, 'BR': 1.98,
    'I': 2.09, 'ZN': 0.6, 'H': 0.0,
}

# SYBYL atom type VdW radii — used when mol2 is provided
SYBYL_VDW = {
    'C.3': 1.88, 'C.2': 1.76, 'C.1': 1.61, 'C.ar': 1.88, 'C.cat': 1.88,
    'N.3': 1.64, 'N.2': 1.64, 'N.1': 1.64, 'N.ar': 1.64, 'N.am':  1.64,
    'N.4': 1.64, 'N.pl3': 1.63,
    'O.3': 1.46, 'O.2': 1.42, 'O.co2': 1.42, 'O.spc': 1.42, 'O.t3p': 1.42,
    'S.3': 1.782, 'S.2': 1.77, 'S.O': 1.77, 'S.O2': 1.77,
    'P': 1.871, 'F': 1.56, 'Cl': 1.735, 'Br': 1.978, 'I': 2.094,
    'ZN': 0.6, 'H': 0.0,
}

# Pharmacophore type → PDB element for output
TYPE_ELEMENT = {'HYD': 'C', 'POL': 'O', 'DON': 'N', 'ACC': 'F', 'POS': 'S', 'NEG': 'P'}


# --- Data structures ---------------------------------------------------------

@dataclass
class Hotspots:
    coords:    np.ndarray        # (N, 3)
    types:     np.ndarray        # (N,) pharmacophore strings
    energies:  np.ndarray        # (N,) kcal/mol


@dataclass
class Pocket:
    coords:    np.ndarray        # (M, 3)
    types:     np.ndarray        # (M,)
    energies:  np.ndarray        # (M,)
    dg:        float = 0.0
    ki:        float = 0.0
    efficiency: float = 0.0
    volume:    float = 0.0
    hyd_count: int = 0
    pol_count: int = 0
    don_count: int = 0
    acc_count: int = 0


# --- Receptor loading --------------------------------------------------------

def _parse_mol2_sybyl(mol2_path: str) -> dict:
    """Returns {atom_serial: sybyl_type} from mol2 ATOM section."""
    atom_flag = False
    id_to_type = {}
    with open(mol2_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('@<TRIPOS>ATOM'):
                atom_flag = True
                continue
            if line.startswith('@<TRIPOS>BOND'):
                break
            if atom_flag and line:
                parts = line.split()
                id_to_type[parts[0]] = parts[5]
    return id_to_type


def read_receptor(pdb_path: str, mol2_path: Optional[str] = None):
    """
    Load receptor coordinates and VdW radii.
    Returns coords (N,3) and vdw_radii (N,) as numpy arrays.
    Uses SYBYL radii from mol2 if provided, else element-based fallback.
    """
    coords = []
    elements = []
    serials = []

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip().upper()
            if not element:
                raw = line[12:16].strip().lstrip('0123456789').upper()
                element = raw[:2] if raw[:2] in ELEMENT_VDW else raw[:1]
            serial = line[6:11].strip()
            coords.append([x, y, z])
            elements.append(element)
            serials.append(serial)

    coords = np.array(coords, dtype=float)

    if mol2_path:
        sybyl_map = _parse_mol2_sybyl(mol2_path)
        vdw_radii = np.array([
            SYBYL_VDW.get(sybyl_map.get(s, ''), ELEMENT_VDW.get(e, 1.7))
            for s, e in zip(serials, elements)
        ])
    else:
        vdw_radii = np.array([ELEMENT_VDW.get(e, 1.7) for e in elements])

    return coords, vdw_radii


# --- Hotspot loading ---------------------------------------------------------

def _detect_format(atom_name: str, res_name: str) -> tuple:
    """
    Detect MDMix hotspot PDB format and return (probe, atom_type).

    Format A (older): atom_name=DG0, res_name=CT/OH/N/O
        -> probe='DG0', atom_type=res_name

    Format B (newer): atom_name=CT/OH/N/O, res_name=ETA/MAM/...
        -> probe=res_name, atom_type=atom_name

    Format C (post-sed): res_name=HYD/POL/DON/ACC
        -> already mapped, return directly
    """
    known_types = {'HYD', 'POL', 'DON', 'ACC', 'POS', 'NEG'}
    if res_name in known_types:
        return None, res_name

    if atom_name == 'DG0':
        return 'DG0', res_name

    if res_name in PROBE_MAP:
        return res_name, atom_name

    return None, None


def read_hotspots(pdb_paths: list) -> Hotspots:
    """
    Load one or more MDMix hotspot PDB files.
    Handles all known MDMix output formats automatically.
    Filters WAT probes.

    Uses fixed PDB column positions for robustness (cols 13-16 atom name,
    18-20 res name, 31-54 coords, 61-66 B-factor as energy).
    """
    all_coords = []
    all_types = []
    all_energies = []
    unknown_probes = set()

    for pdb_path in pdb_paths:
        with open(pdb_path) as f:
            for line in f:
                if not line.startswith('ATOM'):
                    continue
                if len(line.rstrip('\n')) < 66:
                    continue

                try:
                    atom_name = line[12:16].strip()
                    res_name  = line[17:20].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    energy = float(line[60:66])
                except (ValueError, IndexError):
                    continue

                if res_name == 'WAT':
                    continue

                probe, atom_type = _detect_format(atom_name, res_name)

                if atom_type in {'HYD', 'POL', 'DON', 'ACC', 'POS', 'NEG'}:
                    pharm_type = atom_type
                elif probe and probe in PROBE_MAP:
                    pharm_type = PROBE_MAP[probe].get(atom_type)
                    if pharm_type is None:
                        unknown_probes.add(f"{probe}:{atom_type}")
                        continue
                else:
                    unknown_probes.add(f"{res_name}:{atom_name}")
                    continue

                all_coords.append([x, y, z])
                all_types.append(pharm_type)
                all_energies.append(energy)

    if unknown_probes:
        print(f"WARNING: Unknown probe/atom combinations skipped: {unknown_probes}")
        print("  Add them to PROBE_MAP in io.py to include them.")

    return Hotspots(
        coords=np.array(all_coords, dtype=float),
        types=np.array(all_types),
        energies=np.array(all_energies, dtype=float),
    )


# --- Cache -------------------------------------------------------------------

CACHE_FORMAT = 2  # bump when cache structure changes


def save_cache(surface_dm: np.ndarray, removed_ids: np.ndarray,
               meta: dict, path: str):
    """Save surface DM + removed IDs + metadata for cache validation."""
    payload = {
        'format': CACHE_FORMAT,
        'surface_dm': surface_dm,
        'removed_ids': removed_ids,
        'meta': meta,
    }
    with open(path, 'wb') as f:
        pickle.dump(payload, f)


def load_cache(path: str):
    """Load cache. Returns (surface_dm, removed_ids, meta).
    Raises ValueError if cache format is incompatible."""
    with open(path, 'rb') as f:
        data = pickle.load(f)
    if not isinstance(data, dict) or data.get('format') != CACHE_FORMAT:
        raise ValueError(f"incompatible cache format (expected {CACHE_FORMAT})")
    return data['surface_dm'], data['removed_ids'], data.get('meta', {})


# --- Output writers ----------------------------------------------------------

def write_clusts_pdb(pockets: list, path: str, name: str):
    """Write clustering results as PDB with HEAD metadata."""
    ascii_upper = string.ascii_uppercase
    lines_head = []
    lines_atoms = []
    atom_serial = 1

    for i, pocket in enumerate(pockets):
        chain = ascii_upper[i % 26]
        lines_head.append(f"REMARK 999 -----------")
        lines_head.append(f"REMARK 999 Pocket on chain {chain}:")
        lines_head.append(f"REMARK 999   DG: {pocket.dg:.2f} kcal/mol")
        lines_head.append(f"REMARK 999   Ki: {pocket.ki:.4f} nM")
        if pocket.dg < -15:
            lines_head.append(f"REMARK 999   WARNING: DG very low, check clustering parameters")

        for j in range(len(pocket.coords)):
            x, y, z = pocket.coords[j]
            ptype    = pocket.types[j]
            energy   = pocket.energies[j]
            element  = TYPE_ELEMENT.get(ptype, 'C')
            resid    = i + 1
            lines_atoms.append(
                f"ATOM  {atom_serial:5d}{element:>4s}  {ptype:<3s} {chain}{resid:4d}   "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  0.00{energy:6.2f}"
            )
            atom_serial += 1

    with open(path, 'w') as f:
        f.write('\n'.join(lines_head + lines_atoms) + '\n')


def write_summary_csv(pockets: list, path: str, name: str):
    """Write per-pocket metrics to CSV."""
    ascii_upper = string.ascii_uppercase
    rows = []
    for i, pocket in enumerate(pockets):
        rows.append({
            'Pocket':      ascii_upper[i % 26],
            'DG':          round(pocket.dg, 3),
            'Ki_nM':       round(pocket.ki, 4),
            'Efficiency':  round(pocket.efficiency, 3),
            'HYD':         pocket.hyd_count,
            'POL':         pocket.pol_count,
            'DON':         pocket.don_count,
            'ACC':         pocket.acc_count,
            'N_hotspots':  len(pocket.coords),
            'Volume_A3':   round(pocket.volume, 2),
        })
    pd.DataFrame(rows).to_csv(path, index=False)


def write_sgraph_xyz(surface_coords: np.ndarray, path: str):
    """Write surface graph nodes as XYZ for PyMOL visualization."""
    with open(path, 'w') as f:
        f.write(f"{len(surface_coords)}\n\n")
        for x, y, z in surface_coords:
            f.write(f"C {x:.4f} {y:.4f} {z:.4f}\n")


def write_pymol_session(pockets: list, receptor_path: str, path: str, name: str):
    """Write a PyMOL .pml script ready to load and visualize results.
    Paths are relative to the .pml location so the result folder is portable."""
    import os
    ascii_upper = string.ascii_uppercase
    pml_dir = Path(path).parent.resolve()

    # clusts.pdb is always in the same dir as the .pml
    clusts_rel = f"{name}_clusts.pdb"

    # Receptor may live anywhere — try relative, fall back to absolute
    try:
        receptor_rel = os.path.relpath(Path(receptor_path).resolve(), pml_dir)
    except ValueError:
        receptor_rel = str(Path(receptor_path).resolve())

    colors = ['red', 'blue', 'green', 'yellow', 'magenta', 'cyan',
              'orange', 'violet', 'salmon', 'lime']

    lines = [
        f"load {receptor_rel}, receptor",
        f"load {clusts_rel}, pockets",
        "hide everything",
        "show surface, receptor",
        "color grey80, receptor",
        "show spheres, pockets",
        "set sphere_scale, 0.4",
        "",
    ]

    for i, pocket in enumerate(pockets):
        chain = ascii_upper[i % 26]
        color = colors[i % len(colors)]
        ki_str = f"  Ki={pocket.ki:.2f}nM" if pocket.ki < 1e6 else ""
        label = f"pocket_{chain}_DG{pocket.dg:.0f}"
        lines.append(f"select {label}, pockets and chain {chain}")
        lines.append(f"color {color}, {label}")
        lines.append(f"# {label}: N={len(pocket.coords)}{ki_str}")

    lines += [
        "",
        "zoom pockets",
        "set transparency, 0.5, receptor",
        "deselect",
    ]

    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
