import argparse
import os
import pickle
import sys
import time
import yaml
import numpy as np
from pathlib import Path

from surfcluster import __version__
from surfcluster.io import (
    read_receptor, read_hotspots,
    save_cache, load_cache,
    write_clusts_pdb, write_summary_csv,
    write_sgraph_xyz, write_pymol_session,
)
from surfcluster.surface import (
    build_surface_pipeline,
    snap_to_surface, compute_surface_dm, clean_surface_dm,
    filter_roi_sphere, filter_roi_residues,
)
from surfcluster.clustering import cluster_hotspots, force_merge_pockets
from surfcluster.descriptors import annotate_pockets


_EXAMPLES = """
examples:
  # generate a template YAML in the current directory
  surfcluster --init

  # single protein — YAML config (recommended)
  surfcluster myprotein.yml

  # single protein — CLI flags
  surfcluster -r receptor.pdb -p hotspot.pdb -n myprotein -o results/

  # override any parameter on the fly
  surfcluster myprotein.yml --energy-cutoff -0.9 --merge-cut 5.0

  # after viewing in PyMOL: merge chains B and C into one pocket
  surfcluster myprotein.yml --force-merge B,C

  # focus on a region of interest (ROI)
  surfcluster myprotein.yml --roi-center -30 19 19 --roi-radius 15

  # speed up reruns by caching the surface distance matrix
  surfcluster myprotein.yml --cache myprotein_dm.pkl

  # batch — run multiple proteins from one YAML
  surfcluster batch.yml

YAML format (single protein):
  receptor: receptor.pdb
  hotspots:
    - hotspot.pdb
  name: myprotein
  outdir: .                   # output directory (. = same dir as YAML)

  # clustering (all optional — defaults shown)
  energy_cutoff: null         # auto-selected (keeps top 200 hotspots) if not set
  anchor_cutoff: -0.8         # discard pockets with no hotspot below this
  neigh_cut: 6.0              # Euclidean radius for pocket expansion (Å)
  merge_cut: 4.0              # merge adjacent pockets within this surface distance (Å)
  min_hotspots: 3             # minimum hotspots per pocket
  max_hotspots: null          # cap pocket size (no limit if not set)
  force_merge: null           # e.g. [B,C] or [B,C, E,F,G] — merge after viewing in PyMOL
  cache: null                 # path to save/load surface DM (speeds up reruns)

  # region of interest (optional — leave null to use full protein)
  # Tip: change name above (e.g. myprotein_roi) to avoid overwriting full-protein results
  roi_center: null            # [X, Y, Z] sphere center in Å
  roi_radius: 15.0            # sphere radius in Å
  roi_residues: null          # e.g. A:45,A:46 — focus around specific residues
  roi_buffer: 5.0             # buffer around roi_residues in Å

YAML format (batch):
  outdir: results/            # global default for all proteins
  proteins:
    - name: ston2
      receptor: /path/ston2_ref.pdb
      hotspots: [/path/ston2_hot_good.pdb]
    - name: me1
      receptor: /path/me1_ref.pdb
      hotspots: [/path/me1_hot_good.pdb]
      max_hotspots: 20        # per-protein override

output files (per run):
  <name>_clusts.pdb     pockets as PDB chains A, B, C...
  <name>_summary.csv    DG, Ki, efficiency, volume, composition per pocket
  <name>_session.pml    PyMOL script — run with: pymol -r <name>_session.pml
  <name>_sgraph.xyz     surface graph (only with --sgraph flag)
"""


def build_parser():
    p = argparse.ArgumentParser(
        prog='surfcluster',
        description='Surface-aware MDMix hotspot clustering.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_EXAMPLES,
    )

    p.add_argument('config', nargs='?', metavar='config.yml',
                   help='YAML config file (single protein or batch). '
                        'CLI flags override YAML values.')

    # --- Required (can also come from YAML) ---
    req = p.add_argument_group('input (required unless provided in YAML)')
    req.add_argument('-r', '--receptor', metavar='PDB',
                     help='Receptor PDB file.')
    req.add_argument('-p', '--hotspots', nargs='+', metavar='PDB',
                     help='MDMix hotspot PDB file(s). Multiple files are merged.')
    req.add_argument('-n', '--name', metavar='NAME',
                     help='Project name (prefix for all output files).')

    # --- Output ---
    out = p.add_argument_group('output')
    out.add_argument('-o', '--outdir', default=None, metavar='DIR',
                     help='Output directory. (default: .)')
    out.add_argument('-m', '--mol2', metavar='MOL2',
                     help='Receptor mol2 for SYBYL atom-type VdW radii (more accurate).')

    # --- Clustering ---
    cl = p.add_argument_group('clustering parameters')
    cl.add_argument('--energy-cutoff', type=float, default=None, metavar='KCAL',
                    help='Pre-filter: keep hotspots with energy ≤ cutoff. '
                         'Auto-selected (keeps top 200) if not set.')
    cl.add_argument('--anchor-cutoff', type=float, default=None, metavar='KCAL',
                    help='Post-filter: discard pockets with no hotspot below this. '
                         '(default: -0.8)')
    cl.add_argument('--neigh-cut', type=float, default=None, metavar='ANG',
                    help='Euclidean radius for two-hop pocket expansion (Å). '
                         '(default: 6.0 — increase for large binding sites)')
    cl.add_argument('--merge-cut', type=float, default=None, metavar='ANG',
                    help='Merge adjacent pockets within this surface distance (Å). '
                         '(default: 4.0 — increase to 5-6 for kinases)')
    cl.add_argument('--min-hotspots', type=int, default=None, metavar='N',
                    help='Minimum hotspots per pocket. (default: 3)')
    cl.add_argument('--max-hotspots', type=int, default=None, metavar='N',
                    help='Cap pocket size: skips merges that would exceed N, '
                         'then trims. (default: no limit)')
    cl.add_argument('--auto-target', type=int, default=None, metavar='N',
                    help='Auto energy-cutoff target: keep top N hotspots when '
                         'more than N pass the -0.5 kcal/mol threshold. (default: 200)')
    cl.add_argument('--force-merge', action='append', metavar='CHAINS',
                    help='Force-merge chains after viewing in PyMOL. '
                         'E.g. --force-merge B,C  (repeatable)')

    # --- ROI ---
    roi_g = p.add_argument_group('region of interest (ROI)')
    roi = roi_g.add_mutually_exclusive_group()
    roi.add_argument('--roi-center', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                     help='Restrict to sphere of --roi-radius around X Y Z.')
    roi.add_argument('--roi-residues', metavar='CHAIN:RESNUM,...',
                     help='Restrict around residues. Format: A:45,A:46,B:100')
    roi_g.add_argument('--roi-radius', type=float, default=None,
                       help='Sphere radius for --roi-center. (default: 15.0 Å)')
    roi_g.add_argument('--roi-buffer', type=float, default=None,
                       help='Buffer around --roi-residues. (default: 5.0 Å)')

    # --- Surface ---
    surf = p.add_argument_group('surface generation')
    surf.add_argument('--probe', type=float, default=None,
                      help='Probe radius (Å). (default: 1.2)')
    surf.add_argument('--edge-dist', type=float, default=None,
                      help='Surface graph max edge distance (Å). (default: 3.0)')
    surf.add_argument('--n-sphere', type=int, default=None,
                      help='Sphere points per atom for surface sampling. (default: 50)')

    # --- Cache / misc ---
    misc = p.add_argument_group('cache / misc')
    misc.add_argument('--cache', metavar='PATH',
                      help='Save/load surface DM pickle (speeds up reruns).')
    misc.add_argument('--force', action='store_true',
                      help='Ignore existing cache and recompute.')
    misc.add_argument('--sgraph', action='store_true',
                      help='Write surface graph points to <name>_sgraph.xyz.')
    misc.add_argument('--init', nargs='?', const='protein', metavar='NAME',
                      help='Write a template YAML to the current directory and exit. '
                           'Optional name sets the filename and name field. (default: protein)')
    misc.add_argument('--version', action='version',
                      version=f'surfcluster {__version__}')
    misc.add_argument('-v', '--verbose', action='store_true',
                      help='Verbose output.')

    return p


# --- Defaults applied after YAML/CLI merge -----------------------------------

_DEFAULTS = dict(
    outdir='.',
    mol2=None,
    anchor_cutoff=-0.8,
    neigh_cut=6.0,
    merge_cut=4.0,
    min_hotspots=3,
    max_hotspots=None,
    auto_target=200,
    force_merge=None,
    roi_center=None,
    roi_residues=None,
    roi_radius=15.0,
    roi_buffer=5.0,
    probe=1.2,
    edge_dist=3.0,
    n_sphere=50,
    cache=None,
    force=False,
    sgraph=False,
    verbose=False,
    energy_cutoff=None,
)


def _load_yaml(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


_PATH_KEYS = {'receptor', 'hotspots', 'outdir', 'mol2', 'cache'}


def _yaml_to_opts(yaml_data: dict, config_dir: Path, cli_opts) -> argparse.Namespace:
    """Merge YAML values and CLI overrides into a single Namespace.
    CLI values (non-None) always win over YAML. YAML wins over _DEFAULTS."""

    def resolve_path(val):
        """Resolve a relative path string against the YAML file's directory."""
        if isinstance(val, str):
            p = Path(val)
            if not p.is_absolute():
                return str(config_dir / p)
            return val
        if isinstance(val, list):
            return [resolve_path(v) for v in val]
        return val

    merged = dict(_DEFAULTS)

    # Apply YAML values (resolve paths only for known path fields)
    for key, val in yaml_data.items():
        if key == 'proteins':
            continue
        merged[key] = resolve_path(val) if key in _PATH_KEYS else val

    # Apply CLI overrides (argparse stores None for unset flags, False for unset store_true)
    cli_dict = vars(cli_opts)
    for key, val in cli_dict.items():
        if key in ('config', 'init'):
            continue
        if val is not None and val is not False:
            merged[key] = val

    return argparse.Namespace(**merged)


def _cache_meta(opts, hotspot_list: list) -> dict:
    """Build the metadata fingerprint identifying a cache's inputs."""
    def _file_meta(p):
        try:
            return {'path': os.path.abspath(p), 'mtime': os.path.getmtime(p)}
        except OSError:
            return {'path': os.path.abspath(p), 'mtime': None}

    sorted_hot = sorted(hotspot_list)
    return {
        'version':      __version__,
        'receptor':     _file_meta(opts.receptor),
        'hotspots':     [_file_meta(h) for h in sorted_hot],
        'probe':        float(opts.probe),
        'n_sphere':     int(opts.n_sphere),
        'edge_dist':    float(opts.edge_dist),
        'roi_center':   tuple(map(float, opts.roi_center)) if opts.roi_center else None,
        'roi_radius':   float(opts.roi_radius) if opts.roi_center or opts.roi_residues else None,
        'roi_residues': opts.roi_residues,
        'roi_buffer':   float(opts.roi_buffer) if opts.roi_residues else None,
    }


def _meta_mismatch(saved: dict, current: dict) -> list:
    """Return list of fields that differ between saved and current cache meta."""
    keys = set(saved) | set(current)
    return sorted(k for k in keys if saved.get(k) != current.get(k))


def _parse_residues(s: str) -> list:
    pairs = []
    for token in s.split(','):
        token = token.strip()
        if ':' not in token:
            raise ValueError(f"Invalid residue spec '{token}'. Expected CHAIN:RESNUM.")
        chain, resnum = token.split(':', 1)
        pairs.append((chain.strip(), int(resnum.strip())))
    return pairs


# --- Core single-protein run -------------------------------------------------

def run_single(opts: argparse.Namespace):
    """Run clustering for one protein given a fully resolved Namespace."""
    import string as _str
    _alpha = _str.ascii_uppercase

    outdir = Path(opts.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    name = opts.name

    t0 = time.time()

    # --- Load receptor ---
    print(f"[surfcluster] Loading receptor: {opts.receptor}")
    receptor_coords, vdw_radii = read_receptor(opts.receptor, opts.mol2)
    print(f"  {len(receptor_coords)} atoms")

    # --- ROI filtering ---
    if opts.roi_center and opts.roi_residues:
        print("  WARNING: both roi_center and roi_residues specified — using roi_center.")

    roi_center = roi_radius = None
    if opts.roi_center:
        roi_center = np.array(opts.roi_center)
        roi_radius = opts.roi_radius
        print(f"  ROI sphere: center={roi_center}, radius={roi_radius} Å")
        receptor_coords, vdw_radii = filter_roi_sphere(
            receptor_coords, vdw_radii, roi_center, roi_radius)
        print(f"  After ROI filter: {len(receptor_coords)} atoms")
    elif opts.roi_residues:
        res_list = _parse_residues(opts.roi_residues)
        print(f"  ROI residues: {res_list}, buffer={opts.roi_buffer} Å")
        receptor_coords, vdw_radii, roi_center, roi_radius = filter_roi_residues(
            receptor_coords, vdw_radii, opts.receptor, res_list, opts.roi_buffer)
        print(f"  After ROI filter: {len(receptor_coords)} atoms")

    # --- Load hotspots ---
    hotspot_list = opts.hotspots if isinstance(opts.hotspots, list) else [opts.hotspots]
    print(f"[surfcluster] Loading hotspots: {hotspot_list}")
    hotspots = read_hotspots(hotspot_list)
    print(f"  {len(hotspots.coords)} hotspots, "
          f"energy range [{hotspots.energies.min():.2f}, {hotspots.energies.max():.2f}] kcal/mol")

    # Filter hotspots to ROI region
    if roi_center is not None:
        dists = np.linalg.norm(hotspots.coords - roi_center, axis=1)
        roi_mask = dists <= roi_radius
        hotspots = type(hotspots)(
            coords=hotspots.coords[roi_mask],
            types=hotspots.types[roi_mask],
            energies=hotspots.energies[roi_mask],
        )
        print(f"  After ROI filter: {roi_mask.sum()} / {len(roi_mask)} hotspots in region")

    # --- Surface DM: cache or compute ---
    cache_path = opts.cache
    current_meta = _cache_meta(opts, hotspot_list) if cache_path else None

    cache_loaded = False
    if cache_path and Path(cache_path).exists() and not opts.force:
        try:
            surface_dm, removed_ids, saved_meta = load_cache(cache_path)
            mismatches = _meta_mismatch(saved_meta, current_meta)
            if mismatches:
                print(f"  WARNING: cache invalid (changed: {mismatches}). Rebuilding.")
            else:
                print(f"[surfcluster] Loading cached surface DM: {cache_path}")
                print("[surfcluster] Rebuilding surface graph for hotspot snapping...")
                surface_coords, graph, tree = build_surface_pipeline(
                    receptor_coords, vdw_radii,
                    probe=opts.probe, n_sphere=opts.n_sphere, edge_dist=opts.edge_dist)
                node_ids, valid_mask = snap_to_surface(hotspots.coords, tree)
                cache_loaded = True
        except (ValueError, KeyError, EOFError, pickle.UnpicklingError) as e:
            print(f"  WARNING: failed to load cache ({e}). Rebuilding.")

    if not cache_loaded:
        print(f"[surfcluster] Building surface ({opts.n_sphere} pts/atom, edge={opts.edge_dist} Å)...")
        t1 = time.time()
        surface_coords, graph, tree = build_surface_pipeline(
            receptor_coords, vdw_radii,
            probe=opts.probe, n_sphere=opts.n_sphere, edge_dist=opts.edge_dist)
        print(f"  {len(surface_coords)} surface points, "
              f"{graph.number_of_edges()} edges ({time.time()-t1:.1f}s)")

        print("[surfcluster] Snapping hotspots to surface...")
        node_ids, valid_mask = snap_to_surface(hotspots.coords, tree)
        print(f"  {valid_mask.sum()} / {len(hotspots.coords)} hotspots on surface")

        hs_n = node_ids[valid_mask]
        print("[surfcluster] Computing surface distance matrix...")
        t1 = time.time()
        surface_dm = compute_surface_dm(graph, hs_n)
        print(f"  Done in {time.time()-t1:.1f}s")

        surface_dm, removed_ids = clean_surface_dm(surface_dm)

        if cache_path:
            save_cache(surface_dm, removed_ids, current_meta, cache_path)
            print(f"  Cached DM → {cache_path}")

    # Apply valid_mask and removed_ids
    hs_c = hotspots.coords[valid_mask]
    hs_t = hotspots.types[valid_mask]
    hs_e = hotspots.energies[valid_mask]

    if len(removed_ids):
        keep = np.ones(len(hs_c), dtype=bool)
        keep[removed_ids] = False
        hs_c = hs_c[keep]; hs_t = hs_t[keep]; hs_e = hs_e[keep]

    print(f"  {len(hs_c)} hotspots after surface cleaning")

    # --- Auto energy_cutoff --------------------------------------------------
    energy_cutoff = opts.energy_cutoff
    auto_target = int(opts.auto_target)
    if energy_cutoff is None:
        n_pass = (hs_e <= -0.5).sum()
        if n_pass > auto_target:
            sorted_e = np.sort(hs_e)
            energy_cutoff = round(float(sorted_e[auto_target - 1]), 2)
            print(f"  Auto energy_cutoff: {energy_cutoff:.2f} kcal/mol "
                  f"(keeping top {auto_target}/{n_pass} hotspots by energy)")
        else:
            energy_cutoff = -0.5

    # --- Clustering ----------------------------------------------------------
    print(f"[surfcluster] Clustering (eps={opts.neigh_cut}, merge={opts.merge_cut}, "
          f"ecut={energy_cutoff}, acut={opts.anchor_cutoff})...")
    pockets = cluster_hotspots(
        hs_c, hs_t, hs_e, surface_dm,
        energy_cutoff=energy_cutoff,
        anchor_cutoff=opts.anchor_cutoff,
        epsilon=opts.neigh_cut,
        min_samples=opts.min_hotspots,
        merge_cut=opts.merge_cut,
        max_hotspots=opts.max_hotspots,
        verbose=opts.verbose,
    )

    if not pockets:
        print("No pockets found. Exiting.")
        return None

    # --- Force-merge ---------------------------------------------------------
    if opts.force_merge:
        groups = [fm.split(',') for fm in opts.force_merge]
        n_before = len(pockets)
        pockets = force_merge_pockets(pockets, groups)
        print(f"[surfcluster] Force-merge: {n_before} → {len(pockets)} pockets "
              f"(groups: {groups}) — pockets re-sorted by DG, check new chain labels")

    # --- Descriptors ---------------------------------------------------------
    pockets = annotate_pockets(pockets)

    # --- Print summary -------------------------------------------------------
    print(f"\n  Found {len(pockets)} pockets:")
    print(f"  {'#':>3}  {'DG':>8}  {'Ki(nM)':>12}  {'Eff':>6}  {'Vol(Å³)':>8}  {'N':>4}  HYD POL DON ACC")
    flagged = []
    for i, p in enumerate(pockets):
        ki_str = f"{p.ki:.3f}" if p.ki < 1e6 else ">1e6"
        warn = " ← !" if p.dg < -15 else ""
        if p.dg < -15:
            flagged.append((i, p))
        print(f"  {i+1:>3}  {p.dg:>8.2f}  {ki_str:>12}  {p.efficiency:>6.3f}  "
              f"{p.volume:>8.1f}  {len(p.coords):>4}  "
              f"{p.hyd_count:>3} {p.pol_count:>3} {p.don_count:>3} {p.acc_count:>3}{warn}")

    # --- Tuning block --------------------------------------------------------
    print()
    print("  -- Tuning --")
    if flagged:
        print(f"  [!] {len(flagged)} pocket(s) flagged (DG < -15). Possible over-expansion:")
        for idx, p in flagged:
            chain = _alpha[idx % 26]
            n = len(p.coords)
            mx_suggest = max(int(n * 0.6), 5)
            ec_suggest = round(energy_cutoff - 0.1, 2)
            print(f"      Chain {chain} N={n} DG={p.dg:.1f}: "
                  f"try --max-hotspots {mx_suggest}  or  --energy-cutoff {ec_suggest}")
    ec_tighter = round(energy_cutoff - 0.1, 2)
    ec_looser  = round(energy_cutoff + 0.1, 2)
    mx_val = str(opts.max_hotspots) if opts.max_hotspots else 'off'
    print(f"  --energy-cutoff  {energy_cutoff:.2f}  "
          f"(tighten → {ec_tighter} fewer/smaller,  loosen → {ec_looser} more/larger)")
    print(f"  --auto-target    {auto_target:<5}  (auto energy-cutoff keeps top N, e.g. --auto-target 150)")
    print(f"  --max-hotspots   {mx_val:<5}  (cap pocket size, e.g. --max-hotspots 15)")
    print(f"  --merge-cut      {opts.merge_cut:<5.1f}  "
          f"(lower → more separate,  higher → merge adjacent)")
    print(f"  --neigh-cut      {opts.neigh_cut:<5.1f}  "
          f"(lower → smaller pockets,  higher → larger sites)")
    print(f"  --force-merge    B,C   (merge chains after viewing in PyMOL)")

    # --- Write outputs -------------------------------------------------------
    clusts_path = str(outdir / f"{name}_clusts.pdb")
    csv_path    = str(outdir / f"{name}_summary.csv")
    pml_path    = str(outdir / f"{name}_session.pml")
    sgraph_path = str(outdir / f"{name}_sgraph.xyz") if opts.sgraph and surface_coords is not None else None

    write_clusts_pdb(pockets, clusts_path, name)
    write_summary_csv(pockets, csv_path, name)
    write_pymol_session(pockets, opts.receptor, pml_path, name)
    if sgraph_path:
        write_sgraph_xyz(surface_coords, sgraph_path)

    print(f"\n  Output: {clusts_path}")
    print(f"          {csv_path}")
    print(f"          {pml_path}")
    if sgraph_path:
        print(f"          {sgraph_path}")
    print(f"\n  Total time: {time.time()-t0:.1f}s")

    return pockets


# --- Batch mode --------------------------------------------------------------

def run_batch(yaml_data: dict, config_dir: Path, cli_opts: argparse.Namespace):
    proteins = yaml_data.get('proteins', [])
    if not proteins:
        print("ERROR: batch YAML has no 'proteins' list.")
        sys.exit(1)

    global_defaults = {k: v for k, v in yaml_data.items() if k != 'proteins'}
    results = []

    print(f"[surfcluster] Batch mode: {len(proteins)} protein(s)\n")

    for i, prot_cfg in enumerate(proteins):
        name = prot_cfg.get('name', f'protein_{i+1}')
        print(f"{'='*60}")
        print(f"  [{i+1}/{len(proteins)}] {name}")
        print(f"{'='*60}")

        # Merge: global defaults → protein config → CLI overrides
        merged_yaml = {**global_defaults, **prot_cfg}
        opts = _yaml_to_opts(merged_yaml, config_dir, cli_opts)

        if not opts.receptor or not opts.hotspots:
            print(f"  ERROR: missing receptor or hotspots for {name}. Skipping.")
            results.append((name, None))
            continue

        try:
            pockets = run_single(opts)
            results.append((name, pockets))
        except Exception as e:
            print(f"  ERROR processing {name}: {e}")
            results.append((name, None))
        print()

    # Batch summary
    print(f"{'='*60}")
    print(f"  Batch complete: {len(proteins)} protein(s)")
    print(f"  {'Protein':15s}  {'Pockets':>7}  {'Top DG':>8}  {'N':>4}  {'':2}")
    for name, pockets in results:
        if pockets:
            p0 = pockets[0]
            flag = ' !' if p0.dg < -15 else '  '
            print(f"  {name:15s}  {len(pockets):>7}  {p0.dg:>8.1f}  {len(p0.coords):>4}  {flag}")
        else:
            print(f"  {name:15s}  FAILED")


_INIT_TEMPLATE = """\
receptor: receptor.pdb
hotspots:
  - hotspot.pdb
name: myprotein
outdir: .

# Clustering parameters (defaults shown)
energy_cutoff: null       # auto-selected (keeps top 200) if not set
anchor_cutoff: -0.8       # discard pockets with no hotspot below this
neigh_cut: 6.0            # Euclidean radius for pocket expansion (Å)
merge_cut: 4.0            # merge adjacent pockets within this surface distance (Å)
min_hotspots: 3           # minimum hotspots per pocket
max_hotspots: null        # cap pocket size (no limit if not set)
auto_target: 200          # auto energy-cutoff keeps top N hotspots (when >N pass -0.5 threshold)
force_merge: null         # merge chains after viewing in PyMOL, e.g. [B,C] or [B,C, E,F,G]
cache: null               # path to save/load surface DM (speeds up reruns, e.g. myprotein_dm.pkl)

# Region of interest (optional — leave null to use full protein)
# Tip: change name above (e.g. myprotein_roi) to avoid overwriting full-protein results
roi_center: null          # [X, Y, Z] sphere center in Å
roi_radius: 15.0          # sphere radius in Å
roi_residues: null        # e.g. A:45,A:46 — focus around specific residues
roi_buffer: 5.0           # buffer around roi_residues in Å
"""


# --- Entry point -------------------------------------------------------------

def main():
    parser = build_parser()
    cli_opts = parser.parse_args()

    if cli_opts.init is not None:
        init_name = cli_opts.init
        out = Path(f"{init_name}.yml")
        if out.exists():
            print(f"ERROR: {out} already exists. Remove it or rename it first.")
            sys.exit(1)
        content = _INIT_TEMPLATE.replace('name: myprotein', f'name: {init_name}')
        out.write_text(content)
        print(f"Template written to {out}. Edit it and run: surfcluster {out}")
        return

    # Detect YAML config
    if cli_opts.config:
        cfg_path = Path(cli_opts.config)
        if not cfg_path.exists():
            print(f"ERROR: config file not found: {cfg_path}")
            sys.exit(1)
        yaml_data = _load_yaml(str(cfg_path))
        config_dir = cfg_path.parent.resolve()

        # Batch mode: YAML has a 'proteins' key
        if 'proteins' in yaml_data:
            run_batch(yaml_data, config_dir, cli_opts)
            return

        # Single protein from YAML
        opts = _yaml_to_opts(yaml_data, config_dir, cli_opts)
    else:
        # Pure CLI mode — apply defaults manually
        opts = argparse.Namespace(**_DEFAULTS)
        for key, val in vars(cli_opts).items():
            if key in ('config', 'init'):
                continue
            if val is not None:
                setattr(opts, key, val)

    # Validate required fields
    missing = [f for f in ('receptor', 'hotspots', 'name')
               if not getattr(opts, f, None)]
    if missing:
        parser.error(f"Missing required arguments: {', '.join('--' + m.replace('_','-') for m in missing)}")

    run_single(opts)


if __name__ == '__main__':
    main()
