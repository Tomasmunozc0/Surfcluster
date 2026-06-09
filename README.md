# Surfcluster

Surface-aware MDMix hotspot clustering for ligandability analysis.

Clusters MDMix hotspot points onto the receptor solvent-accessible surface using surface-graph distances, producing ranked binding pockets with pharmacophore composition, estimated ΔG, Ki, and volume.

<img width="1774" height="887" alt="image" src="https://github.com/user-attachments/assets/c446f237-cf4b-46aa-ad40-20056f6e9b73" />

---

## Installation and environment setup

```bash
git clone https://github.com/Tomasmunozc0/Surfcluster.git
cd Surfcluster
conda env create -f environment.yml
conda activate surfcluster
```

---

## Tutorial: ULK1

This walks through a complete run using ULK1 as an example. You need two files from your MDMix project:
- `ulk1_ref.pdb` — the receptor structure
- `ulk1_hotspots.pdb` — the MDMix hotspot PDB

### 1. Create a config file

```bash
surfcluster --init ulk1
```

This writes `ulk1.yml` in the current directory. Edit it to point to your files:

```yaml
receptor: /path/to/ulk1_ref.pdb
hotspots:
  - /path/to/ulk1_hotspots.pdb
name: ulk1
outdir: .
energy_cutoff: null
```

### 2. Run

```bash
surfcluster ulk1.yml
```

Output:
```
[surfcluster] Loading receptor: ulk1_ref.pdb
  4397 atoms
[surfcluster] Loading hotspots: [ulk1_hotspots.pdb]
  166 hotspots, energy range [-1.57, -0.84] kcal/mol
[surfcluster] Building surface (50 pts/atom, edge=3.0 Å)...
  2840 surface points, 8574 edges (0.4s)
[surfcluster] Computing surface distance matrix...
  Done in 0.0s

  Found 20 pockets:
    #        DG       Ki(nM)     Eff   Vol(Å³)    N  HYD POL DON ACC
    1    -24.19        0.000  -1.100     731.4   22   11  11   0   0  <- !
    2    -19.81        0.000  -1.101     416.9   18    7  11   0   0  <- !
    3     -9.35      138.289  -1.039      87.0    9    3   6   0   0
    ...

  -- Tuning --
  [!] 2 pocket(s) flagged (DG < -15). Possible over-expansion:
      Chain A N=22 DG=-24.2: try --max-hotspots 13  or  --energy-cutoff -0.6
      Chain B N=18 DG=-19.8: try --max-hotspots 10  or  --energy-cutoff -0.6
  --energy-cutoff  -0.50  (tighten → -0.6 fewer/smaller, loosen → -0.4 more/larger)
```

### 3. Read the output table

| Column | Meaning |
|--------|---------|
| `#` | Pocket rank (1 = most favorable) |
| `DG` | Sum of hotspot energies (kcal/mol). More negative = more druggable |
| `Ki(nM)` | Estimated inhibitor potency derived from ΔG |
| `Eff` | Efficiency = DG / N (energy per hotspot) |
| `Vol(Å³)` | Estimated pocket volume |
| `N` | Number of hotspots in the pocket |
| `HYD POL DON ACC` | Pharmacophore composition (hydrophobic, polar, donor, acceptor) |

A pocket flagged with `← !` has DG < -15 kcal/mol, which usually means over-expansion — too many hotspots were merged into one pocket.

### 4. Tune when pockets are flagged

Pockets A and B are flagged. Tighten `energy_cutoff` in the YAML to only use the most favorable hotspots:

```yaml
energy_cutoff: -1.1   # was null
```

Re-run:

```bash
surfcluster ulk1.yml
```

```
  Found 5 pockets:
    #        DG       Ki(nM)     Eff   Vol(Å³)    N  HYD POL DON ACC
    1     -8.86      316.413  -1.266      75.9    7    4   3   0   0
    2     -3.96         >1e6  -1.320      16.3    3    2   1   0   0
    ...
```

No flags. The top pocket (chain A, DG=-8.86, Ki=316 nM) is the ATP site — 7 hotspots, mix of HYD and POL pharmacophores.

### 5. Visualize in PyMOL

```bash
pymol -r ulk1_session.pml
```

Each pocket is shown as spheres on the receptor surface, colored by chain (A=red, B=blue, C=green, ...).

To merge two pockets after inspecting them in PyMOL, add to the YAML and re-run:

```yaml
force_merge: [A, B]   # merge chains A and B into one pocket
```

> **Note:** after a force merge, pockets are re-sorted by DG — verify the new chain labels in the output.

### 6. Output files

| File | Description |
|------|-------------|
| `ulk1_clusts.pdb` | Pockets as PDB chains A, B, C… |
| `ulk1_summary.csv` | ΔG, Ki, efficiency, volume, composition per pocket |
| `ulk1_session.pml` | PyMOL script — `pymol -r ulk1_session.pml` |

---

## Region of Interest (ROI)

Focus the analysis on a specific part of the protein. All three modes filter the receptor and hotspots to a sphere before surface generation — smaller input, faster run, cleaner results.

### Option 1 — Manual sphere (you know the coordinates)

```yaml
roi_center: [-5.5, -39.3, -15.1]   # X Y Z in Å — read from PyMOL
roi_radius: 12.0
```

### Option 2 — Around specific residues (you know the active site)

```yaml
roi_residues: A:45,A:46,A:168      # chain:resnum, comma-separated
roi_buffer: 5.0                     # Å added around the residue extent
```

### Option 3 — Around a reference ligand (you have a co-crystal structure)

Extract the ligand from your co-crystal PDB as a separate file, then:

```yaml
roi_ligand: ligand.pdb   # center and radius computed automatically
roi_buffer: 5.0          # Å added around the ligand extent
```

Surfcluster computes the ligand's geometric center and sets the radius to the ligand's maximum atom distance from center plus buffer. No need to look up coordinates manually.

**ULK1 example** (inhibitor 3RF from PDB 4WNO):
```
  ROI ligand: ligand.pdb, buffer=5.0 Å
  ROI sphere: center=[-5.47 -39.32 -15.14], radius=12.35 Å
  After ROI filter: 596 / 4397 atoms, 23 / 166 hotspots in region
  Found 1 pocket: DG=-8.86, Ki=316 nM  ← ATP site isolated cleanly
```

> Use ONE option per run. Set the others to null. Tip: change `name` to e.g. `ulk1_roi` so results don't overwrite the full-protein run.

---

## Batch mode

Run multiple proteins from a single YAML:

```yaml
outdir: results/
proteins:
  - name: ulk1
    receptor: /path/ulk1_ref.pdb
    hotspots: [/path/ulk1_hotspots.pdb]
    energy_cutoff: -1.1
  - name: me1
    receptor: /path/me1_ref.pdb
    hotspots: [/path/me1_hotspots.pdb]
    energy_cutoff: -1.0
```

```bash
surfcluster batch.yml
```

Produces one set of output files per protein under `results/`. A summary table is printed at the end.

---

## Caching

For large proteins, surface generation is the slow step. Save the distance matrix with:

```yaml
cache: ulk1_dm.pkl
```

On the first run it is computed and saved. Every subsequent run loads it from disk — reruns for parameter tuning take seconds instead of minutes. The cache is automatically invalidated if the receptor file, hotspot files, or any surface parameter changes.

---

## All parameters

```
surfcluster [-h] [-r PDB] [-p PDB [PDB ...]] [-n NAME] [-o DIR]
            [--energy-cutoff KCAL] [--anchor-cutoff KCAL]
            [--neigh-cut ANG] [--merge-cut ANG]
            [--min-hotspots N] [--max-hotspots N] [--auto-target N]
            [--force-merge CHAINS]
            [--roi-center X Y Z] [--roi-residues CHAIN:RESNUM,...]
            [--roi-ligand PDB] [--roi-radius ANG] [--roi-buffer ANG]
            [--cache PATH] [--force] [--sgraph] [-v]
            [config.yml]
```

### Clustering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `energy_cutoff` | auto | Keep hotspots with energy ≤ this. Auto keeps top `auto_target` if not set |
| `auto_target` | 200 | Number of hotspots kept by auto energy cutoff |
| `anchor_cutoff` | -0.8 | Discard pockets with no hotspot below this |
| `neigh_cut` | 6.0 Å | Euclidean radius for pocket expansion — increase for large binding sites |
| `merge_cut` | 4.0 Å | Merge adjacent pockets within this surface distance — increase for kinases |
| `min_hotspots` | 3 | Minimum hotspots per pocket |
| `max_hotspots` | off | Cap pocket size — trim to N lowest-energy hotspots |
| `force_merge` | null | Merge chains after viewing in PyMOL, e.g. `[B,C]` or `[B,C, E,F,G]` |

### Region of interest

| Parameter | Default | Description |
|-----------|---------|-------------|
| `roi_center` | null | `[X, Y, Z]` — manual sphere center |
| `roi_residues` | null | `A:45,A:46` — sphere around residues |
| `roi_ligand` | null | Path to ligand PDB — sphere centered on ligand |
| `roi_radius` | 15.0 Å | Sphere radius for `roi_center` |
| `roi_buffer` | 5.0 Å | Buffer added around `roi_residues` or `roi_ligand` |

### Surface

| Parameter | Default | Description |
|-----------|---------|-------------|
| `probe` | 1.2 Å | Probe radius for surface generation |
| `n_sphere` | 50 | Points per atom for surface sampling |
| `edge_dist` | 3.0 Å | Max edge distance in surface graph |

### Misc

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cache` | null | Path to save/load surface distance matrix |
| `--force` | off | Ignore existing cache and recompute |
| `--sgraph` | off | Write surface graph nodes to `<name>_sgraph.xyz` |
| `--init NAME` | — | Write a template YAML and exit |
| `--version` | — | Print version |
