# Surfcluster

Surface-aware MDMix hotspot clustering for druggability analysis.

Clusters MDMix hotspot points onto the receptor solvent-accessible surface using surface-graph distances, producing ranked binding pockets with pharmacophore composition, estimated ΔG, Ki, and volume.

## Installation and enviorment setting

```bash
git clone https://github.com/Tomasmunozc0/Surfcluster.git
cd Surfcluster
conda env create -f environment.yml
conda activate surfcluster
```

## Tutorial: ULK1

This walks through a full run using ULK1 as an example. You need two files from your MDMix project:
- `ulk1_ref.pdb` — the receptor structure
- `ulk1_hot_good.pdb` — the MDMix hotspot PDB

### 1. Create a config file

```bash
surfcluster --init ulk1
```

This writes `ulk1.yml` in the current directory. Edit it to point to your files:

```yaml
receptor: /path/to/ulk1_ref.pdb
hotspots:
  - /path/to/ulk1_hot_good.pdb
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
[surfcluster] Loading hotspots: [ulk1_hot_good.pdb]
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

Pockets A and B are flagged above. Tighten `energy_cutoff` in the YAML to only use the most favorable hotspots:

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

No flags. The top pocket (chain A, DG=-8.86, Ki=316 nM) is a well-defined site with 7 hotspots and a healthy mix of HYD and POL pharmacophores.

### 5. Visualize in PyMOL

```bash
pymol -r ulk1_session.pml
```

Each pocket is shown as spheres on the receptor surface, colored by chain (A=red, B=blue, C=green, ...). The YAML `force_merge` parameter lets you merge pockets after inspecting them:

```yaml
force_merge: [A, B]   # merge chains A and B into one pocket
```

### 6. Output files

| File | Description |
|------|-------------|
| `ulk1_clusts.pdb` | Pockets as PDB chains A, B, C… |
| `ulk1_summary.csv` | ΔG, Ki, efficiency, volume, composition per pocket |
| `ulk1_session.pml` | PyMOL script |

---

## Batch mode

Run multiple proteins from a single YAML:

```yaml
outdir: results/
proteins:
  - name: ulk1
    receptor: /path/ulk1_ref.pdb
    hotspots: [/path/ulk1_hot_good.pdb]
    energy_cutoff: -1.1
  - name: me1
    receptor: /path/me1_ref.pdb
    hotspots: [/path/me1_hot_good.pdb]
    energy_cutoff: -1.0
```

```bash
surfcluster batch.yml
```

Produces one set of output files per protein under `results/`.

---

## All parameters

```
surfcluster [-h] [-r PDB] [-p PDB [PDB ...]] [-n NAME] [-o DIR]
            [--energy-cutoff KCAL] [--anchor-cutoff KCAL]
            [--neigh-cut ANG] [--merge-cut ANG]
            [--min-hotspots N] [--max-hotspots N]
            [--force-merge CHAINS]
            [--roi-center X Y Z] [--roi-residues CHAIN:RESNUM,...]
            [--cache PATH] [--force] [--sgraph]
            [config.yml]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `energy_cutoff` | auto | Keep hotspots with energy ≤ this. Auto keeps top 200 if not set |
| `anchor_cutoff` | -0.8 | Discard pockets with no hotspot below this |
| `neigh_cut` | 6.0 Å | Euclidean radius for pocket expansion |
| `merge_cut` | 4.0 Å | Merge adjacent pockets within this surface distance |
| `min_hotspots` | 3 | Minimum hotspots per pocket |
| `max_hotspots` | off | Cap pocket size |
| `force_merge` | null | Merge specific chains after viewing in PyMOL |
| `cache` | null | Save/load surface distance matrix to speed up reruns |
