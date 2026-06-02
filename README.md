# Surfcluster

Surface-aware MDMix hotspot clustering for druggability analysis.

Clusters MDMix hotspot points onto the receptor solvent-accessible surface using surface-graph distances, producing ranked binding pockets with pharmacophore composition, estimated ΔG, Ki, and volume.

## Install

```bash
conda activate surfcluster
pip install -e .
```

## Quick start

```bash
# Generate a template YAML
surfcluster --init myprotein

# Edit myprotein.yml, then run
surfcluster myprotein.yml
```

## Usage

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

YAML config (single protein):
```yaml
receptor: receptor.pdb
hotspots:
  - hotspot.pdb
name: myprotein
outdir: .

energy_cutoff: null     # auto (keeps top 200 hotspots)
anchor_cutoff: -0.8
neigh_cut: 6.0
merge_cut: 4.0
min_hotspots: 3
max_hotspots: null
force_merge: null       # e.g. [B,C] — merge after viewing in PyMOL
cache: null
```

Batch mode (multiple proteins from one YAML):
```yaml
outdir: results/
proteins:
  - name: ston2
    receptor: /path/ston2_ref.pdb
    hotspots: [/path/ston2_hot_good.pdb]
  - name: me1
    receptor: /path/me1_ref.pdb
    hotspots: [/path/me1_hot_good.pdb]
    energy_cutoff: -1.0
```

## Output files

| File | Description |
|------|-------------|
| `<name>_clusts.pdb` | Pockets as PDB chains A, B, C… |
| `<name>_summary.csv` | ΔG, Ki, efficiency, volume, composition per pocket |
| `<name>_session.pml` | PyMOL script — `pymol -r <name>_session.pml` |

## Tuning

Every run prints a tuning block with suggested parameter adjustments:

- **DG < -15 kcal/mol** → over-expansion; try tighter `--energy-cutoff` or `--max-hotspots`
- **Too few pockets** → loosen `--energy-cutoff` or lower `--anchor-cutoff`
- **Pockets too spread** → lower `--merge-cut`

## Environment

```bash
conda env create -f environment.yml
conda activate surfcluster
```
