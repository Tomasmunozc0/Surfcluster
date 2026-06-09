"""Smoke test: run all 6 focus proteins and verify output sanity."""
import subprocess, sys, csv
from pathlib import Path

# Tuple: (name, receptor, hotspots, extra_cli_args)
# extra_cli_args: per-protein tuning flags needed to get sensible results with defaults
PROTEINS = [
    ('ston2',   '/babel/tmunoz/ONR/STON2/mdmix_project/ston2_ref.pdb',     '/babel/tmunoz/ONR/STON2/mdmix_project/ston2_hot_good.pdb',   []),
    ('me1',     '/babel/tmunoz/ONR/ME1/mdmix_project/me1_ref.pdb',          '/babel/tmunoz/ONR/ME1/mdmix_project/me1_hot_good.pdb',        ['--energy-cutoff', '-1.0']),
    ('ppp2r2b', '/babel/tmunoz/ONR/PPP2R2B/mdmix_project/ppp2r2b_ref.pdb', '/babel/tmunoz/ONR/PPP2R2B/mdmix_project/2r2b_hot_good.pdb',  []),
    ('ppp2r2d', '/babel/tmunoz/ONR/PPP2R2D/mdmix_project/ppp2r2d_ref.pdb', '/babel/tmunoz/ONR/PPP2R2D/mdmix_project/2r2d_hot_good.pdb',  ['--energy-cutoff', '-0.96']),
    ('ppp2r1b', '/babel/tmunoz/ONR/PPP2R1B/mdmix_project/ppp2r1b_ref.pdb', '/babel/tmunoz/ONR/PPP2R1B/mdmix_project/2r1b_hot_good.pdb',  []),
    ('ppp2r1a', '/babel/tmunoz/ONR/PPP2R1A/mdmix_project/ppp2r1a_ref.pdb', '/babel/tmunoz/ONR/PPP2R1A/mdmix_project/2r1a_hot_good.pdb',  []),
]

# Sanity bounds for the top pocket
DG_MIN = -15.0   # DG below this → over-expansion → fail
DG_MAX = -4.0    # DG above this → nothing useful found → fail
MIN_N  = 3       # minimum hotspots in top pocket

outdir = Path('/tmp/smoke_test')
outdir.mkdir(exist_ok=True)

passed = failed = 0
print(f"{'Protein':10s}  {'Status':8s}  {'Pockets':>7}  {'Top DG':>8}  {'N':>4}  Notes")
print('-' * 65)

for name, rec, hot, extra in PROTEINS:
    result = subprocess.run(
        ['surfcluster', '-r', rec, '-p', hot, '-n', name, '-o', str(outdir)] + extra,
        capture_output=True, text=True
    )

    clusts = outdir / f'{name}_clusts.pdb'
    csv_f  = outdir / f'{name}_summary.csv'
    pml    = outdir / f'{name}_session.pml'

    if result.returncode != 0:
        print(f"{name:10s}  FAILED    returncode={result.returncode}  "
              f"{result.stderr.strip()[:40]}")
        failed += 1
        continue

    missing = [f.name for f in [clusts, csv_f, pml] if not f.exists()]
    if missing:
        print(f"{name:10s}  FAILED    missing files: {missing}")
        failed += 1
        continue

    # Parse top pocket from stdout
    n_pockets = top_dg = top_n = None
    for line in result.stdout.splitlines():
        if 'Found' in line and 'pockets' in line:
            n_pockets = int(line.split()[1])
        if line.strip().startswith('1 '):
            parts = line.split()
            try:
                top_dg = float(parts[1])
                top_n  = int(parts[5])
            except (IndexError, ValueError):
                pass

    # Verify CSV row count matches pocket count
    issues = []
    with open(csv_f) as f:
        csv_rows = sum(1 for _ in csv.reader(f)) - 1
    if n_pockets is not None and csv_rows != n_pockets:
        issues.append(f"CSV rows={csv_rows} != pockets={n_pockets}")

    # Sanity-check top pocket
    if top_dg is None:
        issues.append("could not parse top DG from stdout")
    else:
        if top_dg < DG_MIN:
            issues.append(f"DG={top_dg:.1f} < {DG_MIN} (over-expansion)")
        if top_dg > DG_MAX:
            issues.append(f"DG={top_dg:.1f} > {DG_MAX} (nothing useful)")

    if top_n is not None and top_n < MIN_N:
        issues.append(f"N={top_n} < {MIN_N} (too few hotspots)")

    if issues:
        print(f"{name:10s}  FAILED    {'; '.join(issues)}")
        failed += 1
    else:
        print(f"{name:10s}  OK        {n_pockets:>7}  {top_dg:>8.1f}  {top_n:>4}")
        passed += 1

print('-' * 65)
print(f"  {passed} passed, {failed} failed")

import shutil; shutil.rmtree(outdir)
sys.exit(0 if failed == 0 else 1)
