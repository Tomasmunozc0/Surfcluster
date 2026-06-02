import math
import numpy as np
from scipy.spatial import ConvexHull

RT = 0.592  # kcal/mol at 298 K


def _sphere_intersection_volume(c1, c2, r: float) -> float:
    """Volume of intersection between two spheres of radius r."""
    d = np.linalg.norm(np.asarray(c1) - np.asarray(c2))
    if d >= 2 * r:
        return 0.0
    return (math.pi / 12.0) * (4 * r + d) * (2 * r - d) ** 2


def pocket_volume(coords: np.ndarray, r: float = 1.2) -> float:
    """
    Pocket volume estimate (Å³).

    Combines sphere-packing volume with convex hull following estuf.py:
      vol = 0.75 * N * (4/3 π r³) - pairwise_overlaps + convex_hull_volume
    """
    n = len(coords)
    if n == 0:
        return 0.0

    sphere_vol = n * (4.0 / 3.0) * math.pi * r ** 3
    overlap = 0.0
    pts = [tuple(c) for c in coords]
    for i in range(n):
        for j in range(i):
            overlap += _sphere_intersection_volume(pts[i], pts[j], r)

    if n >= 4:
        try:
            hull_vol = ConvexHull(coords).volume
        except Exception:
            hull_vol = 0.0
    else:
        hull_vol = 0.0

    return sphere_vol * 0.75 - overlap + hull_vol


def dg_to_ki_nm(dg: float) -> float:
    """Convert ΔG (kcal/mol) to Ki (nM) using ΔG = RT ln(Ki)."""
    if dg >= 0:
        return float('inf')
    return math.exp(dg / RT) * 1e9


def annotate_pockets(pockets: list) -> list:
    """
    Compute and attach descriptors in-place for each Pocket:
      dg, ki, efficiency, volume, hyd/pol/don/acc counts.

    Returns the same list (pockets are modified in place).
    """
    for p in pockets:
        p.dg         = float(p.energies.sum())
        p.ki         = dg_to_ki_nm(p.dg)
        p.efficiency = p.dg / len(p.coords) if len(p.coords) > 0 else 0.0
        p.volume     = pocket_volume(p.coords)

        types_arr = np.asarray(p.types)
        p.hyd_count = int((types_arr == 'HYD').sum())
        p.pol_count = int((types_arr == 'POL').sum())
        p.don_count = int((types_arr == 'DON').sum())
        p.acc_count = int((types_arr == 'ACC').sum())

    return pockets
