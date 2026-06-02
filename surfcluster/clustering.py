import numpy as np
from surfcluster.io import Pocket
from scipy.spatial import cKDTree


# --- Greedy clustering -------------------------------------------------------

def cluster_hotspots(
    coords: np.ndarray,
    types: np.ndarray,
    energies: np.ndarray,
    surface_dm: np.ndarray,
    energy_cutoff: float = -0.5,
    anchor_cutoff: float = -0.8,
    epsilon: float = 6.0,
    min_samples: int = 3,
    merge_cut: float = 4.0,
    max_hotspots: int = None,
    verbose: bool = False,
) -> list:
    """
    Energy-driven greedy clustering with two-hop expansion.

    Pocket formation uses Euclidean distance (spatially compact pockets).
    Merge step uses surface distance (surface-topology-aware adjacency).

    energy_cutoff: pre-filter — discard hotspots weaker than this before clustering.
    anchor_cutoff: post-filter — discard pockets with no hotspot below this value.
    epsilon: Euclidean distance radius for two-hop expansion (Å).
    merge_cut: surface distance threshold for merging adjacent pockets (Å).
    max_hotspots: if set, skip merges that would exceed this size, then trim any
                  remaining oversized pocket to its max_hotspots lowest-energy hotspots.
    """

    # Pre-filter by energy cutoff
    pre_mask = energies <= energy_cutoff
    if verbose:
        print(f"  Pre-filter: {pre_mask.sum()} / {len(energies)} hotspots pass energy_cutoff={energy_cutoff}")

    coords     = coords[pre_mask]
    types      = types[pre_mask]
    energies   = energies[pre_mask]
    surface_dm = surface_dm[np.ix_(pre_mask, pre_mask)]

    n_total = len(coords)
    pool_mask = np.ones(n_total, dtype=bool)
    raw_pockets = []

    # KD-tree for fast Euclidean neighbor lookups during expansion
    pool_tree = cKDTree(coords)

    while pool_mask.sum() > 0:

        pool_ids = np.where(pool_mask)[0]
        pool_energies = energies[pool_ids]
        pool_coords   = coords[pool_ids]

        # Seed: lowest energy hotspot in remaining pool
        seed_local = np.argmin(pool_energies)
        seed_coord = pool_coords[seed_local]

        # Rebuild KD-tree for current pool
        tree = cKDTree(pool_coords)

        # Hop 1 — Euclidean neighbors within epsilon
        hop1_local = np.array(tree.query_ball_point(seed_coord, r=epsilon))

        # Seed has no neighbors — discard
        if len(hop1_local) <= 1:
            pool_mask[pool_ids[seed_local]] = False
            continue

        # Hop 2 — Euclidean neighbors of hop1 neighbors (still within epsilon of each)
        hop2_set = set(hop1_local.tolist())
        for idx in hop1_local:
            nbrs = tree.query_ball_point(pool_coords[idx], r=epsilon)
            hop2_set.update(nbrs)

        pocket_local = np.array(sorted(hop2_set))

        # Map back to original indices
        pocket_orig = pool_ids[pocket_local]

        raw_pockets.append(pocket_orig)
        pool_mask[pocket_orig] = False

    if verbose:
        print(f"  Greedy: {len(raw_pockets)} raw pockets")

    # --- Filter: min size and anchor energy ----------------------------------
    valid_pockets = []
    valid_orig_ids = []

    for orig_ids in raw_pockets:
        p_energies = energies[orig_ids]

        if len(orig_ids) < min_samples:
            continue
        if p_energies.min() > anchor_cutoff:
            continue

        valid_pockets.append(Pocket(
            coords=coords[orig_ids],
            types=types[orig_ids],
            energies=p_energies,
        ))
        valid_orig_ids.append(orig_ids)

    if verbose:
        print(f"  After filter: {len(valid_pockets)} pockets pass size+anchor check")

    if not valid_pockets:
        print("WARNING: No pockets found. Try loosening --energy-cutoff, --anchor-energy or --neigh-cut.")
        return []

    # --- Merge adjacent pockets (surface distance) ---------------------------
    n = len(valid_pockets)
    merge_pairs = []

    for i in range(n):
        for j in range(i + 1, n):
            idx_i = valid_orig_ids[i]
            idx_j = valid_orig_ids[j]
            # Skip merge if it would exceed max_hotspots
            if max_hotspots and len(idx_i) + len(idx_j) > max_hotspots:
                continue
            sub = surface_dm[np.ix_(idx_i, idx_j)]
            if sub.min() <= merge_cut:
                merge_pairs.append((i, j))

    # Union-find to resolve groups
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[py] = px

    for i, j in merge_pairs:
        union(i, j)

    groups = {}
    for i in range(n):
        root = find(i)
        groups.setdefault(root, []).append(i)

    merged_pockets = []
    for root, members in groups.items():
        if len(members) == 1:
            merged_pockets.append(valid_pockets[members[0]])
        else:
            merged_pockets.append(Pocket(
                coords=np.vstack([valid_pockets[m].coords for m in members]),
                types=np.concatenate([valid_pockets[m].types for m in members]),
                energies=np.concatenate([valid_pockets[m].energies for m in members]),
            ))

    if verbose:
        print(f"  After merging: {len(merged_pockets)} final pockets")

    # --- Trim oversized pockets to max_hotspots (keep lowest-energy core) ----
    if max_hotspots:
        trimmed = []
        for p in merged_pockets:
            if len(p.coords) > max_hotspots:
                order = np.argsort(p.energies)[:max_hotspots]
                trimmed.append(Pocket(
                    coords=p.coords[order],
                    types=p.types[order],
                    energies=p.energies[order],
                ))
            else:
                trimmed.append(p)
        merged_pockets = trimmed
        if verbose:
            print(f"  After max_hotspots={max_hotspots} trim: sizes={[len(p.coords) for p in merged_pockets[:8]]}")

    merged_pockets.sort(key=lambda p: p.energies.sum())

    return merged_pockets


# --- Manual merge ------------------------------------------------------------

def force_merge_pockets(pockets: list, groups: list) -> list:
    """
    Merge specified pockets by chain letter (0-indexed A=0, B=1, ...).

    groups: list of lists of chain letters, e.g. [['B', 'C'], ['E', 'F', 'G']]
    Returns new pocket list sorted by DG.
    """
    import string
    alpha = string.ascii_uppercase
    idx_map = {ch: i for i, ch in enumerate(alpha[:len(pockets)])}

    # Collect which pocket indices go into which merge group
    merge_sets = []
    for grp in groups:
        indices = []
        for ch in grp:
            ch = ch.strip().upper()
            if ch not in idx_map:
                print(f"WARNING: Chain '{ch}' not found in output (only {len(pockets)} pockets). Skipping.")
                continue
            indices.append(idx_map[ch])
        if len(indices) >= 2:
            merge_sets.append(indices)

    # Union-find across all merge groups
    parent = list(range(len(pockets)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    def union(a, b):
        pa, pb = find(a), find(b)
        if pa != pb:
            parent[pb] = pa

    for indices in merge_sets:
        for k in range(1, len(indices)):
            union(indices[0], indices[k])

    groups_out = {}
    for i in range(len(pockets)):
        groups_out.setdefault(find(i), []).append(i)

    result = []
    for root, members in groups_out.items():
        if len(members) == 1:
            result.append(pockets[members[0]])
        else:
            result.append(Pocket(
                coords=np.vstack([pockets[m].coords for m in members]),
                types=np.concatenate([pockets[m].types for m in members]),
                energies=np.concatenate([pockets[m].energies for m in members]),
            ))

    result.sort(key=lambda p: p.energies.sum())
    return result
