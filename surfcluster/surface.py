import numpy as np
import networkx as nx
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix
import pickle


# --- Sphere point generation -------------------------------------------------

def _points_on_sphere(n: int) -> np.ndarray:
    """Fibonacci spiral — evenly distributed points on unit sphere."""
    n = float(n)
    pts = []
    inc = np.pi * (3.0 - np.sqrt(5.0))
    off = 2.0 / n
    for k in range(int(n)):
        y = k * off - 1.0 + (off / 2.0)
        r = np.sqrt(max(1.0 - y * y, 0.0))
        phi = k * inc
        pts.append([np.cos(phi) * r, y, np.sin(phi) * r])
    return np.array(pts)


# --- Surface point generation ------------------------------------------------

def generate_surface_points(
    coords: np.ndarray,
    vdw_radii: np.ndarray,
    probe: float = 1.2,
    n_sphere: int = 50,
) -> np.ndarray:
    """
    Vectorized Shrake-Rupley surface point generation.

    For each atom, places n_sphere points on its surface (VdW + probe radius)
    and keeps only those not occluded by any neighboring atom.

    Returns array of accessible surface point coordinates (M, 3).
    """
    radii = vdw_radii + probe
    sphere_pts = _points_on_sphere(n_sphere)

    # Neighbor lookup: only check atoms within 8 Angstrom
    tree = cKDTree(coords)
    accessible_points = []

    for i in range(len(coords)):
        r_i = radii[i]
        if r_i == probe:
            continue

        # Test points on this atom's surface
        test_pts = sphere_pts * r_i + coords[i]

        # Find nearby atoms (within max possible occlusion distance)
        neigh_ids = tree.query_ball_point(coords[i], r=8.0)
        neigh_ids = [j for j in neigh_ids if j != i]

        if not neigh_ids:
            accessible_points.append(test_pts)
            continue

        neigh_coords = coords[neigh_ids]
        neigh_radii  = radii[neigh_ids]

        # For each test point, check if inside any neighbor sphere
        # diff shape: (n_sphere, n_neighbors, 3)
        diff = test_pts[:, np.newaxis, :] - neigh_coords[np.newaxis, :, :]
        sq_dist = (diff ** 2).sum(axis=2)                 # (n_sphere, n_neighbors)
        sq_radii = (neigh_radii ** 2)[np.newaxis, :]      # (1, n_neighbors)

        # Point is accessible if outside ALL neighbor spheres
        accessible_mask = np.all(sq_dist > sq_radii, axis=1)

        if accessible_mask.any():
            accessible_points.append(test_pts[accessible_mask])

    if not accessible_points:
        return np.empty((0, 3))

    return np.vstack(accessible_points)


# --- Deduplication -----------------------------------------------------------

def deduplicate_points(surface_coords: np.ndarray, min_dist: float = 1.5) -> np.ndarray:
    """
    Remove surface points closer than min_dist to each other.
    Keeps the first encountered point in each cluster.
    """
    if len(surface_coords) == 0:
        return surface_coords

    tree = cKDTree(surface_coords)
    keep = np.ones(len(surface_coords), dtype=bool)

    for i in range(len(surface_coords)):
        if not keep[i]:
            continue
        neighbors = tree.query_ball_point(surface_coords[i], r=min_dist)
        for n in neighbors:
            if n != i:
                keep[n] = False

    return surface_coords[keep]


# --- Surface graph -----------------------------------------------------------

def build_surface_graph(
    surface_coords: np.ndarray,
    edge_dist: float = 2.5,
) -> tuple:
    """
    Build a networkx graph from surface points.
    Nodes are point indices, edges connect points within edge_dist.
    Edge weights are Euclidean distances.

    Returns (graph, kd_tree).
    """
    tree = cKDTree(surface_coords)
    g = nx.Graph()
    g.add_nodes_from(range(len(surface_coords)))

    pairs = tree.query_pairs(r=edge_dist)
    for i, j in pairs:
        d = np.linalg.norm(surface_coords[i] - surface_coords[j])
        g.add_edge(i, j, weight=d)

    # Isolated nodes are left disconnected — they produce inf surface distances
    # and are removed by clean_surface_dm. Bridging them would create artificial
    # long edges that silently distort distances across disjoint regions.

    return g, tree


# --- Hotspot snapping and surface distance matrix ----------------------------

def snap_to_surface(
    hotspot_coords: np.ndarray,
    surface_tree: cKDTree,
    max_dist: float = 4.0,
) -> tuple:
    """
    Snap each hotspot to its nearest surface node.

    Returns:
        node_ids    — surface node index for each hotspot (N,)
        valid_mask  — boolean mask of hotspots close enough to surface (N,)
    """
    dists, node_ids = surface_tree.query(hotspot_coords, k=1)
    valid_mask = dists <= max_dist
    return node_ids, valid_mask


def compute_surface_dm(
    graph: nx.Graph,
    node_ids: np.ndarray,
) -> np.ndarray:
    """
    Compute pairwise surface distances between hotspot-snapped surface nodes.

    Uses scipy dijkstra with all hotspot nodes as sources simultaneously.
    Returns NxN distance matrix (N = number of valid hotspots).
    """
    n_nodes = graph.number_of_nodes()
    unique_nodes = np.unique(node_ids)

    # Build sparse adjacency matrix
    rows, cols, weights = [], [], []
    for u, v, data in graph.edges(data=True):
        w = data.get('weight', 1.0)
        rows += [u, v]
        cols += [v, u]
        weights += [w, w]

    sparse = csr_matrix(
        (weights, (rows, cols)),
        shape=(n_nodes, n_nodes)
    )

    # Run dijkstra from all hotspot-snapped nodes at once
    dist_matrix = dijkstra(
        sparse,
        indices=unique_nodes,
        directed=False,
    )

    # Map back to hotspot indices — vectorized
    node_to_row = {node: i for i, node in enumerate(unique_nodes)}
    row_indices = np.array([node_to_row[n] for n in node_ids])
    surface_dm = dist_matrix[np.ix_(row_indices, node_ids)]

    return surface_dm


# --- ROI filtering -----------------------------------------------------------

def filter_roi_sphere(
    coords: np.ndarray,
    vdw_radii: np.ndarray,
    center: np.ndarray,
    radius: float,
) -> tuple:
    """
    Filter receptor atoms to those within radius of center.
    Returns filtered coords and vdw_radii.
    """
    dists = np.linalg.norm(coords - center, axis=1)
    mask = dists <= radius
    return coords[mask], vdw_radii[mask]


def filter_roi_residues(
    coords: np.ndarray,
    vdw_radii: np.ndarray,
    pdb_path: str,
    residue_list: list,
    buffer: float = 5.0,
) -> tuple:
    """
    Filter receptor atoms to those within buffer of specified residues.
    Returns (coords, vdw_radii, center, radius) so the same sphere can be
    applied to hotspots.
    """
    res_coords = []

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            chain  = line[21].strip()
            resnum = int(line[22:26].strip())
            if (chain, resnum) in residue_list:
                res_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])

    if not res_coords:
        raise ValueError(f"No atoms found for residues: {residue_list}")

    res_coords_arr = np.array(res_coords)
    center = np.mean(res_coords_arr, axis=0)
    radius = np.max(np.linalg.norm(res_coords_arr - center, axis=1)) + buffer

    mask = np.linalg.norm(coords - center, axis=1) <= radius
    return coords[mask], vdw_radii[mask], center, radius


def filter_roi_ligand(
    coords: np.ndarray,
    vdw_radii: np.ndarray,
    ligand_path: str,
    buffer: float = 5.0,
) -> tuple:
    """
    Filter receptor atoms to those within buffer of a reference ligand.
    Center is the geometric center of the ligand; radius is the max distance
    from center to any ligand atom plus buffer.
    Returns (coords, vdw_radii, center, radius) — same signature as
    filter_roi_residues so the same sphere can be applied to hotspots.
    """
    lig_coords = []

    with open(ligand_path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            try:
                lig_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue

    if not lig_coords:
        raise ValueError(f"No ATOM/HETATM records found in ligand file: {ligand_path}")

    lig_arr = np.array(lig_coords)
    center  = np.mean(lig_arr, axis=0)
    radius  = np.max(np.linalg.norm(lig_arr - center, axis=1)) + buffer

    mask = np.linalg.norm(coords - center, axis=1) <= radius
    return coords[mask], vdw_radii[mask], center, radius


# --- DM cleaning -------------------------------------------------------------

def clean_surface_dm(surface_dm: np.ndarray) -> tuple:
    """
    Remove hotspots that are disconnected from the main surface graph.
    A hotspot is disconnected if it has more inf values than the minimum
    observed across all hotspots — same logic as original surf_cluster.py.

    Returns (cleaned_dm, removed_ids).
    """
    inf_mask = np.isinf(surface_dm)
    inf_count = inf_mask.sum(axis=0)
    min_inf = inf_count.min()
    removed_ids = np.where(inf_count > min_inf)[0]

    if len(removed_ids) > 0:
        cleaned = np.delete(surface_dm, removed_ids, axis=0)
        cleaned = np.delete(cleaned, removed_ids, axis=1)
    else:
        cleaned = surface_dm

    return cleaned, removed_ids


# --- Full pipeline helper ----------------------------------------------------

def build_surface_pipeline(
    receptor_coords: np.ndarray,
    vdw_radii: np.ndarray,
    probe: float = 1.2,
    n_sphere: int = 50,
    edge_dist: float = 3.0,
    min_dedup_dist: float = 1.5,
) -> tuple:
    """
    Full surface pipeline: generate points → deduplicate → build graph.
    Returns (surface_coords, graph, kd_tree).
    """
    surface_coords = generate_surface_points(
        receptor_coords, vdw_radii, probe=probe, n_sphere=n_sphere
    )
    surface_coords = deduplicate_points(surface_coords, min_dist=min_dedup_dist)
    graph, tree = build_surface_graph(surface_coords, edge_dist=edge_dist)
    return surface_coords, graph, tree
