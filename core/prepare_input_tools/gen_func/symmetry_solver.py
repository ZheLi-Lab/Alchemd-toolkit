import numpy as np
from rdkit.Chem import AllChem

def find_best_fitting_plane(points):
    centroid = np.mean(points, axis=0)
    centered = points - centroid
    cov_matrix = np.cov(centered, rowvar=False)
    eigen_values, eigen_vectors = np.linalg.eigh(cov_matrix)
    normal = eigen_vectors[:, 0]  # The eigenvector corresponding to the smallest eigenvalue
    # Adjust the normal vector to point towards the positive z-axis
    reference_direction = np.array([0, 0, 1])
    if np.dot(normal, reference_direction) < 0:
        normal = -normal
    # Construct an orthonormal basis
    if np.allclose(normal, [1, 0, 0], atol=1e-6):
        vec = np.array([0, 1, 0])
    else:
        vec = np.array([1, 0, 0])
    u = np.cross(normal, vec)
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    v /= np.linalg.norm(v)
    return centroid, normal, u, v

def project_points(points, centroid, u, v):
    proj_coords = []
    for point in points:
        vec = point - centroid
        x = np.dot(vec, u)
        y = np.dot(vec, v)
        proj_coords.append([x, y])
    return np.array(proj_coords)

def sort_ring_indices(proj_coords, ring_atom_indices, start_idx):
    centroid_2d = np.mean(proj_coords, axis=0)
    angles = []
    for x, y in proj_coords:
        dx = x - centroid_2d[0]
        dy = y - centroid_2d[1]
        angle = np.arctan2(dy, dx)
        angles.append(angle)
    angles = (np.array(angles) + 2 * np.pi) % (2 * np.pi)
    sorted_indices = np.argsort(-angles)  # Descending order for clockwise
    sorted_atom_indices = [ring_atom_indices[i] for i in sorted_indices]
    # Adjust the starting point
    try:
        start_pos = sorted_atom_indices.index(start_idx)
    except ValueError:
        raise ValueError(f"Starting atom index {start_idx} is not in the atom list of this ring")
    ordered = sorted_atom_indices[start_pos:] + sorted_atom_indices[:start_pos]
    return ordered

def process_rings(mol1, mol2, ring1_atoms, ring2_atoms, start_idx1, start_idx2):
    # Get coordinates of the first ring
    coords_ring1 = np.array([mol1.GetConformer().GetAtomPosition(idx) for idx in ring1_atoms])
    # Calculate the best fitting plane
    centroid, normal, u, v = find_best_fitting_plane(coords_ring1)
    # Get coordinates of both rings and project
    coords_ring1 = np.array([mol1.GetConformer().GetAtomPosition(idx) for idx in ring1_atoms])
    coords_ring2 = np.array([mol2.GetConformer().GetAtomPosition(idx) for idx in ring2_atoms])
    proj_ring1 = project_points(coords_ring1, centroid, u, v)
    proj_ring2 = project_points(coords_ring2, centroid, u, v)
    # Sort both rings
    sorted_ring1 = sort_ring_indices(proj_ring1, ring1_atoms, start_idx1)
    sorted_ring2 = sort_ring_indices(proj_ring2, ring2_atoms, start_idx2)
    return sorted_ring1, sorted_ring2

def symmetry_mapping_patch_in_ring(mol1, mol2, exist_mapping) -> dict:
        """
        Decide the symmetry mapping between mol1 and mol2
        Only consider rings where all atoms are in exist_mapping
        Return: patch mapping from mol1 to mol2
        """
        mol1_rings = mol1.GetRingInfo().AtomRings()
        mol2_rings = mol2.GetRingInfo().AtomRings()
        return_patch = {}

        mol1_to_check_rings = []
        mol2_to_check_rings = []
        for ring in mol1_rings:
            full_ring = True
            for idx in ring:
                if idx not in exist_mapping.keys():
                    full_ring = False
                    break
            if full_ring:
                mol1_to_check_rings.append(ring)

        for ring in mol2_rings:
            full_ring = True
            for idx in ring:
                if idx not in exist_mapping.values():
                    full_ring = False
                    break
            if full_ring:
                mol2_to_check_rings.append(ring)

        _ring_map = {}
        for ring1 in mol1_to_check_rings:
            for ring2 in mol2_to_check_rings:
                if len(set([exist_mapping[idx] for idx in ring1]) & set(ring2)) == len(ring1):
                    _ring_map[ring1] = ring2
                    break
        # Get the ring mapping

        for ring1, ring2 in _ring_map.items():
            using_ring1_idx = None
            is_odd_ring = len(ring1) % 2 == 1
            ring1_direct_connect_idx = []
            for idx1 in ring1:
                _direct_connect_heavy_idx = [i.GetIdx() for i in mol1.GetAtomWithIdx(idx1).GetNeighbors()
                                             if i.GetSymbol() != 'H' and (i.GetIdx() not in ring1) and (i.GetIdx() in exist_mapping.keys())]
                # Directly connected heavy atoms that are in the MCS but not in this ring
                if len(_direct_connect_heavy_idx) > 0:
                    # Junction between the scaffold and the ring
                    ring1_direct_connect_idx.extend(_direct_connect_heavy_idx)
                    using_ring1_idx = idx1
    
            if len(ring1_direct_connect_idx) == 0:
                # This only occurs when the entire MCS is a single ring
                continue
            if len(ring1_direct_connect_idx) != 1:
                # This function only handles the case of a single substituent
                continue
            
            using_ring2_idx = exist_mapping[using_ring1_idx]

            clockwise_ring1_order, clockwise_ring2_order = process_rings(mol1, mol2, ring1, ring2, using_ring1_idx, using_ring2_idx)
            all_match = True
            for i in range(len(clockwise_ring1_order)):
                if exist_mapping[clockwise_ring1_order[i]] != clockwise_ring2_order[i]:
                    all_match = False
                    break
            if all_match:
                continue

            for i, j in zip(clockwise_ring1_order, clockwise_ring2_order):
                return_patch[i] = j
                i_H = [ii for ii in mol1.GetAtomWithIdx(i).GetNeighbors() if ii.GetSymbol() == 'H']
                j_H = [jj for jj in mol2.GetAtomWithIdx(j).GetNeighbors() if jj.GetSymbol() == 'H']
                if len(i_H) > 0 and len(j_H) > 0:
                    return_patch[i_H[0].GetIdx()] = j_H[0].GetIdx()
        return return_patch

if __name__ == "__main__":
    pass
