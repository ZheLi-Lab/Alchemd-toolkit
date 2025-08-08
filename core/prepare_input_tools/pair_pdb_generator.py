import copy
import logging
import warnings
import tempfile
import re
import time

import numpy as np
from typing import List, Tuple, Union
from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdDepictor
from rdkit.Chem import AllChem
from rdkit.Geometry import rdGeometry
from rdkit.Chem.Draw import MolsToGridImage
import shutil
import subprocess
import os
from contextlib import contextmanager
try:
    from .pdbqt_generator import MolToPDBQTBlock, MolCoreToPDBQTBlock
    from .gen_func.symmetry_solver import symmetry_mapping_patch_in_ring
    from .gen_func.pdb_func import calculate_formal_charge, pdb_formatter, pdb_info_clean, protein_pdb_reader
except:
    from pdbqt_generator import MolToPDBQTBlock, MolCoreToPDBQTBlock
    from gen_func.symmetry_solver import symmetry_mapping_patch_in_ring
    from gen_func.pdb_func import calculate_formal_charge, pdb_formatter, pdb_info_clean, protein_pdb_reader


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# logger.setLevel(logging.DEBUG)

config_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'configs.toml')
WATVINA_PATH = None
with open(config_file, 'r') as f:
    for line in f:
        if line.startswith('Watvina_Path'):
            WATVINA_PATH = line.strip().split('=')[1].strip().strip('"')
            WATVINA_PATH = os.path.join(WATVINA_PATH, 'watvina')


class Executable_checker:
    def __init__(self, ):
        pass

    @staticmethod
    def is_executable_available(program_name):
        """Checks if the program is available in the system PATH and executable."""
        return shutil.which(program_name) is not None


@contextmanager
def working_dir(mkdir_,):
    pwd = os.getcwd()
    if os.path.exists(mkdir_) and os.path.isdir(mkdir_):
        pass
    else:
        os.makedirs(mkdir_)
    os.chdir(mkdir_)
    yield
    os.chdir(pwd)


class Mol3dTransformMethod:

    """
    A utility class that provides static methods for 3d coordinate transformations
    Using rdkit.Chem classes.
    """

    @staticmethod
    def get_atom_xyz(conformer: Chem.Conformer, atom_idx: int) -> Tuple[float, float, float]:
        p3d = conformer.GetAtomPosition(atom_idx)
        return p3d.x, p3d.y, p3d.z

    @staticmethod
    def get_vector(p1: rdGeometry.Point3D, p2: rdGeometry.Point3D) -> Tuple[float, float, float]:
        return p2.x - p1.x, p2.y - p1.y, p2.z - p1.z

    @staticmethod
    def _get_mol_cal_matrix(_mol: Union[Chem.Mol, Chem.RWMol], atom_idx_list: list) -> np.ndarray:
        _mol_conformer = _mol.GetConformer()
        tmp_x = []
        tmp_y = []
        tmp_z = []
        tmp_t = [1 for _ in range(len(atom_idx_list))]
        for atom_idx in atom_idx_list:
            x, y, z = Mol3dTransformMethod.get_atom_xyz(_mol_conformer, atom_idx)
            tmp_x.append(x)
            tmp_y.append(y)
            tmp_z.append(z)

        return np.array([tmp_x, tmp_y, tmp_z, tmp_t])

    @staticmethod
    def get_relative_matrix(_conformer: Chem.Conformer, ref_idx: int, other_idxes: Union[set, list]) -> np.ndarray:
        tmp_x = [0]
        tmp_y = [0]
        tmp_z = [0]
        tmp_t = [1 for _ in range(len(other_idxes)+1)]
        for idx in other_idxes:
            p_o = _conformer.GetAtomPosition(idx)
            p_r = _conformer.GetAtomPosition(ref_idx)
            x, y, z = Mol3dTransformMethod.get_vector(p_r, p_o)
            tmp_x.append(x)
            tmp_y.append(y)
            tmp_z.append(z)

        return np.array([tmp_x, tmp_y, tmp_z, tmp_t])

    @staticmethod
    def get_relative_matrix3(_conformer: Chem.Conformer, ref_idx: int, other_idxes: Union[set, list]) -> np.ndarray:
        tmp_x = []
        tmp_y = []
        tmp_z = []
        for idx in other_idxes:
            p_o = _conformer.GetAtomPosition(idx)
            p_r = _conformer.GetAtomPosition(ref_idx)
            x, y, z = Mol3dTransformMethod.get_vector(p_r, p_o)
            tmp_x.append(x)
            tmp_y.append(y)
            tmp_z.append(z)

        return np.array([tmp_x, tmp_y, tmp_z])

    @staticmethod
    def Kabsch_rotation_matrix(vectors_a: np.ndarray, vectors_b: np.ndarray) -> np.ndarray:
        U, s, Vt = np.linalg.svd(vectors_a.T.dot(vectors_b))
        R = U.dot(Vt).T
        if np.linalg.det(R) < 0:
            Vt[2, :] *= -1
            R = U.dot(Vt).T

        return R

    @staticmethod
    def Rodrigues_rotation_matrix(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """
        np.dot(R,a) = b
        :param a: vector A
        :param b: vector B
        :return: Rotation_matrix
        """
        # Normalize the vectors
        a = a / np.linalg.norm(a)
        b = b / np.linalg.norm(b)

        # Compute the cross product and sine of the angle
        v = np.cross(a, b, axis=0)
        sin_theta = np.linalg.norm(v)
        cos_theta = np.dot(a, b)

        # If the vectors are parallel or anti-parallel
        if sin_theta == 0:
            if cos_theta > 0:
                return np.eye(3)  # Same direction
            else:
                # 180-degree rotation: choose any perpendicular vector for axis
                perp_vec = np.array([1, 0, 0]) if np.abs(a[0]) < 0.9 else np.array([0, 1, 0])
                v = np.cross(a, perp_vec)
                v = v / np.linalg.norm(v)
                sin_theta = 1
                cos_theta = -1

        # Skew-symmetric cross-product matrix of v
        K = np.array([[0, -v[2], v[1]],
                      [v[2], 0, -v[0]],
                      [-v[1], v[0], 0]])

        # Rodrigues' rotation formula
        R = np.eye(3) + K * sin_theta + np.dot(K, K) * (1 - cos_theta)

        return R

    @staticmethod
    def get_transform_matrix(mol1: Union[Chem.Mol, Chem.RWMol], mol2: Union[Chem.Mol, Chem.RWMol],
                             using_atoms: dict) -> np.ndarray:
        """
        Return a homogeneous transform matrix from mol1 to mol2
        :param mol1:
        :param mol2:
        :param using_atoms: mapping atom id between mol1 and mol2, at least 4 len
        :return: a homogeneous transform matrix in np.ndarray
        """

        mol1_used_atom_idx = []
        mol2_used_atom_idx = []
        for mol1_atom_idx, mol2_atom_idx in using_atoms.items():
            mol1_used_atom_idx.append(mol1_atom_idx)
            mol2_used_atom_idx.append(mol2_atom_idx)

        mol1_cal_matrix = Mol3dTransformMethod._get_mol_cal_matrix(mol1, mol1_used_atom_idx)
        mol2_cal_matrix = Mol3dTransformMethod._get_mol_cal_matrix(mol2, mol2_used_atom_idx)

        return np.dot(mol2_cal_matrix, np.linalg.inv(mol1_cal_matrix))

    @staticmethod
    def get_transformed_xyz(xyz: Tuple[float, float, float], transformation_matrix: np.ndarray) -> Tuple[float, float, float]:
        """
        Return a new xyz coordinates after applying transformation.
        :param xyz:
        :param transformation_matrix:
        :return:
        """
        return_xyz_array = np.array([[xyz[0]], [xyz[1]], [xyz[2]], [1]])
        return return_xyz_array[0][0], return_xyz_array[1][0], return_xyz_array[2][0]


class UniGroupDecider:
    @dataclass(frozen=True)
    class uni_group:
        edge: int
        group: set

        def __eq__(self, other):
            return self.edge == other.edge and self.group == other.group

        def is_empty(self) -> bool:
            return len(self.group) == 0

        def __str__(self):
            return f"{self.edge}: {self.group}"

    @dataclass(frozen=True)
    class mol_with_uni_group:
        mol: Chem.Mol
        uni_groups: List['UniGroupDecider.uni_group']

        def __post_init__(self):
            object.__setattr__(self, 'edges', [i.edge for i in self.uni_groups])

    @staticmethod
    def _get_alchem_groups(mol, mcs_mol, match_idxes = None) -> List['UniGroupDecider.uni_group']:
        if match_idxes is None:
            match_idxes = mol.GetSubstructMatch(mcs_mol)

        alchem_atom_idxes = set()
        for atom_idx in range(mol.GetNumAtoms()):
            if atom_idx not in match_idxes:
                alchem_atom_idxes.add(atom_idx)

        alchem_heavy_atom_idxes = set()
        alchem_unique_groups = []
        for atom_idx in alchem_atom_idxes:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() != 'H':
                alchem_heavy_atom_idxes.add(atom_idx)
                boo = list(atom.GetNeighbors())
                boo = [i.GetIdx() for i in boo if i.GetSymbol() != 'H']
                boo.append(atom_idx)
                alchem_unique_groups.append(set(boo))

        _modded = True
        while _modded:
            _modded = False
            for i, group in enumerate(alchem_unique_groups):
                for j, other_group in enumerate(alchem_unique_groups[i + 1:]):
                    for atom_idx in group:
                        if atom_idx in other_group and atom_idx in alchem_heavy_atom_idxes:
                            _modded = True
                            g2 = alchem_unique_groups.pop(i + j + 1)
                            g1 = alchem_unique_groups.pop(i)
                            g3 = set.union(g1, g2)
                            alchem_unique_groups.append(g3)
                            break
                    if _modded:
                        break
                if _modded:
                    break
        # Merge the directly connected groups into one set

        group_list = []
        for group in alchem_unique_groups:
            tmp_in_set = set()
            for atom_idx in group:
                if atom_idx in alchem_heavy_atom_idxes:
                    tmp_in_set.add(atom_idx)
            edges = sorted(list(group - tmp_in_set))

            group_edge = edges[0]
            group_group = tmp_in_set

            tmp_uni_group = UniGroupDecider.uni_group(group_edge, group_group)
            group_list.append(tmp_uni_group)

        return group_list
    
    @staticmethod
    def is_edge_split_bone(mol: Chem.Mol, edge_idx: int, bone_idxes: list) -> bool:
        edge_bonds_with_bone = []
        edge_atom = mol.GetAtomWithIdx(edge_idx)
        for neighbor_atom in edge_atom.GetNeighbors():
            if neighbor_atom.GetIdx() in bone_idxes:
                edge_bonds_with_bone.append(mol.GetBondBetweenAtoms(edge_idx, neighbor_atom.GetIdx()).GetIdx())
        tmp_frags= Chem.FragmentOnBonds(mol, edge_bonds_with_bone, addDummies=False)
        tmp_frags_list=list(Chem.GetMolFrags(tmp_frags))
        return len(tmp_frags_list) > 1
    
    @staticmethod
    def get_unique_mapping(mol1: Chem.Mol, mol2: Chem.Mol, mcs_mol, rebuild_3d: bool = False) -> Tuple[list, list]:
        mol1_matches = mol1.GetSubstructMatches(mcs_mol)
        mol2_matches = mol2.GetSubstructMatches(mcs_mol)
        if len(mol1_matches) == 0 or len(mol2_matches) == 0 or mcs_mol.GetNumAtoms() == 0:
            raise ValueError("No matches found for mol1 and mol2 in MCS mol.")

        # To avoid a bug when mol1 and mol2 have different orientation and 3D coordinates far away.
        # We need to align the mol1 and mol2 before calculating the RMSD by regenerating their 3D conformer.
        tmp_mol1 = copy.deepcopy(mol1)
        tmp_mol2 = copy.deepcopy(mol2)
        if rebuild_3d and MolCoordModer.has_3d_conformer(tmp_mol1) and MolCoordModer.has_3d_conformer(tmp_mol2):
            tmp_mol1.RemoveAllConformers()
            tmp_mol2.RemoveAllConformers()
            tmp_mol1 = Chem.AddHs(tmp_mol1)
            tmp_mol2 = Chem.AddHs(tmp_mol2)
            AllChem.EmbedMolecule(tmp_mol1, randomSeed=42)
            AllChem.EmbedMolecule(tmp_mol2, randomSeed=42)
            rdForceFieldHelpers.MMFFOptimizeMolecule(tmp_mol1, mmffVariant='MMFF94s')
            rdForceFieldHelpers.MMFFOptimizeMolecule(tmp_mol2, mmffVariant='MMFF94s')

        mol1_conformer = tmp_mol1.GetConformer()
        mol2_conformer = tmp_mol2.GetConformer()

        rmsd_dict = {}

        for i, mol1_match in enumerate(mol1_matches):
            for j, mol2_match in enumerate(mol2_matches):
                mapping_dict = {mol1_match[atm_count]: mol2_match[atm_count] for atm_count in range(mcs_mol.GetNumAtoms())}
                distance_square_sum = 0
                for mol1_atom_idx, mol2_atom_idx in mapping_dict.items():
                    mol1_atom_position = mol1_conformer.GetAtomPosition(mol1_atom_idx)
                    mol2_atom_position = mol2_conformer.GetAtomPosition(mol2_atom_idx)
                    distance_square_sum += (mol1_atom_position.Distance(mol2_atom_position))** 2
                rmsd = np.sqrt(distance_square_sum / mcs_mol.GetNumAtoms())
                rmsd_dict[f"{i}-{j}"] = rmsd

        sorted_rmsd_dict = sorted(rmsd_dict.items(), key=lambda x: x[1])
        logger.debug(f"sorted_rmsd_dict: {sorted_rmsd_dict}")
        using_mol1_match = mol1_matches[int(sorted_rmsd_dict[0][0].split("-")[0])]
        using_mol2_match = mol2_matches[int(sorted_rmsd_dict[0][0].split("-")[1])]

        tmp_mapping = {i: j for i, j in zip(using_mol1_match, using_mol2_match)}
        mapping_patch = symmetry_mapping_patch_in_ring(mol1, mol2, tmp_mapping)
        modded_mol2_match = []
        if len(mapping_patch) > 0:
            logger.debug(f"have symmetry mapping patch: {mapping_patch}")
            for i, j in mapping_patch.items():
                tmp_mapping[i] = j
            for idx1 in using_mol1_match:
                modded_mol2_match.append(tmp_mapping[idx1])
            using_mol2_match = modded_mol2_match

        return using_mol1_match, using_mol2_match

    @staticmethod
    def find_uni_group(mol1: Chem.Mol, mol2: Chem.Mol, enforce_3d_distance: bool = True) -> Tuple[mol_with_uni_group, mol_with_uni_group]:
        """
        find uni_group in mol1 and mol2, base on MCS.
        :param mol1:
        :param mol2:
        :return: mol1 with uni_groups info and mol2 with uni_groups info
        """
        mcs = MolCoordModer.customized_FindMCS([mol1, mol2], enforce_3d_distance=enforce_3d_distance)

        if enforce_3d_distance:
            group_list1_match, group_list2_match = UniGroupDecider.get_unique_mapping(mol1, mol2, mcs_mol=mcs.queryMol)
            group_list1 = UniGroupDecider._get_alchem_groups(mol1, mcs.queryMol, group_list1_match)
            group_list2 = UniGroupDecider._get_alchem_groups(mol2, mcs.queryMol, group_list2_match)
            mol1_mcs_map = group_list1_match
            mol2_mcs_map = group_list2_match
        else:
            logger.debug(f"MCS is unique, using MCS to find uni groups.")
            group_list1 = UniGroupDecider._get_alchem_groups(mol1, mcs.queryMol)
            group_list2 = UniGroupDecider._get_alchem_groups(mol2, mcs.queryMol)
            mol1_mcs_map = mol1.GetSubstructMatch(mcs.queryMol)
            mol2_mcs_map = mol2.GetSubstructMatch(mcs.queryMol)

        mol1_mol2_map = {mol1_mcs_map[i]: mol2_mcs_map[i] for i in range(mcs.queryMol.GetNumAtoms())}
        mol2_mol1_map = {mol2_mcs_map[i]: mol1_mcs_map[i] for i in range(mcs.queryMol.GetNumAtoms())}

        group_index_map = {}
        visited_index_list = []
        for i, group1 in enumerate(group_list1):
            for j, group2 in enumerate(group_list2):
                if mol1_mol2_map[group1.edge] == group2.edge and j not in visited_index_list:
                    group_index_map[i] = j
                    visited_index_list.append(j)
                    break
                    # Find correspond uni group.

        no_match_group_indexes1 = set(range(len(group_list1))) - set(group_index_map.keys())
        no_match_group_indexes2 = set(range(len(group_list2))) - set(group_index_map.values())
        return_group_list1 = []
        return_group_list2 = []
        for i, j in group_index_map.items():
            return_group_list1.append(group_list1[i])
            return_group_list2.append(group_list2[j])
        for i in no_match_group_indexes1:
            return_group_list1.append(group_list1[i])
            return_group_list2.append(UniGroupDecider.uni_group(mol1_mol2_map[group_list1[i].edge], set()))
        for i in no_match_group_indexes2:
            return_group_list1.append(UniGroupDecider.uni_group(mol2_mol1_map[group_list2[i].edge], set()))
            return_group_list2.append(group_list2[i])

        logger.debug(f"mol1 uni groups: {[i.edge for i in return_group_list1]}")
        logger.debug(f"mol2 uni groups: {[i.edge for i in return_group_list2]}")

        return (UniGroupDecider.mol_with_uni_group(mol1, return_group_list1),
                UniGroupDecider.mol_with_uni_group(mol2, return_group_list2))

    @staticmethod
    def get_atom_ring_sizes(mol, atom_idx):
        ring_info = mol.GetRingInfo()

        atom_ring_sizes = ring_info.AtomRingSizes(atom_idx)

        return sorted(atom_ring_sizes)

    @staticmethod
    def get_core_atom_list(mol_with_uni_group1: 'UniGroupDecider.mol_with_uni_group',
                           mol_with_uni_group2: 'UniGroupDecider.mol_with_uni_group',
                           enforce_3d_distance: bool = True) -> set:
        """
        returns core atom list of mol1
        :param mol_with_uni_group1:
        :param mol_with_uni_group2:
        :return:
        """
        uni_groups1 = mol_with_uni_group1.uni_groups
        uni_groups2 = mol_with_uni_group2.uni_groups
        mol1 = mol_with_uni_group1.mol
        mol2 = mol_with_uni_group2.mol

        alchem_atoms = set()
        for uni_group in uni_groups1:
            alchem_atoms.update(uni_group.group)

        core_atoms = set()
        for atom in mol1.GetAtoms():
            idx = atom.GetIdx()
            if idx not in alchem_atoms and atom.GetSymbol() != "H":
                core_atoms.add(idx)

        if enforce_3d_distance:
            max_distance = 1.0
        else:
            max_distance = None


        for i, group in enumerate(uni_groups1):
            corr_group = uni_groups2[i]
            if MolCoordModer.atom_compare(mol1=mol1, atom1=group.edge, mol2=mol2, atom2=corr_group.edge,
                                          max_distance=max_distance):
                continue
            core_atoms.discard(group.edge)

        return core_atoms


class _UniqueFragmentTransplant:

    def __init__(self, old_mol: Chem.Mol, bone_mol: Chem.RWMol, ref_mol: Chem.Mol,
                 enforce_3d_distance: bool = True):
        """
        Consider this class as a pirate class.
        Use to transplant groups of atoms to a bone_mol
        :param old_mol: the mol
        :param bone_mol: the bone part
        :param ref_mol: the reference mol
        """

        self.bone_mol = bone_mol
        self.old_mol = old_mol
        self.ref_mol = ref_mol

        self.enforce_3d_distance = enforce_3d_distance

        old_mol_with_uni_group, ref_mol_with_uni_group = UniGroupDecider.find_uni_group(self.old_mol, self.ref_mol,
                                                                                        enforce_3d_distance=self.enforce_3d_distance)

        self.old_unique_group = old_mol_with_uni_group.uni_groups
        self.ref_unique_group = ref_mol_with_uni_group.uni_groups

        logger.debug(f"uni groups: {len(self.old_unique_group)}")
        logger.debug(f"edges: {[i.edge for i in self.old_unique_group]}")

    @staticmethod
    def _set_group_xyz(old_conformer: Chem.Conformer, ref_atom_idx: int, group_idxes: Union[list, set],
                       mod_conformer: Chem.Conformer, idx_map: dict):
        """
        TODO: No longer using this method. Just for docstring.
        This function is used to set the atoms' 3D coordinates properly in mod conformer.
        Principle:
        1. Use old_conformer's conformer to strict the relative position of the atoms in group,
        use the ref_atom_idx (usually the edge atom of the group) to anchor.
        2. Calculate the 3D vectors from the ref_atom to the atoms in group.
        3. Try to solve the 3D coordinates of the atoms in group in mod conformer,
        by using the calculated relative vectors.

        :param old_conformer: the conformer contains the 3D coordinates of the atoms in group and the ref atom.
        Use to get relative position between the atoms in group and the ref atom.
        :param ref_atom_idx: the index of the ref atom in old_conformer
        :param group_idxes: the indexes of the atoms in group in old_conformer
        :param mod_conformer: the conformer to be modified
        :param idx_map: an index map between atoms in the old_conformer and the mod_conformer.
        Type: dict. {old_conformer_idx: mod_conformer_idx}
        :return: None
        """
        raise NotImplementedError('This method has been abandoned.')

    @staticmethod
    def _transform_set_group_xyz(old_mol: Chem.Mol, mod_mol: Union[Chem.Mol, Chem.RWMol],
                                 ref_atom_idx: int, group_idxes: Union[list, set],
                                 not_trusted_atom_idxes: Union[list, set],
                                 idx_map: dict,
                                 ):
        """
        This method using Kabsch algorithm to find a good transformation between the minimized state and the bone
        the transform matrix is calculated by closest 4 atom including edge
        :param old_conformer:
        :param ref_atom_idx:
        :param group_idxes:
        :param mod_conformer:
        :param idx_map:
        :param old_mol:
        :return:
        """
        if len(group_idxes) == 0:
            return False

        old_conformer = old_mol.GetConformer()
        mod_conformer = mod_mol.GetConformer()

        visit_queue = [ref_atom_idx]
        visited = set()
        while len(visited) < 4:
            now_atom_idx = visit_queue.pop(0)
            visited.add(now_atom_idx)
            for neighbor_atom in old_mol.GetAtomWithIdx(now_atom_idx).GetNeighbors():
                idx = neighbor_atom.GetIdx()
                if neighbor_atom.GetSymbol() != 'H':
                    if idx not in visited and idx not in group_idxes and idx not in not_trusted_atom_idxes:
                        visit_queue.append(idx)
                        visited.add(idx)
                        if len(visited) >= 4:
                            break
        visited.remove(ref_atom_idx)
        closest_3_atoms = list(visited)

        logger.debug(f"Using {closest_3_atoms} to align edge atom {ref_atom_idx} for atoms {group_idxes}")

        m_old = Mol3dTransformMethod.get_relative_matrix3(old_conformer, ref_atom_idx, closest_3_atoms)
        m_new = Mol3dTransformMethod.get_relative_matrix3(mod_conformer, idx_map[ref_atom_idx],
                                                           [idx_map[i] for i in closest_3_atoms])
        r_t = Mol3dTransformMethod.Kabsch_rotation_matrix(m_old.T, m_new.T)

        ref_p_old = old_conformer.GetAtomPosition(ref_atom_idx)
        ref_p_new = mod_conformer.GetAtomPosition(idx_map[ref_atom_idx])

        length_factor = 0.9
        for group_idx in group_idxes:
            p = old_conformer.GetAtomPosition(group_idx)
            old_v = Mol3dTransformMethod.get_vector(ref_p_old, p)
            old_v = np.array([[old_v[0], old_v[1], old_v[2]]])
            new_v = np.dot(old_v, r_t)
            new_xyz = (ref_p_new.x + new_v[0][0] * length_factor,
                       ref_p_new.y + new_v[0][1] * length_factor,
                       ref_p_new.z + new_v[0][2] * length_factor)
            mod_conformer.SetAtomPosition(idx_map[group_idx], new_xyz)

        return True

    @staticmethod
    def _set_edge_H_xyz(mod_mol: Union[Chem.Mol, Chem.RWMol], old_mol: Chem.Mol, ref_mol: Chem.Mol,
                        not_trusted_atom_idxes: Union[list, set],
                        old_edge_atom_idx: int, ref_edge_atom_idx: int):
        # Notice that the atom index in old mol is the same as the atom index in mod mol.
        # Note: atom indices in old mol and mod mol are the same
        _mod_conformer = mod_mol.GetConformer()
        _ref_conformer = ref_mol.GetConformer()

        _mod_edge_nighbor_H_idxes = []
        for neighbor_atom in mod_mol.GetAtomWithIdx(old_edge_atom_idx).GetNeighbors():
            if neighbor_atom.GetSymbol() == 'H':
                _mod_edge_nighbor_H_idxes.append(neighbor_atom.GetIdx())

        _ref_edge_nighbor_H_idxes = []
        for neighbor_atom in ref_mol.GetAtomWithIdx(ref_edge_atom_idx).GetNeighbors():
            if neighbor_atom.GetSymbol() == 'H':
                _ref_edge_nighbor_H_idxes.append(neighbor_atom.GetIdx())

        if len(_mod_edge_nighbor_H_idxes) == 0 and len(_ref_edge_nighbor_H_idxes) == 0:
            return False

        _UniqueFragmentTransplant._transform_set_group_xyz(old_mol=old_mol, mod_mol=mod_mol,
                                                           ref_atom_idx=old_edge_atom_idx,
                                                           group_idxes=_mod_edge_nighbor_H_idxes,
                                                           not_trusted_atom_idxes=not_trusted_atom_idxes,
                                                           idx_map={i: i for i in range(mod_mol.GetNumAtoms())},
                                                           )

        ref_binded_H_idxes = []
        mod_binded_H_idxes = []
        H_distance_threshold = 0.7
        is_H_num_equal = len(_mod_edge_nighbor_H_idxes) == len(_ref_edge_nighbor_H_idxes)
        for H_idx in _mod_edge_nighbor_H_idxes:
            mod_H_position = _mod_conformer.GetAtomPosition(H_idx)
            for ref_H_idx in _ref_edge_nighbor_H_idxes: 
                ref_H_position = _ref_conformer.GetAtomPosition(ref_H_idx)
                if mod_H_position.Distance(ref_H_position) < H_distance_threshold and ref_H_idx not in ref_binded_H_idxes:
                    ref_binded_H_idxes.append(ref_H_idx)
                    mod_binded_H_idxes.append(H_idx)
                    _mod_conformer.SetAtomPosition(H_idx, ref_H_position)
                    break

        if is_H_num_equal and (len(ref_binded_H_idxes) != len(_ref_edge_nighbor_H_idxes)):
            # H num equal and not all H atoms are binded.
            for H_idx in _mod_edge_nighbor_H_idxes:
                if H_idx in mod_binded_H_idxes:
                    continue
                for ref_H_idx in _ref_edge_nighbor_H_idxes:
                    if ref_H_idx in ref_binded_H_idxes:
                        continue
                    ref_H_position = _ref_conformer.GetAtomPosition(ref_H_idx)
                    ref_binded_H_idxes.append(ref_H_idx)
                    mod_binded_H_idxes.append(H_idx)
                    _mod_conformer.SetAtomPosition(H_idx, ref_H_position)
                    break

        return True

    @staticmethod
    def _set_alchem_H_xyz(mod_mol: Union[Chem.Mol, Chem.RWMol], old_mol: Chem.Mol, alchem_heavy_atom_idx: int):
        neighbor_H_idxes = set()
        for neighbor_atom in old_mol.GetAtomWithIdx(alchem_heavy_atom_idx).GetNeighbors():
            if neighbor_atom.GetSymbol() == 'H':
                neighbor_H_idxes.add(neighbor_atom.GetIdx())

        _UniqueFragmentTransplant._transform_set_group_xyz(old_mol=old_mol, mod_mol=mod_mol,
                                                           ref_atom_idx=alchem_heavy_atom_idx,
                                                           group_idxes=neighbor_H_idxes,
                                                           not_trusted_atom_idxes=[],
                                                           idx_map={i: i for i in range(mod_mol.GetNumAtoms())},
                                                           )

    @staticmethod
    def hypothetical_addHs(mol: Chem.Mol):
        added_mol = copy.deepcopy(mol)
        added_mol.UpdatePropertyCache()
        Chem.GetSSSR(added_mol)
        added_mol.GetRingInfo()

        b4_neighbors = {}
        for atom in added_mol.GetAtoms():
            if atom.GetSymbol() == "H":
                continue
            b4_neighbors[atom.GetIdx()] = len(atom.GetNeighbors())

        added_mol = Chem.AddHs(added_mol)

        after_neighbors = {}
        for atom in added_mol.GetAtoms():
            if atom.GetSymbol() == "H":
                continue
            after_neighbors[atom.GetIdx()] = len(atom.GetNeighbors())

        missingH_atom_idxes = set()
        for idx, neighbors in b4_neighbors.items():
            if neighbors != after_neighbors[idx]:
                missingH_atom_idxes.add(idx)

        logger.debug(f"Add Hs to atoms: {missingH_atom_idxes}")

        return Chem.AddHs(mol, addCoords=True, onlyOnAtoms=missingH_atom_idxes), missingH_atom_idxes

    def transplant(self) -> Chem.Mol:
        mod_mol = copy.deepcopy(self.old_mol)
        _mcs = MolCoordModer.customized_FindMCS([self.old_mol, self.ref_mol], enforce_3d_distance=True)
        old_mol_match, ref_mol_match = UniGroupDecider.get_unique_mapping(self.old_mol, self.ref_mol, _mcs.queryMol)
        old2ref_map = {old_mol_match[i]: ref_mol_match[i] for i in range(len(old_mol_match))}
        old2mod_map = {i: i for i in range(mod_mol.GetNumAtoms())}

        alchem_heavy_atom_idxes = set()
        allow_optimize_H_edge_idxes = set()
        for group in self.old_unique_group:
            alchem_heavy_atom_idxes.update(group.group)
        not_trusted_atom_idxes = copy.deepcopy(alchem_heavy_atom_idxes)

        logger.debug(f"old2ref_map: {old2ref_map}")

        ref_conformer = self.ref_mol.GetConformer()
        mod_conformer = mod_mol.GetConformer()

        for old_idx, ref_idx in old2ref_map.items():
            if self.old_mol.GetAtomWithIdx(old_idx).GetSymbol() == 'H' or self.ref_mol.GetAtomWithIdx(ref_idx).GetSymbol() == 'H':
                # Due to distance constraints in MCS, do not trust H atoms found by MCS search
                continue
            mod_conformer.SetAtomPosition(old_idx, ref_conformer.GetAtomPosition(ref_idx))
            neighbor_H_idxes = []
            for neighbor_atom in mod_mol.GetAtomWithIdx(old_idx).GetNeighbors():
                if neighbor_atom.GetSymbol() == 'H':
                    neighbor_H_idxes.append(neighbor_atom.GetIdx())
            if len(neighbor_H_idxes) == 0:
                # For scaffold atoms with no H atoms, H atom mapping is not needed
                continue

            corr_H_idxes = []
            for neighbor_atom in self.ref_mol.GetAtomWithIdx(ref_idx).GetNeighbors():
                if neighbor_atom.GetSymbol() == 'H':
                    corr_H_idxes.append(neighbor_atom.GetIdx())
            if len(neighbor_H_idxes) == len(corr_H_idxes):
                # If H atom counts are the same, it's not a boundary, H atom mapping in function already handles equal H count cases
                self._set_edge_H_xyz(mod_mol=mod_mol, old_mol=self.old_mol, ref_mol=self.ref_mol,
                                     old_edge_atom_idx=old_idx,
                                     not_trusted_atom_idxes=not_trusted_atom_idxes,
                                     ref_edge_atom_idx=ref_idx)
                # Different cases are handled at edge

        for i, uni_group in enumerate(self.old_unique_group):
            # for each group, find ref atom and optimize atoms' coordinate in group separately.
            logger.debug(f'edge: {uni_group.edge} - {uni_group.group}')
            self._set_edge_H_xyz(mod_mol=mod_mol, old_mol=self.old_mol, ref_mol=self.ref_mol,
                                 old_edge_atom_idx=uni_group.edge,
                                 not_trusted_atom_idxes=not_trusted_atom_idxes,
                                 ref_edge_atom_idx=self.ref_unique_group[i].edge)
            
            edge_H_nums = len([i for i in mod_mol.GetAtomWithIdx(uni_group.edge).GetNeighbors() if i.GetSymbol() == 'H'])
            ref_edge_H_nums = len([i for i in self.ref_mol.GetAtomWithIdx(self.ref_unique_group[i].edge).GetNeighbors() if i.GetSymbol() == 'H'])
            if edge_H_nums != ref_edge_H_nums:
                allow_optimize_H_edge_idxes.add(uni_group.edge)

            if len(uni_group.group) == 0:
                # A group in other mol transform to H in this mol.
                continue

            self._transform_set_group_xyz(old_mol=self.old_mol, mod_mol=mod_mol, 
                                          not_trusted_atom_idxes=not_trusted_atom_idxes,
                                          ref_atom_idx=uni_group.edge,
                                          group_idxes=uni_group.group,
                                          idx_map=old2mod_map)
            log_info = 'Kabsch method'
            not_trusted_atom_idxes.difference_update(uni_group.group)

            logger.debug(f'group{i} edge: {uni_group.edge} - {log_info}')

        for alchem_heavy_atom_idx in alchem_heavy_atom_idxes:
            self._set_alchem_H_xyz(mod_mol=mod_mol, old_mol=self.old_mol, alchem_heavy_atom_idx=alchem_heavy_atom_idx)

        b4_opt_conf = copy.deepcopy(mod_mol).GetConformer()
        fixed_atom_idxes = []
        mmff_p = rdForceFieldHelpers.MMFFGetMoleculeProperties(mod_mol, mmffVariant='MMFF94s')
        mmff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mod_mol, mmff_p)

        for atom_id in range(mod_mol.GetNumAtoms()):
            if atom_id not in alchem_heavy_atom_idxes:
                # Set fixed for bone atoms.
                _dealing_atom = mod_mol.GetAtomWithIdx(atom_id)
                if _dealing_atom.GetSymbol() == 'H':
                    if (_dealing_atom.GetNeighbors()[0].GetIdx() in alchem_heavy_atom_idxes or
                        _dealing_atom.GetNeighbors()[0].GetIdx() in allow_optimize_H_edge_idxes):
                        # Allow the H atom connecting edge to be minimized.
                        # logger.debug(f'Allow H atom to be minimized: {atom_id}')
                        continue
                mmff.AddFixedPoint(atom_id)
                fixed_atom_idxes.append(atom_id)

        allow_optimize_atom_idxes = set(range(mod_mol.GetNumAtoms())) - set(fixed_atom_idxes)
        logger.debug(f"allow_optimize_atom_idxes: {allow_optimize_atom_idxes}")

        if logger.level == logging.DEBUG:
            Chem.MolToMolFile(mod_mol, 'mol_b4_opt.mol')

        if len(allow_optimize_atom_idxes) == 0:
            warnings.warn("No atoms are allowed to be optimized."
                          "This usually occurs when a molecule is a subset of the other molecule.")
            return mod_mol

        try:
            rdForceFieldHelpers.OptimizeMoleculeConfs(mod_mol, ff=mmff)
            after_opt_conf = mod_mol.GetConformer()
            rmsd = 0
            for i in range(mod_mol.GetNumAtoms()):
                if i not in allow_optimize_atom_idxes:
                    continue
                rmsd += (b4_opt_conf.GetAtomPosition(i).Distance(after_opt_conf.GetAtomPosition(i))) ** 2
            rmsd = np.sqrt(rmsd / len(allow_optimize_atom_idxes))
            logger.debug(f"Alchem group RMSD after opt: {rmsd}")
        except RuntimeError:
            # All atoms are fixed, so optimization is failed.
            warnings.warn("Couldn't optimize conf.")

        return mod_mol


class WatVina_api:
    def __init__(self, receptor='rec.pdbqt', ligand='lig.pdbqt', _core_sdf='core.sdf', _lig_sdf='lig.sdf',
                 _protein_pdb='protein.pdb', watvina_path=WATVINA_PATH):
        """

        :param receptor: the name of OUTPUT receptor file name
        :param ligand: the name of OUTPUT ligand file name
        :param _core_sdf:
        :param _lig_sdf:
        :param _protein_pdb:
        :param watvina_path:
        """

        self.receptor = receptor
        self.ligand = ligand
        self.center_x = 0
        self.center_y = 0
        self.center_z = 0
        self.size_x = 0
        self.size_y = 0
        self.size_z = 0
        self._core_sdf = os.path.abspath(_core_sdf)
        self._lig_sdf = os.path.abspath(_lig_sdf)
        self._protein_pdb = os.path.abspath(_protein_pdb)
        self._watvina_excute_path = watvina_path

    def write_conf(self, filename='conf.txt'):
        with open(filename, 'w') as f:
            for attr, value in self.__dict__.items():
                if not attr.startswith('_'):
                    f.write(f'{attr} = {value}\n')

    def write_protein_pdbqt(self, ):
        receptor_mol = protein_pdb_reader(self._protein_pdb)
        receptor_lines = MolToPDBQTBlock(receptor_mol, False, False, True)
        with open(self.receptor, 'w') as f:
            f.writelines(receptor_lines)

    def write_ligand_pdbqt(self, ):
        """
        No longer using this method.
        """
        raise NotImplementedError("This method is not implemented.")

    def find_box(self, extend_length=5.0):
        """
        Find the box information according to a bound ligand.
        :param extend_length: the extend length of the box size, default is 5.0 Angstrom.
        """
        mol_obj = Chem.MolFromMolFile(self._lig_sdf, removeHs=False)
        # Ensure the molecule has 3D coordinates
        if not mol_obj.GetNumConformers():
            raise ValueError("The molecule has no 3D coordinates.")

        conformer = mol_obj.GetConformer()
        coords = np.array([conformer.GetAtomPosition(i) for i in range(mol_obj.GetNumAtoms())])

        # Calculate the centroid
        centroid = coords.mean(axis=0)

        # Calculate size in each dimension
        min_coords = coords.min(axis=0)
        max_coords = coords.max(axis=0)
        size = max_coords - min_coords + extend_length

        self.center_x = centroid[0]
        self.center_y = centroid[1]
        self.center_z = centroid[2]
        self.size_x = size[0]
        self.size_y = size[1]
        self.size_z = size[2]

    @staticmethod
    def perform_docking(watvina_path, conf='conf.txt', ifconstrain_dock=True, max_retries=3):
        if not Executable_checker.is_executable_available(watvina_path):
            raise FileNotFoundError(f"{watvina_path} not found."
                                    f"Please check if the program is installed and available in the system PATH.")
        else:
            logger.debug(f"{watvina_path} is available in the system PATH and executable.")
        if ifconstrain_dock:
            command = [f"{watvina_path}", "--config", conf, "--tramplitude", '0']
        else:
            command = [f"{watvina_path}", "--config", conf]

        retry_count = 0
        while retry_count < max_retries:
            try:
                with open("dock.log", "w") as dock_log:
                    process = subprocess.Popen(" ".join(command), shell=True, 
                                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                              universal_newlines=True, bufsize=1)

                    for line in process.stdout:
                        dock_log.write(line)
                        dock_log.flush()

                    stderr = process.stderr.read()
                    if stderr:
                        dock_log.write("-"*21+"ERROR"+"-"*21+"\n")
                        dock_log.write(stderr)
                        dock_log.flush()

                    return_code = process.wait()
                    if return_code != 0:
                        raise subprocess.CalledProcessError(return_code, " ".join(command))

                    break

            except subprocess.CalledProcessError as e:
                retry_count += 1
                logger.warning(f"Docking attempt {retry_count} failed: {e}")
                # if retry_count >= max_retries:
                #     logger.warning(f"Failed to run watvina after {max_retries} attempts: {e}")

                time.sleep(1)

    @staticmethod
    def convert_pdbqt2pdb(pdbqt_file, pdb_file, res_name):
        pdb_lines = []
        with open(pdbqt_file, 'r') as pdbqt_f:
            pdbqt_lines = pdbqt_f.readlines()
        for line in pdbqt_lines:
            if line.startswith('ENDMDL'):
                pdb_lines.append(line)
                break
            if line.startswith("ATOM") or line.startswith("HETATM"):
                line.replace('UNL', res_name)
                pdb_lines.append(line[:68]+'\n')
        with open(pdb_file, 'w') as pdb_f:
            pdb_f.writelines(pdb_lines)

    @staticmethod
    def single_dock(core_sdf, lig_sdf, protein_pdb, output_pdb, _aligned_mol=None, spec_lig_pdbqt=None) -> Union[Chem.Mol, bool]:
        logger.debug("Docking...")
        output_pdb_filename = os.path.basename(output_pdb)
        watvina_dock = WatVina_api(_core_sdf=core_sdf, _lig_sdf=lig_sdf, _protein_pdb=protein_pdb)
        watvina_dock.docking_stream(working_path='dock', ligand_pdb_name=output_pdb_filename,
                                    specified_lig_pdbqt=spec_lig_pdbqt)

        shutil.copy2(os.path.join('./dock', output_pdb_filename), output_pdb)
        if _aligned_mol is not None:
            docked_mol = MolCoordModer.update_mol_coordinates(_aligned_mol, output_pdb)
            return docked_mol

    def docking_stream(self, working_path, extend_length=5.0, ligand_pdb_name='lig.pdb', ligand_res_name='LIG',
                       specified_lig_pdbqt=None) -> bool:
        with working_dir(working_path):
            self.find_box(extend_length)
            self.write_conf('conf.txt')
            self.write_protein_pdbqt()
            if specified_lig_pdbqt is not None:
                shutil.copy2(specified_lig_pdbqt, self.ligand)
            else:
                self.write_ligand_pdbqt()
            WatVina_api.perform_docking(self._watvina_excute_path, 'conf.txt', True)
            result_pdbqt = self.ligand.split('.')[0]+'_out.pdbqt'
            if not os.path.exists(result_pdbqt):
                retry_count = 0
                while retry_count < 3:
                    WatVina_api.perform_docking(self._watvina_excute_path, 'conf.txt', True)
                    if os.path.exists(result_pdbqt):
                        break
                    retry_count += 1
                if retry_count >= 3:
                    logger.warning("Failed to dock after 3*3 attempts. Please check this pair manually.")
                    shutil.copy2(self.ligand, result_pdbqt)
                    WatVina_api.convert_pdbqt2pdb(result_pdbqt, ligand_pdb_name, ligand_res_name)
                    return False
            WatVina_api.convert_pdbqt2pdb(result_pdbqt, ligand_pdb_name, ligand_res_name)
            return True


class MolCoordModer:
    def __init__(self, bone_mol: Union[Chem.Mol, Chem.RWMol]):
        """
        Using to align molecule MCS's coordinate to a bone molecule.
        :param bone_mol: MUST embed the 3D coordinates of both explicit and implicit H atoms !!!
        """
        if isinstance(bone_mol, Chem.Mol):
            bone_mol = Chem.RWMol(bone_mol)
        self.bone_mol = bone_mol

    @staticmethod
    def has_3d_conformer(mol: Chem.Mol):
        for conformer in mol.GetConformers():
            if conformer.Is3D():
                return True
        return False
    
    @staticmethod
    def align(other_mol: Chem.Mol, ref_mol: Chem.Mol, bone_mol: Union[Chem.Mol, Chem.RWMol],
              enforce_3d_distance: bool = True) -> Tuple[_UniqueFragmentTransplant, Chem.Mol]:
        """
        Returns a mol that has the same coordinates as the atoms in the same structure between itself and bone_mol.
        :param other_mol:
        :param ref_mol:
        :return: 
        """
        _other_mol = copy.deepcopy(other_mol)

        if _other_mol.GetNumAtoms() != Chem.AddHs(_other_mol).GetNumAtoms():
            _other_mol = Chem.RemoveAllHs(_other_mol)
            _other_mol = Chem.AddHs(_other_mol)
            # In case the target mol doesn't have implicit H.

        if not MolCoordModer.has_3d_conformer(other_mol):
            enforce_3d_distance = False
            rdDistGeom.EmbedMolecule(_other_mol)
            rdForceFieldHelpers.UFFOptimizeMolecule(_other_mol)
            if rdForceFieldHelpers.MMFFHasAllMoleculeParams(_other_mol):
                rdForceFieldHelpers.MMFFOptimizeMolecule(_other_mol)
        # Add 3D coordinates in case the target mol doesn't have it.

        _transplant_tool = _UniqueFragmentTransplant(old_mol=_other_mol,
                                                     bone_mol=bone_mol,
                                                     ref_mol=ref_mol,
                                                     enforce_3d_distance=enforce_3d_distance)
        return _transplant_tool, _transplant_tool.transplant()

    @staticmethod
    def align_with_MCS(ref_mol: Chem.Mol, other_mol: Chem.Mol, enforce_3d_distance: bool = True) -> Tuple[Chem.Mol, Chem.Mol, _UniqueFragmentTransplant]:
        """
        Align target molecule's MCS part's 3D coordinates with the MCS between the reference and target molecule
        :param other_mol: the mol to be aligned.
        :param ref_mol: the mol provide 3d coordinates. Must embed 3d coordinates before aligning!!!
        :param MCS_parameters: the MCS parameters for rdFMCS.FindMCS.
                            If not provided, defaults to a little changes from default MCS_parameters
        :return: Modified copy of target molecule and MCS mol
        """

        bone_mol = MolCoordModer.extract_bone_mol(ref_mol=ref_mol,
                                                  other_mol=other_mol,
                                                  enforce_3d_distance=enforce_3d_distance)
        # Extract bone mol from ref mol.

        return_tool, return_mol = MolCoordModer.align(other_mol=other_mol,
                                                      ref_mol=ref_mol,
                                                      bone_mol=bone_mol,
                                                      enforce_3d_distance=enforce_3d_distance)
        return_mol.SetProp('_Name', other_mol.GetProp('_Name'))

        return return_mol, bone_mol, return_tool

    @staticmethod
    def deal_amibiguous_aromatic_H(bone_mol: Chem.RWMol) -> Chem.RWMol:
        need_fix_H = False
        edit_mol = copy.deepcopy(bone_mol)

        for atom in edit_mol.GetAtoms():
            if atom.GetSymbol() == 'H':
                continue
            actual_bonds = sorted([i.GetBondTypeAsDouble() for i in atom.GetBonds()])
            expect_bonds = sorted([float(i) for i in edit_mol.GetAtomWithIdx(atom.GetIdx()).GetProp('expect_bonds').split(',')])
            if expect_bonds != actual_bonds:
                edit_mol.GetAtomWithIdx(atom.GetIdx()).SetNumExplicitHs(int(sum(expect_bonds)-sum(actual_bonds)))
                need_fix_H = True

        if need_fix_H:
            rwmol = Chem.RWMol(edit_mol)
            for atom in rwmol.GetAtoms():
                atom.SetIsAromatic(False)
            for bond in rwmol.GetBonds():
                bond.SetIsAromatic(False)
            Chem.SanitizeMol(rwmol, sanitizeOps=Chem.SANITIZE_ALL^Chem.SANITIZE_SETAROMATICITY)
            Chem.Kekulize(rwmol, clearAromaticFlags=True)
            Chem.SetAromaticity(rwmol)
            rwmol, _ = _UniqueFragmentTransplant.hypothetical_addHs(rwmol)
            edit_mol = Chem.RWMol(rwmol)
            # Add explicit H atoms to some ambiguous aromatic situation.

        return edit_mol

    @staticmethod
    def extract_bone_mol(ref_mol: Chem.Mol, other_mol: Chem.Mol, enforce_3d_distance: bool = True) -> Chem.RWMol:
        """
        Extract bone molecule from reference molecule.
        :param ref_mol: the reference molecule. Should embed 3D Coordinates.
        :param mcs_query_mol: the bone molecule that generate by rdFMCS.FindMCS
        :return: a bone mol with 3D Coordinates.
        """
        bone_mol = Chem.RWMol(ref_mol)

        for atom in ref_mol.GetAtoms():
            if atom.GetSymbol() == 'H':
                continue
            neighbors = [i.GetBondTypeAsDouble() for i in atom.GetBonds()]
            expect_bonds = ','.join(sorted([str(i) for i in neighbors]))
            bone_mol.GetAtomWithIdx(atom.GetIdx()).SetProp('expect_bonds', expect_bonds)

        ref_info, other_info = UniGroupDecider.find_uni_group(ref_mol, other_mol, enforce_3d_distance=enforce_3d_distance)
        core_atom_idx_set = UniGroupDecider.get_core_atom_list(ref_info, other_info, enforce_3d_distance=enforce_3d_distance)
        tmp_remove_atom_list = []

        for atom_id in range(bone_mol.GetNumAtoms()):
            if atom_id not in core_atom_idx_set:
                if bone_mol.GetAtomWithIdx(atom_id).GetSymbol() == 'H':
                    if bone_mol.GetAtomWithIdx(atom_id).GetNeighbors()[0].GetIdx() in core_atom_idx_set:
                        core_atom_idx_set.add(atom_id)
                        continue
                tmp_remove_atom_list.append(atom_id)

        # Some symmetic situation here.
        tmp_remain_atom_idx_set = set(range(bone_mol.GetNumAtoms())) - set(tmp_remove_atom_list)
        # tmp_ringinfo = bone_mol.GetRingInfo()
        h_to_remove = set()
        for atom_id in tmp_remain_atom_idx_set:
            atom = bone_mol.GetAtomWithIdx(atom_id)
            if atom.GetSymbol() == 'H':
                h_to_remove.add(atom_id)
        tmp_remove_atom_list.extend(list(h_to_remove))

        edge_atom_idx_set = set([grp.edge for grp in ref_info.uni_groups])
        for edge_atom_idx in edge_atom_idx_set:
            for neighbor_atom in ref_mol.GetAtomWithIdx(edge_atom_idx).GetNeighbors():
                if neighbor_atom.GetSymbol() == 'H' and neighbor_atom.GetIdx() not in tmp_remove_atom_list:
                    logger.debug(f"Clean H {neighbor_atom.GetIdx()} on {edge_atom_idx}")
                    tmp_remove_atom_list.append(neighbor_atom.GetIdx())
        # Here we clean H atoms on edge atom to avoid colliding.

        tmp_remove_atom_list.sort(reverse=True)
        for remove_id in tmp_remove_atom_list:
            bone_mol.RemoveAtom(remove_id)
        # Remove not included atoms.

        for atom in bone_mol.GetAtoms():
            atom.UpdatePropertyCache()
            if not atom.IsInRing():
                atom.SetIsAromatic(False)
        # Remove Aromatic mark.

        return bone_mol

    @staticmethod
    def atom_compare(mol1, atom1, mol2, atom2, p=None, self=None, max_distance: Union[float, None] = None):
        a1 = mol1.GetAtomWithIdx(atom1)
        a2 = mol2.GetAtomWithIdx(atom2)
        if a1.GetAtomicNum() != a2.GetAtomicNum():
            return False
        if p is not None and self is not None:
            if p.MatchChiralTag and not self.CheckAtomChirality(p, mol1, atom1, mol2, atom2):
                return False
        if max_distance is not None:
            c1 = mol1.GetConformer().GetAtomPosition(atom1)
            c2 = mol2.GetConformer().GetAtomPosition(atom2)
            if c1.Distance(c2) > max_distance:
                return False

        if a1.GetHybridization() == a2.GetHybridization():
            if a1.IsInRing() and a2.IsInRing():
                a1_ring_sizes = UniGroupDecider.get_atom_ring_sizes(mol1, atom1)
                a2_ring_sizes = UniGroupDecider.get_atom_ring_sizes(mol2, atom2)
                if a1_ring_sizes == a2_ring_sizes:
                    return True
            elif not a1.IsInRing() and not a2.IsInRing():
                return True
        return False

    class CustomCompareElements(rdFMCS.MCSAtomCompare):
        def __init__(self, max_distance: Union[float, None] = None):
            super().__init__()
            self.max_distance = max_distance

        def __call__(self, p, mol1, atom1, mol2, atom2):
            return MolCoordModer.atom_compare(mol1, atom1, mol2, atom2, p=p, self=self,
                                              max_distance=self.max_distance)

    @staticmethod
    def customized_FindMCS(mol_list: list, MCS_parameters: rdFMCS.MCSParameters = None, enforce_3d_distance: bool = False) -> rdFMCS.MCSResult:
        if MCS_parameters is None:
            params = rdFMCS.MCSParameters()
            if enforce_3d_distance:
                cce = MolCoordModer.CustomCompareElements(max_distance=99999.0)
                # Change to 5.0 means consider 3d distance no mare
            else:
                cce = MolCoordModer.CustomCompareElements()
            params.AtomTyper = cce
            params.BondTyper = rdFMCS.BondCompare.CompareOrderExact
            params.BondCompareParameters.RingMatchesRingOnly = True
            params.BondCompareParameters.CompleteRingsOnly = False
            params.AtomCompareParameters.MatchChiralTag = False
            params.ringMatchesRingOnly = True
            params.Timeout = 300
        else:
            params = MCS_parameters

        mcs = rdFMCS.FindMCS(mols=mol_list, parameters=params)
        return mcs

    @staticmethod
    def write_sdf_file_from_mol(mol: Chem.Mol, file_path: str):
        """
        Write a sdf file from a Chem.Mol.
        :param mol: the Chem.Mol object to be write as sdf file.
        :param file_path: the str for the path sdf file.
        """
        with Chem.SDWriter(file_path) as sdf_file:
            sdf_file.write(mol)
    
    @staticmethod
    def _get_atom_idx_dict(pdb_file):
        return_dict = {}
        counter = 0
        with open(pdb_file, 'r') as f:
            for line in f.readlines():
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # idx = int(line[7:11].strip()) - 1
                    atom_name = line[12:16].strip()
                    return_dict[atom_name] = counter
                    counter += 1
        return return_dict

    @staticmethod
    def update_mol_coordinates(mol: Chem.Mol, pdb_file: str):
        """
        Update the 3D coordinates of the mol atoms with the dock result.
        :param mol: the Chem.Mol object to be updated 3D coordinate
        :param pdb_file: the docked pdb file, contain dock result.
        :return:
        """

        dock_xyz_result = {}
        docked_mol = Chem.MolFromPDBFile(pdb_file, sanitize=False, removeHs=False)
        docked_conformer = docked_mol.GetConformer()
        for atom_idx in range(docked_conformer.GetNumAtoms()):
            dock_xyz_result[atom_idx] = docked_conformer.GetAtomPosition(atom_idx)

        if docked_mol.GetNumAtoms() != mol.GetNumAtoms():
            raise ValueError("mol and pdb file must have the same number of atoms.")

        with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=True) as temp_pdb:
            Chem.MolToPDBFile(mol, temp_pdb.name)
            mol_atom_dict = MolCoordModer._get_atom_idx_dict(temp_pdb.name)
            
        docked_atom_dict = MolCoordModer._get_atom_idx_dict(pdb_file)
        idx_map = {}
        for atom_name, atom_idx in mol_atom_dict.items():
            idx_map[atom_idx] = docked_atom_dict[atom_name]

        mod_mol_conformer = mol.GetConformer()
        with open("corr.txt", 'w') as f:
            for a, b in idx_map.items():
                f.write(f"{a} : {b}\n")
            f.write("\n")
            for a, b in mol_atom_dict.items():
                f.write(f"{a} : {b}\n")
            f.write("\n")
            for a, b in docked_atom_dict.items():
                f.write(f"{a} : {b}\n")

        for atom_idx, docked_atom_idx in idx_map.items():
            mod_mol_conformer.SetAtomPosition(atom_idx, dock_xyz_result[docked_atom_idx])

        return_mol = copy.deepcopy(mol)
        return_mol.RemoveAllConformers()
        return_mol.AddConformer(mod_mol_conformer)

        return return_mol

    @staticmethod
    def remove_mol_info(mol):
        for atom in mol.GetAtoms():
            atom.UpdatePropertyCache()
            atom.SetIsAromatic(False)


class PairProcess:  # Temp Name
    def __init__(self, ref_mol: Chem.Mol, other_mol: Chem.Mol, protein_pdb: str,
                 trust_input_conformation: bool = True):
        self.ref_mol = ref_mol
        self.other_mol = other_mol
        self.trust_input_conformation = trust_input_conformation

        if ref_mol.HasProp('_Name') and ref_mol.GetProp('_Name') != 'other_mol':
            self.ref_mol_name = self.ref_mol.GetProp("_Name")
        else:
            self.ref_mol_name = 'ref_mol'
            self.ref_mol.SetProp('_Name', self.ref_mol_name)

        if other_mol.HasProp('_Name'):
            self.other_mol_name = self.other_mol.GetProp("_Name")
        else:
            self.other_mol_name = 'other_mol'
            self.other_mol.SetProp('_Name', self.other_mol_name)

        self.protein_pdb = protein_pdb

    def process_stream(self, _if_dock: bool = True, _working_dir: Union[str, None] = 'preprocess'):
        """
        Process the molecule
        :param _if_dock:
        :param _working_dir:
        :return: path to generated other_mol.pdb and other_mol.mol
        """
        ori_cwd = os.getcwd()

        err_happen = False

        if _working_dir is not None:
            os.makedirs(_working_dir, exist_ok=True)
            os.chdir(_working_dir)

        os.makedirs('./gen_files', exist_ok=True)
        gen_files_folder = os.path.abspath('./gen_files')

        collect_aligned_mol_pdb = f'{self.other_mol_name}.pdb'
        collect_aligned_mol_mol = f'{self.other_mol_name}.mol'
        collect_ref_mol_pdb = f'{self.ref_mol_name}.pdb'

        _ref_mol_pdb = os.path.join(gen_files_folder, f'{self.ref_mol_name}.pdb')
        _aligned_mol_pdb = os.path.join(gen_files_folder, f'{self.other_mol_name}.pdb')
        _aligned_mol_sdf = os.path.join(gen_files_folder, f'{self.other_mol_name}.sdf')
        _bone_sdf = os.path.join(gen_files_folder, 'bone.sdf')

        mol_2d = [Chem.MolFromSmiles(Chem.MolToSmiles(i, allHsExplicit=True, allBondsExplicit=True))
                  for i in [self.ref_mol, self.other_mol]]

        old_mol_with_group_info, ref_mol_with_group_info = UniGroupDecider.find_uni_group(mol_2d[1], mol_2d[0],
                                                                                          enforce_3d_distance=False)
        core_atms_2d = UniGroupDecider.get_core_atom_list(old_mol_with_group_info,
                                                          ref_mol_with_group_info,
                                                          enforce_3d_distance=False)
        highlight_ref_atoms = set()
        for i in ref_mol_with_group_info.uni_groups:
            highlight_ref_atoms.add(i.edge)

        mol_2d_mcs = MolCoordModer.customized_FindMCS(mol_2d).queryMol
        rdDepictor.Compute2DCoords(mol_2d_mcs)
        for _mol in mol_2d:
            rdDepictor.GenerateDepictionMatching2DStructure(_mol, mol_2d_mcs)

        pic = MolsToGridImage(mol_2d, subImgSize=(500, 500),
                              legends=['ref', 'other', 'mcs'],
                              highlightAtomLists=[highlight_ref_atoms, core_atms_2d, []],
                              addAtomIndices=True,
                              )
        pic.save("mark.png")

        Chem.MolToPDBFile(self.ref_mol, _ref_mol_pdb)
        _aligned_mol, _bone_mol, _align_tool = MolCoordModer.align_with_MCS(ref_mol=self.ref_mol,
                                                                            other_mol=self.other_mol,
                                                                            enforce_3d_distance=self.trust_input_conformation)

        MolCoordModer.remove_mol_info(_bone_mol)

        Chem.MolToPDBFile(_aligned_mol, _aligned_mol_pdb)
        if logger.level == logging.DEBUG:
            try:
                # The generated bone sdf file is not longer needed. Only used for debug.
                MolCoordModer.write_sdf_file_from_mol(_bone_mol, _bone_sdf)
            except Chem.rdchem.KekulizeException as e:
                _bone_mol = MolCoordModer.deal_amibiguous_aromatic_H(_bone_mol)
                MolCoordModer.write_sdf_file_from_mol(_bone_mol, _bone_sdf)
        MolCoordModer.write_sdf_file_from_mol(_aligned_mol, _aligned_mol_sdf)

        if _if_dock:
            output_pdb = os.path.join(os.path.abspath(os.getcwd()), collect_aligned_mol_pdb)

            fixed_atom_list = set([i.GetIdx() for i in _aligned_mol.GetAtoms() if i.GetSymbol() != 'H'])
            for group in _align_tool.old_unique_group:
                if not UniGroupDecider.is_edge_split_bone(_aligned_mol, group.edge, fixed_atom_list):
                    logger.debug(f"edge {group.edge} would split bone")
                    fixed_atom_list.discard(group.edge)
                fixed_atom_list = fixed_atom_list.symmetric_difference([i for i in group.group])
            
            try:
                mol_3d = [self.ref_mol, _aligned_mol, self.other_mol , _bone_mol, _aligned_mol]
                pic = MolsToGridImage(mol_3d, subImgSize=(500, 500),
                                    legends=['ref', 'aligned','old', 'bone', 'fixed'],
                                    highlightAtomLists=[[], [], [], [], fixed_atom_list],
                                    addAtomIndices=True,
                                    )
                pic.save("Summary.png")
            except Exception as e:
                logger.warning(f"Failed to generate summary image: {e}")
            logger.debug(f"fixed_atom_list: {fixed_atom_list}")

            pdbqtlines = MolCoreToPDBQTBlock(_aligned_mol, fixed_atom_list,
                                             True, False, True)
            dock_lig_pdbqt = os.path.join(gen_files_folder, f'dock_lig.pdbqt')
            with open(dock_lig_pdbqt, 'w') as f:
                f.write(pdbqtlines)

            docked_mol = WatVina_api.single_dock(core_sdf=_bone_sdf, lig_sdf=_aligned_mol_sdf,
                                                 protein_pdb=self.protein_pdb, output_pdb=output_pdb,
                                                 _aligned_mol=_aligned_mol, spec_lig_pdbqt=dock_lig_pdbqt)
            Chem.MolToMolFile(docked_mol, collect_aligned_mol_mol)
        else:
            shutil.copy2(_aligned_mol_pdb, collect_aligned_mol_pdb)
            Chem.MolToMolFile(_aligned_mol, collect_aligned_mol_mol)

        shutil.copy2(_ref_mol_pdb, collect_ref_mol_pdb)

        pdb_info_clean(collect_aligned_mol_pdb, res_name='LGO')
        pdb_info_clean(collect_ref_mol_pdb, res_name='LGR')

        pdb_formatter(collect_aligned_mol_pdb, 'title')
        pdb_formatter(collect_ref_mol_pdb, 'title')

        os.chdir(ori_cwd)

        return (os.path.abspath(collect_ref_mol_pdb),
                os.path.abspath(collect_aligned_mol_pdb),
                os.path.abspath(collect_aligned_mol_mol))


def read_mol_from_file(filename) -> Chem.Mol:
    if os.path.exists(filename) and os.path.isfile(filename):
        # Assume that input is a path to file.
        if filename.endswith('.mol'):
            return Chem.MolFromMolFile(filename, removeHs=False)
        elif filename.endswith('.pdb'):
            return Chem.MolFromPDBFile(filename, removeHs=False)
        elif filename.endswith('.sdf'):
            sd_reader = Chem.SDMolSupplier(filename, removeHs=False)
            for mol in sd_reader:
                return mol
        elif filename.endswith('.mol2'):
            return Chem.MolFromMol2File(filename, removeHs=False)
        else:
            # Assume that input is a text file containing smiles expression
            with open(filename, 'r') as f:
                smi = f.readline().strip()
                try:
                    return Chem.MolFromSmiles(smi)
                except Exception:
                    raise ValueError(f'Invalid SMILES in file {filename}')
    else:
        # Directly read smiles from CLI input
        smi = filename
        try:
            return Chem.MolFromSmiles(smi)
        except Exception:
            raise ValueError(f'Invalid SMILES expression {smi}')


def GeneratePairPDB(ref_mol: Union[str, Chem.Mol], other_mol: Union[str, Chem.Mol],
                    protein_pdb: str, _output_dir: str, _if_dock=True,
                    trust_input_conformation: bool = True) -> Tuple[str, str, str]:

    _output_dir = os.path.abspath(_output_dir)
    os.makedirs(_output_dir, exist_ok=True)

    if isinstance(ref_mol, str):
        ref_Mol = read_mol_from_file(ref_mol)
    elif isinstance(ref_mol, Chem.Mol):
        ref_Mol = ref_mol
    else:
        raise ValueError(f'ref_mol must be str or Chem.Mol')

    if isinstance(other_mol, str):
        other_Mol = read_mol_from_file(other_mol)
    elif isinstance(other_mol, Chem.Mol):
        other_Mol = other_mol
    else:
        raise ValueError(f'other_mol must be str or Chem.Mol')

    ref_Mol_atom_nums = ref_Mol.GetNumAtoms()
    if Chem.AddHs(ref_Mol).GetNumAtoms() != ref_Mol_atom_nums:
        raise ValueError(f'The given reference molecule does not have all the explicit H atoms.'
                         f'{Chem.AddHs(ref_Mol).GetNumAtoms()} > {ref_Mol_atom_nums}')
        # Molecule lost some explicit Hs.
    other_Mol = Chem.AddHs(other_Mol)

    origin_cwd = os.getcwd()

    gen_file_dir = rf'./gen_file/{ref_Mol.GetProp("_Name")}-{other_Mol.GetProp("_Name")}'
    if os.path.exists(gen_file_dir):
        shutil.rmtree(gen_file_dir)
    os.makedirs(gen_file_dir, exist_ok=True)

    os.chdir(gen_file_dir)

    pair_preprocessor = PairProcess(ref_mol=ref_Mol, other_mol=other_Mol, protein_pdb=protein_pdb,
                                    trust_input_conformation=trust_input_conformation)
    ref_mol_pdb, other_mol_pdb, other_mol_mol = pair_preprocessor.process_stream(_if_dock=_if_dock, _working_dir=None)
    del pair_preprocessor

    shutil.copy2(ref_mol_pdb, _output_dir)
    shutil.copy2(other_mol_pdb, _output_dir)
    shutil.copy2(other_mol_mol, _output_dir)

    os.chdir(origin_cwd)

    return (os.path.join(_output_dir, ref_mol_pdb),
            os.path.join(_output_dir, other_mol_pdb),
            os.path.join(_output_dir, other_mol_mol))


def CLI_main():

    import argparse
    logger.setLevel(logging.INFO)

    arg_parser = argparse.ArgumentParser(description='Automatically generate pdb file for ref_mol and other_mol,'
                                                     'a mol file contains bond information for other_mol.'
                                                     'Make a ./preprocess dir on cwd.',
                                         epilog='Default: use Watvina to dock after alignment')
    arg_parser.add_argument('-r', '--receptor_pdb', type=str, help='the protein pdb file.\nSupported file format: .pdb')
    arg_parser.add_argument('-a', '--A_ref_mol', type=str,
                            help='the reference mol file.\nBond information, 3D coordinates are required, '
                                 'both implicit and explicit H atom must be recorded in file.'
                                 '\nSupported file format: .pdb, .mol, .mol2, .sdf')
    arg_parser.add_argument('-b', '--B_other_mol', type=str,
                            help='the mol file to be aligned.\nBond information is required'
                                 '\nSupported file format: .pdb, .mol, .mol2, .sdf, smiles')
    arg_parser.add_argument('-nd', '--NO_DOCK', nargs='?', const=True, default=False,
                            help='whether dock after align.')

    args = arg_parser.parse_args()
    ref_mol = read_mol_from_file(args.A_ref_mol)
    other_mol = read_mol_from_file(args.B_other_mol)
    receptor_pdb = os.path.abspath(args.receptor_pdb)
    if_dock = not args.NO_DOCK

    ref_mol_atom_nums = ref_mol.GetNumAtoms()
    if Chem.AddHs(ref_mol).GetNumAtoms() != ref_mol_atom_nums:
        raise ValueError(f'The given reference molecule does not have all the explicit H atoms.'
                         f'{Chem.AddHs(ref_mol).GetNumAtoms()} != {ref_mol_atom_nums}')
        # Molecule lost some explicit Hs.

    other_mol = Chem.AddHs(other_mol)

    logger.info(f'{"-"*21} Input Arguments {"-"*21}')
    logger.info(f"Receptor pdb: {receptor_pdb}")
    logger.info(f'If DOCK: {if_dock}')
    logger.info('-'*59)

    pair_preprocessor = PairProcess(ref_mol=ref_mol, other_mol=other_mol, protein_pdb=receptor_pdb,
                                    trust_input_conformation=True)
    pair_preprocessor.process_stream(_if_dock=if_dock)


if __name__ == '__main__':
    CLI_main()
    pass
