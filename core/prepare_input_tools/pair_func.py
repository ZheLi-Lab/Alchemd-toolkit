import lomap
import os, sys, shutil
import argparse
from typing import Set, List, Union, Iterable
import random

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFMCS
import numpy as np


class MyMol:
    def __init__(self, mol):
        self.mol = mol
        self.name = mol.GetProp('_Name')
        self.tpsa = rdMolDescriptors.CalcTPSA(mol)
        self.atom_nums = len([a.GetSymbol() for a in mol.GetAtoms() if a.GetSymbol() != "H"])
        self.hbd = rdMolDescriptors.CalcNumHBD(mol)
        self.hba = rdMolDescriptors.CalcNumHBA(mol)
        self.rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        self.mw = rdMolDescriptors.CalcExactMolWt(mol)
        self.vdw_volume = rdMolDescriptors.DoubleCubicLatticeVolume(mol=mol).GetVDWVolume()
        self.mcs_gap = 0  # Will be updated in subsequent calculations

    @staticmethod
    def calculate_mcs_gap(mol_list):
        """Calculate MCS for all molecules and update each molecule's mcs_gap"""
        params = rdFMCS.MCSParameters()
        params.BondCompareParameters.RingMatchesRingOnly = False
        params.BondCompareParameters.CompleteRingsOnly = False
        params.AtomCompareParameters.MatchChiralTag = False
        params.Timeout = 10

        try:
            mcs = rdFMCS.FindMCS([m.mol for m in mol_list], params)
            mcs_mol = mcs.queryMol
            if mcs_mol is None:
                return
            
            mcs_atoms = mcs_mol.GetNumAtoms()
            for mol in mol_list:
                mol.mcs_gap = abs(mol.atom_nums - mcs_atoms)
        except Exception as e:
            print(f"Error calculating MCS: {e}")


class Pair:
    def __init__(self, node1, node2):
        self.node1 = str(node1)
        self.node2 = str(node2)
        self.name = f"{self.node1}-{self.node2}"

    def __str__(self):
        return f"{self.node1}-{self.node2}"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if self.node1 == other.node1 and self.node2 == other.node2:
            return True
        elif self.node1 == other.node2 and self.node2 == other.node1:
            return True
        return False

    def __hash__(self):
        return hash(frozenset({self.node1, self.node2}))

    @staticmethod
    def from_string(string):
        return Pair(string.split('-')[0], string.split('-')[1])
    
    def has_node(self, node: str):
        if self.node1 == node or self.node2 == node:
            return True
        return False
    
    def flip(self):
        boo = self.node1
        self.node1 = self.node2
        self.node2 = boo


class SimilarityCalculator:
    @staticmethod
    def calculate_psa_similarity(mol1: MyMol, mol2: MyMol):
        max_psa = max(mol1.tpsa, mol2.tpsa)
        if max_psa == 0:
            return 1.0
        return 1 - abs(mol1.tpsa - mol2.tpsa) / max_psa

    @staticmethod
    def calculate_mcs_similarity(mol1: MyMol, mol2: MyMol):
        max_gap = max(mol1.mcs_gap, mol2.mcs_gap)
        if max_gap == 0:
            return 1.0
        return 1 - abs(mol1.mcs_gap - mol2.mcs_gap) / max_gap

    @staticmethod
    def calculate_atom_nums_similarity(mol1: MyMol, mol2: MyMol):
        diff = abs(mol1.atom_nums - mol2.atom_nums)
        if diff < 5:
            return (5-diff)/5
        return 1 / diff

    @staticmethod
    def calculate_hbd_hba_similarity(mol1: MyMol, mol2: MyMol):
        hbd_diff = abs(mol1.hbd - mol2.hbd)
        hba_diff = abs(mol1.hba - mol2.hba)
        
        max_hbd = max(mol1.hbd, mol2.hbd)
        max_hba = max(mol1.hba, mol2.hba)
        
        hbd_sim = 1 - (hbd_diff / max_hbd if max_hbd > 0 else 0)
        hba_sim = 1 - (hba_diff / max_hba if max_hba > 0 else 0)
        
        return (hbd_sim + hba_sim) / 2

    @staticmethod
    def calculate_rotatable_bonds_similarity(mol1: MyMol, mol2: MyMol):
        max_rot = max(mol1.rotatable_bonds, mol2.rotatable_bonds)
        if max_rot == 0:
            return 1.0
        
        return 1 - abs(mol1.rotatable_bonds - mol2.rotatable_bonds) / max_rot

    @staticmethod
    def calculate_mw_similarity(mol1: MyMol, mol2: MyMol):
        max_mw = max(mol1.mw, mol2.mw)
        return 1 - abs(mol1.mw - mol2.mw) / max_mw

    @staticmethod
    def calculate_mol_vdw_volume_similarity(mol1: MyMol, mol2: MyMol):
        max_vdw_v = max(mol1.vdw_volume, mol2.vdw_volume)
        return 1 - abs(mol1.vdw_volume - mol2.vdw_volume) / max_vdw_v

    @staticmethod
    def my_score(mol1: MyMol, mol2: MyMol):
        mcs_sim = SimilarityCalculator.calculate_mcs_similarity(mol1, mol2)
        hbd_hba_sim = SimilarityCalculator.calculate_hbd_hba_similarity(mol1, mol2)
        rot_bonds_sim = SimilarityCalculator.calculate_rotatable_bonds_similarity(mol1, mol2)
        mw_sim = SimilarityCalculator.calculate_mw_similarity(mol1, mol2)
        an_sim = SimilarityCalculator.calculate_atom_nums_similarity(mol1, mol2)
        vdw_sim = SimilarityCalculator.calculate_mol_vdw_volume_similarity(mol1, mol2)

        weights = [0.3, 0.1, 0.1, 0.1, 0.1, 0.3]
        scores = [mcs_sim, hbd_hba_sim, rot_bonds_sim, mw_sim, an_sim, vdw_sim]

        return np.average(scores, weights=weights)


def check_connected_loops(pair_list: Set[Pair]) -> Set[set]:

    def find_connected_components_dfs(pairs):
        graph = {}
        for node1, node2 in pairs:
            if node1 not in graph:
                graph[node1] = set()
            if node2 not in graph:
                graph[node2] = set()
            graph[node1].add(node2)
            graph[node2].add(node1)
        
        def dfs(node, visited, component):
            visited.add(node)
            component.add(node)
            for neighbor in graph[node]:
                if neighbor not in visited:
                    dfs(neighbor, visited, component)

        visited = set()
        components = set()

        for node in graph:
            if node not in visited:
                current_component = set()
                dfs(node, visited, current_component)
                components.add(frozenset(current_component))

        return components

    connected_components = []
    for pair in pair_list:
        connected_components.append((pair.node1, pair.node2))

    return find_connected_components_dfs(connected_components)

def get_scores_between_groups(group1, group2, mol_dict):
    return_results = []
    for i in group1:
        for j in group2:
            mol1 = mol_dict[i]
            mol2 = mol_dict[j]
            p = Pair(i, j)
            s = SimilarityCalculator.my_score(mol1=mol1, mol2=mol2)
            return_results.append((p, s))
            # print(f"Pair {p} score {s}")
    return sorted(return_results, key=lambda x:x[1], reverse=True)

def arg_parser():
    parser = argparse.ArgumentParser(description='Generate supplementary pairs')
    parser.add_argument('-a', '--all-pairs', type=str, help='Path to the pair list file', dest="all_pairs")
    parser.add_argument('-e', '--err-pairs', type=str, help='Path to the error pair list file', dest='err_pairs')
    parser.add_argument('-p', '--moles-path', type=str, dest='moles_path')
    parser.add_argument('-mpe', '--min-pairs-per-edge', type=int, dest='min_pairs_per_edge', default=3)
    return parser.parse_args()

def load_pairs(file_path) -> Set[Pair]:
    with open(file_path, 'r') as f:
        pairs = f.readlines()
    return set([Pair.from_string(p.strip()) for p in pairs])

def convert2pair(iter_boo: Union[List, Set]) -> Iterable[Pair]:
    if type(iter_boo) is list:
        return_boo = list()
        for boo in iter_boo:
            return_boo.append(Pair.from_string(boo))
    elif type(iter_boo) is set:
        return_boo = set()
        for boo in iter_boo:
            return_boo.add(Pair.from_string(boo))
    else:
        raise TypeError

    return return_boo


def main():
    args = arg_parser()
    all_pairs = load_pairs(os.path.abspath(args.all_pairs))
    err_pairs = load_pairs(os.path.abspath(args.err_pairs))
    moles_path = os.path.abspath(args.moles_path)
    min_pair_an_edge = args.min_pairs_per_edge
    max_pair_a_group_cnt_factor = 3
    threshold = 0.3

    mol_dict = {}
    my_mol_list = []
    for mol_sdf in os.listdir(moles_path):
        mol = Chem.SDMolSupplier(os.path.join(moles_path, mol_sdf), removeHs=False)[0]
        my_mol = MyMol(mol)
        mol_dict[mol.GetProp("_Name")] = my_mol
        my_mol_list.append(my_mol)
    MyMol.calculate_mcs_gap(my_mol_list)

    all_moles = lomap.DBMolecules(directory=moles_path, parallel=4, cutoff=0.2, max=50, )
    all_moles.read_molecule_files()
    all_moles.build_matrices()
    s = all_moles.strict_mtx.to_numpy_2D_array()
    lomap_avail_pairs = set()
    rows, cols = s.shape
    for i in range(rows):
        for j in range(cols):
            if s[i, j] >= threshold:
                node1 = all_moles[i].getMolecule().GetProp('_Name')
                node2 = all_moles[j].getMolecule().GetProp('_Name')
                lomap_avail_pairs.add((Pair(node1, node2), s[i, j]))  #(Pair, float)

    success_pairs = all_pairs - err_pairs

    groups = check_connected_loops(success_pairs)
    
    for avail_pair, score in lomap_avail_pairs:
        cnt_group_1 = None
        cnt_group_2 = None
        for i, group in enumerate(groups):
            if avail_pair.node1 in group:
                cnt_group_1 = i
            if avail_pair.node2 in group:
                cnt_group_2 = i
        if cnt_group_1 is None or cnt_group_2 is None:
            if cnt_group_1 is None:
                groups.add(frozenset([avail_pair.node1]))
            if cnt_group_2 is None:
                groups.add(frozenset([avail_pair.node2]))

    for i,group in enumerate(groups):
        print(f"group {i}: {group}")
    group_cnt_nums = {str(i): 0 for i in range(len(groups))}  # This dict records how many outgoing pairs each group has

    _to_be_discard_pairs = set()
    for pair, score in lomap_avail_pairs:
        if pair in err_pairs:
            _to_be_discard_pairs.add((pair, score))
    lomap_avail_pairs = lomap_avail_pairs - _to_be_discard_pairs
    print(f"lomap provides {len(lomap_avail_pairs)} pairs")

    cnt_pair_dict = {}
    for avail_pair, score in lomap_avail_pairs:
        cnt_group_1 = None
        cnt_group_2 = None
        for i, group in enumerate(groups):
            if avail_pair.node1 in group:
                cnt_group_1 = i
            if avail_pair.node2 in group:
                cnt_group_2 = i
        cnt_relationship = Pair(cnt_group_1, cnt_group_2)
        if cnt_relationship not in cnt_pair_dict:
            cnt_pair_dict[cnt_relationship] = set()
        cnt_pair_dict[cnt_relationship].add((avail_pair, score))

    keys_to_remove = [key for key in cnt_pair_dict.keys() if key.node1 == key.node2]
    for key in keys_to_remove:
        del cnt_pair_dict[key]

    groups_num = len(groups)

    if len(cnt_pair_dict.keys()) == 0:  # all lomap pairs connect in their own group.
        lomap_output_cnt = [[]]
    else:
        lomap_output_cnt = check_connected_loops(set(cnt_pair_dict.keys()))
    _lomap_cnt = list(lomap_output_cnt)[0]
    sup_pairs = set()

    for _cnt, _pairs in cnt_pair_dict.items():
        _sup_pairs = [i[0] for i in set(sorted(_pairs, key=lambda x: x[1], reverse=True)[:min_pair_an_edge])]
        sup_pairs.update(_sup_pairs)
        group_cnt_nums[_cnt.node1] += len(_sup_pairs)
        group_cnt_nums[_cnt.node2] += len(_sup_pairs)
        print(f"Added {_sup_pairs} connecting {_cnt} [Lomap]")
    print(f"after lomap: {group_cnt_nums}")

    if len(lomap_output_cnt) == 1 and len(_lomap_cnt) == groups_num:
        # All supplementary pairs used to connect groups can be found in lomap
        pass
    else:
        sorted_groups = sorted([(group, i) for i, group in enumerate(groups)], key=lambda x: len(x[0]), reverse=True)
        for i, (group, idx) in enumerate(sorted_groups):
            if len(group) == 1:
                break
            for j in range(i+1, groups_num):
                other_group, other_idx = sorted_groups[j]
                now_cnt = Pair(idx, other_idx)
                if now_cnt in cnt_pair_dict.keys() and len(cnt_pair_dict[now_cnt]) >= min_pair_an_edge:
                    continue
                if group_cnt_nums[now_cnt.node1] >= max_pair_a_group_cnt_factor * len(group):
                    break
                if group_cnt_nums[now_cnt.node2] >= max_pair_a_group_cnt_factor * len(other_group):
                    continue
                my_score_results = get_scores_between_groups(group, other_group, mol_dict=mol_dict)
                _sup_pairs = set()
                for _pair, _score in my_score_results:
                    if len(_sup_pairs) >= min_pair_an_edge:
                        break
                    if _pair not in err_pairs:
                        _sup_pairs.add(_pair)
                        print(_pair, _score)
                sup_pairs.update(_sup_pairs)
                group_cnt_nums[now_cnt.node1] += len(_sup_pairs)
                group_cnt_nums[now_cnt.node2] += len(_sup_pairs)
                print(f"Added {_sup_pairs} connecting {now_cnt} [Custom]")

    all_nodes = set()
    min_pair_per_node = 3
    for pair in success_pairs:
        all_nodes.add(pair.node1)
        all_nodes.add(pair.node2)
    for pair in sup_pairs:
        all_nodes.add(pair.node1)
        all_nodes.add(pair.node2)

    node_pair_count = {node: 0 for node in all_nodes}
    for pair in success_pairs:
        node_pair_count[pair.node1] += 1
        node_pair_count[pair.node2] += 1
    for pair in sup_pairs:
        node_pair_count[pair.node1] += 1
        node_pair_count[pair.node2] += 1

    for node, count in node_pair_count.items():
        if count >= min_pair_per_node:
            continue

        other_nodes = list(all_nodes - {node})

        candidate_pairs = []
        for other_node in other_nodes:
            pair = Pair(node, other_node)
            if (pair not in err_pairs) and (pair not in success_pairs) and (pair not in sup_pairs):
                mol1 = mol_dict[node]
                mol2 = mol_dict[other_node]
                score = SimilarityCalculator.my_score(mol1, mol2)
                candidate_pairs.append((pair, score))
        for avail_pair, score in lomap_avail_pairs:
            if avail_pair.has_node(node):
                candidate_pairs.append((avail_pair, 1))  # Use lomap2 pair first

        candidate_pairs.sort(key=lambda x: x[1], reverse=True)
        needed_pairs = min_pair_per_node - count
        for pair, _ in candidate_pairs[:needed_pairs]:
            if pair not in sup_pairs and pair not in success_pairs:
                sup_pairs.add(pair)
                node_pair_count[pair.node1] +=1
                node_pair_count[pair.node2] +=1
                print(f"Added cyc pair {pair} for node {node} [cyc]")

    with open(os.path.join(os.path.dirname(os.path.abspath(args.all_pairs)), "supplement_pair.lst"),"w") as f:
        for _sup_pair in sup_pairs:
            f.write(f"{str(_sup_pair)}\n")

def generate_legal_visit_queue(pairs: Set[Pair], start_nodes: Set[str], allow_sep: bool = False) -> List[Pair]:
    """
    Generate a legal visit queue with minimal pair flips
    
    Args:
        pairs: Set of all pairs to be visited
        start_nodes: Set of starting nodes that are already visited
        
    Returns:
        list: Queue of pairs in visit order
    """
    # Check connectivity first
    def check_connectivity(pairs: Set[Pair]) -> bool:
        # Build adjacency list
        graph = {}
        all_nodes = set()
        for pair in pairs:
            all_nodes.add(pair.node1)
            all_nodes.add(pair.node2)
            if pair.node1 not in graph:
                graph[pair.node1] = set()
            if pair.node2 not in graph:
                graph[pair.node2] = set()
            graph[pair.node1].add(pair.node2)
            graph[pair.node2].add(pair.node1)
        
        # Start DFS from any node
        if not all_nodes:
            return True
            
        visited = set()
        stack = [next(iter(all_nodes))]
        
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            for neighbor in graph.get(node, set()):
                if neighbor not in visited:
                    stack.append(neighbor)
        
        # If visited nodes count equals total nodes count, the graph is connected
        return len(visited) == len(all_nodes)
    
    # Check connectivity
    if not check_connectivity(pairs):
        if allow_sep:
            return list(pairs)
        raise ValueError("Given pairs cannot form a connected graph")
    
    visited = set(start_nodes)
    queue = []
    remaining_pairs = set(pairs)
    
    while remaining_pairs:
        # Find all immediately accessible pairs
        available_pairs = set()
        for pair in remaining_pairs:
            if pair.node1 in visited and pair.node2 not in visited:
                available_pairs.add((pair, False))  # (pair, needs_flip)
            elif pair.node2 in visited and pair.node1 not in visited:
                available_pairs.add((pair, True))   # (pair, needs_flip)
            elif pair.node1 in visited and pair.node2 in visited:
                available_pairs.add((pair, False))
        
        if not available_pairs:
            raise ValueError("Cannot generate legal visit queue: unreachable nodes exist")
        
        # Select the first available pair
        pair, needs_flip = next(iter(available_pairs))
        if needs_flip:
            # If flip is needed, create new pair
            new_pair = Pair(pair.node2, pair.node1)
            queue.append(new_pair)
            visited.add(pair.node1)
        else:
            queue.append(pair)
            visited.add(pair.node2)
        
        remaining_pairs.remove(pair)
    
    return queue

if __name__ == "__main__":
    main()