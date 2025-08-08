import os
import sys
import pickle
import copy
import time
import shutil
import subprocess
import logging
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict, Set

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from prepare_input_tools.pdb_preparer import PrepareProteinPDB
from prepare_input_tools.pair_pdb_generator import GeneratePairPDB, MolCoordModer, WatVina_api
from prepare_input_tools.pair_pdb_generator import pdb_formatter, pdb_info_clean, working_dir, calculate_formal_charge
from prepare_input_tools.lomap_api import GeneratePairMap
from analyze_tools.SettingManager import get_default_path
from prepare_input_tools.gen_func.collision_checker import CollisionChecker
from prepare_input_tools.pair_func import Pair, generate_legal_visit_queue, load_pairs, convert2pair

# Constants
UNCOMMON_AMINO_ACIDS = ['TLA', 'ZNA', 'TPO']
# GENAMBRBFE_PATH = os.path.join(get_default_path('GenambRBFE_Path'), 'genambRBFE.sh')
GENAMBRBFE_BASE_PATH = get_default_path('GenambRBFE_Path')
GENAMBRBFE_PATH = os.path.join(GENAMBRBFE_BASE_PATH, 'genambRBFE')

def setup_logger(logger_name: str, log_level: int = logging.INFO) -> logging.Logger:
    """Setup standardized logger for CADD toolkit."""
    logger = logging.getLogger(logger_name)
    logger.setLevel(log_level)

    if logger.handlers:
        logger.handlers.clear()

    is_wsl = _detect_wsl_environment()
    
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)

    if is_wsl:
        formatter = logging.Formatter('%(levelname)s: %(message)s')
    else:
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    return logger

def _detect_wsl_environment() -> bool:
    """Detect if running in WSL environment."""
    try:
        with open('/proc/version', 'r') as f:
            return 'microsoft' in f.read().lower()
    except (OSError, IOError):
        return False

logger = setup_logger(__name__, logging.INFO)


class PairMap:
    """Graph-based pair mapping for molecular transformations."""
    
    def __init__(self, pair_list: List[Pair]):
        self.pair_list = pair_list
        self.connect_dict: Dict[str, Set[str]] = {}
        self._build_connection_graph()
    
    def _build_connection_graph(self) -> None:
        """Build connection graph from pair list."""
        for pair in self.pair_list:
            self._add_connection(pair.node1, pair.node2)
            self._add_connection(pair.node2, pair.node1)
    
    def _add_connection(self, node_a: str, node_b: str) -> None:
        """Add bidirectional connection between nodes."""
        if node_a not in self.connect_dict:
            self.connect_dict[node_a] = set()
        self.connect_dict[node_a].add(node_b)

    @staticmethod
    def _is_same_pair(pair1: Pair, pair2: Pair) -> bool:
        """Check if two pairs represent the same connection."""
        return pair1 == pair2

    @staticmethod
    def can_generate_from_nodes(nodes: List[str], pairs: List[Pair], 
                               strict_direction: bool = False) -> Tuple[bool, List[Pair]]:
        """Check if all pairs can be generated from starting nodes."""
        reached_nodes = set(nodes)
        valid_pairs = set()
        remaining_pairs = pairs.copy()
        
        found_new = True
        while found_new:
            found_new = False
            for pair in remaining_pairs:
                if PairMap._can_reach_pair(pair, reached_nodes, strict_direction):
                    valid_pairs.add(pair)
                    remaining_pairs.remove(pair)
                    reached_nodes.update([pair.node1, pair.node2])
                    found_new = True
                    break
        
        return len(remaining_pairs) == 0, list(valid_pairs)
    
    @staticmethod
    def _can_reach_pair(pair: Pair, reached_nodes: Set[str], strict_direction: bool) -> bool:
        """Check if a pair can be reached from current nodes."""
        if strict_direction:
            return pair.node1 in reached_nodes
        return pair.node1 in reached_nodes or pair.node2 in reached_nodes

    def BFS_from_node(self, node) -> Tuple[list, dict]:
        if node not in self.connect_dict.keys():
            raise ValueError(f'Node {node} is not connected to any other node')

        visit_queue = [node]
        visited = set()
        BFS_pair_queue = []

        BFS_edge_layer = {}
        layer_counter = 0

        while visit_queue:
            now = visit_queue.pop(0)
            visited.add(now)
            layer_counter += 1
            for neighbor in self.connect_dict[now]:
                if neighbor not in visited:
                    if neighbor not in visit_queue:
                        visit_queue.append(neighbor)
                    BFS_pair_queue.append(f"{now}-{neighbor}")
                    BFS_edge_layer[f"{now}-{neighbor}"] = layer_counter
            # print(f"visited: {visited} visit_queue: {visit_queue}")

        return BFS_pair_queue, BFS_edge_layer

    def BFS_from_nodes(self, node_list: List[str]) -> Tuple[list, dict]:
        min_layers = {}
        for node in node_list:
            if node not in self.connect_dict.keys():
                continue
            _, edge_layer = self.BFS_from_node(node)
            for edge, layer in edge_layer.items():

                if edge not in min_layers.keys():
                    # No edge record in min layer.
                    min_layers[edge] = layer
                    continue
                if layer < min_layers[edge]:
                    min_layers[edge] = layer

        return_pair_queue = [i[0] for i in sorted(min_layers.items(), key=lambda item: item[1])]

        remove_pairs = []
        for i in range(len(return_pair_queue)):
            if return_pair_queue[i] in remove_pairs:
                continue
            for j in range(i+1, len(return_pair_queue)):
                if self._is_same_pair(return_pair_queue[i], return_pair_queue[j]):
                    remove_pairs.append(return_pair_queue[j])
                    break

        for pair in remove_pairs:
            return_pair_queue.remove(pair)
            min_layers.pop(pair)

        return_pair_queue = convert2pair(return_pair_queue)
        return return_pair_queue, min_layers

    def visualize(self, output_png, nodes: list):
        try:
            import networkx as nx
            import matplotlib.pyplot as plt
            G = nx.Graph([(i.node1, i.node2) for i in self.pair_list])

            edge_queue, edge_layer = self.BFS_from_nodes(nodes)
            logger.debug(edge_queue)

            for u, v, d in G.edges(data=True):
                if f'{u}-{v}' == str(edge_queue):
                    d['order'] = edge_layer[f'{u}-{v}']
                elif f'{v}-{u}' == str(edge_queue):
                    d['order'] = edge_layer[f'{v}-{u}']

            edges, orders = zip(*nx.get_edge_attributes(G, 'order').items())

            fig = plt.figure(figsize=(17, 17))
            nx.draw(G, edgelist=edges, edge_color=orders, width=6,
                    node_color='#0e1730', node_size=1000, with_labels=True,
                    font_color='w', font_weight='bold', font_size='20',
                    edge_cmap=plt.cm.Greys)
            fig.set_facecolor('#9abbc2')
            plt.savefig(output_png, dpi=500)
            plt.close()

        except ImportError:
            logger.warning('Please install networkx package to use this function.')


@dataclass
class MolBase:
    ref_moles: List[Chem.Mol]
    ac_moles: List[Chem.Mol]
    protein_pdb: str
    nsaa_params_dir: str = None

    def __post_init__(self):
        self.all_moles = self.ref_moles + self.ac_moles
        self.all_moles_dict = {mol.GetProp('_Name'): mol for mol in self.all_moles}
        if self.nsaa_params_dir is not None:
            self.nsaa_params_dir = os.path.abspath(self.nsaa_params_dir)
    
    def get_molecule(self, name: str) -> Optional[Chem.Mol]:
        """Get molecule by name."""
        return self.all_moles_dict.get(name)


@dataclass
class PairInfo:
    pair_name: str

    lig_a_pdb: str
    lig_b_pdb: str
    protein_pdb: str

    lig_a_chg: int = 0
    lig_b_chg: int = 0


class RawFileManager:
    def __init__(self, mol_base: MolBase, output_dir, pair_queue: List):

        self.mol_base = mol_base
        self._output_dir = os.path.abspath(output_dir)
        self._preprocess_dir = os.path.join(self._output_dir, 'preprocess')
        self.prepare_dir = os.path.join(self._output_dir, 'prepare')
        os.makedirs(self._preprocess_dir, exist_ok=True)

        self.pair_generate_queue = pair_queue
        self.run_list = []

        self.pair_info_dict = {}
        self.sbond_lst = None
        self.dock_pdb = None
        self.prepared_protein_pdb = os.path.join(self._preprocess_dir, 'protein.pdb')
        self.collision_checker = None

    def prepare_protein_pdb(self, if_prepare: bool = True):
        prepare_dir = os.path.join(self._preprocess_dir, 'rec_prepare')
        os.makedirs(prepare_dir, exist_ok=True)
        if not if_prepare:
            shutil.copy2(self.mol_base.protein_pdb, self.prepared_protein_pdb)
            self.dock_pdb = self.prepared_protein_pdb
            self.sbond_lst = None
            return False

        try:
            self.dock_pdb, self.sbond_lst = PrepareProteinPDB(
                input_pdb=self.mol_base.protein_pdb,
                output_pdb=self.prepared_protein_pdb,
                prepare_dir=prepare_dir
            )
        except ValueError as e:
            logger.warning(f"PDB preparation failed: {e}. Using original PDB.")
            shutil.copy2(self.mol_base.protein_pdb, self.prepared_protein_pdb)
            self.dock_pdb = self.prepared_protein_pdb
            self.sbond_lst = None
        logger.debug(f"Prepared protein PDB: {self.dock_pdb}")

    def _prepare_folders(self):
        gen_mol_dir = os.path.join(self._preprocess_dir, 'gen_mol')
        if os.path.exists(gen_mol_dir):
            shutil.rmtree(gen_mol_dir)
        os.makedirs(gen_mol_dir, exist_ok=True)

        pair_pdb_collect_dir = os.path.join(self._preprocess_dir, 'pair')
        if os.path.exists(pair_pdb_collect_dir):
            shutil.rmtree(pair_pdb_collect_dir)
        os.makedirs(pair_pdb_collect_dir, exist_ok=True)

        gen_file_dir = os.path.join(self._preprocess_dir, 'gen_file')
        if os.path.exists(gen_file_dir):
            shutil.rmtree(gen_file_dir)
        os.makedirs(gen_file_dir, exist_ok=True)
        
        return gen_mol_dir, pair_pdb_collect_dir, gen_file_dir

    @staticmethod
    def _find_special_amino_acids(pdb_file: str) -> Optional[str]:
        """Find special amino acids in PDB file."""
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    residue_name = line[17:20]
                    if residue_name in UNCOMMON_AMINO_ACIDS:
                        return residue_name
        except (OSError, IOError) as e:
            logger.error(f"Failed to read PDB file {pdb_file}: {e}")
        return None

    def prepare_pair_files(self, if_cs_align=True, if_dock=True, if_auto_run=False,
                           **kwargs):
        origin_cwd = os.getcwd()
        os.chdir(self._preprocess_dir)

        gen_mol_dir, pair_pdb_collect_dir, gen_file_dir = self._prepare_folders()
        given_mol_dir = os.path.join(self._preprocess_dir, 'given_mol')
        if os.path.exists(given_mol_dir):
            shutil.rmtree(given_mol_dir)
        os.makedirs(given_mol_dir, exist_ok=True)

        def get_ref_mol_file(node_name, info=False):
            if os.path.exists(os.path.join(given_mol_dir, f"{node_name}.mol")):
                if info:
                    logger.debug(f'Using given mol\'s conformer: {node_name}.mol')
                return os.path.join(given_mol_dir, f"{node_name}.mol")
            return os.path.join(gen_mol_dir, f"{node_name}.mol")

        mol_base = self.mol_base
        if 'mol_base' in kwargs.keys():
            mol_base = kwargs['mol_base']
        if if_cs_align:
            collision_ref_moles = []
            for ref_mol in mol_base.ref_moles:
                Chem.MolToMolFile(ref_mol, os.path.join(given_mol_dir, f"{ref_mol.GetProp('_Name')}.mol"))
                if self.collision_checker is not None:
                    collision_idxes = self.collision_checker.get_collision_idxes(ref_mol)
                    if len(collision_idxes) > 0:
                        collision_ref_moles.append(ref_mol.GetProp('_Name'))
            if len(collision_ref_moles) > 0:
                collision_ref_moles.sort()
                logger.warning(f"Protein collision detected on following reference moles: {collision_ref_moles}")
        # This part is used in generate
        sections = 1
        if if_auto_run:
            sections += 1

        if if_cs_align:
            generated_moles = []

        for i, pair in enumerate(self.pair_generate_queue):
            logger.info(f"Generating pair {pair} ({i+1}/{len(self.pair_generate_queue)}) Section [1/{sections}]")
            ref_node = pair.node1
            other_node = pair.node2
            pair_pdb_dir = os.path.join(pair_pdb_collect_dir, pair.name)
            os.makedirs(pair_pdb_dir, exist_ok=True)

            ref_mol_pdb = os.path.join(pair_pdb_dir, f'{ref_node}.pdb')
            other_mol_pdb = os.path.join(pair_pdb_dir, f'{other_node}.pdb')
            self.pair_info_dict[pair] = PairInfo(pair_name=pair.name,
                                                 lig_a_pdb=os.path.abspath(ref_mol_pdb),
                                                 lig_b_pdb=os.path.abspath(other_mol_pdb),
                                                 protein_pdb=os.path.abspath(self.prepared_protein_pdb),
                                                 lig_a_chg=calculate_formal_charge(mol_base.all_moles_dict[ref_node]),
                                                 lig_b_chg=calculate_formal_charge(mol_base.all_moles_dict[other_node]))

            if if_cs_align:
                if not os.path.exists(get_ref_mol_file(ref_node)):
                    raise FileNotFoundError(f"Pair {pair.name} is not accessible in this order."
                                            f"Ref mol file: {get_ref_mol_file(ref_node)} does not exist.")

                ref_mol_names = [mol.GetProp('_Name') for mol in mol_base.ref_moles]
                trust_input_conformation = other_node in ref_mol_names

                ref_mol = Chem.MolFromMolFile(get_ref_mol_file(ref_node, info=True), removeHs=False)
                other_mol = mol_base.all_moles_dict[other_node]
                ref_mol_pdb, other_mol_pdb, other_mol_mol = GeneratePairPDB(ref_mol=ref_mol,
                                                                            other_mol=other_mol,
                                                                            protein_pdb=self.dock_pdb,
                                                                            _output_dir=pair_pdb_dir,
                                                                            _if_dock=if_dock,
                                                                            trust_input_conformation=trust_input_conformation)

                if not os.path.exists(get_ref_mol_file(other_node)):
                    shutil.copy2(other_mol_mol, get_ref_mol_file(other_node))

                tmp_mol = Chem.MolFromMolFile(other_mol_mol, removeHs=False)
                tmp_mol.SetProp('_Name', pair.name)
                generated_moles.append(tmp_mol)
                del tmp_mol

                os.remove(os.path.join(pair_pdb_dir, os.path.basename(other_mol_mol)))

            else:
                ref_mol = mol_base.all_moles_dict[ref_node]
                other_mol = mol_base.all_moles_dict[other_node]

                Chem.MolToPDBFile(ref_mol, ref_mol_pdb)
                Chem.MolToPDBFile(other_mol, other_mol_pdb)

            pdb_info_clean(other_mol_pdb, res_name='LGO')
            pdb_info_clean(ref_mol_pdb, res_name='LGR')

            pdb_formatter(other_mol_pdb, atom_name_formatter='title')
            pdb_formatter(ref_mol_pdb, atom_name_formatter='title')

            sys.stdout.flush()
            time.sleep(1)

        if os.path.exists(self.prepare_dir):
            shutil.rmtree(self.prepare_dir)
        shutil.copytree(pair_pdb_collect_dir, self.prepare_dir)

        if if_cs_align:
            all_generate_mol_sdf = Chem.SDWriter(os.path.join(self._preprocess_dir, 'all_generate_mol.sdf'))
            sorted_generated_moles = sorted(generated_moles, key=lambda x: x.GetProp('_Name').split('-')[1], reverse=False)
            collision_generated_moles = []
            for mol in sorted_generated_moles:
                all_generate_mol_sdf.write(mol)
                if self.collision_checker is not None:
                    collision_idxes = self.collision_checker.get_collision_idxes(mol)
                    if len(collision_idxes) > 0:
                        collision_generated_moles.append(mol.GetProp('_Name'))
            if len(collision_generated_moles) > 0:
                logger.warning(f"Protein collision detected on following pairs: {collision_generated_moles}")
            all_generate_mol_sdf.close()

        os.chdir(origin_cwd)

    def load_prepare_from_custom_folder(self, folder):
        """
        Assume the folder stores a list of pair folder.
        -folder - 'A-B'
                - 'C-D'
        :param folder: 
        :return: 
        """
        tmp_prepare_dir = os.path.abspath(folder)
        pair_list = [pair for pair in os.listdir(tmp_prepare_dir) if pair.find('-') != -1]
        for pair in pair_list:
            pair_pdb_dir = os.path.join(tmp_prepare_dir, pair)
            ref_node = pair.split('-')[0]
            other_node = pair.split('-')[1]
            _pair_Pair = Pair.from_string(pair)
            os.makedirs(os.path.join(self.prepare_dir, pair), exist_ok=True)

            ref_mol_pdb = os.path.join(pair_pdb_dir, f'{ref_node}.pdb')
            other_mol_pdb = os.path.join(pair_pdb_dir, f'{other_node}.pdb')

            copyed_ref_mol_pdb = os.path.join(self.prepare_dir, pair, f'{ref_node}.pdb')
            copyed_other_mol_pdb = os.path.join(self.prepare_dir, pair, f'{other_node}.pdb')

            shutil.copy2(ref_mol_pdb, copyed_ref_mol_pdb)
            shutil.copy2(other_mol_pdb, copyed_other_mol_pdb)
            pdb_formatter(copyed_ref_mol_pdb, atom_name_formatter='title')
            pdb_formatter(copyed_other_mol_pdb, atom_name_formatter='title')

            ref_mol = Chem.MolFromPDBFile(copyed_ref_mol_pdb, removeHs=False)
            other_mol = Chem.MolFromPDBFile(copyed_other_mol_pdb, removeHs=False)
            lig_a_chg = calculate_formal_charge(ref_mol)
            lig_b_chg = calculate_formal_charge(other_mol)
            if ref_node not in self.mol_base.all_moles_dict.keys():
                self.mol_base.all_moles_dict[ref_node] = ref_mol
            if other_node not in self.mol_base.all_moles_dict.keys():
                self.mol_base.all_moles_dict[other_node] = other_mol
            self.pair_info_dict[_pair_Pair] = PairInfo(pair_name=pair,  # Notice here pair type is !!!str!!!
                                                        lig_a_pdb=os.path.abspath(ref_mol_pdb),
                                                        lig_b_pdb=os.path.abspath(other_mol_pdb),
                                                        protein_pdb=os.path.abspath(self.prepared_protein_pdb),
                                                        lig_a_chg=lig_a_chg,
                                                        lig_b_chg=lig_b_chg
                                                        )

    def dump_infos_to_pickle(self):
        pair_list = list(self.pair_info_dict.keys())
        if len(self.mol_base.all_moles_dict) > 0:
            all_moles = copy.deepcopy(self.mol_base.all_moles_dict)
        else:
            all_moles = None
        info = {'pair_list': pair_list,
                'all_moles': all_moles}
        info_pickle_file = os.path.join(self._output_dir, 'info.pickle')
        with open(info_pickle_file, 'wb') as f:
            pickle.dump(info, f)

    def initial_prepare_folder(self) -> List[Tuple[str, str]]:
        prepare_queue_lst = os.path.join(self._output_dir, 'prepare_queue.lst')
        special_AA = self._find_special_amino_acids(self.prepared_protein_pdb)

        if self.sbond_lst is not None and os.path.exists(self.sbond_lst):
            sbond_arg = f'-s {self.sbond_lst}'
        else:
            sbond_arg = ''
        if special_AA is not None:
            special_AA_arg = f'-A "{special_AA}"'
        else:
            special_AA_arg = ''

        if self.mol_base.nsaa_params_dir is not None:
            if not os.path.exists(self.mol_base.nsaa_params_dir):
                raise ValueError(f"NSAA parameters directory {self.mol_base.nsaa_params_dir} does not exist")

            subfolders = [f for f in os.listdir(self.mol_base.nsaa_params_dir) if os.path.isdir(os.path.join(self.mol_base.nsaa_params_dir, f))]
            if not subfolders:
                raise ValueError(f"No subfolders found in {self.mol_base.nsaa_params_dir}")

            for subfolder in subfolders:
                src = os.path.join(self.mol_base.nsaa_params_dir, subfolder)
                dst = os.path.join(GENAMBRBFE_BASE_PATH, subfolder)
                if os.path.exists(dst):
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)

            extra_args = ' '.join(subfolders)
            if special_AA_arg:
                special_AA_arg += f' {extra_args}'
            else:
                special_AA_arg = f'-A "{extra_args}"'

        run_list = []
        collision_checker_valid = False
        try:
            tmp_collision_checker = CollisionChecker.from_pdb_file(self.mol_base.protein_pdb)
            collision_checker_valid = True
        except:
            pass
        # Must be prepare_protein_pdb here.

        prepare_queue_writer = open(prepare_queue_lst, 'w')

        for pair, pair_info in self.pair_info_dict.items():
            submit_bash = os.path.join(self.prepare_dir, pair.name, 'submit.sh')
            lig_a_basename = os.path.basename(pair_info.lig_a_pdb)
            lig_b_basename = os.path.basename(pair_info.lig_b_pdb)

            lig_a_mol = Chem.MolFromPDBFile(pair_info.lig_a_pdb, removeHs=False, sanitize=False)
            lig_b_mol = Chem.MolFromPDBFile(pair_info.lig_b_pdb, removeHs=False, sanitize=False)
            pair_rec_pdb_file = os.path.join(self.prepare_dir, pair.name, 'protein_del_collision_water.pdb')
            using_rec_pdb_file = pair_info.protein_pdb
            if collision_checker_valid:
                if tmp_collision_checker.get_water_remove_pdb(moles=[lig_a_mol, lig_b_mol], output_pdb_file=pair_rec_pdb_file,
                                                              if_reprepare=True):
                    using_rec_pdb_file = "protein_del_collision_water.pdb"
                    logger.info(f'Detect water collision in complex, using pdb deleted collision water. ({pair})')
            with open(submit_bash, 'w') as f:
                f.write('#!/bin/bash\n')
                f.write('printf "BEGIN AT:%s \\n" "$(date -R)" > prepare.log\n')
                f.write(f'{GENAMBRBFE_PATH} -r {using_rec_pdb_file} -a {lig_a_basename} -b {lig_b_basename} '
                        f'-c {pair_info.lig_a_chg} -C {pair_info.lig_b_chg} {special_AA_arg} {sbond_arg} '
                        f'>> prepare.log 2>&1\n')
                f.write('printf "END AT:%s \\n" "$(date -R)" >> prepare.log\n')
            prepare_queue_writer.write(f"{os.path.join(self.prepare_dir, pair.name)} {submit_bash}\n")
            run_list.append((os.path.join(self.prepare_dir, pair.name), submit_bash))

        prepare_queue_writer.close()
        
        return run_list
        
    @staticmethod
    def auto_run_prepare(run_list):
        for i, (run_dir, submit_bash) in enumerate(run_list):
            logger.info(f'Preparing pair {os.path.split(run_dir)[1]} ({i + 1}/{len(run_list)}) Section [2/2]')
            os.chdir(run_dir)
            p = subprocess.Popen(['/bin/bash', submit_bash])
            p.wait()


def Initial_ligands(ligand_file, ref_ligand_list) -> Tuple[List[Chem.Mol], List[Chem.Mol]]:
    from rdkit.Chem import rdForceFieldHelpers
    from rdkit.Chem import rdDistGeom
    if ligand_file.endswith('.sdf'):
        ligands = []
        mol_supplier = Chem.SDMolSupplier(ligand_file, removeHs=False)
        if mol_supplier is None:
            raise Exception('Something wrong in given ref sdf file.')
        for mol in mol_supplier:
            if mol.GetNumAtoms() != Chem.AddHs(mol).GetNumAtoms():
                raise ValueError('Given molecule does not have all implicit H.')
            ligands.append(mol)
    else:
        raise ValueError(f'Unsupported file of type {ligand_file} for ligands.')

    _ref_ligand_list = []
    if ref_ligand_list is not None:
        with open(ref_ligand_list, 'r') as f:
            for line in f:
                _ref_ligand_list.append(line.strip())

    ref_moles = []
    ac_moles_un_init = []
    ac_moles = []
    if len(_ref_ligand_list) > 0:
        for mol in ligands:
            if mol.GetProp('_Name') in _ref_ligand_list:
                ref_moles.append(mol)
            else:
                ac_moles_un_init.append(mol)
    else:
        # print('all ligands as reference ligands')
        ref_moles = ligands

    if len(ac_moles_un_init) > 0:
        for mol_un_init in ac_moles_un_init:
            if not MolCoordModer.has_3d_conformer(mol_un_init):
                mol_init = Chem.AddHs(mol_un_init, addCoords=True)
                rdDistGeom.EmbedMolecule(mol_init)
                rdForceFieldHelpers.MMFFOptimizeMolecule(mol_init)
            else:
                mol_init = Chem.AddHs(mol_un_init, addCoords=True)

            ac_moles.append(mol_init)
            
    return ref_moles, ac_moles


def get_mol_vdw_volume(_mol: Chem.Mol):
    _boo = rdMolDescriptors.DoubleCubicLatticeVolume(mol=_mol)
    return _boo.GetVDWVolume()


def CLI_main():
    import argparse
    import textwrap
    logger.setLevel(logging.INFO)

    arg_parser = argparse.ArgumentParser(description='Auto generate mol pair map then generate pdb files. '
                                                     'Make ./prepare dir and ./preprocess dir under cwd.',
                                         formatter_class=argparse.RawTextHelpFormatter)
    arg_parser.add_argument('-l', '--ligands', type=str, dest='ligands', default=None,
                            help=textwrap.dedent('''\
                                path to ligands, like ./ref_mol.sdf
                                Bond info are REQUIRED.
                                The implicit H and explicit H are both REQUIRED.
                                [Supported formats]: .sdf'''))
    arg_parser.add_argument('-p', '--protein', type=str, dest='protein',
                            help=textwrap.dedent('''\
                                path to protein pdb, like ./protein.pdb
                                [Supported formats]: .pdb'''))
    arg_parser.add_argument('-af', '--advance-folder', type=str, dest='advance_folder', default=None,
                            help=textwrap.dedent('''\
                                path to a well organized prepare folder.
                                This folder should hold all the folders of pairs and save their molecular pdb files
                                in these subfolders.
                                [Example] A subfolder named A-B, where A.pdb is saved along with B.pdb'''))
    arg_parser.add_argument('-r', '--ref-ligand-list', type=str, default=None, dest='ref_ligand_list',
                            help=textwrap.dedent('''\
                                file contains the name of reference ligands, like ./ref_ligand.lst
                                [Supported formats]: .sdf'''))
    arg_parser.add_argument('-pl', '--pair-list', type=str, default=None, dest='pair_list',
                            help=textwrap.dedent('''\
                                file contains pairs info, like ./pairs.lst\
                                '''))
    arg_parser.add_argument('-o', '--output-dir', type=str, default=None, dest='output_dir',
                            help=textwrap.dedent('''\
                                where you want to save the processed files.
                                defaults to a new ./prepare folder in the current directory.\
                                '''))
    arg_parser.add_argument('-sp', '--strict-pair-list', dest='strict_pair_list',
                            nargs='?', const=True, default=False,
                            help='enforce strict adherence to user-defined pair list.'
                                 'Only avail when --cs-align and -p are given.')
    arg_parser.add_argument('--cs-align', dest='cs_align', nargs='?', const=True, default=False,
                            help='whether using coordinate system alignment to generate other mol from ref_mol in same pair')
    arg_parser.add_argument('-d', '--dock', nargs='?', const=True, default=False,
                            help='whether docking molecules')
    arg_parser.add_argument('-f', '--flip', nargs='?', const=True, default=False, dest='flip_by_vdw',
                            help='whether flip pair by vdw volume follow bigger to smaller')
    arg_parser.add_argument('-as', '--allow-sep', nargs='?', const=True, default=False, dest='allow_sep',
                            help='whether allow given mutiple unconnect pairs')
    arg_parser.add_argument('-npp', '--not-prepare-protein', dest='not_prepare_protein',
                            nargs='?', const=True, default=False,
                            help='whether preparing protein')
    arg_parser.add_argument('-auto', dest='auto_run',
                            nargs='?', const=True, default=False,
                            help='auto prepare top and crd files after preparing pair files.')
    arg_parser.add_argument('-npd', '--nsaa-params-dir', type=str, dest='nsaa_params_dir', default=None,
                            help=textwrap.dedent('''\
                                path to a folder containing extra parameters for NSAA.'''))

    args = arg_parser.parse_args()

    protein_pdb = os.path.abspath(args.protein) if args.protein is not None else None
    advance_folder = os.path.abspath(args.advance_folder) if args.advance_folder is not None else None
    pair_list = args.pair_list
    if_dock = args.dock
    if_cs_align = args.cs_align
    if_prepare_protein = not args.not_prepare_protein
    if_strict_pair_list = (pair_list is not None and if_cs_align) and args.strict_pair_list
    if_auto_run = args.auto_run
    if_flip_by_vdw = args.flip_by_vdw
    if_allow_sep = args.allow_sep

    if protein_pdb is None or (not os.path.isfile(protein_pdb)):
        raise ValueError('Protein pdb file required.')

    if args.output_dir is not None:
        output_dir = os.path.abspath(args.output_dir)
        os.makedirs(output_dir, exist_ok=True)
        os.chdir(output_dir)
    else:
        raise ValueError(f'Output directory is required.')

    preprocess_dir = os.path.join(output_dir, 'preprocess')

    run_mode = 'ligands' if args.ligands is not None else 'folder'

    logger.info("-"*42)
    logger.info(f"run_mode: {run_mode}")
    logger.info(f"if_dock: {if_dock}")
    logger.info(f'if_cs_align: {if_cs_align}')
    logger.info(f"if_prepare_protein: {if_prepare_protein}")
    logger.info(f"if_strict_pair_list: {if_strict_pair_list}")
    logger.info("-"*42)

    if run_mode == 'folder':
        if advance_folder is None:
            raise ValueError("No ligands provided and no advance_folder provided")
        # Folder mode
        mol_base = MolBase([], [], protein_pdb=protein_pdb, nsaa_params_dir=args.nsaa_params_dir)
        pair_queue = None
    elif run_mode == 'ligands':
        # Prepare from ligands
        ref_moles, ac_moles = Initial_ligands(args.ligands, args.ref_ligand_list)

        mol_base = MolBase(ref_moles=ref_moles, ac_moles=ac_moles, protein_pdb=protein_pdb, nsaa_params_dir=args.nsaa_params_dir)

        if pair_list is None:
            logger.info(f"Generate pair map using lomap2.")
            pair_list = GeneratePairMap(mol_base.all_moles, os.path.join(preprocess_dir, 'pairs.lst'))
        else:
            pair_list = load_pairs(pair_list)
        
        pair_list = set(pair_list)

        _pair_map = PairMap(pair_list)
        start_nodes = [i.GetProp('_Name') for i in mol_base.ref_moles]

        not_using_start_nodes = set()
        for node in start_nodes:
            show_in_pair_list = False
            for pair in pair_list:
                if pair.has_node(node) != -1:
                    show_in_pair_list = True
                    break
            if not show_in_pair_list:
                not_using_start_nodes.add(node)

        if len(not_using_start_nodes) > 0:
            logger.warning(f"Warning: {not_using_start_nodes} are not in the pair list. This maybe caused by lomap2.")
            start_nodes = [node for node in start_nodes if node not in not_using_start_nodes]
        # Sometimes the pair list generated by lomap2 is not covering all the nodes, which means some nodes are abandoned in calculation.
        # We will delete these nodes from the start nodes.

        if if_strict_pair_list:
            avail_gen, avail_queue = PairMap.can_generate_from_nodes(start_nodes, pair_list, strict_direction=True)
            if avail_gen:
                logger.info('Strictly follow the preparation order in the given list.')
                pair_queue = avail_queue
            else:
                logger.info('Automatically arrange the preparation order in the given list.')
                pair_queue, _ = _pair_map.BFS_from_nodes(start_nodes)
        else:
            if if_flip_by_vdw:
                vdw_volume_cache = {}
                for pair in pair_list:
                    if pair.node1 not in vdw_volume_cache:
                        vdw_volume_cache[pair.node1] = get_mol_vdw_volume(_mol=mol_base.getMol(pair.node1))
                    if pair.node2 not in vdw_volume_cache:
                        vdw_volume_cache[pair.node2] = get_mol_vdw_volume(_mol=mol_base.getMol(pair.node2))

                    _mol1_v = vdw_volume_cache[pair.node1]
                    _mol2_v = vdw_volume_cache[pair.node2]

                    if _mol1_v < _mol2_v:
                        pair.flip()

                logger.debug(f"VDW volume: {vdw_volume_cache}")
            
            pair_queue = generate_legal_visit_queue(pair_list, set(start_nodes), allow_sep=if_allow_sep)

    else:
        raise ValueError("Shouldn't be able to show this.")

    r = RawFileManager(mol_base=mol_base, output_dir=output_dir, pair_queue=pair_queue)

    r.prepare_protein_pdb(if_prepare=if_prepare_protein)
    # We need to update the protein info in the FileManager, so we do the 'if_prepare' check in the function.
    try:
        r.collision_checker = CollisionChecker.from_pdb_file(protein_pdb)
    except:
        logger.warning("Failed to setup collision checker from pdb file, ignore this if you don't need it.")

    if run_mode == 'folder':
        r.load_prepare_from_custom_folder(advance_folder)
    elif run_mode == 'ligands':
        r.prepare_pair_files(if_cs_align=if_cs_align, if_dock=if_dock,
                             if_auto_run=if_auto_run)

    run_list = r.initial_prepare_folder()
    r.dump_infos_to_pickle()

    if if_auto_run:
        r.auto_run_prepare(run_list)


if __name__ == '__main__':
    CLI_main()
    pass
