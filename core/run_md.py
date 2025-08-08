import os
import copy
import shutil
import logging
import sys
import warnings
import re
import argparse
from dataclasses import dataclass, field
from typing import List, Tuple, Union, Dict, Optional

import parmed as pmd

from top_mod_tools.top_atomtype_fix import top_atomtype_fix
from top_mod_tools.top_nonbonded_exclusion_fix import check_and_fix_exclusions_from_file
from analyze_tools.SettingManager import SectionSettings, AllSettings, get_default_path

# Constants
TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), 'templates')
DEFAULT_EDGES = ['lM2A', 'lM2B', 'cM2A', 'cM2B']

@dataclass
class SimulationConfig:
    """Configuration for MD simulation."""
    edges: List[str] = field(default_factory=lambda: DEFAULT_EDGES.copy())
    work_section_numbers: int = 1
    local_run: bool = False

def setup_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """Setup logger for run_md module."""
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    if logger.handlers:
        logger.handlers.clear()
    
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    
    # Detect WSL environment
    try:
        with open('/proc/version', 'r') as f:
            is_wsl = 'microsoft' in f.read().lower()
    except (OSError, IOError):
        is_wsl = False
    
    if is_wsl:
        formatter = logging.Formatter('%(levelname)s: %(message)s')
    else:
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger

logger = setup_logger(__name__)


AUTOFILL_TOKEN = '<AUTOFILL>'


class _DefaultProperties:
    """Default properties for MD simulations."""
    props = {
        # Json Part
        'lambda_nums': 50,
        'vdw_json_dict': {'lambda_restraints': (0, 0),
                          'lambda_electrostatics': (0, 0),
                          'lambda_sterics': (1, 0)},
        'charge_json_dict': {'lambda_restraints': (0, 0),
                             'lambda_electrostatics': (1, 0),
                             'lambda_sterics': (1, 1)},

        'normal_json_dict': {'lambda_restraints': (0, 0),
                             'lambda_electrostatics': (1, 0),
                             'lambda_sterics': (1, 0)},

        # amber json part. 0 to 1.
        'amber_vdw_json_dict': {'lambda_restraints': (0, 0),
                                'lambda_electrostatics': (0, 1),
                                'lambda_sterics': (0, 1)},
        'amber_charge_json_dict': {'lambda_restraints': (0, 0),
                                   'lambda_electrostatics': (0, 1),
                                   'lambda_sterics': (0, 1)},
        'amber_normal_json_dict': {'lambda_restraints': (0, 0),
                                   'lambda_electrostatics': (0, 1),
                                   'lambda_sterics': (0, 1)},

        # openmm json part. 1 to 0.
        'openmm_vdw_json_dict': {'lambda_restraints': (0, 0),
                                 'lambda_electrostatics': (0, 0),
                                 'lambda_sterics': (1, 0)},
        'openmm_charge_json_dict': {'lambda_restraints': (0, 0),
                                    'lambda_electrostatics': (1, 0),
                                    'lambda_sterics': (1, 1)},
        'json_scheme_name': 'default',
        # Bash Part
        'bash_filename': 'submit.sh',
        'python_script': '',
        'python_args': {},
        # Input Part
        'input_filename': 'input.txt',
        'segment_input_filename': 'segment_input.txt',
        'final_state_filename': 'alc_final_state.xml',
        # Static File Part
        'default_static_ref_dict': {
            'top_file': {
                'lM2A': 'FEP_openmm/ligands/lM2A/lM2A.prmtop',
                'lM2B': 'FEP_openmm/ligands/lM2B/lM2B.prmtop',
                'cM2A': 'FEP_openmm/complex/cM2A/cM2A.prmtop',
                'cM2B': 'FEP_openmm/complex/cM2B/cM2B.prmtop'
            },
            'crd_file': {
                'lM2A': 'FEP_openmm/ligands/lM2A/lM2A.prmcrd',
                'lM2B': 'FEP_openmm/ligands/lM2B/lM2B.prmcrd',
                'cM2A': 'FEP_openmm/complex/cM2A/cM2A.prmcrd',
                'cM2B': 'FEP_openmm/complex/cM2B/cM2B.prmcrd'
            },
            'pdbx_file': {
                'lM2A': 'FEP_openmm/ligands/lM2A/ligAM.pdbx',
                'lM2B': 'FEP_openmm/ligands/lM2B/ligBM.pdbx',
                'cM2A': 'FEP_openmm/complex/cM2A/ligAM.pdbx',
                'cM2B': 'FEP_openmm/complex/cM2B/ligBM.pdbx'
            }
        },
        'amber_static_ref_dict': {
            'top_file': {
                'lM2A': 'FEP_amber/ligands/lM2A/lM2A.prmtop',
                'lM2B': 'FEP_amber/ligands/lM2B/lM2B.prmtop',
                'cM2A': 'FEP_amber/complex/cM2A/cM2A.prmtop',
                'cM2B': 'FEP_amber/complex/cM2B/cM2B.prmtop'
            },
            'crd_file': {
                'lM2A': 'FEP_amber/ligands/lM2A/lM2A.prmcrd',
                'lM2B': 'FEP_amber/ligands/lM2B/lM2B.prmcrd',
                'cM2A': 'FEP_amber/complex/cM2A/cM2A.prmcrd',
                'cM2B': 'FEP_amber/complex/cM2B/cM2B.prmcrd'
            },
            'pdbx_file': {
                # pdbx file is unnecessary in amber.
                'lM2A': 'FEP_openmm/ligands/lM2A/ligAM.pdbx',
                'lM2B': 'FEP_openmm/ligands/lM2B/ligBM.pdbx',
                'cM2A': 'FEP_openmm/complex/cM2A/ligAM.pdbx',
                'cM2B': 'FEP_openmm/complex/cM2B/ligBM.pdbx'
            }
        },
        'openmm_static_ref_dict': {
            'top_file': {
                'lM2A': 'FEP_openmm/ligands/lM2A/lM2A.prmtop',
                'lM2B': 'FEP_openmm/ligands/lM2B/lM2B.prmtop',
                'cM2A': 'FEP_openmm/complex/cM2A/cM2A.prmtop',
                'cM2B': 'FEP_openmm/complex/cM2B/cM2B.prmtop'
            },
            'crd_file': {
                'lM2A': 'FEP_openmm/ligands/lM2A/lM2A.prmcrd',
                'lM2B': 'FEP_openmm/ligands/lM2B/lM2B.prmcrd',
                'cM2A': 'FEP_openmm/complex/cM2A/cM2A.prmcrd',
                'cM2B': 'FEP_openmm/complex/cM2B/cM2B.prmcrd'
            },
            'pdbx_file': {
                'lM2A': 'FEP_openmm/ligands/lM2A/ligAM.pdbx',
                'lM2B': 'FEP_openmm/ligands/lM2B/ligBM.pdbx',
                'cM2A': 'FEP_openmm/complex/cM2A/ligAM.pdbx',
                'cM2B': 'FEP_openmm/complex/cM2B/ligBM.pdbx'
            }
        },
        # Submit Job Part
        'to_run_queue_lst': r''
    }

    def __init__(self):
        self.__dict__.update(_DefaultProperties.props)


_ActiveProps = _DefaultProperties()


def _pdbx_charge_mod(ori_pdbx_file, new_pdbx_file) -> str:
    """
    change the charge in pdbx file to keep the same charge between start and end.
    :param ori_pdbx_file: Abs path of origin pdbx file
    :param new_pdbx_file: Abs path of new pdbx file
    :return: the name of new pdbx file
    """
    new_pdbx_file_name = os.path.split(new_pdbx_file)[1]
    new_content = ''
    with open(ori_pdbx_file, 'r') as f:
        for line in f:
            # bg_charge = line[58:69]
            ed_charge = line[73:84]
            new_line = line[:58] + ed_charge + '    ' + ed_charge + line[84:]
            new_content += new_line
    with open(new_pdbx_file, 'w') as f:
        f.write(new_content)
    return new_pdbx_file_name


def _pdbx_mark(ori_pdbx_file) -> Union[str, None]:
    def _extract_pdbx_info(aline):
        tmp_line = copy.deepcopy(aline)
        while tmp_line.find('  ') != -1:
            tmp_line = tmp_line.replace("  ", ' ')
        tmp_line = tmp_line.strip()
        r_info = tmp_line.split(' ')
        return r_info

    if not (os.path.exists(ori_pdbx_file) and ori_pdbx_file.endswith('.pdbx')):
        return None

    new_pdbx_content = ''
    with open(ori_pdbx_file, 'r') as f:
        for line in f:
            info = _extract_pdbx_info(line)
            atom = info[2]
            res_name = info[3]
            atom_group = int(info[10])
            bf_chg = info[8]
            af_chg = info[9]

            if bf_chg != af_chg and bf_chg == '0.00000000':
                new_pdbx_content += line[:-2]+'2\n'
            elif atom_group == 1:
                new_pdbx_content += line
            else:
                new_pdbx_content += line[:17]+"COM"+line[20:-2]+'3\n'

    new_pdbx_file = ori_pdbx_file.replace('.pdbx', '_n.pdbx')
    with open(new_pdbx_file, 'w') as f:
        f.write(new_pdbx_content)

    return new_pdbx_file


def _pdbx_id_match(pdbx_file, top_file, output_pdbx_file: str = None) -> Union[str, None]:
    if pdbx_file is None:
        return None

    pdbx_file = os.path.abspath(pdbx_file)

    if not (os.path.exists(pdbx_file) and pdbx_file.endswith('.pdbx')):
        return None

    if not (os.path.exists(top_file) and top_file.endswith('.prmtop')):
        return None

    atom_idx_dict = {}
    register_atom_name = []
    top = pmd.amber.AmberParm(top_file)

    if output_pdbx_file is None:
        output_pdbx_file = pdbx_file

    for residue in top.residues:
        if residue.name == 'LAM' or residue.name == 'LBM':
            for atom in residue.atoms:
                if atom.name not in register_atom_name:
                    atom_idx_dict[f'{residue.name}_{atom.name.upper()}'] = atom.idx
                    register_atom_name.append(f'{residue.name}_{atom.name.upper()}')
                else:
                    logger.debug(f'{atom.name} is already in list')
                    warnings.warn(f'Duplicate atom name in Ligand {residue.name}')

    new_pdb_content = ''
    ori_resname_map = {'LAM': 'LAM', 'LBM': 'LBM', 'COM': None}
    with open(pdbx_file, 'r') as f:
        id_map = {'LAM': [], 'LBM': [], 'COM': []}
        for line in f:
            tmp_line = line.strip()
            while tmp_line.find('  ') != -1:
                tmp_line = tmp_line.replace('  ', ' ')
            info = tmp_line.split(' ')
            residue_name = info[3]
            group_id = int(info[-1])
            id_map[residue_name].append(group_id)
        if id_map['COM'][0]-1 in id_map['LAM']:
            ori_resname_map['COM'] = 'LAM'
        else:
            ori_resname_map['COM'] = 'LBM'

    with open(pdbx_file, 'r') as f:
        for line in f:
            tmp_line = line.strip()
            while tmp_line.find('  ') != -1:
                tmp_line = tmp_line.replace('  ', ' ')
            info = tmp_line.split(' ')
            atom_name = info[2]
            residue_name = info[3]
            # group_id = info[-1]
            new_id = atom_idx_dict[f'{ori_resname_map[residue_name]}_{atom_name}']
            new_id = f'{new_id:>5}'
            new_line = line[:6] + new_id + line[11:]
            new_pdb_content += new_line

    with open(output_pdbx_file, 'w') as f:
        f.write(new_pdb_content)

    return output_pdbx_file


def _opposite_edge(this_edge):
    if this_edge == 'lM2A':
        return 'lM2B'
    if this_edge == 'lM2B':
        return 'lM2A'
    if this_edge == 'cM2A':
        return 'cM2B'
    if this_edge == 'cM2B':
        return 'cM2A'


def _fix_13_16_dislocation(pdbx):
    """
    To fix 13-16 dislocation after genambRBFE's pdbx
    Parameters
    ----------
    pdbx

    Returns
    -------

    """
    new_pdbx_content = ""

    def _get_char(_str):
        return ''.join(c for c in _str if c.isalpha())

    def _get_num(_str):
        return ''.join(c for c in _str if c.isnumeric())

    with open(pdbx, 'r') as f:
        for line in f.readlines():
            tmp_line = line
            while tmp_line.find("  ") != -1:
                tmp_line = tmp_line.replace("  ", ' ')
            atom_text = line[12:17].strip()
            atom_name = _get_char(atom_text)
            atom_id = _get_num(atom_text)
            line = line[:12] + f'{atom_name.upper():>2}{atom_id:<2}' + " " + line[17:]

            new_pdbx_content += line

    with open(pdbx, 'w') as f:
        f.write(new_pdbx_content)

    return os.path.abspath(pdbx)


class _pair_static_file_path:
    default_ref_dict = {
        'top_file': {
            'lM2A': 'FEP_openmm/ligands/lM2A/lM2A.prmtop',
            'lM2B': 'FEP_openmm/ligands/lM2B/lM2B.prmtop',
            'cM2A': 'FEP_openmm/complex/cM2A/cM2A.prmtop',
            'cM2B': 'FEP_openmm/complex/cM2B/cM2B.prmtop'
        },
        'crd_file': {
            'lM2A': 'FEP_openmm/ligands/lM2A/lM2A.prmcrd',
            'lM2B': 'FEP_openmm/ligands/lM2B/lM2B.prmcrd',
            'cM2A': 'FEP_openmm/complex/cM2A/cM2A.prmcrd',
            'cM2B': 'FEP_openmm/complex/cM2B/cM2B.prmcrd'
        },
        'pdbx_file': {
            'lM2A': 'FEP_openmm/ligands/lM2A/ligAM.pdbx',
            'lM2B': 'FEP_openmm/ligands/lM2B/ligBM.pdbx',
            'cM2A': 'FEP_openmm/complex/cM2A/ligAM.pdbx',
            'cM2B': 'FEP_openmm/complex/cM2B/ligBM.pdbx'
        }
    }

    def __init__(self, pair_dir, ref_dict: dict):
        self.pair_dir = os.path.abspath(pair_dir)
        all_edge_abspath_dict = {}
        for file_type, file_path_dict in ref_dict.items():
            tmp_abs_dict = {}
            for edge_name, file_path in file_path_dict.items():
                tmp_abs_dict[edge_name] = os.path.join(self.pair_dir, file_path)
            all_edge_abspath_dict[file_type] = tmp_abs_dict
        self.__dict__.update(all_edge_abspath_dict)


class BashGenerator:
    def __init__(self, run_python_file: str = None, run_args: dict = None, log_file: str = 'run.log'):
        """
        :param log_file: the name of the log file. File name, Ref path or Abs path are both ok.
        """
        self.log_file = log_file
        self.log_file_name = os.path.splitext(self.log_file)[1]
        if local_run:
            self.gpu_setup_command = '#!/bin/bash\n'
        else:
            self.gpu_setup_command = '#!/bin/bash\nexport CUDA_DEVICE_ORDER=PCI_BUS_ID\nexport CUDA_VISIBLE_DEVICES=$1\n'
        self.begin_timestamp_command = f'printf "BEGIN AT:%s \\n" "$(date -R)" > {log_file}\n'
        self.end_timestamp_command = f'printf "END AT:%s \\n" "$(date -R)" >> {log_file}\n'
        
        self.before_run_command = ''
        self.after_run_command = ''
        
        self.run_python_file = run_python_file
        self.run_args = {}
        if run_args is not None:
            self.run_args.update(run_args)

    def replace_args_placeholder(self, placeholder, replace):
        for key, value in self.run_args.items():
            if value.find(placeholder) != -1:
                self.run_args[key] = self.run_args[key].replace(placeholder, replace)
        
    def gen_file(self, output_file: str, python_file: str = None, args: dict = None):
        """
        Generate a bash
        :param output_file: The output file name of the bash.
        :param python_file: which python script to run.
        :param args: arguments to pass to the python script, save as a dict. like {'-i': 'input.txt'}
        :return: None
        """
        if python_file is None:
            python_file = self.run_python_file
        if args is None:
            args = self.run_args
        
        if python_file is None:
            raise Exception(f'No python file is specified in this bash generator')
        run_command = f'python {python_file} '
        
        if len(args) > 0:
            for key, value in args.items():
                run_command += f'{key} {value} '
                
        run_command += f'>> {self.log_file} 2>&1\n'
        
        with open(output_file, 'w') as f:
            f.write(self.gpu_setup_command)
            f.write(self.before_run_command)
            f.write(self.begin_timestamp_command)
            f.write(run_command)
            f.write(self.end_timestamp_command)
            f.write(self.after_run_command)

    def valid(self) -> bool:
        if self.run_python_file is None:
            raise Exception(f'No python file is specified in this bash')
        return True


class JsonGenerator:
    example_lambda_range_dict = {'lambda_restraints': (0, 0),
                                 'lambda_electrostatics': (1, 0),
                                 'lambda_sterics': (1, 0),
                                 'lambda_electrostatics_env': (1, 1),
                                 'lambda_electrostatics_chgprod_square': (1, 1)}

    def __init__(self, lambda_range_dict: dict, window_num: int, scheme_name: str) -> None:
        self.lambdas = {}
        for lambda_name, lambda_range in lambda_range_dict.items():
            tmp_list = []
            bg_lambda, el_lambda = lambda_range
            if bg_lambda == el_lambda:
                tmp_list = [float(f'{bg_lambda:.4f}') for _ in range(window_num)]
            else:
                if bg_lambda < el_lambda:
                    step = 1
                else:
                    step = -1
                for i in range(bg_lambda*window_num, el_lambda*window_num, step):
                    tmp_list.append(float(f'{i / window_num:.4f}'))
            tmp_list.append(float(f'{el_lambda:.4f}'))
            self.lambdas[lambda_name] = tmp_list
        self.lambdas_num = len(self.lambdas)
        self.scheme_name = scheme_name

    @classmethod
    def gen_json_from_dict(cls, lambda_range_dict: dict, window_num: int, scheme_name: str, output_file: str) -> str:
        tmp_jg = cls(lambda_range_dict=lambda_range_dict, window_num=window_num, scheme_name=scheme_name)
        tmp_jg.gen_json(output_file=output_file, lambdas=tmp_jg.lambdas, scheme_name=tmp_jg.scheme_name)
        return output_file

    def interpolate_env_lambda(self):
        new_lambdas = {}
        for lambda_name, lambda_list in self.lambdas.items():
            tmp_list = []
            for lambda_value in lambda_list:
                tmp_list.append(lambda_value)
                tmp_list.append(lambda_value)
            new_lambdas[lambda_name] = tmp_list
        lambda_electrostatics_env_list = []
        lambda_electrostatics_chgprod_square_list = []
        for i in range(len(new_lambdas['lambda_electrostatics'])):
            if i % 2 == 0:
                lambda_electrostatics_env_list.append(1)
            else:
                lambda_electrostatics_env_list.append(0)
            lambda_electrostatics_chgprod_square_list.append(float(f"{new_lambdas['lambda_electrostatics'][i] ** 2:.4f}"))
        new_lambdas['lambda_electrostatics_env'] = lambda_electrostatics_env_list
        new_lambdas['lambda_electrostatics_chgprod_square'] = lambda_electrostatics_chgprod_square_list
        return new_lambdas

    @staticmethod
    def gen_json(output_file: str, lambdas: dict = None, scheme_name: str = 'default'):
        active_lambda = lambdas
        active_lambda_num = len(active_lambda)
        with open(output_file, "w") as fo:
            fo.write("{\n\t\"" + scheme_name + "\":\n\t{\n")
            name_max_len = -1
            for name in lambdas.keys():
                if name_max_len < len(name):
                    name_max_len = len(name)
                # Get max length in lambdas' name.

            for lambda_counter, (lambda_name, lambda_list) in enumerate(active_lambda.items()):
                align_space = ' ' * (name_max_len - len(lambda_name))
                # out_str = '\t\t"' + lambda_name + '"' + align_space + ':'+str(lambda_list)
                out_str = f'\t\t"{lambda_name}"{align_space}:{lambda_list}'
                if lambda_counter < active_lambda_num-1:
                    out_str += ',\n'
                else:
                    out_str += '\n'
                fo.write(out_str)
            fo.write('\t}\n}')


class SingleProcessFileManager:
    def __init__(self, top_file: str, crd_file: str, pdbx_file: str, json_file: str, new_dir: str,
                 additional_stc_files: dict = None):
        """
        Using Abspath is better.
        :param top_file: path to top file
        :param crd_file: path to crd file
        :param pdbx_file: path to pdbx file
        :param new_dir: where to save the files, also the work directory.

        :param additional_stc_files: a dict containing additional static files. should be: {filename: path, }
        """

        self.__new_dir = new_dir
        
        self.__top_file = top_file
        self.__crd_file = crd_file
        self.__pdbx_file = pdbx_file
        self.__json_file = json_file

        self.__additional_stc_files_dict = additional_stc_files

        self.__json_name = os.path.basename(self.__json_file)
        self.__new_json_file = os.path.join(self.__new_dir, self.__json_name)

        self.__input_file_name = _ActiveProps.input_filename
        self.__input_file = os.path.join(self.__new_dir, self.__input_file_name)
        # 5 Necessary file to run

        self.settings = None  # Type: AllSettings or SectionSetting, use to generate self.__input_file. (input.txt)
        self.__submit_bash_name = 'submit.sh'
        self.__submit_bash_file = os.path.join(self.__new_dir, self.__submit_bash_name)

        self.submit_bash_generator = None  # Type: BashGenerator. Has its own config to gen diff bash.

    @property
    def input_file(self):
        return self.__input_file

    @input_file.setter
    def input_file(self, name):
        self.__input_file_name = name
        self.__input_file = os.path.join(self.__new_dir, self.__input_file_name)

    @property
    def submit_bash(self):
        return self.__submit_bash_file

    @submit_bash.setter
    def submit_bash(self, name):
        self.__submit_bash_name = name
        self.__submit_bash_file = os.path.join(self.__new_dir, self.__submit_bash_name)
    
    def import_settings(self, input_setting: Union[SectionSettings, AllSettings]):
        if isinstance(input_setting, SectionSettings) or isinstance(input_setting, AllSettings):
            self.settings = input_setting
        else:
            raise TypeError('Cannot import settings from given obj. Should be class: SectionSettings or AllSettings')

    def _generate_variable_files(self):
        if self.settings is not None:
            if os.path.exists(self.__input_file):
                os.remove(self.__input_file)
            self.settings.gen_file(self.__input_file)
        else:
            raise ValueError('Cannot generate input file. self.settings is not imported.')
            # Generate settings. (input.txt)

        if self.submit_bash_generator is not None:
            self.submit_bash_generator.gen_file(self.__submit_bash_file)
        else:
            raise ValueError('Cannot generate bash file. self.submit_bash_generator is not set.')
            # Generate bash script (submit.sh)

    def prepare_run_files(self, queue_lst_file: str):
        self._generate_variable_files()
        # Generate files. (input.txt, submit.sh, submit_json.json)

        self._copy_static_files()
        # Copy files. (.prmtop, .prmcrd, .pdbx)

        with open(queue_lst_file, 'a') as f:
            f.write(f'{self.__new_dir} {self.__submit_bash_file}\n')
            # Add submit command to the file storages queue list.

    def _copy_static_files(self):
        if not os.path.exists(self.__new_dir):
            os.mkdir(self.__new_dir)

        def new_file_path(x):
            return os.path.join(self.__new_dir, os.path.basename(x))

        shutil.copy(self.__top_file, new_file_path(self.__top_file))
        shutil.copy(self.__crd_file, new_file_path(self.__crd_file))
        shutil.copy(self.__pdbx_file, new_file_path(self.__pdbx_file))

        shutil.copy(self.__json_file, new_file_path(self.__json_file))

        if self.__additional_stc_files_dict is not None:
            if len(self.__additional_stc_files_dict) > 0:
                for filename, filepath in self.__additional_stc_files_dict.items():
                    shutil.copy(filepath, new_file_path(filename))


class EdgeFileManager(object):
    def __init__(self, settings: AllSettings, edge_name: str, edge_dir: str,
                 top_file: str, crd_file: str, group_pdbx_file: str,
                 lambda_nums: int, **kwargs):
        """
        File manager on an edge.
        :param settings: A setting to same edge. Only diff on current_group_nb and input_state. Type: AllSettings
        :param edge_dir: the FEP dir of Edge.
        :param top_file: Edge's top file. lM2A.prmtop
        :param crd_file: Edge's crd file. lM2A.prmcrd
        :param pdbx_file: Edge's pdbx file. lM2A.pdbx
        :param lambda_nums: json window size
        """
        self.settings = settings
        self._edge_dir = edge_dir
        self._top_file = top_file
        self._crd_file = crd_file
        self._group_pdbx_file = group_pdbx_file  # The pdbx ONLY contains group 1, 2, 3.
        # Static files and setting (input.txt)

        self._lambda_nums = lambda_nums
        self.overwrite_setting_dict = {
            'alchemical': {'lambdas_group': _ActiveProps.json_scheme_name,
                           'simulation_lambdas_name': _ActiveProps.json_scheme_name}
        }
        self.settings.apply_settings_with_dict(self.overwrite_setting_dict)
        # Use to generate individual settings for each single process. (input.txt)
        # Only fill lambdas_json, pdbx, current_group_nb and input_state in this layer.

        self.vdw_bash_if_check_end = kwargs['vdw_bash_if_check_end']

        self._edge_name = edge_name
        self._tmp_files = []

        if not os.path.exists(self._edge_dir):
            os.makedirs(self._edge_dir)

    def _fill_prm_files_in_bash_generator(self, bash_generator):
        bash_generator.replace_args_placeholder('<top_file>', os.path.basename(self._top_file))
        bash_generator.replace_args_placeholder('<crd_file>', os.path.basename(self._crd_file))

    def _clean_tmp_files(self):
        for file in self._tmp_files:
            os.remove(file)
        self._tmp_files.clear()

    def _prepare_single_process(self, json_file, pdbx_file, bash_generator, queue_list_file,
                                apply_setting_dict, process_dir, next_process_dir=None):
        if not os.path.exists(process_dir):
            os.makedirs(process_dir)

        tmp_singe_file_manager = SingleProcessFileManager(top_file=self._top_file,
                                                          crd_file=self._crd_file,
                                                          json_file=json_file,
                                                          pdbx_file=pdbx_file,
                                                          new_dir=process_dir)

        tmp_singe_file_manager.submit_bash_generator = copy.deepcopy(bash_generator)
        self._fill_prm_files_in_bash_generator(tmp_singe_file_manager.submit_bash_generator)
        if next_process_dir is not None:
            cp_file_shell = f'cp ./{_ActiveProps.final_state_filename} {next_process_dir}/{_ActiveProps.final_state_filename}\n'
            tmp_singe_file_manager.submit_bash_generator.after_run_command += cp_file_shell
        # Replace args things like {edge}.prmtop to lM2A.prmtop

        tmp_singe_file_manager.settings = copy.deepcopy(self.settings)
        tmp_singe_file_manager.settings.apply_settings_with_dict(apply_setting_dict)

        tmp_singe_file_manager.input_file = _ActiveProps.input_filename
        tmp_singe_file_manager.submit_bash = _ActiveProps.bash_filename

        tmp_singe_file_manager.prepare_run_files(queue_lst_file=queue_list_file)

    def prepare_files(self, bash_generator: BashGenerator, queue_list_file, **kwargs):
        _normal_json_file = JsonGenerator.gen_json_from_dict(lambda_range_dict=_ActiveProps.normal_json_dict,
                                                             window_num=self._lambda_nums,
                                                             scheme_name=_ActiveProps.json_scheme_name,
                                                             output_file=os.path.join(os.getcwd(), 'normal.json'))
        self._tmp_files.append(_normal_json_file)

        _normal_pdbx_file = self._group_pdbx_file
        normal_dir = os.path.join(self._edge_dir)
        tmp_normal_setting_dict = {'alchemical': {'lambdas_json': os.path.basename(_normal_json_file),
                                                  'pdbx': os.path.basename(_normal_pdbx_file),
                                                  'current_group_nb': 1,
                                                  'input_state': None
                                                  }
                                   }
        self._prepare_single_process(json_file=_normal_json_file,
                                     pdbx_file=_normal_pdbx_file,
                                     bash_generator=bash_generator,
                                     queue_list_file=queue_list_file,
                                     apply_setting_dict=tmp_normal_setting_dict,
                                     process_dir=normal_dir,
                                     next_process_dir=None)

        self._clean_tmp_files()


class MainFileManager:
    def __init__(self, setting_file: str, lambda_nums: int, bash_generator: BashGenerator,
                 path_to_pairs: str, pairs_list: Optional[List] = None, work_dir: Optional[str] = None,
                 edge_file_manager_class=EdgeFileManager,
                 vdw_bash_if_check_end: bool = True,
                 config: Optional[SimulationConfig] = None):
        """Initialize MainFileManager.
        
        Args:
            setting_file: Path to settings file
            lambda_nums: Number of lambda values
            bash_generator: Bash script generator
            path_to_pairs: Path to pairs directory
            pairs_list: List of pair names
            work_dir: Working directory
            edge_file_manager_class: Edge file manager class
            vdw_bash_if_check_end: Check VDW bash end
            config: Simulation configuration
        """
        self.__path_to_pairs = path_to_pairs
        self.config = config or SimulationConfig()
        
        if pairs_list is None:
            pairs_list = []
            for folder_name in os.listdir(self.__path_to_pairs):
                if os.path.isdir(os.path.join(self.__path_to_pairs, folder_name)) and folder_name.find('-') != -1:
                    pairs_list.append(folder_name)
        self.pairs_list = pairs_list

        self.bash_generator = bash_generator
        self.settings = AllSettings.from_file(setting_file)
        self.lambda_nums = lambda_nums
        self.edge_file_manager_class = edge_file_manager_class
        self.vdw_bash_if_check_end = vdw_bash_if_check_end

        # Edge File Manager Part
        self.edge_file_manager_class = edge_file_manager_class

        self.vdw_bash_if_check_end = vdw_bash_if_check_end

        # Where Work.
        if work_dir is None:
            work_dir = os.path.join(os.getcwd(), 'run')
        self.work_dir = os.path.abspath(work_dir)
        if not os.path.exists(self.work_dir):
            os.makedirs(self.work_dir)

        self.settings.apply_settings_with_dict({'alchemical': {'lambdas_per_state': 1}})

    @staticmethod
    def _optimize_queue_list_file(queue_list_file, pattern='group_[0-9]+'):
        """
        optimize the queue list file in order to run efficiently.
        :param queue_list_file:
        :return:
        """
        new_queue_list = []
        with open(queue_list_file, 'r') as f:
            for line in f.readlines():
                foo = re.search(pattern, line).group()
                group_num = int(foo.split('_')[1])
                new_queue_list.append((line, group_num))

        new_queue_list = sorted(new_queue_list, key=lambda x: x[1])

        with open(queue_list_file, 'w') as f:
            for line, group_num in new_queue_list:
                f.write(line)

    def prepare_files(self, if_amber=False):
        self.bash_generator.valid()

        # Queue list Part
        outer_work_dir = os.path.dirname(self.work_dir)
        vdw_queue_file = os.path.join(outer_work_dir, 'run_vdw_queue.lst')
        charge_queue_file = os.path.join(outer_work_dir, 'run_charge_queue.lst')
        all_queue_file = os.path.join(outer_work_dir, 'run_all_queue.lst')

        open(vdw_queue_file, 'w').close()
        open(charge_queue_file, 'w').close()
        # Make clean queue files.
        success_pair_count = 0

        err_pairs = set()

        for i, pair in enumerate(self.pairs_list):
            logger.info(f'Preparing MD files... {pair} ({i+1}/{len(self.pairs_list)}) Section [1/{self.config.work_section_numbers}]')
            origin_pair_path = os.path.join(self.__path_to_pairs, pair)
            pair_static_files = _pair_static_file_path(pair_dir=origin_pair_path,
                                                       ref_dict=_ActiveProps.default_static_ref_dict)

            work_pair_path = os.path.join(self.work_dir, pair)
            if not os.path.exists(work_pair_path):
                os.makedirs(work_pair_path)

            error_happened = False

            for edge in DEFAULT_EDGES:
                if error_happened:
                    continue
                    # If one edge error, abandon the whole pair.

                edge_work_path = os.path.join(work_pair_path, edge)
                if not os.path.exists(edge_work_path):
                    os.makedirs(edge_work_path)
                
                mutant_lst_file = os.path.join(origin_pair_path, 'mutant.lst')

                if not os.path.exists(pair_static_files.top_file[edge]):
                    error_happened = True
                    logger.warning(f'Skipping {pair} {edge}: no top file found')
                    break

                if_nonbonded_exclusion_fix = check_and_fix_exclusions_from_file(top_file=pair_static_files.top_file[edge],
                                                                                mutant_list_file=mutant_lst_file,
                                                                                output_file=pair_static_files.top_file[edge].replace(".prmtop", "_fixed.prmtop"))
                if if_nonbonded_exclusion_fix:
                    logger.info(f"Fixed nonbonded exclusions in {pair}-{edge}")
                    shutil.copy2(pair_static_files.top_file[edge], f"{pair_static_files.top_file[edge]}.bck")
                    shutil.copy2(pair_static_files.top_file[edge].replace(".prmtop", "_fixed.prmtop"), pair_static_files.top_file[edge])

                if not if_amber:
                    fixed_pdbx = _fix_13_16_dislocation(pair_static_files.pdbx_file[edge])
                    marked_pdbx = _pdbx_mark(fixed_pdbx)

                    shutil.copy2(pair_static_files.top_file[edge], f"{pair_static_files.top_file[edge]}.bck")
                    top_atomtype_fix(pair_static_files.top_file[edge])
                else:
                    marked_pdbx = pair_static_files.pdbx_file[edge]
                
                if marked_pdbx is None:
                    error_happened = True
                    logger.warning(f'Skipping {pair} {edge}: no pdbx file found')
                    break
                # Generate marked pdbx to define groups

                tmp_EdgeFileManager = self.edge_file_manager_class(settings=self.settings,
                                                                   edge_name=edge,
                                                                   edge_dir=edge_work_path,
                                                                   top_file=pair_static_files.top_file[edge],
                                                                   crd_file=pair_static_files.crd_file[edge],
                                                                   group_pdbx_file=marked_pdbx,
                                                                   lambda_nums=self.lambda_nums,
                                                                   vdw_bash_if_check_end=self.vdw_bash_if_check_end)
                try:
                    tmp_EdgeFileManager.prepare_files(bash_generator=self.bash_generator,
                                                      queue_list_file=charge_queue_file,
                                                      vdw_queue_list_file=vdw_queue_file)
                except AttributeError:
                    error_happened = True
                    logger.error(f"{edge}'s pdbx has error")
                except (FileNotFoundError, ImportError) as e:
                    error_happened = True
                    logger.error(f"{edge}'s preparation has error")

                if not if_amber:
                    os.remove(marked_pdbx)

            if error_happened:
                logger.error('Preparation error occurred')
                err_pairs.add(pair)
            else:
                success_pair_count += 1

            if error_happened:
                shutil.rmtree(work_pair_path)

        logger.info(f'Preparation finished. Success: {success_pair_count} / {len(self.pairs_list)}')
        if success_pair_count != len(self.pairs_list):
            logger.warning(f"Error pairs: {err_pairs}")

        self._optimize_queue_list_file(vdw_queue_file)

        with open(charge_queue_file, 'r') as charge_queue:
            with open(vdw_queue_file, 'r') as vdw_queue:
                with open(all_queue_file, 'w') as all_queue:
                    all_queue.write(charge_queue.read() + vdw_queue.read())
                    # Merge two queue.

        os.remove(charge_queue_file)
        os.remove(vdw_queue_file)


class MainFileManager_CAR(MainFileManager):

    default_amber_input_file = os.path.join(TEMPLATE_DIR, 'amber_input.json')
    only_charge_input_file = os.path.join(TEMPLATE_DIR, 'amber_charge_input.json')
    only_vdw_input_file = os.path.join(TEMPLATE_DIR, 'amber_vdw_input.json')

    def __init__(self, setting_file: str, lambda_nums: int, bash_generator: BashGenerator,
                 segment_run_settings: AllSettings,
                 path_to_pairs: str, pairs_list: Optional[List] = None, work_dir: Optional[str] = None,
                 if_AMBER: bool = False, config: Optional[SimulationConfig] = None, **kwargs):
        super().__init__(setting_file, lambda_nums, bash_generator, path_to_pairs, pairs_list, work_dir, config=config, **kwargs)
        self.if_AMBER = if_AMBER

        self.segment_run_settings = segment_run_settings
        _default_CAR_run_setting_dict = {'normal_alc_md': {'prod_md_time': 10,
                                                           'input_file': _ActiveProps.input_filename}}
        self.segment_run_settings.apply_settings_with_dict(_default_CAR_run_setting_dict)

        _openmm_run_dict = {'normal_alc_md': {'simulation_software': 'openmm'}}
        self.segment_run_settings.apply_settings_with_dict(_openmm_run_dict)

        if self.if_AMBER:
            _ActiveProps.input_filename = 'input.json'
            tmp_amber_setting_dict = {'normal_alc_md': {'simulation_software': 'amber',
                                                        'input_file': 'input.json'}}
            self.segment_run_settings.apply_settings_with_dict(tmp_amber_setting_dict)

            _ActiveProps.default_static_ref_dict = _ActiveProps.amber_static_ref_dict
            _ActiveProps.vdw_json_dict = _ActiveProps.amber_vdw_json_dict
            _ActiveProps.charge_json_dict = _ActiveProps.amber_charge_json_dict

        lM2A_setting_dict = {'normal_alc_md': {'coordinate_file': 'lM2A.prmcrd', 'topology_file': 'lM2A.prmtop'}}
        lM2B_setting_dict = {'normal_alc_md': {'coordinate_file': 'lM2B.prmcrd', 'topology_file': 'lM2B.prmtop'}}
        cM2A_setting_dict = {'normal_alc_md': {'coordinate_file': 'cM2A.prmcrd', 'topology_file': 'cM2A.prmtop'}}
        cM2B_setting_dict = {'normal_alc_md': {'coordinate_file': 'cM2B.prmcrd', 'topology_file': 'cM2B.prmtop'}}

        edge_setting_apply_dict = {'lM2A': lM2A_setting_dict, 'lM2B': lM2B_setting_dict,
                                   'cM2A': cM2A_setting_dict, 'cM2B': cM2B_setting_dict}

        self.edge_setting_dict = {}
        for edge, apply_dict in edge_setting_apply_dict.items():
            tmp_setting = copy.deepcopy(self.segment_run_settings)
            tmp_setting.apply_settings_with_dict(apply_dict)
            self.edge_setting_dict[edge] = tmp_setting

    def prepare_files(self):

        super().prepare_files(if_amber=self.if_AMBER)

        for pair in self.pairs_list:
            work_pair_path = os.path.join(self.work_dir, pair)
            
            if not os.path.exists(work_pair_path):
                continue

            for edge in DEFAULT_EDGES:
                edge_work_path = os.path.join(work_pair_path, edge)

                normal_path = edge_work_path
                normal_input_file = os.path.join(normal_path, _ActiveProps.segment_input_filename)
                open(normal_input_file, 'w').close()
                normal_apply_dict = {'normal_alc_md': {'mbar_lambda_dict_file': 'normal.json'}}
                self.edge_setting_dict[edge].gen_file(normal_input_file, normal_apply_dict)
                if self.if_AMBER:
                    shutil.copy(MainFileManager_CAR.default_amber_input_file,
                                os.path.join(normal_path, _ActiveProps.input_filename))


def auto_run_MD(to_run_queue_lst, all_queue_lst, _config: SimulationConfig):
    import time
    import shutil
    import subprocess
    shutil.copy2(all_queue_lst, to_run_queue_lst)

    with open(to_run_queue_lst, 'r') as f:
        all_jobs = len(f.readlines())

    while True:
        with open(to_run_queue_lst, 'r') as f:
            line = f.readline()
            line = line.strip()
            if len(line) == 0:
                break
            path = line.split(" ")[0]
            submit_bash = line.split(" ")[1]
            remain_jobs = f.readlines()

        with open(to_run_queue_lst, 'w') as f:
            for line in remain_jobs:
                f.write(line)

        os.chdir(path)
        now_process = all_jobs - len(remain_jobs)
        logger.info(f'Running MD {submit_bash} ({now_process}/{all_jobs}) Section [2/{_config.work_section_numbers}]')
        p = subprocess.Popen(['/bin/bash', submit_bash])
        p.wait()
        time.sleep(10)


if __name__ == '__main__':

    ALCHEMD_PATH = get_default_path('Alchemd_Path')
    CAR_PATH = get_default_path('CAR_Path')

    if ALCHEMD_PATH is None or CAR_PATH is None:
        warnings.warn('ALCHEMD_PATH or CAR_PATH are not defined')

    arg_parser = argparse.ArgumentParser(description='Copy files after preparation. '
                                                     'Make a work dir on cwd.',
                                         epilog='With no flag input means MD platform is Openmm and no CAR')
    arg_parser.add_argument('-p', '--preparation_dir', dest='preparation_dir', type=str,
                            help='path to preparation directory')
    arg_parser.add_argument('-ln', '--lambda_nums', type=int, default=50, dest='lambda_nums',
                            help='number PLUS 1 equal to the windows in json.\nInput with 10 means lambda gap is 0.1.\n'
                                 'Default is 50, which means lambda gap is 0.02.')
    arg_parser.add_argument('-m', '--mode', default='Amber', dest='mode',
                            help='AMBER.')
    arg_parser.add_argument('-ci', '--CAR_input', dest='car_input', default=None,
                            help='customize the car input file.')
    arg_parser.add_argument('-auto', nargs='?', dest='auto_run', const=True, default=False,
                            help='auto run MD after the preparation is finished.')


    args = arg_parser.parse_args()

    pair_dir = os.path.abspath(args.preparation_dir)
    _lambda_nums = args.lambda_nums
    _if_CAR = True
    _if_auto_run = args.auto_run
    if args.mode.lower() == 'amber' and _if_CAR:
        _if_AMBER = True
    else:
        _if_AMBER = False
    car_input = args.car_input

    # Create configuration object
    config = SimulationConfig(
        work_section_numbers=1 if not args.auto_run else 2
    )

    WORK_PATH = os.path.dirname(pair_dir)
    ORIGINAL_PATH = os.getcwd()
    os.chdir(WORK_PATH)

    active_edge_file_manager = EdgeFileManager

    if _if_AMBER:
        _ActiveProps.normal_json_dict = _ActiveProps.amber_normal_json_dict

    if _if_auto_run:
        config.work_section_numbers = 2
        local_run = True
        _ActiveProps.to_run_queue_lst = os.path.join(WORK_PATH, 'to_run_queue.lst')

    _py_script = os.path.join(ALCHEMD_PATH, 'openmm-FEP-run.py')
    _py_args = {'-i': _ActiveProps.input_filename,
                '-p': '<top_file>',
                '-c': '<crd_file>'}
    default_bash_generator = BashGenerator(run_python_file=_py_script, run_args=_py_args)

    _CAR_py_script = os.path.join(CAR_PATH, 'segmented_converge_control.py')
    _CAR_py_args = {'-i': _ActiveProps.segment_input_filename}
    segment_bash_generator = BashGenerator(run_python_file=_CAR_py_script, run_args=_CAR_py_args)

    default_input_template = os.path.join(TEMPLATE_DIR, 'input_template.txt')
    if car_input is not None and os.path.exists(car_input) and os.path.isfile(car_input):
        default_CAR_input_template = car_input
    else:
        default_CAR_input_template = os.path.join(TEMPLATE_DIR, 'segment_input_template.txt')
    shutil.copy(default_input_template, os.path.join(os.getcwd(), 'input.txt'))
    default_input_file = os.path.join(os.getcwd(), 'input.txt')

    main_manager = MainFileManager(setting_file=default_input_file,
                                   lambda_nums=_lambda_nums,
                                   bash_generator=default_bash_generator,
                                   path_to_pairs=pair_dir,
                                   edge_file_manager_class=active_edge_file_manager,
                                   config=config)

    CAR_manager = MainFileManager_CAR(setting_file=default_input_file,
                                      lambda_nums=_lambda_nums,
                                      bash_generator=segment_bash_generator,
                                      path_to_pairs=pair_dir,
                                      segment_run_settings=AllSettings.from_file(default_CAR_input_template),
                                      if_AMBER=_if_AMBER,
                                      vdw_bash_if_check_end=False,
                                      config=config)
    CAR_manager.edge_file_manager_class = active_edge_file_manager

    if _if_CAR:
        CAR_manager.prepare_files()
    else:
        main_manager.prepare_files()
    os.remove(default_input_file)

    if _if_auto_run:
        auto_run_MD(_ActiveProps.to_run_queue_lst, all_queue_lst=os.path.join(WORK_PATH, 'run_all_queue.lst'),
                    _config=config)

    os.chdir(ORIGINAL_PATH)
