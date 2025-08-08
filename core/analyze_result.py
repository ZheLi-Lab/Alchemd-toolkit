import os
import warnings
import shutil
import glob
import copy
import subprocess
import multiprocessing
import logging
import sys
import re
import argparse
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional, Union
from pathlib import Path
from multiprocessing import Pool
from concurrent import futures
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

from analyze_tools.SettingManager import AllSettings, get_default_path

# Setup logger
def setup_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """Setup logger for analyze_result module."""
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

# Configuration constants

# Adjust this value to balance speed and reliability on your device.
# Excessive parallelism may cause undefined errors, while too low a value will slow performance.
# This controls the maximum number of worker processes for both file preparation and calculation tasks.
MAX_WORKER_COUNT = min(multiprocessing.cpu_count(), 4)
MAX_WORKER_CAL = min(MAX_WORKER_COUNT, 2)

MAX_RETRIES = 3
LOG_CHECK_LINES = 2

# Pyarmor compatibility mode - use os.system instead of subprocess for better compatibility
USE_SYSTEM_CALL = True

script_dir = os.path.dirname(os.path.realpath(__file__))

ALCHEMCONVTOOLS_PATH = get_default_path('AlchemdConvTools_Path')
calculate_py_script = os.path.join(ALCHEMCONVTOOLS_PATH, 'one_end_fe_aly.py')
default_input_template = os.path.join(script_dir, 'templates', 'calculate_input_template.txt')

ene_csv_file_prefix = "\'state_g<n>_s\'"
calculate_input_filename = 'calculate_input.txt'
calculate_submit_bash_filename = 'calculate_submit.sh'


def _detect_lambda_indexes(csv_file: Union[str, Path]) -> Tuple[List[int], List[str]]:
    """Detect lambda parameter indexes from CSV file.
    
    Args:
        csv_file: Path to CSV file
        
    Returns:
        Tuple of index list and parameter names
        
    Raises:
        FileNotFoundError: If CSV file doesn't exist
        ValueError: If no lambda parameters found
    """
    try:
        with open(csv_file, 'r') as f:
            first_line = f.readline().strip()
            elements = first_line.split('|')
            
            # Find boundary where lambda parameters end
            lambda_count = 0
            for i, element in enumerate(elements):
                if re.search(r'(.+,.+)', element):
                    lambda_count = i
                    break
            
            if lambda_count == 0:
                raise ValueError(f"No lambda parameters found in {csv_file}")
                
            parameter_names = elements[:lambda_count]
            index_list = list(range(lambda_count))
            
            logger.debug(f"Detected {lambda_count} lambda parameters in {csv_file}")
            return index_list, parameter_names
            
    except FileNotFoundError:
        logger.error(f"CSV file not found: {csv_file}")
        raise
    except Exception as e:
        logger.error(f"Failed to detect indexes in {csv_file}: {e}")
        raise ValueError(f"Index detection failed: {e}")


def process_csv_environment(input_csv: Union[str, Path], output_csv: Union[str, Path]) -> None:
    """Process CSV file to handle environment contributions.
    
    Args:
        input_csv: Input CSV file path
        output_csv: Output CSV file path
    """
    try:
        index_columns, parameter_names = _detect_lambda_indexes(input_csv)
        env_control_lambda = 'lambda_electrostatics_env'
        
        if env_control_lambda in parameter_names:
            logger.debug(f"Found environment control parameter: {env_control_lambda}")
            lambda_parameters = parameter_names.copy()
            lambda_parameters.pop(0)  # Remove first parameter
            # the first item in list is 'time(ps)', remove it.
            env_control_index = parameter_names.index(env_control_lambda)
        else:
            # warnings.warn(f'Could not find env control lambda in {os.path.basename(input_csv)}, use the origin csv.')
            shutil.copy(input_csv, output_csv)
            return False
    except:
        pass

    df = pd.read_csv(input_csv, index_col=index_columns, delimiter='|')
    # index_list = ['times(ps)', 'lambda_restraints', 'lambda_electrostatics',
    #               'lambda_sterics','lambda_electrostatics_env', 'lambda_electrostatics_chgprod_square']
    for i in df.index.names:
        if i not in parameter_names:
            df = df.reset_index(level=i, drop=True)
    headers = df.columns.tolist()
    tuples = []
    for h in headers:
        if h.startswith('('):
            single_list = []
            for ll in h.strip('()').split(', '):
                single_list.append(float(ll))
            tuples.append(tuple(single_list))
    df.columns = tuples
    intra_mol_tuples = [t for t in tuples if t[env_control_index] == 0.0]
    all_mol_tuples = [t for t in tuples if t[env_control_index] == 1.0]
    intra_mol_ene = df[intra_mol_tuples]
    all_mol_ene = df[all_mol_tuples]
    mol_env_ene_data = {}
    for intra_t in intra_mol_tuples:
        for all_t in all_mol_tuples:
            if intra_t[:3] == all_t[:3]:
                new_col = tuple(intra_t[:3])
                mol_env_ene_data[new_col] = all_mol_ene[all_t]
                break
    mol_env_ene = pd.DataFrame(mol_env_ene_data, index=df.index, columns=list(mol_env_ene_data.keys()))
    mol_env_ene.to_csv(output_csv, sep='|')

    del2_df = pd.read_csv(output_csv, sep='|', index_col=None)
    remove_columns = ['lambda_electrostatics_env', 'lambda_electrostatics_chgprod_square']
    time = [f'{int(i) / 10:.1f}' for i in range(len(del2_df))]
    time_df = pd.DataFrame({'times(ps)': time})
    for remove_column in remove_columns:
        if remove_column in del2_df.columns:
            del2_df.drop(remove_column, axis=1, inplace=True)

    del2_df = pd.concat([time_df, del2_df], axis=1)
    del2_df.to_csv(output_csv, sep='|', index=False)


def md_success_check(md_dir, log_name='run.log'):
    log_path = os.path.join(md_dir, log_name)
    if not os.path.exists(log_path):
        return False
    with open(log_path, 'r') as f:
        lines = f.readlines()[-2:]
        for line in lines:
            if 'ERROR' in line or 'error' in line or 'Error' in line:
                return False
    return True


class SingleProcessCalculateFileManager:
    def __init__(self, _setting: AllSettings, _path_to_csvs: str):
        self._path_to_csvs = os.path.abspath(_path_to_csvs)
        self.setting = _setting

    def prepare_calculate_files(self, output_dir, group_num: int = 1):
        tmp_apply_dict = {'Basic_settings': {'file_prefix': ene_csv_file_prefix.replace('<n>', str(group_num))}}
        self.setting.apply_settings_with_dict(tmp_apply_dict)
        open(os.path.join(output_dir, calculate_input_filename), 'w').close()
        self.setting.gen_file(os.path.join(output_dir, calculate_input_filename))

        with open(os.path.join(output_dir, calculate_submit_bash_filename), 'w') as f:
            f.write(f'cd {output_dir}\n')
            f.write(f'python {calculate_py_script} -i {calculate_input_filename} > calculate.log 2>&1')
        # Copy files here.
        tmp_csv_pattern = ene_csv_file_prefix.replace('<n>', str(group_num)).replace("'", "")
        tmp_csv_pattern += "*.csv"
        csvs = glob.glob(os.path.join(self._path_to_csvs, tmp_csv_pattern))
        if len(csvs) == 0:
            logger.warning(f'No csv files found in {"_".join(self._path_to_csvs.split("/")[-3:])}')
            shutil.rmtree(output_dir)
            return False

        for csv in csvs:
            csv = os.path.join(os.path.join(self._path_to_csvs, csv))
            process_csv_environment(csv, os.path.join(output_dir, os.path.basename(csv)))

        return f'{output_dir} {os.path.join(output_dir, calculate_submit_bash_filename)}\n'


@dataclass
class _path_with_deduct_mark:
    path: str
    need_deduct: bool = False

    def __str__(self):
        return self.path


@dataclass
class ConnectMethod:
    name: str
    _path_connector_pattern: str
    need_deduct: bool = False

    def __call__(self, _root_path):
        return_list = []
        for path in glob.glob(os.path.join(_root_path, self._path_connector_pattern)):
            return_list.append(_path_with_deduct_mark(path=path, need_deduct=self.need_deduct))
        return return_list


class _simpleConnectMethod(ConnectMethod):
    def __call__(self, _root_path):
        return [_path_with_deduct_mark(_root_path, self.need_deduct)]


class PathCollector:
    def __init__(self, _connect_method_list: list, elements: list = None):
        self.elements = elements
        if self.elements is None:
            self.elements = ['lM2A', 'lM2B', 'cM2A', 'cM2B']
        self._connect_method_list = _connect_method_list
        pass

    def __call__(self, _root_path):
        if len(self._connect_method_list) == 0:
            raise Exception('No connect method defined for PathCollector')
        return_path_list = []
        for connect_method in self._connect_method_list:
            for element in self.elements:
                for connected_path_tuple in connect_method(os.path.join(_root_path, element)):
                    return_path_list.append(connected_path_tuple)

        return return_path_list


class MainCalculateFileManager:
    def __init__(self, _setting_template_file: str, _path_to_pairs: str, _path_collector: PathCollector,
                 _pairs_list: list = None):
        self._setting = AllSettings.from_file(_setting_template_file)
        self._path_to_pairs = _path_to_pairs
        self._path_collector = _path_collector
        self.calculate_work_queue_lst_file = None
        self._pairs_list = []
        if _pairs_list is None:
            tmp_file_list = os.listdir(self._path_to_pairs)
            for tmp_file in tmp_file_list:
                if os.path.isdir(os.path.join(self._path_to_pairs, tmp_file)) and tmp_file.find('-') != -1:
                    self._pairs_list.append(tmp_file)
        else:
            self._pairs_list = _pairs_list

    def single_pair_prepare(self, pair, _output_root_dir: str = None, _is_seg: bool = False,
                            if_bygrp: bool = False, if_amber: bool = False):
        origin_path = os.path.join(self._path_to_pairs, pair)
        process_paths = self._path_collector(origin_path)
        return_process = []
        for process_path in process_paths:

            err_happened = False
            _path = process_path.path
            _path_deduct_mark = process_path.need_deduct

            destination_path = _path.replace(self._path_to_pairs, _output_root_dir)
            os.makedirs(destination_path, exist_ok=True)
            now_group_num = 0
            method = destination_path.split('/')[-1]
            if method.find('group_') != -1:
                now_group_num = int(method.replace('group_', ''))

            if not if_bygrp:
                now_group_num = 1

            if if_amber:
                lambda_descending = False
            else:
                lambda_descending = True

            if not md_success_check(_path):
                logger.error(f'MD failed in [{pair}] {_path}, skipping this pair')
                shutil.rmtree(os.path.join(_output_root_dir, pair))
                return []

            if _is_seg:
                err_happened = _no_deduct_segment_preprocess(single_process_path=_path,
                                                             group_num=now_group_num,
                                                             lambda_descending=lambda_descending,
                                                             )

            if not err_happened:
                tmp_single_process_file_manager = SingleProcessCalculateFileManager(self._setting, _path)
                pcmd = tmp_single_process_file_manager.prepare_calculate_files(output_dir=destination_path,
                                                                               group_num=now_group_num
                                                                               )
                return_process.append(pcmd)

        return return_process

    def single_prepare_for_multiprocessing(self, args):
        now_id = args['id']
        job_nums = args['job_nums']
        pair = args['pair']
        _output_root_dir = args['_output_root_dir']
        _is_seg = args['_is_seg']
        if_bygrp = args['if_bygrp']
        if_amber = args['if_amber']
        logger.info(f"Multiprocessing file preprocess started: ({now_id}/{job_nums}) Section [1/2]")
        return self.single_pair_prepare(pair=pair,
                                        _output_root_dir=_output_root_dir,
                                        _is_seg=_is_seg,
                                        if_bygrp=if_bygrp,
                                        if_amber=if_amber,
                                        )

    def prepare_files(self, _output_root_dir: str = None, _is_seg: bool = False,
                      if_bygrp: bool = False, if_amber: bool = True):
        if _output_root_dir is None:
            _output_root_dir = os.path.join(os.getcwd(), 'calculate')
            os.makedirs(_output_root_dir, exist_ok=True)
        _output_root_dir = os.path.abspath(_output_root_dir)

        calculate_work_queue_lst_file = os.path.join(os.path.dirname(_output_root_dir), 'calculate_queue.lst')
        err_pair_lst_file = os.path.join(os.path.dirname(_output_root_dir), 'err_pairs.lst')
        open(calculate_work_queue_lst_file, 'w').close()

        job_pool = []

        job_nums = len(self._pairs_list)
        for i, pair in enumerate(self._pairs_list):
            tmp_args = {
                'pair': pair,
                '_output_root_dir': _output_root_dir,
                '_is_seg': _is_seg,
                'if_bygrp': if_bygrp,
                'id': i+1,
                'job_nums': job_nums,
                'if_amber': if_amber,
            }
            job_pool.append(tmp_args)
        with Pool(processes=MAX_WORKER_COUNT) as p:
            return_work_queue_lst = p.map(self.single_prepare_for_multiprocessing, job_pool)

        actual_work_queue_lst = []
        error_pairs = []
        for i, lst in enumerate(return_work_queue_lst):
            if len(lst) == 0:
                error_pairs.append(self._pairs_list[i])
            actual_work_queue_lst.extend(lst)

        if error_pairs:
            with open(err_pair_lst_file, 'w') as f:
                for pair in error_pairs:
                    f.write(f"{pair}\n")

        logger.info(f"File preparation finished: {job_nums - len(error_pairs)}/{job_nums} pairs successful")

        with open(calculate_work_queue_lst_file, 'w') as f:
            for line in actual_work_queue_lst:
                f.write(line)

        self.calculate_work_queue_lst_file = calculate_work_queue_lst_file


def _no_deduct_segment_preprocess(single_process_path, group_num: int, lambda_descending: bool = True,
                                  prefix_pattern: str = 'lambda', suffix_pattern: str = '.csv', sep: str = '_'):
    single_process_path = os.path.abspath(single_process_path)
    no_deduct_result_folder = 'ana_used_data'
    data_path = os.path.join(single_process_path, no_deduct_result_folder)
    if not os.path.exists(data_path):
        return True
    file_compare_list = []
    for csv_file in os.listdir(data_path):
        numbers = csv_file.replace(prefix_pattern,'').replace(suffix_pattern, '')
        sum_up_value = sum([float(number) for number in numbers.split(sep)])
        file_compare_list.append((csv_file, sum_up_value))

    ranked_file_list = sorted(file_compare_list, key=lambda x: x[1], reverse=lambda_descending)

    for i, (file, value) in enumerate(ranked_file_list):
        new_csv_name = f'state_g{group_num}_s{i}.csv'
        shutil.copy(os.path.join(data_path, file), os.path.join(single_process_path, new_csv_name))


def run_submit_bash(all_args):
    now_job, job_nums = all_args[2]
    command = all_args[1]
    work_dir = all_args[0]
    logger.info(f"Calculating job ({now_job}/{job_nums})")
    
    try:
        if USE_SYSTEM_CALL and len(command) >= 2 and command[0] == '/bin/bash':
            # Pyarmor compatibility mode: use subprocess.Popen with shell=True for complete isolation
            bash_script = command[1]
            
            # Use subprocess.Popen with shell=True for complete process isolation
            process = subprocess.Popen(
                f'bash "{bash_script}"',
                cwd=work_dir,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=dict(os.environ)  # Pass clean environment
            )
            
            stdout, stderr = process.communicate(timeout=600)
            exit_code = process.returncode
            
            if stderr and stderr.strip():
                logger.warning(f"Job {now_job} stderr: {stderr.strip()}")
            
            if exit_code != 0:
                logger.error(f"Job {now_job} failed with exit code {exit_code}. Path: {work_dir}")
                if stderr and stderr.strip():
                    logger.error(f"Job {now_job} error details: {stderr.strip()}")
                return 1
            else:
                logger.info(f"Job {now_job} completed successfully")
                return 0
        else:
            # Standard subprocess mode for better error handling
            if len(command) >= 2 and command[0] == '/bin/bash':
                bash_script = command[1]
                shell_command = f'bash "{bash_script}"'
            else:
                shell_command = ' '.join(f'"{arg}"' for arg in command)
            
            result = subprocess.run(
                shell_command,
                cwd=work_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=600,
                shell=True
            )
            
            if result.stderr and result.stderr.strip():
                logger.warning(f"Job {now_job} stderr: {result.stderr.strip()}")

            if result.returncode != 0:
                logger.error(f"Job {now_job} failed with return code {result.returncode}")
                if result.stderr and result.stderr.strip():
                    logger.error(f"Job {now_job} error details: {result.stderr.strip()}")
                return 1
            
            logger.info(f"Job {now_job} completed successfully")
            return 0
        
    except subprocess.TimeoutExpired:
        logger.error(f"Job {now_job} timed out after 600 seconds")
        return 1
    except Exception as e:
        logger.error(f"Job {now_job} failed with exception: {str(e)}")
        return 1


def auto_cal_and_stat(cal_queue, cal_dir, work_dir):
    all_cal_jobs = []
    with open(cal_queue, 'r') as f:
        lines = f.readlines()
        total_lines = len(lines)
        for i, line in enumerate(lines):
            elements = line.strip().split(" ")
            tmp_dir = elements[0]
            submit_bash = elements[1]
            all_cal_jobs.append((tmp_dir, ['/bin/bash', submit_bash], (i+1, total_lines)))

    warnings.filterwarnings("ignore")

    rerun_count = 0
    failed_jobs = {}
    
    while rerun_count < 3:
        jobs_to_run = all_cal_jobs if rerun_count == 0 else list(failed_jobs.values())
        
        if not jobs_to_run:
            break
            
        logger.info(f"Running {'initial' if rerun_count == 0 else 'failed'} jobs in parallel (attempt {rerun_count+1}/{MAX_RETRIES})")
        failed_jobs = {}
        
        with futures.ProcessPoolExecutor(max_workers=MAX_WORKER_CAL) as executor:
            future_to_job = {executor.submit(run_submit_bash, job): job for job in jobs_to_run}
            
            for future in futures.as_completed(future_to_job):
                job = future_to_job[future]
                job_id = job[2][0]
                try:
                    result = future.result()
                    if result != 0:
                        failed_jobs[job_id] = job
                except Exception as exc:
                    logger.error(f"Job {job_id} generated an exception: {exc}")
                    failed_jobs[job_id] = job
        
        if not failed_jobs:
            logger.info("All jobs completed successfully!")
            break
            
        logger.warning(f"{len(failed_jobs)} jobs failed. Retrying...")
        rerun_count += 1
    
    if failed_jobs:
        logger.warning(f"{len(failed_jobs)} jobs still failed after {MAX_RETRIES} attempts: {', '.join(map(str, failed_jobs.keys()))}")

    stat_script = os.path.join(script_dir, 'stat_result.py')
    analyze_args = ['python', stat_script, '-W', work_dir, '-C', cal_dir]
    logger.info('Beginning statistical analysis...')
    p = subprocess.Popen(analyze_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         encoding='utf-8', errors='replace',
                         universal_newlines=True, bufsize=1)
    for line in iter(p.stdout.readline, ''):
        output = line
        if output == '' and p.poll() is not None:
            break
        if output:
            output = output.strip()
            logger.info(output)
    for line in iter(p.stderr.readline, ''):
        output = line
        if output == '' and p.poll() is not None:
            break
        if output:
            output = output.strip()
            logger.info(output)
    p.wait()


if __name__ == '__main__':
    ori_dir = os.getcwd()

    arg_parser = argparse.ArgumentParser(description='Copy calculate files after working. Make a calculate dir on cwd.')
    arg_parser.add_argument('-w', dest='work_dir', type=str, default='./work',
                            help='path to working directory, like ./work/')
    arg_parser.add_argument('-m', '--mode', default='Amber', dest='mode',
                            help='The MD engine.')
    arg_parser.add_argument('-auto', nargs='?', dest='auto_run', const=True, default=False,
                            help='auto calculate energy and stat consequently after the preparation.')

    args = arg_parser.parse_args()

    path_to_pairs = os.path.abspath(args.work_dir)
    if_CAR = True
    if_auto_run = args.auto_run
    if args.mode.lower() == 'amber':
        if_amber = True
    else:
        if_amber = False

    os.chdir(os.path.dirname(path_to_pairs))

    simple_path_collector = _simpleConnectMethod('no_bygrp', '', need_deduct=True)

    logger.info(f'Configuration - path_to_pairs: {path_to_pairs}')

    output_root_dir = os.path.join(os.path.dirname(path_to_pairs), 'calculate')
    if os.path.exists(output_root_dir):
        shutil.rmtree(output_root_dir)
    os.makedirs(output_root_dir, exist_ok=True)

    active_path_collector = PathCollector([simple_path_collector])

    a = MainCalculateFileManager(_setting_template_file=default_input_template,
                                 _path_to_pairs=path_to_pairs,
                                 _path_collector=active_path_collector)
    a.prepare_files(_output_root_dir=output_root_dir, _is_seg=if_CAR)

    if if_auto_run:
        auto_cal_and_stat(a.calculate_work_queue_lst_file, output_root_dir, path_to_pairs)

    os.chdir(ori_dir)

    pass
