import copy, pickle, io
import os, argparse, subprocess, shutil
import sys
import re

from typing import List, Any, Union, Tuple

from PIL import Image, ImageFont
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PyPDF2 import PdfReader, PdfWriter
from rdkit.Chem import AllChem

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak
from reportlab.platypus import Image as reportImg
from reportlab.lib.units import inch
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfgen import canvas

from analyze_tools.static_pdf_page import PDF_TITLE_STYLE, TABLE_STYLE, static_info_page, cover_page, delete_pics

from analyze_tools.DateParser import DateParser
from analyze_tools.fast_pdf import PDFGenerator, PDFFormatter
from analyze_tools.SettingManager import get_default_path
from plotting_tools.plotting_tools import PLOTTING
from prepare_input_tools.pair_pdb_generator import MolCoordModer

WCC_PATH = os.path.join(get_default_path('Wcc_Path'), 'wcc_main.py')

SUMMARY_PDF_FILENAME = 'summary.pdf'

edges = ['lM2A', 'lM2B', 'cM2A', 'cM2B']

collect_font = canvas.Canvas("font_check.pdf", pagesize=letter)
avail_fonts = collect_font.getAvailableFonts()
PDF_font = 'Helvetica'
if PDF_font not in avail_fonts:
    PDF_font = 'Courier'
PDF_bold_font = f"{PDF_font}-Bold"

default_convergence_files = {'forward_df_com.csv', 'reverse_df_com.csv', 'moving_df_com.csv'}
default_convergence_folder = 'fe_cal_out/time_serial_check'

default_edges = {'lM2A', 'lM2B', 'cM2A', 'cM2B'}
default_convergence_methods = {'forward', 'reverse', 'moving'}


def __pair_ene_add_decider(csv_path):
    if csv_path.find('lM2A') != -1:
        return True
    elif csv_path.find('cM2A') != -1:
        return False
    elif csv_path.find('lM2B') != -1:
        return False
    elif csv_path.find('cM2B') != -1:
        return True
    else:
        raise Exception(f'Could not find edge value in {csv_path}')


def formalize_edge_name(input_str) -> str:
    if input_str.find('lM2A') != -1:
        return input_str.replace('lM2A', 'lig_CStoA')
    elif input_str.find('lM2B') != -1:
        return input_str.replace('lM2B', 'lig_CStoB')
    elif input_str.find('cM2A') != -1:
        return input_str.replace('cM2A', 'com_CStoA')
    elif input_str.find('cM2B') != -1:
        return input_str.replace('cM2B', 'com_CStoB')
    return input_str


def merge_convergence_csvs(csv_list, output_csv, std_step_ratio: float = 0.1, if_pair: bool = False) -> bool:
    """
    Merge csv in csv_list into output_csv (Sum up all convergence data)
    :param csv_list: a list containing csvs to merge
    :param output_csv: output csv file
    :param std_step_ratio: the ratio use to rolling standard deviation
    :param if_pair: if merging a pair's 4 edges file
    :return: True if merge succeed, else False
    """
    free_energy_column_title = 'free_energy(kcal/mol)'
    rolling_std_column_title = 'Rolling_STD(kcal/mol)'

    if len(csv_list) == 0:
        print('No files to merge.')
        return False

    # print(f'Merging {csv_list}')

    sum_df = None
    for csv in csv_list:  # Allow using one existed csv to set up the sum_df format
        try:
            sum_df = pd.read_csv(csv, delimiter='|')
            sum_df.drop(columns=rolling_std_column_title, inplace=True)
            sum_df[free_energy_column_title] = 0
            break
        except FileNotFoundError:
            print(f'not csv was found.')

    if sum_df is None:
        return False

    std_step = int(np.floor(len(sum_df)*std_step_ratio))


    for csv in csv_list:
        if not os.path.exists(csv):
            continue

        tmp_df = pd.read_csv(csv, delimiter='|')
        if if_pair:
            if __pair_ene_add_decider(csv):
                sum_df[free_energy_column_title] = sum_df[free_energy_column_title] + tmp_df[free_energy_column_title]
                # print(f'Add {csv}')
            else:
                # print(f'deduct {csv}')
                sum_df[free_energy_column_title] = sum_df[free_energy_column_title] - tmp_df[free_energy_column_title]
        else:
            sum_df[free_energy_column_title] = sum_df[free_energy_column_title] + tmp_df[free_energy_column_title]

    to_cal_df = pd.DataFrame(sum_df.loc[:, free_energy_column_title])
    to_cal_df.columns = [rolling_std_column_title]

    std_df = to_cal_df.rolling(std_step).std()

    nan_replace_value = std_df.loc[std_step-1, rolling_std_column_title]
    for i in range(std_step):
        std_df.loc[i, rolling_std_column_title] = nan_replace_value

    all_df = pd.concat([sum_df, std_df], axis=1)
    all_df.to_csv(output_csv, index=False, encoding='utf-8', header=True, sep='|')


def get_normal_processes_folders(path) -> list:
    return [os.path.abspath(path)]


def _connect_two_tables(title: List[str], table1: List[List], table2: Union[List[List], None]) -> List[List[str]]:
    if table2 is None:
        return [title]+table1
    full_title = [*title, ' ', *title]  # Sep
    len_table2 = len(table2)
    len_normal_line = len(title)
    full_line = []
    for i, line in enumerate(table1):
        if i >= len_table2:
            extend_table2_line = [" " for _ in range(len_normal_line)]
        else:
            extend_table2_line = table2[i]
        full_line.append([*line, ' ', *extend_table2_line])

    return [full_title]+full_line


def formalize_df(df: pd.DataFrame, dec_format=None, specialize_title=None) -> list:
    """
    ['ligand1', 'ligand2', 'ΔΔG_ori',  'ΔΔG_corrected'] (0, 0, 3, 3)
    ['ligand_ID', 'ΔG'] (0, 3)
    :param df:
    :param dec_format:
    :param specialize_title:
    :return:
    """
    using_dec_format = (0, 0, 3, 3)
    if dec_format is not None:
        using_dec_format = dec_format

    title_list = df.columns.tolist()
    if specialize_title is not None:
        title_list = specialize_title
    formalized_title_list = []
    for title in title_list:
        if title.lower().find('unnamed') != -1:
            title = ' '
        title = title.split('.')[0]
        formalized_title_list.append(title)

    contents = df.values.tolist()
    formalized_contents = []
    for line in contents:
        formalized_line = []
        for i, ele in enumerate(line):
            if str(ele) == 'nan' or ele is None:
                formalized_line.append(' ')
                continue
            if using_dec_format[i] == 0:
                formalized_line.append(ele)
            else:
                if isinstance(ele, str):
                    ele = float(ele)
                boo = f" {ele:.{using_dec_format[i]}f}" if ele >= 0 else f"{ele:.{using_dec_format[i]}f}"
                formalized_line.append(boo)
        formalized_contents.append(formalized_line)

    # Fix pagination logic
    line_per_page = 32
    cut_contents = []
    for i in range(0, len(formalized_contents), line_per_page):
        page = formalized_contents[i:i+line_per_page]
        if page:  # Avoid empty pages
            cut_contents.append(page)

    pages = []
    # while len(cut_contents) > 2:
    #     table1 = cut_contents.pop(0)
    #     table2 = cut_contents.pop(0)
    #     merge_table = _connect_two_tables(formalized_title_list, table1, table2)
    #     pages.append(merge_table)
    # if len(cut_contents) > 0:
    #     pages.append(_connect_two_tables(formalized_title_list, cut_contents[0], None))
    for table in cut_contents:
        pages.append(_connect_two_tables(formalized_title_list, table, None))

    return pages


def df_to_pages(df: pd.DataFrame, dec_format: tuple, title=None) -> list:
    pages = formalize_df(df, dec_format=dec_format)
    all_pages = []
    for page in pages:
        table = Table(page)
        page_title_line = list(page[0])
        page_merged = True if page_title_line.count(page_title_line[0]) > 1 else False
        line_len = len(page_title_line)

        table_style = TABLE_STYLE

        table.setStyle(TableStyle(table_style))

        if title is not None:
            page_title = Paragraph(title, style=PDF_TITLE_STYLE)
            page_ele = [page_title, table, PageBreak()]
        else:
            page_ele = [table, PageBreak()]
        all_pages.extend(page_ele)
    return all_pages


def moles_to_pages(mol_list: list, title: str = None, page_matrix=(4, 5)):
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import rdDepictor

    mol_per_row = page_matrix[0]
    row_per_page = page_matrix[1]
    mol_per_page = mol_per_row * row_per_page

    mol_list_2d = []
    for mol in mol_list:
        tmp_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol, allHsExplicit=True, allBondsExplicit=True))
        AllChem.Compute2DCoords(tmp_mol)
        _name = mol.GetProp("_Name") if mol.HasProp("_Name") else ''
        tmp_mol.SetProp("_Name", _name)
        mol_list_2d.append(tmp_mol)
    MCS_bone = MolCoordModer.customized_FindMCS(mol_list_2d)
    MCS_bone_mol = MCS_bone.queryMol
    rdDepictor.Compute2DCoords(MCS_bone_mol)

    highlight_atoms_list = []

    for _mol in mol_list_2d:
        match = _mol.GetSubstructMatch(MCS_bone_mol)
        if match:
            rdDepictor.GenerateDepictionMatching2DStructure(_mol, MCS_bone_mol)
            highlight_atoms_list.append(list(match))
        else:
            rdDepictor.Compute2DCoords(_mol)

    draw_options = Draw.MolDrawOptions()

    draw_options.atomLabelFontSize = 12
    draw_options.legendFontSize = 22
    draw_options.bondLineWidth = 2
    draw_options.highlightRadius = 0.3
    draw_options.highlightBondWidthMultiplier = 10
    draw_options.variableAtomRadius = 0.2
    draw_options.setHighlightColour((170 / 255, 220 / 255, 255 / 255, 1))

    page_moles = []
    page_legends = []
    page_highlight_atoms = []
    while len(mol_list_2d) > mol_per_page:
        slice_moles = mol_list_2d[:mol_per_page]
        page_moles.append(slice_moles)
        page_legends.append([_mol.GetProp('_Name') if _mol.HasProp('_Name') else '' for _mol in slice_moles])
        page_highlight_atoms.append(highlight_atoms_list[:mol_per_page])

        mol_list_2d = mol_list_2d[mol_per_page:]
        highlight_atoms_list = highlight_atoms_list[mol_per_page:]

    if len(mol_list_2d) > 0:
        page_moles.append(mol_list_2d)
        page_legends.append([_mol.GetProp('_Name') if _mol.HasProp('_Name') else '' for _mol in mol_list_2d])
        page_highlight_atoms.append(highlight_atoms_list)

    all_pages = []
    for i in range(len(page_moles)):
        page_matrix_mol = []
        page_matrix_legend = []
        page_matrix_highlight_atoms = []

        now_page_moles = copy.deepcopy(page_moles[i])
        now_page_legends = copy.deepcopy(page_legends[i])
        now_page_highlight_atoms = copy.deepcopy(page_highlight_atoms[i])

        gap_to_full_page = mol_per_page - len(now_page_moles)
        now_page_moles.extend([None for _ in range(gap_to_full_page)])
        now_page_legends.extend([" " for _ in range(gap_to_full_page)])
        now_page_highlight_atoms.extend([[] for _ in range(gap_to_full_page)])

        for j in range(row_per_page):
            page_matrix_mol.append(now_page_moles[mol_per_row * j: mol_per_row * (j + 1)])
            page_matrix_legend.append(now_page_legends[mol_per_row * j: mol_per_row * (j + 1)])
            page_matrix_highlight_atoms.append(now_page_highlight_atoms[mol_per_row * j: mol_per_row * (j + 1)])

        img = Draw.MolsMatrixToGridImage(
            page_matrix_mol,
            legendsMatrix=page_matrix_legend,
            subImgSize=(300, 300),
            highlightAtomListsMatrix=page_matrix_highlight_atoms,
            returnPNG=False,
            useSVG=False,
            drawOptions=draw_options,
        )

        buffer = io.BytesIO()
        img.save(buffer, format='PNG')
        buffer.seek(0)
        print_img = reportImg(buffer, width=444, height=555, kind='proportional')
        if title is not None:
            page_title = Paragraph(title, style=PDF_TITLE_STYLE)
            page_ele = [page_title, print_img, PageBreak()]
        else:
            page_ele = [print_img, PageBreak()]
        all_pages.extend(page_ele)
        # img.save(f'test{i}.png')
    return all_pages


class EdgeConvergenceMerger:

    def __init__(self, _path_to_edge, _processes_paths: list,
                 _convergence_files: set = default_convergence_files,
                 _convergence_folder: str = default_convergence_folder):
        self._path_to_edge = os.path.abspath(_path_to_edge)
        self._processes_paths = _processes_paths

        self._convergence_files = _convergence_files
        self._convergence_folder = _convergence_folder

    def merge(self, _output_path: str) -> None:
        """
        Merge an edge's diff convergence data into diff section csv. (csv in self._convergence_files)
        :param _output_path: output path of the merged csv files.
        :return: None
        """
        processes_paths = self._processes_paths
        for file in self._convergence_files:
            section_csv_list = []
            for process_path in processes_paths:
                file_abs_path = os.path.join(process_path, self._convergence_folder, file)
                section_csv_list.append(file_abs_path)
            # print(f"List {section_csv_list}")
            merge_convergence_csvs(section_csv_list, os.path.join(_output_path, file))


class MainConvergenceMerger:

    def __init__(self, _path_to_pairs: str, pair_list: list = None,
                 _convergence_files: set = default_convergence_files,
                 _convergence_folder: str = default_convergence_folder,
                 _convergence_method: set = default_convergence_methods,
                 _edges: set = default_edges
                 ):
        self._path_to_pairs = _path_to_pairs

        self._convergence_files = _convergence_files
        self._convergence_folder = _convergence_folder
        self.__convergence_methods = _convergence_method
        self._edges = _edges

        self._pairs_list = []
        if pair_list is None:
            tmp_file_list = os.listdir(self._path_to_pairs)
            for tmp_file in tmp_file_list:
                if os.path.isdir(os.path.join(self._path_to_pairs, tmp_file)) and tmp_file.find('-') != -1:
                    self._pairs_list.append(tmp_file)
        else:
            self._pairs_list = pair_list

    @staticmethod
    def __formalize_df(_df: pd.DataFrame) -> pd.DataFrame:
        _df.rename(index=_df.iloc[:, 0], inplace=True)
        _df.drop(columns='time_ratio', inplace=True)
        return _df

    def merge(self, _output_path: str = None, get_folder_method=get_normal_processes_folders) -> None:
        if _output_path is None:
            _output_path = os.path.join(os.getcwd(), 'results')
        else:
            _output_path = os.path.abspath(_output_path)
        os.makedirs(_output_path, exist_ok=True)

        plot_tool = PLOTTING()

        # for pair in tqdm(self._pairs_list, desc='Merging energy', unit='pair'):
        for i, pair in enumerate(self._pairs_list):
            print(f'Merging {pair} ({i+1}/{len(self._pairs_list)}) Section [2/4]', flush=True)
            pair_output_path = os.path.join(_output_path, pair)
            os.makedirs(pair_output_path, exist_ok=True)

            pair_cal_path = os.path.join(self._path_to_pairs, pair)
            pair_convergence_files_dict = {filename: [] for filename in self._convergence_files}

            for edge in self._edges:
                edge_output_path = os.path.join(pair_output_path, edge)
                os.makedirs(edge_output_path, exist_ok=True)

                edge_cal_path = os.path.join(pair_cal_path, edge)
                processes_paths = get_folder_method(edge_cal_path)
                # print(processes_paths)
                tmp_edge_merger = EdgeConvergenceMerger(edge_cal_path, processes_paths,
                                                        _convergence_files=self._convergence_files,
                                                        _convergence_folder=self._convergence_folder)
                tmp_edge_merger.merge(edge_output_path)

                std_plot_dict = {}
                for method in self.__convergence_methods:
                    accord_csv = os.path.join(edge_output_path, f'{method}_df_com.csv')
                    accord_df = self.__formalize_df(pd.read_csv(accord_csv, sep='|'))
                    std_plot_dict[method] = accord_df
                    # print(f"{method}\n{accord_df}")

                plot_png = os.path.join(edge_output_path, 'edge_sum.png')
                fe = std_plot_dict['forward'].iloc[len(std_plot_dict['forward'])-1, 0]
                fe_std = std_plot_dict['forward'].iloc[len(std_plot_dict['forward'])-1, 0]
                plot_tool.plot_fe_time_serial(plot_png,
                                              moving=std_plot_dict['moving'],
                                              forward=std_plot_dict['forward'],
                                              reverse=std_plot_dict['reverse'],
                                              fe=fe,
                                              fe_std=fe_std
                                              )
                plt.close()
                # Draw plot here.

                for file in self._convergence_files:
                    pair_convergence_files_dict[file].append(os.path.join(edge_output_path, file))
                # Append same convergence file in each edge to a list in order to merge the whole pair's convergence

            for filename, filelist in pair_convergence_files_dict.items():
                merged_file = os.path.join(pair_output_path, filename)
                merge_convergence_csvs(filelist, merged_file, if_pair=True)
                # Merge each section in edges to pair folder.

            std_plot_dict = {}
            for method in self.__convergence_methods:
                accord_csv = os.path.join(pair_output_path, f'{method}_df_com.csv')
                accord_df = self.__formalize_df(pd.read_csv(accord_csv, sep='|'))
                std_plot_dict[method] = accord_df

            plot_png = os.path.join(pair_output_path, 'pair_sum.png')
            fe = std_plot_dict['forward'].iloc[len(std_plot_dict['forward'])-1, 0]
            fe_std = std_plot_dict['forward'].iloc[len(std_plot_dict['forward'])-1, 0]
            plot_tool.plot_fe_time_serial(plot_png,
                                          moving=std_plot_dict['moving'],
                                          forward=std_plot_dict['forward'],
                                          reverse=std_plot_dict['reverse'],
                                          fe=fe, fe_std=fe_std
                                          )
            plt.close()
            # Draw plot here.


class SingleStater:
    def __init__(self, _name, _ene_result_folder, _cal_log=None):
        self._ene_csv = os.path.join(_ene_result_folder, 'free_ene.csv')
        self._convergence_pic = os.path.join(_ene_result_folder, 'time_serial_check',
                                             'bar_time_series_forward_reverse_moving_sum_all.png')
        self.wall_time = 0  # minutes, type = int
        self.energy = 0.0  # float
        self._name = _name
        self.simu_time = 0

        if os.path.exists(self._ene_csv):
            with open(self._ene_csv, 'r') as f:
                last_line = f.readlines()[-1]
                energy = float(last_line.split('|')[1])
                self.energy = energy

        if _cal_log is not None:
            if os.path.exists(_cal_log):
                try:
                    with open(_cal_log, 'r') as f:
                        lines = f.readlines()
                        first_line = lines[0]
                        last_line = lines[-1]
                        begin_time = DateParser.str2datetime(first_line.replace("BEGIN AT:", "").strip())
                        end_time = DateParser.str2datetime(last_line.replace("END AT:", "").strip())
                        used_time_text = DateParser.get_time_gap_str(begin_time, end_time)
                        self.wall_time = DateParser.convert2min(used_time_text)

                        # Count all simul_time: xx ps
                        simu_time_sum = 0
                        for line in lines:
                            match = re.search(r'simul_time:\s*(\d+)ps', line)
                            if match:
                                simu_time_sum += int(match.group(1))
                        self.simu_time = simu_time_sum
                except:
                    self.wall_time = -9999 # May choose a better value

        if not os.path.exists(self._convergence_pic):
            self._convergence_pic = Image.new('RGBA',(1200, 900),color=(0, 0, 0, 0))

    @property
    def plot(self):
        # description_text = ''
        # if self.wall_time is not None:
        #     description_text = f'tc: {DateParser.convert_from_min(self.wall_time)}'
        # if self.energy is not None:
        #     description_text += f' ene: {self.energy:.2f}'
        description_text = None

        return self._convergence_pic, formalize_edge_name(self._name), description_text


class EdgeStater:
    def __init__(self, _name, _edge_cal_folder, _edge_work_folder, _log_name):
        self._name = _name
        self._edge_cal_folder = os.path.abspath(_edge_cal_folder)
        self._edge_work_folder = os.path.abspath(_edge_work_folder)
        self._log_name = _log_name

        self.charge_ene = 0.0
        self.vdw_ene = 0.0
        self.wall_time = 0
        self.simu_time = 0

        self.detail_ene_dict = {}
        self._print_content = []
        self.detail_ene_sorted_list = []

        self.total_ene = 0

        self.stat()

    @property
    def content(self):
        """
         -> (List[(Image.Image or str, str, str),], str, str)
        :return:
        """
        # description_text = ''
        # if self.wall_time is not None:
        #     description_text = f'tc: {DateParser.convert_from_min(self.wall_time)}'
        # description_text += f' ene: {self.total_ene:.2f}'
        description_text = None
        return self._print_content, formalize_edge_name(self._name), description_text

    def stat(self):
        _result_folder = os.path.join(self._edge_cal_folder, 'fe_cal_out')
        tmp_cal_log = os.path.join(self._edge_work_folder, self._log_name)
        tmp_single_stater = SingleStater(self._name, _result_folder, tmp_cal_log)

        self._print_content.append(tmp_single_stater.plot)
        self.total_ene = tmp_single_stater.energy
        self.wall_time = tmp_single_stater.wall_time
        self.simu_time = tmp_single_stater.simu_time

        self.detail_ene_dict[self._name] = tmp_single_stater.energy
        
        self.detail_ene_sorted_list = sorted(self.detail_ene_dict.items(), key=lambda item: item[0], reverse=False)


class MainStater:
    def __init__(self, _calculate_path: str, _work_path: str, _log_name: str,
                 pair_list=None, edge_stater=EdgeStater, results_path=None):
        self._calculate_path = _calculate_path
        self._work_path = _work_path
        self.pair_list = pair_list
        self._cal_pair_path_dict = {}
        self._work_pair_path_list = {}
        if results_path is not None:
            self.__stat_path = results_path
        else:
            self.__stat_path = os.path.join(os.path.dirname(self._calculate_path), 'results')
        self.ene_csv = os.path.join(self.__stat_path, 'results.csv')
        self.__log_name = _log_name

        self.__edge_stater = edge_stater

        self.page_content_dict = {}
        self._edge_description = {}
        self._pair_description = {}

        self._edge_info_dict = {}
        self._pair_info_dict = {}

        self.total_detail_ene_dict = {}

        self.print_content_dict = {}

        self.connent_group_num = 0
        self.simu_time = 0

        if self.pair_list is None:
            self.pair_list = []
            for folder in os.listdir(self._calculate_path):
                if os.path.isdir(os.path.join(self._calculate_path, folder)) and folder.find('-') != -1:
                    self.pair_list.append(folder)

        for pair in self.pair_list:
            self._cal_pair_path_dict[pair] = os.path.join(self._calculate_path, pair)
            self._work_pair_path_list[pair] = os.path.join(self._work_path, pair)

        self.clean_file()

    def get_pair_and_ligand_nums(self):
        pair_num = len(self.pair_list)
        ligand_set = set()
        for pair in self.pair_list:
            ligand_set.add(pair.split('-')[0])
            ligand_set.add(pair.split('-')[1])
        ligand_num = len(ligand_set)
        return pair_num, ligand_num

    @staticmethod
    def check_connected_loops(pair_list):

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
            node1 = pair.split('-')[0]
            node2 = pair.split('-')[1]
            connected_components.append((node1, node2))

        return find_connected_components_dfs(connected_components)

    def stat(self, output_path=None):
        if output_path is None:
            stat_path = self.__stat_path
        else:
            stat_path = output_path
        os.makedirs(stat_path, exist_ok=True)
        ene_csv = os.path.join(stat_path, 'results.csv')
        self.ene_csv = ene_csv
        node_content_lines = {}

        sus_pairs = set()
        self.simu_time = 0

        with open(ene_csv, 'w', encoding='utf-8-sig') as f:
            f.write('pair,ΔG_lig_CS_to_A,ΔG_lig_CS_to_B,ΔG_com_CS_to_A,ΔG_com_CS_to_B,total\n')

        # for pair in tqdm(self.pair_list, desc='Stating pairs', unit='pair'):
        for i, pair in enumerate(self.pair_list):
            print(f"Stating {pair} ({i + 1}/{len(self.pair_list)}) Section [1/4]", flush=True)
            pair_wall_time = 0
            pair_energy = {}
            pair_page_content = []
            pair_simu_time = 0

            node1 = pair.split('-')[0]
            node2 = pair.split('-')[1]
            if node1 not in node_content_lines:
                node_content_lines[node1] = {i}
            else:
                node_content_lines[node1].add(i)
            if node2 not in node_content_lines:
                node_content_lines[node2] = {i}
            else:
                node_content_lines[node2].add(i)

            self._edge_description[pair] = {}
            self._edge_info_dict[pair] = {}
            self.total_detail_ene_dict[pair] = {}

            for edge in edges:
                # print(edge)
                edge_cal_path = os.path.join(self._cal_pair_path_dict[pair], edge)
                edge_work_path = os.path.join(self._work_pair_path_list[pair], edge)
                tmp_edge_stater = self.__edge_stater(edge, edge_cal_path, edge_work_path,
                                                     _log_name=self.__log_name)
                pair_wall_time += tmp_edge_stater.wall_time
                pair_simu_time += tmp_edge_stater.simu_time
                # print(tmp_edge_stater.simu_time)
                # print(pair_wall_time)
                pair_energy[edge] = tmp_edge_stater.total_ene

                edge_print_content, edge_title, edge_description = tmp_edge_stater.content
                pair_page_content.extend(edge_print_content)

                self._edge_description[pair][edge] = edge_description
                self._edge_info_dict[pair][edge] = (tmp_edge_stater.total_ene,
                                                    f"{tmp_edge_stater.simu_time / 1000 :.2f} ns")
                self.total_detail_ene_dict[pair][edge] = tmp_edge_stater.detail_ene_sorted_list
                # pair_description += f'|<b>{edge}</b>|{edge_description}'

            total_energy = pair_energy['cM2B'] - pair_energy['cM2A'] - (pair_energy['lM2B'] - pair_energy['lM2A'])
            # pair_description = f'<b>Total tc:</b> {DateParser.convert_from_min(pair_wall_time)} <b>Total ene:</b> {total_energy:.2f}'
            pair_description = None
            # self._pair_info_dict[pair] = (total_energy, DateParser.convert_from_min(pair_wall_time))
            self._pair_info_dict[pair] = (total_energy, f"{pair_simu_time / 1000 :.2f} ns")

            self._pair_description[pair] = pair_description
            self.page_content_dict[pair] = {'content': pair_page_content, 'title': pair, 'description': pair_description}
            _csv_line = f"{pair},{pair_energy['lM2A']},{pair_energy['lM2B']},{pair_energy['cM2A']},{pair_energy['cM2B']},{total_energy:.3f}\n"
            self.print_content_dict[i] = _csv_line

            if abs(total_energy) > 4:
                sus_pairs.add(pair)

            with open(ene_csv, 'a') as f:
                f.write(_csv_line)
            
            self.simu_time += pair_simu_time
        
        groups = self.check_connected_loops(self.pair_list)
        self.connent_group_num = len(groups)
        if len(groups) > 1:
            print(f"Found {len(groups)} seperate groups in result.", flush=True)
            for group_counter, group in enumerate(groups):
                print(f"group {group_counter+1} : {group}")
                group_line_idxes = set()
                for node in group:
                    group_line_idxes.update(node_content_lines[node])
                with open(f'{ene_csv.split(".")[0]}_group{group_counter+1}.csv', 'w') as f:
                    f.write('pair,ΔG_lig_CS_to_A,ΔG_lig_CS_to_B,ΔG_com_CS_to_A,ΔG_com_CS_to_B,total\n')
                    for line_idx in group_line_idxes:
                        f.write(self.print_content_dict[line_idx])
        
        if len(sus_pairs) > 0:
            sus_csv = os.path.join(os.path.dirname(stat_path), 'sus_pairs.lst')
            with open(sus_csv, 'w') as f:
                for sus_pair in sus_pairs:
                    f.write(f"{sus_pair}\n")


    def get_pdf_conv_tables(self):
        all_data_tables = []
        sorted_pair_list = sorted(self._edge_info_dict.keys(), reverse=False)
        for pair in sorted_pair_list:
            now_info_dict = self._edge_info_dict[pair]
            temp_table = [
                ['Edge', 'Free\nenergy\n(kcal/mol)', 'Simu.\ntime'],
                ['lig_CStoA', f'{now_info_dict["lM2A"][0]:.3f}', f'{now_info_dict["lM2A"][1]}'],
                ['lig_CStoB', f'{now_info_dict["lM2B"][0]:.3f}', f'{now_info_dict["lM2B"][1]}'],
                ['com_CStoA', f'{now_info_dict["cM2A"][0]:.3f}', f'{now_info_dict["cM2A"][1]}'],
                ['com_CStoB', f'{now_info_dict["cM2B"][0]:.3f}', f'{now_info_dict["cM2B"][1]}'],
                ['Total', f'{self._pair_info_dict[pair][0]:.3f}', f'{self._pair_info_dict[pair][1]}'],
            ]
            all_data_tables.append(temp_table)
        return all_data_tables

    def append_sum_convergence_plots(self, if_edge_sum: bool = True):
        stat_path = self.__stat_path
        for pair in self.pair_list:
            pair_stat_path = os.path.join(stat_path, pair)
            pair_convergence_plots = []

            if if_edge_sum:
                for edge in edges:
                    edge_sum_plot_img = os.path.join(pair_stat_path, edge, 'edge_sum.png')
                    pair_convergence_plots.append((edge_sum_plot_img, f'<b>{edge}</b>', self._edge_description[pair][edge]))

            pair_sum_plot_img = os.path.join(pair_stat_path, 'pair_sum.png')
            pair_convergence_plots.append((pair_sum_plot_img, f'<b>{pair}</b>', self._pair_description[pair]))

            self.page_content_dict[pair]['content'].extend(pair_convergence_plots)

    def gen_detail_ene_csv(self, output_csv=None):
        def get_thermo_ene(_ene_dict):
            return _ene_dict['cM2B'] - _ene_dict['cM2A'] - _ene_dict['lM2B'] + _ene_dict['lM2A']
        if output_csv is None:
            output_csv = os.path.join(self.__stat_path, 'ene_detail.csv')

        detail_df = pd.DataFrame(columns=['pair'])

        sep_detail_df = pd.DataFrame(columns=['pair_name', 'total_chg', 'total_vdw'])

        # detail_merge_df = pd.DataFrame(columns=['pair'])
        line_count = 0

        for pair, pair_ene_detail_dict in self.total_detail_ene_dict.items():
            _pair_each_process_ene_dict = {}
            _pair_chg_sep_ene_dict = {}
            _pair_vdw_sep_ene_dict = {}

            for edge, edge_ene_detail_tuple in pair_ene_detail_dict.items():
                now_pair_name = f"{pair}_{formalize_edge_name(edge)}"
                detail_df.loc[now_pair_name] = ''
                detail_df.loc[now_pair_name, 'pair'] = now_pair_name

                _edge_process_dict = {}
                _edge_chg_ene = 0
                _edge_vdw_ene = 0

                for process_name, process_ene in edge_ene_detail_tuple:
                    if process_name not in detail_df.columns.values:
                        detail_df.loc[:, process_name] = None

                    detail_df.loc[now_pair_name, process_name] = process_ene

                    _edge_process_dict[process_name] = process_ene

                    if process_name.find("charge") != -1:
                        _edge_chg_ene += process_ene
                    elif process_name.find("vdw") != -1:
                        _edge_vdw_ene += process_ene

                _pair_each_process_ene_dict[edge] = _edge_process_dict
                _pair_vdw_sep_ene_dict[edge] = _edge_vdw_ene
                _pair_chg_sep_ene_dict[edge] = _edge_chg_ene

            pair_chg_ene = get_thermo_ene(_pair_chg_sep_ene_dict)
            pair_vdw_ene = get_thermo_ene(_pair_vdw_sep_ene_dict)
            sep_detail_df.loc[line_count] = [pair, pair_chg_ene, pair_vdw_ene]
            # print(line_count)

            line_count += 1

        detail_df.loc[:, None] = None
        
        original_index = detail_df.index
        
        detail_df = detail_df.reset_index(drop=True)
        sep_detail_df = sep_detail_df.reset_index(drop=True)
        new_detail_df = pd.concat([detail_df, sep_detail_df], axis=1)
        
        new_detail_df.index = original_index
        new_detail_df.to_csv(output_csv, index=False, encoding='utf-8-sig', header=True, sep=',')

    def analyze_ene_csv(self):
        if self.connent_group_num > 1:
            return_pages = []
            for i in range(self.connent_group_num):
                group_csv = os.path.join(self.__stat_path, f'{self.ene_csv.split(".")[0]}_group{i+1}.csv')
                return_pages.extend(WCC_API.analyze(group_csv, i+1)[1])
            return None, return_pages
        return WCC_API.analyze(self.ene_csv)

    def gen_pdf(self, output_pdf=None):
        if output_pdf is None:
            output_pdf = os.path.join(self.__stat_path, SUMMARY_PDF_FILENAME)

        output_tiny_pdf = output_pdf.replace('.pdf', '_abstract.pdf')

        # if os.path.exists(output_tiny_pdf):
        #     os.remove(output_tiny_pdf)
        # if os.path.exists(output_pdf):
        #     os.remove(output_pdf)
        specify_fonts = ['FreeMonoBold.ttf', 'NotoSans-VF.ttf', 'arial.ttf', 'DejaVuSansMono.ttf']
        specify_font = 'FreeMonoBold.ttf'
        for font in specify_fonts:
            try:
                ImageFont.truetype(font)
                specify_font = font
                break
            except OSError:
                pass

        pdf_format = PDFFormatter.get_default_format(specify_font)
        pdf_format.plot_title_font_args = {'name': specify_font, 'fontsize': 60, 'color': (0, 0, 0, 255)}
        pdf_format.plot_description_font_args = {'name': specify_font, 'fontsize': 40, 'color': (0, 0, 0, 255)}
        pdf_format.page_description_font_args = {'name': specify_font, 'fontsize': 80, 'color': (0, 0, 0, 255)}
        pdf_generator = PDFGenerator.setup_from_Formatter(pdf_format)
        tiny_pdf_generator = PDFGenerator.setup_from_Formatter(pdf_format)

        print('Generating PDF...', flush=True)

        reader = PdfReader(output_pdf)
        exist_pages_num = len(reader.pages)
        # print(pdf_format)

        page_queue = sorted(list(self.page_content_dict.keys()), reverse=False)

        for i, page in enumerate(page_queue):
            print(f"Generating detail PDF page ({i + 1}/{len(page_queue)}) Section [3/4]", flush=True)
            pair_page_content = self.page_content_dict[page]
            full_page = pdf_generator.fast_build_page(pair_page_content['content'],
                                                      title=pair_page_content['title'],
                                                      description=pair_page_content['description'],
                                                      resize=(1200, 1200),
                                                      restrict_page_size=False)
            pdf_generator.append_page_to_pdf(full_page,
                                             pdf_file=output_pdf,
                                             )
            del full_page
        del pdf_generator

        for i, page in enumerate(page_queue):
            print(f"Generating abstract PDF page ({i + 1}/{len(page_queue)}) Section [4/4]", flush=True)
            pair_page_content = self.page_content_dict[page]
            abstract_page = tiny_pdf_generator.fast_build_page([pair_page_content['content'][-1]],
                                                               title=pair_page_content['title'],
                                                               description=pair_page_content['description'],
                                                               resize=(1200, 1200),
                                                               restrict_page_size=False)
            tiny_pdf_generator.append_page_to_pdf(abstract_page,
                                                  pdf_file=output_tiny_pdf,
                                                  )
            del abstract_page
        del tiny_pdf_generator

        conv_data_tables = self.get_pdf_conv_tables()
        print_tables = [None for _ in range(exist_pages_num)]
        print_tables.extend(conv_data_tables)
        # print(print_tables, flush=True)

        PDFGenerator.add_header_and_pagenum2pdf(input_pdf=output_pdf, output_pdf=output_pdf,
                                                header_text='Alchemd')
        PDFGenerator.add_header_and_pagenum2pdf(input_pdf=output_tiny_pdf, output_pdf=output_tiny_pdf,
                                                header_text='Alchemd')
        PDFGenerator.add_tables2pdf(input_pdf=output_pdf, output_pdf=output_pdf, table_datas=print_tables,
                                    style=TABLE_STYLE,
                                    x=360, y=130, width=400, height=180,)
        # Avoid OOM.

    def clean_file(self):
        if os.path.exists(self.__stat_path):
            for file in os.listdir(self.__stat_path):
                if file.endswith('.pdf'):
                    os.remove(os.path.join(self.__stat_path, file))


class WCC_API:
    def __init__(self):
        pass

    @staticmethod
    def auto_wcc(csv) -> Tuple[pd.DataFrame, pd.DataFrame]:
        tmp_wcc_csv = os.path.join(os.path.dirname(csv), f'{os.path.basename(csv).replace(".csv", "")}_wcc_tmp')

        cmd = f'python {WCC_PATH} -f {csv} -p yes > {tmp_wcc_csv}'
        try:
            p = subprocess.Popen(cmd, shell=True)
            p.wait()
        except IndexError:
            print('Something wrong with given csv.')
        wcc_ene_list = []
        node_name_list = []
        node_ene_list = []
        with open(tmp_wcc_csv, 'r') as f:
            f.readline()
            f.readline()
            # Skip first 2 lines
            for line in f:
                line = line.strip()
                if line.startswith('***'):
                    break
                while line.find('  ') != -1:
                    line = line.replace('  ',' ')
                _ene = line.split(' ')[1]
                wcc_ene_list.append(_ene)

            # nodes after.
            f.readline()
            for line in f:
                line = line.strip()
                if len(line) == 0:
                    break
                while line.find('  ') != -1:
                    line = line.replace('  ',' ')
                _ene = line.split(' ')[1]
                _node = line.split(' ')[0]
                node_ene_list.append(_ene)
                node_name_list.append(str(_node))

        wcc_df = pd.DataFrame({'cal_wcc': wcc_ene_list})
        wcc_node_df = pd.DataFrame({'node': node_name_list, 'cal_wcc': node_ene_list})
        return wcc_df, wcc_node_df

    @staticmethod
    def standardize_csv(csv, output_delimiter='\t') -> str:
        """
        Output format: 23\t55\tene\n , use to input wcc
        :param csv:
        :param output_delimiter:
        :return:
        """
        possible_delimiter_list = [',', '|', '\t']
        csv_df = None
        for delimiter in possible_delimiter_list:
            tmp_df = pd.read_csv(csv, delimiter=delimiter)
            if len(tmp_df.columns) > 1:
                csv_df = tmp_df
        ene_column = csv_df.iloc[:, -1]
        pair_column = csv_df.iloc[:, 0]
        new_df = pd.concat([pair_column.str.split('-').str[0], pair_column.str.split('-').str[1],
                            ene_column], axis=1)
        new_csv = str(os.path.abspath(csv).replace('.csv', '_std.csv'))
        new_df.to_csv(new_csv, sep=output_delimiter, encoding='utf-8', index=False, header=False)
        return str(new_csv)

    @staticmethod
    def analyze(csv: str, group_idx: int = None) -> Tuple[str, list]:
        if group_idx is None:
            group_str = ''
        else:
            group_str = f'Group {group_idx}'
        csv = os.path.abspath(csv)
        tmp_dir = os.path.join(os.path.dirname(csv), 'tmp')
        os.makedirs(tmp_dir, exist_ok=True)
        shutil.copy(csv, tmp_dir)
        tmp_csv = os.path.join(tmp_dir, os.path.basename(csv))
        std_csv = WCC_API.standardize_csv(tmp_csv)
        analyze_df = pd.read_csv(std_csv, sep='\t', names=['node1', 'node2', 'cal'])
        wcc, wcc_node_pd = WCC_API.auto_wcc(std_csv)
        analyze_df = pd.concat([analyze_df, wcc], axis=1)
        pair_ana_df = copy.deepcopy(analyze_df)
        node_ana_df = copy.deepcopy(wcc_node_pd)

        analyze_df[''] = None
        analyze_df = pd.concat([analyze_df, wcc_node_pd], axis=1)
        output_columns = ['ligand1', 'ligand2', 'ΔΔG_ori (kcal/mol)',  'ΔΔG_corrected (kcal/mol)',
                          '   ', 'ligand', 'ΔG (kcal/mol)']
        analyze_df.columns = output_columns
        pair_ana_df.columns = output_columns[:4]
        node_ana_df.columns = output_columns[5:]

        # pdf_writer = SimpleDocTemplate(pdf_file, pagesize=letter)
        pair_ana_pages = df_to_pages(pair_ana_df, dec_format=(0, 0, 3, 3), title=f"{group_str} Pair Results")
        node_ana_pages = df_to_pages(node_ana_df, dec_format=(0, 3), title=f"{group_str} Ligand Results")
        all_ana_pages = [*pair_ana_pages, *node_ana_pages]
        # pdf_writer.build(all_ana_pages)

        output_csv = os.path.join(os.path.dirname(csv), csv.replace('.csv', '_analyze.csv'))
        analyze_df.to_csv(output_csv, sep=',', encoding='utf-8-sig', index=False)
        shutil.rmtree(tmp_dir)
        return output_csv, all_ana_pages


def CLIMain(**kwargs):
    ori_dir = os.getcwd()
    default_log_filename = 'run.log'
    parser = argparse.ArgumentParser()
    parser.add_argument('-C', '--calculate_path', type=str, help='Where the calculations are done.')
    parser.add_argument('-W', '--work_path', type=str, help='Where the works are done.')
    parser.add_argument('-L', '--log_name', type=str, help='Name of the log file.', default=default_log_filename)
    args = parser.parse_args()

    calculate_path = os.path.abspath(args.calculate_path)
    work_path = os.path.abspath(args.work_path)
    log_name = args.log_name

    pair_list = None
    get_folder_method = get_normal_processes_folders
    edge_stater = EdgeStater

    base_path = os.path.dirname(calculate_path)
    os.chdir(base_path)
    results_path = os.path.join(base_path, 'results')
    summary_pdf = os.path.join(results_path, SUMMARY_PDF_FILENAME)

    info_pickle_file = os.path.join(base_path, 'info.pickle')
    pdf_mol_content = []
    if os.path.exists(info_pickle_file):
        with open(info_pickle_file, 'rb') as f:
            infos = pickle.load(f)
        ligands = infos['all_moles']
        ligand_list = []
        for name, ligand in ligands.items():
            ligand.SetProp('_Name', name)
            ligand_list.append(ligand)
        ligand_list = sorted(ligand_list, key=lambda x: x.GetProp('_Name'), reverse=False)
        pairs = infos['pair_list']
        pdf_mol_content = moles_to_pages(ligand_list, title='Ligands')

    a = MainStater(calculate_path, work_path, pair_list=pair_list, edge_stater=edge_stater, _log_name=log_name,
                   results_path=results_path)
    a.stat()

    pair_num, ligand_num = a.get_pair_and_ligand_nums()
    cover_pdf_page = cover_page(pair_num=pair_num, ligand_num=ligand_num)
    info_pdf_page = static_info_page()
    _, pdf_ene_content = a.analyze_ene_csv()

    all_pdf_content = [*cover_pdf_page, *info_pdf_page, *pdf_mol_content, *pdf_ene_content]
    c = SimpleDocTemplate(summary_pdf, pagesize=letter)
    c.build(all_pdf_content)
    delete_pics()

    a.gen_detail_ene_csv()
    b = MainConvergenceMerger(calculate_path, pair_list=pair_list)
    b.merge(get_folder_method=get_folder_method, _output_path=results_path)
    a.append_sum_convergence_plots(if_edge_sum=False)
    a.gen_pdf(output_pdf=summary_pdf)

    os.chdir(ori_dir)


if __name__ == '__main__':
    CLIMain()
    pass
