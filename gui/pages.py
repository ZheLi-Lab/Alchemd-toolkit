import tkinter as tk
from tkinter import messagebox
import os, toml
from typing import Union, Tuple
from enum import Enum, unique
from typing import Any, Callable
from textwrap import dedent
from tkinter import *
from tkinter import ttk
from tkinter import filedialog

from dataclasses import dataclass

AVAIL_FONT = 'Arial'
BIG_FONT = (AVAIL_FONT, 13)  # For Titles
MIDDLE_FONT = (AVAIL_FONT, 10)  # For Widgets content
SMALL_FONT = (AVAIL_FONT, 8)  # For Descriptions

FRAME_FORMAT = {'sticky': 'ew', 'padx': 5, 'ipady': 3}


def fontsize_decide(font: tuple) -> tuple:
    return font, (font[0], int(font[1]-2)), (font[0], int(font[1]-3))


@unique
class ArgType(Enum):
    File = 'file'
    Directory = 'dir'
    Bool = 'bool'
    DoubleBool = 'doublebool'
    Int = 'int'
    Str = 'str'
    ComboBox = 'combobox'


@dataclass
class SubmitInfo:
    state: bool  # Ready to run.
    work_dir: Union[str, None]  # The work path, the directory to be modded
    log_dir: str = None  # Where the log should store
    begin_text: str = None  # The text show before run
    end_text: str = None  # The text show after run
    page_type: str = 'info'  # 'info' or 'work'
    run_cmd: str = None  # The command to run
    extra_info: str = None  # The info given in specify situation

    def __str__(self):
        return dedent(f'''
        Ready to submit: {self.state}
        Working directory: {self.work_dir}
        run_cmd: {self.run_cmd}
        extra_info: {self.extra_info}
        ''')


@dataclass
class arg_info:
    arg_flag: Union[str, None]
    arg_type: ArgType
    arg_desc: str  # The title of arg in GUI
    arg_name: str = None  # The arg name using in this script
    arg_note: str = None  # The info of arg in GUI
    arg_value: Any = None
    is_explicit: bool = True
    widget: tk.Widget = None
    widget_args: dict = None
    must_fill: bool = False

    def set_value(self, arg_value):
        self.arg_value = arg_value
        if isinstance(self.widget, PathWidget):
            self.widget.update_text(self.arg_value)
        if isinstance(self.widget, ComboBoxWidget):
            self.widget.set_value(self.arg_value)

    def get_value(self):
        if self.arg_type is ArgType.File or self.arg_type is ArgType.Directory:
            return self.arg_value
        return self.arg_value

    def is_empty(self) -> bool:
        if self.arg_value is None:
            return True
        if self.arg_type is ArgType.File or self.arg_type is ArgType.Directory:
            if self.arg_value == '':
                return True
        return False


class ScrollableFrame(tk.Frame):
    def __init__(self, root, *args, **kwargs):
        super().__init__(root, *args, **kwargs)
        self.configure()

        self.allow_scroll = False

        self.canvas = tk.Canvas(self, relief='flat', bd=0)
        self.scrollbar = tk.Scrollbar(self, orient="vertical", command=self.canvas.yview)

        self.scrollable_frame = tk.Frame(self.canvas, relief='flat', padx=0)

        self.canvas_frame = self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        self.scrollable_frame.bind("<Configure>", self.on_frame_configure)
        self.canvas.bind("<Configure>", self.on_canvas_configure)

        self.canvas.grid(row=0, column=0, sticky="nsew")
        self.scrollbar.grid(row=0, column=1, sticky="ns")

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.bind_mouse_wheel()

        self.scrollable_frame.grid_columnconfigure(0, weight=1)

    def on_frame_configure(self, event):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        self.update_scrollbar_visibility()

    def on_canvas_configure(self, event):
        self.canvas.itemconfig(self.canvas_frame, width=event.width)
        self.update_scrollbar_visibility()

    def update_scrollbar_visibility(self):
        canvas_height = self.canvas.winfo_height()
        frame_height = self.scrollable_frame.winfo_reqheight()

        if frame_height > canvas_height:
            self.allow_scroll = True
            self.scrollbar.grid()
        else:
            self.allow_scroll = False
            self.scrollbar.grid_remove()

        self.bind_mouse_wheel()

    def switch_update(self):
        self.update_scrollbar_visibility()
        self.bind_mouse_wheel()

    def bind_mouse_wheel(self):
        def on_mouse_wheel(event):
            self.canvas.yview_scroll(-1 * (event.delta // 120), "units")

        def nothing(event):
            pass

        if self.allow_scroll:
            self.canvas.bind_all("<MouseWheel>", on_mouse_wheel)
        else:
            self.canvas.bind_all("<MouseWheel>", nothing)


class LockFrame(tk.LabelFrame):
    def __init__(self, root, button_text, default_toggle: bool = False, toggle_is_lock: bool = False, mode='normal',
                 text2=None,
                 font=BIG_FONT, *args, **kwargs):
        title_font, widget_font, description_font = fontsize_decide(font)
        super().__init__(root, font=title_font, *args, **kwargs)

        self.is_toggled = False
        self.toggle_is_lock = toggle_is_lock
        self.is_locked = False
        self.is_locked = self.lock_update()
        self.all_exclusive_frames = []

        self.toggle_var = tk.BooleanVar()
        self.toggle_var.trace("w", self.toggle_switch)

        # self.lock_button = tk.Checkbutton(self, text=button_text, variable=self.toggle_var, font=widget_font)
        if mode == 'normal':
            self.lock_button = tk.Checkbutton(self, text=button_text, variable=self.toggle_var, font=widget_font)
        else:
            self.lock_button = DoubleBoolWidget(self, text1=button_text, text2=text2, report_method=None,
                                                same_line=False, tk_var=self.toggle_var)
        self.lock_button.grid(row=0, column=0, sticky='wn', padx=5)

        self.grid_configure(sticky='ew', padx=5, ipady=3)

        if default_toggle:
            self.toggle_var.set(True)

    def lock_update(self) -> bool:
        if self.toggle_is_lock == self.is_toggled:
            return True
        return False

    def toggle_switch(self, *args):
        self.is_toggled = self.toggle_var.get()
        self.is_locked = self.lock_update()
        if not self.is_locked:
            for frame in self.all_exclusive_frames:
                if frame != self:
                    frame.toggle_var.set(False)
        else:
            for frame in self.all_exclusive_frames:
                if frame != self:
                    frame.toggle_var.set(True)
                    break
        try:
            self.update_widgets_state()
        except:
            pass

    def force_update_widgets_state(self, state=tk.DISABLED or tk.NORMAL):
        outer_lock_state = True if state == tk.DISABLED else False
        if not (outer_lock_state or self.is_locked):
            inner_state = tk.NORMAL
        else:
            inner_state = tk.DISABLED

        for widget in self.winfo_children():
            if widget != self.lock_button:
                try:
                    if isinstance(widget, PathWidget):
                        widget.switch_toggle_without_button(inner_state)
                    elif isinstance(widget, LockFrame):
                        widget.force_update_widgets_state(inner_state)
                    elif isinstance(widget, BoolWidget):
                        widget.switch_lock(inner_state)
                    else:
                        widget.config(state=inner_state)
                except:
                    pass
            else:
                if isinstance(widget, DoubleBoolWidget):
                    widget.switch_lock(state)
                else:
                    widget.config(state=state)

    def update_widgets_state(self):
        for widget in self.winfo_children():
            if widget != self.lock_button:
                now_state = tk.DISABLED if self.is_locked else tk.NORMAL
                try:
                    if isinstance(widget, PathWidget):
                        widget.switch_toggle_without_button(now_state)
                    elif isinstance(widget, LockFrame):
                        widget.force_update_widgets_state(now_state)
                    elif isinstance(widget, BoolWidget):
                        widget.switch_lock(now_state)
                    else:
                        widget.config(now_state)
                except:
                    pass


class PathWidget(tk.LabelFrame):
    def __init__(self, root, text, font=BIG_FONT, note=None,
                 lock_button=False, lock_button_text='lock', lock_status_default=True,
                 report_method=None,
                 ask_method: Callable = tk.filedialog.askdirectory,
                 title: Union[str, None] = None,):

        title_font, widget_font, description_font = fontsize_decide(font)

        super().__init__(root)
        if lock_button or title is not None:
            using_relief = tk.GROOVE
        else:
            using_relief = tk.FLAT
        self.config(relief=using_relief, bd=2, text=title, font=title_font)
        self.root = root

        self.path = tk.StringVar()
        self.lock_state = tk.IntVar()
        self.report_method = report_method
        self.text = text

        if lock_button:
            self.lock_radio_button = tk.Checkbutton(self, text=lock_button_text, variable=self.lock_state,
                                                    command=self.switch_lock_status, font=widget_font)
            self.lock_radio_button.grid(row=0, column=0, sticky='sw', pady=0, columnspan=3)

        self.label = tk.Label(self, text=f"{text}:", font=widget_font)
        self.label.grid(row=1, column=0, sticky='w')

        self.text_box = tk.Entry(self, font=widget_font, textvariable=self.path)
        self.path.trace('w', self.update_input)
        self.text_box.grid(row=1, column=1, padx=5, sticky='ew')

        self.select_button = tk.Button(self, text='...', font=(AVAIL_FONT, int(font[1]*0.8), 'bold'),
                                       command=self.select_path, height=1)
        self.select_button.grid(row=1, column=2, sticky='w')

        if note is not None:
            tk.Label(self, text=note, font=description_font,
                     relief='groove', justify='left').grid(row=2, column=0, columnspan=3, sticky='nw')

        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(1, pad=5)
        # self.grid_configure(ipady=5)

        self.ask_method = ask_method

        if lock_status_default and lock_button:
            self.lock_radio_button.toggle()
            self.switch_lock_status()

    def switch_lock_status(self):
        v = self.lock_state.get()
        self.switch_toggle_without_button(v)

    def switch_toggle_without_button(self, v):
        if v == 1 or v == tk.DISABLED:
            # being locked
            self.text_box.config(state=tk.DISABLED)
            self.select_button.config(state=tk.DISABLED)
        else:
            self.text_box.config(state=tk.NORMAL)
            self.select_button.config(state=tk.NORMAL)

    def select_path(self):
        path_ = self.ask_method(title=self.text)
        if path_ != '':
            self.update_text(path_)

    def report(self):
        self.report_method(self.text_box.get())

    def update_text(self, text):
        self.path.set(text)
        self.text_box.delete(0, tk.END)
        self.text_box.insert(tk.END, text)

    def get_text(self):
        if self.path.get() != self.text_box.get():
            return self.text_box.get()
        return self.path.get()

    def update_input(self, *args, **kwargs):
        self.update_text(self.text_box.get())
        try:
            self.report()
        except:
            pass


class BoolWidget(tk.Frame):
    def __init__(self, root, report_method, if_toggle=False, note=None, font=BIG_FONT,
                 *args, **kwargs):
        self.report_value = tk.BooleanVar()
        title_font, widget_font, description_font = fontsize_decide(font)

        super().__init__(root)
        self.checkbutton = tk.Checkbutton(self, variable=self.report_value, font=widget_font, *args, **kwargs)
        self.checkbutton.grid(row=0, column=0, sticky='w', pady=0)
        if note is not None:
            tk.Label(self, text=note, font=description_font,
                     relief='groove').grid(row=1, column=0, columnspan=3, sticky='nw')

        self.report_method = report_method
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, pad=5)

        if if_toggle:
            self.checkbutton.toggle()

    def report(self):
        self.report_method(self.report_value.get())

    def switch_lock(self, state=tk.DISABLED or tk.NORMAL):
        self.checkbutton.config(state=state)


class DoubleBoolWidget(tk.Frame):
    def __init__(self, root, report_method, text1, text2, font=BIG_FONT, same_line: bool = True,
                 tk_var: tk.Variable = None,
                 *args, **kwargs):
        super().__init__(root, *args, **kwargs)
        if tk_var is None:
            self.report_value = tk.BooleanVar()
        else:
            self.report_value = tk_var

        self.other_var = tk.BooleanVar()

        title_font, widget_font, description_font = fontsize_decide(font)
        self.checkbutton1 = tk.Checkbutton(self, variable=self.report_value, font=widget_font,
                                           text=text1, command=self.checkbutton1_click)
        self.checkbutton2 = tk.Checkbutton(self, variable=self.other_var, font=widget_font, text=text2,
                                           command=self.checkbutton2_click)

        if same_line:
            self.checkbutton1.grid(row=0, column=0, sticky='w')
            tk.Label(self, text=" Or ", font=widget_font, padx=5).grid(row=0, column=1, sticky='w')
            self.checkbutton2.grid(row=0, column=2, sticky='w')
        else:
            self.checkbutton1.grid(row=0, column=0, sticky='w', padx=5)
            tmp_frame = tk.Frame(self)
            tk.Label(tmp_frame, text=" Or ", font=widget_font, justify='left').grid(row=0, column=0, sticky='w')
            self.checkbutton2 = tk.Checkbutton(tmp_frame, variable=self.other_var, font=widget_font, text=text2,
                                               command=self.checkbutton2_click)
            self.checkbutton2.grid(row=0, column=1, sticky='w')
            tmp_frame.grid(row=1, column=0, sticky='ew')

        self.report_method = report_method
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, pad=5)

        self.report_value.set(True)

    def checkbutton1_click(self, *args):
        if self.report_value.get():
            self.other_var.set(False)
        else:
            self.other_var.set(True)

    def checkbutton2_click(self, *args):
        if self.other_var.get():
            self.report_value.set(False)
        else:
            self.report_value.set(True)

    def report(self):
        self.report_method(self.report_value.get())

    def switch_lock(self, state=tk.DISABLED or tk.NORMAL):
        self.checkbutton1.config(state=state)
        self.checkbutton2.config(state=state)


class ComboBoxWidget(tk.Frame):
    def __init__(self, root, report_method, note=None, title: str = None, default: int = 0, font: tuple = BIG_FONT,
                 *args, **kwargs):

        title_font, widget_font, description_font = fontsize_decide(font)

        super().__init__(root)
        self.root = root
        self.report_value = tk.StringVar()
        self.report_method = report_method

        self.label = tk.Label(self, text=f'{title}:', justify='left', font=widget_font)
        self.label.grid(row=0, column=0, sticky='w', pady=5)

        self.ComboBox = ttk.Combobox(self, textvariable=self.report_value, *args, **kwargs, font=widget_font)
        self.report_value.trace('w', self.report)
        self.ComboBox.grid(row=0, column=1, sticky='ew', padx=5)
        self.ComboBox.current(default)
        self.ComboBox['state'] = 'readonly'

        if note is not None:
            tk.Label(self, text=note, font=description_font,
                     relief='groove').grid(row=1, column=0, columnspan=3, sticky='nw')

        self.grid_columnconfigure(1, weight=1)

        self.ComboBox.bind('<<ComboboxSelected>>', self.prevent_selection)
        self.ComboBox.bind('<FocusIn>', self.prevent_selection)

    def prevent_selection(self, event):
        # 移除焦点
        self.root.focus_set()
        self.ComboBox.selection_clear()

    def report(self, *args):
        self.report_method(self.report_value.get())

    def switch_lock(self, state=tk.DISABLED or tk.NORMAL):
        self.ComboBox.config(state=state)

    def set_value(self, value):
        if value not in self.ComboBox["values"]:
            return False
        idx = list(self.ComboBox["values"]).index(value)
        self.ComboBox.current(idx)


class IntWidget(tk.Frame):
    def __init__(self, root, report_method, desc, note=None, default=None, font=BIG_FONT):
        self.report_value = tk.IntVar()
        title_font, widget_font, description_font = fontsize_decide(font)

        super().__init__(root)

        self.label = tk.Label(self, text=f"{desc}:", justify='left', font=widget_font)
        self.label.grid(row=0, column=0, sticky='w', pady=5)

        self.entry_box = tk.Entry(self, width=8, validate="key", font=widget_font,
                                  validatecommand=(root.register(self.validate_int), "%P"))
        self.entry_box.grid(row=0, column=1, sticky='ws', padx=5, pady=5)
        if default is not None:
            self.entry_box.insert(0, str(default))

        if note is not None:
            tk.Label(self, text=note, font=widget_font,
                     relief='groove').grid(row=1, column=0, columnspan=3, sticky='nw')

        self.report_method = report_method

    @staticmethod
    def validate_int(value):
        try:
            if len(value) > 0:
                int(value)
            return True
        except ValueError:
            return False

    def report(self):
        self.report_method(int(self.entry_box.get()))


class ABC_ConfigPage(ScrollableFrame):
    def __init__(self, root, script):
        super().__init__(root)
        self.script = script
        self.script_base = ''

        self.submit_info = SubmitInfo(False, None)
        self.submit_info.page_type = 'work'
        self.submit_info.begin_text = 'The script begins.'
        self.submit_info.state = False

        self.not_filled_arg_list = []
        self.active_arg_dict = {}

        self.toml_dict = {}
        self.toml_file = 'job.toml'

        self.grid_columnconfigure(0, weight=1)
        self.scrollable_frame.grid_rowconfigure(0, weight=1)

    def _check_arg(self, arg: arg_info):
        arg.widget.report()
        if arg.is_empty() and arg.must_fill:
            self.not_filled_arg_list.append(arg.arg_desc)
        else:
            if arg.arg_type is ArgType.Bool or arg.arg_type is ArgType.DoubleBool:
                if arg.arg_value:
                    self.active_arg_dict[arg.arg_flag] = None
            elif arg.arg_type is ArgType.File or arg.arg_type is ArgType.Directory:
                path = arg.get_value()
                if os.path.exists(path):
                    if arg.arg_type is ArgType.File and os.path.isfile(path):
                        self.active_arg_dict[arg.arg_flag] = path
                    elif arg.arg_type is ArgType.Directory and os.path.isdir(path):
                        self.active_arg_dict[arg.arg_flag] = path
                else:
                    self.not_filled_arg_list.append(arg.arg_desc)
            else:
                self.active_arg_dict[arg.arg_flag] = arg.get_value()

    def _read_toml(self, toml_file):
        if os.path.exists(toml_file) and os.path.isfile(toml_file):
            with open(toml_file, 'r') as f:
                self.toml_dict = toml.load(f)

    def sync_from_toml(self, arg: arg_info, toml_key: str):
        if toml_key in self.toml_dict.keys():
            arg.set_value(self.toml_dict[toml_key])

    def sync_to_toml(self, arg: arg_info, toml_key: str):
        if not arg.is_empty():
            self.toml_dict[toml_key] = arg.get_value()

    def _dump_toml(self, toml_file):
        with open(toml_file, 'w') as f:
            toml.dump(self.toml_dict, f)

    def report_args(self) -> Tuple[list, dict]:
        """
        Report args according to the gui page.
        :return: the list of args active but not filled, the dict of active and filled args.
        """
        return_tuple = (self.not_filled_arg_list, self.active_arg_dict)
        raise NotImplementedError

    def sync_args(self, sync_method: str = 'r', dump: bool = False, from_dict: dict = None):
        """
        Sync args according to the toml file.
        Use this function to sync some autofill args between pages
        - sync_method: read ('r') or write ('w')
        - dump: whether to dump to a toml file, default is False
        - from_dict: read from a dict instead of reading from toml file
        :return:
        """
        if self.toml_file is None:
            return None

        if sync_method == "r":
            if from_dict is None:
                self._read_toml(self.toml_file)
            else:
                self.toml_dict = from_dict
            self._sync_args_read()
        elif sync_method == "w":
            self._sync_args_write()
            if dump:
                self._dump_toml(self.toml_file)
        else:
            raise ValueError(f"Unknown sync method: {sync_method}")

    def _sync_args_read(self):
        # Do something after read
        pass

    def _sync_args_write(self):
        # Do something before write
        pass

    def report_submit_info(self) -> SubmitInfo:
        """
        Report submit_info according to the args.
        The check of argument filling state is done by this method, the child class must implement this method further.
        :return:
        """
        self.report_args()

        if len(self.not_filled_arg_list) > 0:
            # Not filled all active args.
            self.submit_info.state = False
            self.submit_info.extra_info = self.not_filled_arg_list
            return self.submit_info

        # All active and must_filled args filled.
        self.submit_info.state = True
        self.submit_info.extra_info = None

    def gen_cmd(self) -> SubmitInfo:
        """
        Generate a submit command according to the args
        :return:
        """
        py_script = os.path.join(self.script_base, 'core', self.script)
        run_cmd = f"python {py_script} "
        for key, value in self.active_arg_dict.items():
            if value is None:
                run_cmd += f"{key} "
            else:
                run_cmd += f"{key} {value} "
        self.submit_info.run_cmd = run_cmd

        self.sync_args('w', dump=True)

        return self.submit_info


def quick_create_arg(_arg_info: arg_info, _root):
    tmp_arg = _arg_info
    if tmp_arg.arg_type is ArgType.File or tmp_arg.arg_type is ArgType.Directory:
        ask_method = tk.filedialog.askopenfilename if tmp_arg.arg_type is ArgType.File else tk.filedialog.askdirectory
        tmp_arg.widget = PathWidget(_root,
                                    text=tmp_arg.arg_desc,
                                    note=tmp_arg.arg_note,
                                    report_method=tmp_arg.set_value,
                                    ask_method=ask_method)
    elif tmp_arg.arg_type is ArgType.Bool:
        if_toggle = False
        if tmp_arg.widget_args is not None and 'if_toggle' in tmp_arg.widget_args.keys():
            if_toggle = tmp_arg.widget_args['if_toggle']
        tmp_arg.widget = BoolWidget(_root, report_method=tmp_arg.set_value, if_toggle=if_toggle,
                                    text=tmp_arg.arg_desc, note=tmp_arg.arg_note)
    elif tmp_arg.arg_type is ArgType.ComboBox:
        tmp_arg.widget = ComboBoxWidget(_root, report_method=tmp_arg.set_value, title=tmp_arg.arg_desc,
                                        note=tmp_arg.arg_note, **tmp_arg.widget_args)
    elif tmp_arg.arg_type is ArgType.DoubleBool:
        tmp_arg.widget = DoubleBoolWidget(_root, report_method=tmp_arg.set_value, text1=tmp_arg.arg_desc,
                                          **tmp_arg.widget_args)

    return tmp_arg


class PreparePage(ABC_ConfigPage):
    def __init__(self, root):
        super().__init__(root, 'prepare_file.py')

        self._prepare_mode = 'ligand'
        self.clean_last_toml_record()

        arg_dict = {'ProteinFile': arg_info(arg_flag='-p', arg_type=ArgType.File, arg_desc='Protein file',
                                            arg_note='Supported format: pdb  ',
                                            must_fill=True),
                    'LigandFile': arg_info(arg_flag='-l', arg_type=ArgType.File, arg_desc='Ligands file',
                                           arg_note='Supported format: sdf  ',
                                           must_fill=True
                                           ),
                    'RefLigandFile': arg_info(arg_flag='-r', arg_type=ArgType.File, arg_desc='Reference ligand list',
                                              # arg_note='Customize reference ligands  ',
                                              must_fill=True
                                              ),
                    'PairFile': arg_info(arg_flag='-pl', arg_type=ArgType.File, arg_desc='Pair list',
                                         # arg_note='Customize pair list  ',
                                         must_fill=True
                                         ),
                    'AdvanceFolder': arg_info(arg_flag='-af', arg_type=ArgType.Directory, arg_desc='Designated directory',
                                              # arg_note='Path to a designated directory  ',
                                              must_fill=True),
                    'OutputDirectory': arg_info(arg_flag='-o', arg_type=ArgType.Directory, arg_desc='Output directory',
                                                arg_name='output_dir',
                                                arg_note='Specify the root directory for output',
                                                must_fill=True),
                    # 'StrictPairList': arg_info(arg_flag='-sp', arg_type=ArgType.Bool, arg_desc='Strict pair list',
                    #                            arg_note='Enforce strict adherence to user-defined pair list  ',
                    #                            widget_args={'if_toggle': False}),
                    'NSAAParamsDir': arg_info(arg_flag='-npd', arg_type=ArgType.Directory,
                                                     arg_desc='NSAA parameters directory',
                                                     arg_note='Directory containing parameters for non-standard amino acids',
                                                     must_fill=True,
                                                     )
                    # TODO: fix this.
                    }

        protein_select_frame = tk.LabelFrame(self.scrollable_frame, text='Protein', relief='groove', bd=4, font=BIG_FONT)
        self.ArgProteinFile = quick_create_arg(arg_dict['ProteinFile'], protein_select_frame)
        self.ArgProteinFile.widget.grid(row=0, column=0, sticky='ew', padx=5)
        self.LockFrameNSAAParams = LockFrame(protein_select_frame,
                                                       button_text='Parameters for non-standard residues', relief=tk.FLAT,
                                                       default_toggle=False, toggle_is_lock=False,
                                                       )
        self.ArgNSAAParams = quick_create_arg(arg_dict['NSAAParamsDir'],
                                                        self.LockFrameNSAAParams)
        self.ArgNSAAParams.widget.grid(row=1, column=0, sticky='ew', padx=5)
        self.LockFrameNSAAParams.grid(row=1, column=0, sticky='ew')
        protein_select_frame.grid(row=0, column=0, sticky='ew', padx=5, pady=5)
        protein_select_frame.columnconfigure(0, weight=1)

        ligand_label_frame = LabelFrame(self.scrollable_frame, text='Ligands', relief=tk.GROOVE, bd=4, font=BIG_FONT)
        self.LockFrameLigand = LockFrame(ligand_label_frame,
                                         button_text='Automatic preparation of ligand pairs',
                                         default_toggle=True,
                                         text='Automatic ligand preparation', relief='groove', bd=3)
        self.ArgLigandFile = quick_create_arg(arg_dict['LigandFile'], self.LockFrameLigand)
        self.ArgLigandFile.widget.grid(row=1, column=0, sticky='ew', padx=5)

        self.LockFrameRefLigandList = LockFrame(self.LockFrameLigand, button_text='Consider all ligands as reference',
                                                default_toggle=True, toggle_is_lock=True,
                                                text='Reference ligands', relief='groove',
                                                mode='double_options', text2='Customize reference ligands'
                                                )
        self.ArgRefLigandFile = quick_create_arg(arg_dict['RefLigandFile'], self.LockFrameRefLigandList)

        self.LockFramePairList = LockFrame(self.LockFrameLigand, button_text='Auto generate pair list',
                                           default_toggle=True, toggle_is_lock=True,
                                           text='Pair list', relief='groove',
                                           mode='double_options', text2='Customize pairs'
                                           )
        self.ArgPairFile = quick_create_arg(arg_dict['PairFile'], self.LockFramePairList)
        # self.ArgStrictPairList = quick_create_arg(arg_dict['StrictPairList'], self.LockFramePairList)

        self.ArgPairFile.widget.grid(row=1, column=0, sticky='ew', padx=5)
        # self.ArgStrictPairList.widget.grid(row=2, column=0, sticky='ew', padx=5)
        self.ArgRefLigandFile.widget.grid(row=1, column=0, sticky='ew', padx=5)

        self.LockFrameRefLigandList.grid(row=2, column=0, sticky='ew', padx=5, ipady=3)
        self.LockFramePairList.grid(row=3, column=0, ipady=3)
        self.LockFrameLigand.grid(row=0, column=0, sticky='ew', padx=5, pady=5)
        ligand_label_frame.grid(row=1, column=0, sticky='ew', padx=5, ipady=5)

        self.LockFrameAdvanceFolder = LockFrame(ligand_label_frame, text='Manual ligand preparation',
                                                default_toggle=False,
                                                button_text='Import manually prepared ligand pairs',
                                                bd=3)
        self.ArgAdvanceFolder = quick_create_arg(arg_dict['AdvanceFolder'], self.LockFrameAdvanceFolder)
        self.ArgAdvanceFolder.widget.grid(row=1, column=0, sticky='ew', padx=5)
        self.LockFrameAdvanceFolder.grid(row=1, column=0, sticky='ew', padx=5, ipady=3)

        self.ArgOutputDirectory = quick_create_arg(arg_dict['OutputDirectory'], self.scrollable_frame)
        self.ArgOutputDirectory.widget.grid(row=3, column=0, sticky='ew', padx=5)

        self.LockFrameAdvanceFolder.update_widgets_state()
        self.LockFrameLigand.update_widgets_state()
        self.LockFrameRefLigandList.update_widgets_state()
        self.LockFramePairList.update_widgets_state()
        self.LockFrameNSAAParams.update_widgets_state()
        self.LockFramePairList.columnconfigure(0, weight=1)
        ligand_label_frame.columnconfigure(0, weight=1)
        self.LockFrameAdvanceFolder.columnconfigure(0, weight=1)
        self.LockFrameLigand.columnconfigure(0, weight=1)
        self.LockFrameRefLigandList.columnconfigure(0, weight=1)
        self.LockFrameNSAAParams.columnconfigure(0, weight=1)

        all_exclusive_frames = [self.LockFrameLigand,
                                self.LockFrameAdvanceFolder,
                                ]
        self.LockFrameLigand.all_exclusive_frames = all_exclusive_frames
        self.LockFrameAdvanceFolder.all_exclusive_frames = all_exclusive_frames

        # self.sync_args('r')

    def report_args(self) -> Tuple[list, dict]:
        self.not_filled_arg_list = []
        self.active_arg_dict = {'-auto': None}

        self._check_arg(self.ArgProteinFile)
        if not self.LockFrameNSAAParams.is_locked:
            self._check_arg(self.ArgNSAAParams)

        if not self.LockFrameLigand.is_locked:
            # Prepare from ligand
            self._prepare_mode = 'ligand'
            self.active_arg_dict['--cs-align -d'] = None

            self._check_arg(self.ArgLigandFile)
            if not self.LockFrameRefLigandList.is_locked:
                self._check_arg(self.ArgRefLigandFile)
            if not self.LockFramePairList.is_locked:
                self.active_arg_dict['-sp'] = None
                self._check_arg(self.ArgPairFile)
                # self._check_arg(self.ArgStrictPairList)
        else:
            # Prepare from folder
            self._prepare_mode = 'folder'

            self._check_arg(self.ArgAdvanceFolder)

        self._check_arg(self.ArgOutputDirectory)

        return self.not_filled_arg_list, self.active_arg_dict

    def report_submit_info(self) -> SubmitInfo:
        super().report_submit_info()

        if not self.submit_info.state:
            # No completely filled.
            return self.submit_info

        if self._prepare_mode == 'ligand':
            self.submit_info.work_dir = self.ArgOutputDirectory.get_value()
        else:
            self.submit_info.work_dir = os.path.dirname(self.ArgAdvanceFolder.get_value())

        self.submit_info.log_dir = self.submit_info.work_dir
        self.submit_info.end_text = f"Collect results in {self.submit_info.work_dir}\n"
        self.submit_info.extra_info = f"Prepare in mode {self._prepare_mode}\n"

        self.gen_cmd()

        return self.submit_info

    def _sync_args_write(self):
        self.toml_dict['PrepareMode'] = self._prepare_mode
        if self.ArgOutputDirectory.is_empty():
            return False

        if self._prepare_mode == 'ligand':
            work_base_dir = self.ArgOutputDirectory.get_value()
        else:
            work_base_dir = os.path.dirname(self.ArgAdvanceFolder.get_value())
        self.toml_dict['PrepareDir'] = os.path.join(work_base_dir, 'prepare')
        self.toml_dict['MD_Dir'] = os.path.join(work_base_dir, 'run')

    def clean_last_toml_record(self):
        open(self.toml_file, 'w').close()


class RunPage(ABC_ConfigPage):
    def __init__(self, root):
        super().__init__(root, 'run_md.py')

        arg_dict = {'PrepareDirectory': arg_info(arg_flag='-p', arg_type=ArgType.Directory, arg_desc='Prepared directory',
                                                 arg_name='prepare_dir', must_fill=True,
                                                 arg_note='Path to molecule preparation directory '),
                    'MD_Program': arg_info(arg_flag='-m', arg_type=ArgType.ComboBox, arg_desc='MD simulation program',
                                           arg_name='md_program', widget_args={'values': ('Amber',
                                                                                          # 'Openmm'
                                                                                          ),
                                                                               }
                                           ),
                    'LocalRun': arg_info(arg_flag='-auto', arg_type=ArgType.DoubleBool,
                                         arg_desc='Run simulation locally',
                                         arg_name='local_run', widget_args={'text2': "Create a job list only",
                                                                            'same_line': False,
                                                                            }
                                         ),
                    }

        prepare_dir_frame = LabelFrame(self.scrollable_frame, text='MD Input', relief='groove', bd=4,
                                       font=BIG_FONT)
        self.ArgPrepareDirectory = quick_create_arg(arg_dict['PrepareDirectory'], prepare_dir_frame)
        self.ArgPrepareDirectory.widget.grid(row=0, column=0, sticky='ew', padx=5)
        prepare_dir_frame.grid(row=0, column=0, pady=5, **FRAME_FORMAT)
        prepare_dir_frame.columnconfigure(0, weight=1)

        run_frame = LabelFrame(self.scrollable_frame, text='Run with CAR', relief='groove', bd=4, font=BIG_FONT)
        self.ArgMD_Program = quick_create_arg(arg_dict['MD_Program'], run_frame)
        self.ArgLocalRun = quick_create_arg(arg_dict['LocalRun'], run_frame)
        self.ArgMD_Program.widget.grid(row=0, column=0, sticky='ew', padx=5)
        self.ArgLocalRun.widget.grid(row=1, column=0, sticky='ew')
        run_frame.grid(row=1, column=0, **FRAME_FORMAT)
        run_frame.columnconfigure(0, weight=1)

    def report_args(self) -> Tuple[list, dict]:
        self.not_filled_arg_list = []
        self.active_arg_dict = {'-ln': 50,}

        self._check_arg(self.ArgPrepareDirectory)
        self._check_arg(self.ArgMD_Program)
        self._check_arg(self.ArgLocalRun)

        return self.not_filled_arg_list, self.active_arg_dict

    def report_submit_info(self) -> SubmitInfo:
        super().report_submit_info()

        if not self.submit_info.state:
            # Not completely filled
            return self.submit_info

        self.submit_info.work_dir = os.path.join(os.path.dirname(self.ArgPrepareDirectory.get_value()), 'run')
        self.submit_info.end_text = f"Collect results in {self.submit_info.work_dir}\n"
        self.submit_info.log_dir = os.path.dirname(self.ArgPrepareDirectory.get_value())

        self.gen_cmd()
        return self.submit_info

    def _sync_args_read(self):
        self.sync_from_toml(self.ArgPrepareDirectory, 'PrepareDir')

    def _sync_args_write(self):
        self.sync_to_toml(self.ArgMD_Program, 'MD_Program')
        if self.ArgPrepareDirectory.is_empty():
            return False
        self.toml_dict['MD_Dir'] = os.path.join(os.path.dirname(self.ArgPrepareDirectory.get_value()), 'run')


class AnalyzePage(ABC_ConfigPage):
    def __init__(self, root):
        super().__init__(root, 'analyze_result.py')

        self.result_dir = None
        self.EnergyRawFile = None
        self.EnergyWCCFile = None
        self.ConvergenceFile = None

        title_font, widget_font, description_font = fontsize_decide(BIG_FONT)

        arg_dict = {'MD_Directory': arg_info(arg_flag='-w', arg_type=ArgType.Directory, arg_desc='MD directory',
                                             arg_name='MD_dir', must_fill=True,
                                             arg_note='Path to MD directory '),
                    'MD_Program': arg_info(arg_flag='-m', arg_type=ArgType.ComboBox, arg_desc='MD simulation program',
                                           arg_name='md_program', widget_args={'values': ('Amber',
                                                                                          # 'Openmm'
                                                                                          ),
                                                                               }
                                           ),
                    }

        inputFrame = LabelFrame(self.scrollable_frame, text="Analyze Input", relief='groove', bd=4, font=BIG_FONT)
        self.ArgMD_Directory = quick_create_arg(arg_dict['MD_Directory'], inputFrame)
        self.ArgMD_Program = quick_create_arg(arg_dict['MD_Program'], inputFrame)
        self.ArgMD_Directory.widget.grid(row=0, column=0, sticky='ew', padx=5)
        self.ArgMD_Program.widget.grid(row=1, column=0, sticky='ew', padx=5)
        inputFrame.grid(row=0, column=0, pady=5, **FRAME_FORMAT)
        inputFrame.columnconfigure(0, weight=1)

        def quick_output_path_widget(row, title) -> tk.Text:
            _label = tk.Label(outputFrame, text=title, font=widget_font)
            _widget = tk.Text(outputFrame, height=1, state=tk.DISABLED, font=widget_font, padx=5, wrap=tk.NONE)
            _place_holder = tk.Label(outputFrame, text=" ", font=widget_font)
            _label.grid(row=row, column=0, sticky='w', padx=5)
            _widget.grid(row=row, column=1, sticky='ew', padx=5, pady=5)
            _place_holder.grid(row=row, column=2, sticky='e')
            return _widget

        outputFrame = LabelFrame(self.scrollable_frame, text="Analyze Results", relief='groove', bd=4, font=BIG_FONT)

        self.raw_ene_path_widget = quick_output_path_widget(0, 'Results file:')
        self.ana_ene_path_widget = quick_output_path_widget(1, "Analyzed results file:")
        self.conv_path_widget = quick_output_path_widget(2, 'Summary file:')
        self.folder_path_widget = quick_output_path_widget(3, 'Result directory:')

        open_button = tk.Button(outputFrame, text="Open result directory", command=self.open_result_folder,
                                font=widget_font, justify='center')
        open_button.grid(row=4, column=0, sticky='we', padx=5, pady=5, columnspan=2)

        outputFrame.grid(row=1, column=0, **FRAME_FORMAT)
        outputFrame.columnconfigure(1, weight=1)

    def open_result_folder(self, *args):
        self.update_result_path()
        if self.result_dir is not None and os.path.isdir(self.result_dir) and os.path.exists(self.result_dir):
            if os.name == 'nt':  # Windows
                os.system(f"start {self.result_dir}")
            else:
                os.system(f"nautilus {self.result_dir}")
        messagebox.showwarning('Warning', f"Result directory does not exist.")

    def update_result_path(self, *args):
        if self.ArgMD_Directory.is_empty():
            return False
        self.result_dir = os.path.join(os.path.dirname(self.ArgMD_Directory.get_value()), 'results')

        self.EnergyRawFile = os.path.join(self.result_dir, 'results.csv')
        self.EnergyWCCFile = os.path.join(self.result_dir, 'results_analyze.csv')
        self.ConvergenceFile = os.path.join(self.result_dir, 'summary.pdf')

        def update_text(widget, text):
            widget.config(state='normal')
            widget.delete(1.0, tk.END)
            widget.insert(tk.END, text)
            widget.config(state='disabled')

        update_text(self.raw_ene_path_widget, self.EnergyRawFile)
        update_text(self.ana_ene_path_widget, self.EnergyWCCFile)
        update_text(self.conv_path_widget, self.ConvergenceFile)
        update_text(self.folder_path_widget, self.result_dir)

    def _sync_args_read(self):
        self.sync_from_toml(self.ArgMD_Directory, 'MD_Dir')
        self.sync_from_toml(self.ArgMD_Program, 'MD_Program')
        self.update_result_path()

    def _sync_args_write(self):
        if self.ArgMD_Directory.is_empty():
            return False
        self.update_result_path()
        self.toml_dict['ResultDir'] = self.result_dir
        self.toml_dict['EnergyRawFile'] = self.EnergyRawFile
        self.toml_dict['EnergyWCCFile'] = self.EnergyWCCFile
        self.toml_dict['ConvergenceFile'] = self.ConvergenceFile

    def report_args(self) -> Tuple[list, dict]:
        self.not_filled_arg_list = []
        self.active_arg_dict = {'-auto': None,
                                }

        self._check_arg(self.ArgMD_Directory)
        self._check_arg(self.ArgMD_Program)

        return self.not_filled_arg_list, self.active_arg_dict

    def report_submit_info(self) -> SubmitInfo:
        super().report_submit_info()

        if not self.submit_info.state:
            # Not completely filled
            return self.submit_info

        self.submit_info.work_dir = os.path.join(os.path.dirname(self.ArgMD_Directory.get_value()), 'results')
        self.submit_info.end_text = f"Collect results in {self.submit_info.work_dir}\n"
        self.submit_info.log_dir = os.path.dirname(self.ArgMD_Directory.get_value())

        self.gen_cmd()
        return self.submit_info
