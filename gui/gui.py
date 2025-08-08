# Pyinstaller -F --clean -w .\gui.py
import os.path
from tkinter import messagebox
import sys, textwrap, subprocess, signal, base64
import datetime as dt
from datetime import datetime
import toml
import threading
import signal

from pages import *


def get_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except AttributeError:
        base_path = os.path.abspath(".")

    return os.path.normpath(os.path.join(base_path, relative_path))


def get_application_path():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    else:
        return os.path.dirname(os.path.abspath(__file__))


ALCHEMD_TOOLKIT_DIR = os.path.dirname(get_application_path())
ALCHEMD_CORE_DIR = os.path.join(ALCHEMD_TOOLKIT_DIR, "core")
AVAIL_FONT = 'Arial'

autofill_token = '<AUTOFILL>'


def timedelta2str(delta: dt.timedelta):
    seconds = delta.seconds % 60
    minutes = delta.seconds % 3600 // 60
    hours = delta.seconds // 3600
    days = delta.days
    return_str = ''
    if days > 0:
        return_str += f"{days}day "
    return_str += f"{hours}h{minutes}m{seconds}s"
    return return_str


class InfoPanel(tk.LabelFrame):
    break_line = '-'*42+"\n"

    def __init__(self, root, request_command_method):

        self.process = None
        self.running = False

        super().__init__(root)
        self.config(relief='groove', borderwidth=2, text='Log', font=(AVAIL_FONT, 12))

        self.request_command_method = request_command_method

        self.log_frame = tk.Frame(self, relief='flat')
        self.log_frame.pack(fill=tk.BOTH, side='top', expand=False)
        self.log_window = tk.Text(self.log_frame, height=12, relief='flat', state='disabled')
        self.scrollbar = tk.Scrollbar(self.log_frame, orient=tk.VERTICAL, command=self.log_window.yview)
        self.log_window.configure(yscrollcommand=self.scrollbar.set)

        self.log_window.grid(row=0, column=0, sticky='nsew')
        self.scrollbar.grid(row=0, column=1, sticky='ns')
        self.log_frame.grid_rowconfigure(0, weight=1)
        self.log_frame.grid_columnconfigure(0, weight=1)

        style = ttk.Style()
        style.theme_use('default')
        style.configure("Custom.Horizontal.TProgressbar",
                        # troughcolor='#ecf1f7',
                        background='#0869a0',
                        )
        self.progress_bar = ttk.Progressbar(self, style='Custom.Horizontal.TProgressbar',
                                            orient='horizontal',
                                            mode="determinate", length=200)
        self.progress_bar.pack(expand=False, fill='x', side='bottom')

        self.go_button = tk.Button(self, text="GO", justify='center',
                                   font=(AVAIL_FONT, 18, 'bold'), fg='#0869a0', bg='#ecf1f7',
                                   activebackground='#8da9ba',
                                   activeforeground='#0869a0',
                                   command=lambda: self.go())
        self.go_button.pack(expand=False, fill='x', side='bottom')

        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=0)
        self.grid_rowconfigure(2, weight=0)

    def add_log(self, text):
        self.log_window.config(state='normal')
        lines = set(text.split('\n'))
        for line in lines:
            if line == '':
                continue
            self.log_window.insert(tk.END, line+"\n")
            self.log_window.see(tk.END)
        self.log_window.config(state='disabled')
        self.log_window.see(tk.END)

    def go(self, **kwargs):
        if self.running and self.process is not None:
            if os.name == 'nt':  # Windows
                subprocess.run(['taskkill', '/F', '/T', '/PID', str(self.process.pid)])
            else:  # Linux/Unix
                os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
                
            self.add_log(f'Terminated process group {self.process.pid}')
            self.go_button.config(text='Go', state='normal', relief='raised')
            self.progress_bar['value'] = 0
            self.running = False
            self.process = None
            return 0

        submit_info = self.request_command_method()

        self.log_window.config(state='normal')
        self.log_window.delete(1.0, tk.END)
        self.log_window.config(state='disabled')
        # Clean log

        if submit_info.page_type == 'info':
            self.add_log(submit_info.extra_info+"\n")
            return 0
        if not submit_info.state:
            not_filled_arg_list = submit_info.extra_info
            self.add_log(f"These arguments are not filled, invalid or the path does not exist:\n{not_filled_arg_list}")
            return 0

        command = submit_info.run_cmd
        work_dir = submit_info.log_dir
        self.finish_info = submit_info.end_text
        self.add_log(command+"\n")
        self.add_log(submit_info.begin_text+"\n")

        self.add_log('\nGO.\n')
        self.progress_bar['value'] = 0
        self.go_button.config(text='Stop')
        self.running = True
        begin_time = datetime.now()
        self._last_log = ''

        self.update()

        # Set environment variables to ensure unbuffered output
        env = os.environ.copy()
        env['PYTHONUNBUFFERED'] = '1'
        env['PYTHONIOENCODING'] = 'utf-8'
        
        # Prepare process creation arguments
        popen_args = {
            'stdout': subprocess.PIPE,
            'stderr': subprocess.STDOUT,  # Merge stderr to stdout
            'encoding': 'utf-8',
            'errors': 'replace',
            'universal_newlines': True,
            'bufsize': 0,  # Unbuffered
            'shell': True,
            'env': env
        }
        
        # Add platform-specific arguments
        if os.name != 'nt':  # Unix/Linux
            popen_args['preexec_fn'] = os.setsid
        else:  # Windows
            popen_args['creationflags'] = subprocess.CREATE_NEW_PROCESS_GROUP
        
        self.process = subprocess.Popen(command, **popen_args)

        def check_process():
            if not self.running:
                return

            while self.running and self.process and self.process.poll() is None:
                try:
                    output = self.process.stdout.readline()
                    if not output:
                        continue
                        
                    output = output.strip()
                    if not output:
                        continue
                        
                    log, progress_now = self.parse_child_process_log(output)
                    if log and log != self._last_log:
                        def update_log(log_text=log):
                            self.add_log(log_text + '\n')
                        self.after(0, update_log)
                        self._last_log = log
                    if progress_now is not None:
                        def update_progress(progress=progress_now):
                            self.progress_bar.configure(value=progress)
                        self.after(0, update_progress)
                except (UnicodeDecodeError, IOError) as e:
                    # Handle specific encoding or I/O errors
                    error_msg = f"Output reading error: {str(e)}"
                    self.after(0, lambda msg=error_msg: self.add_log(f"{msg}\n"))
                    break
                except Exception as e:
                    # Handle unexpected errors
                    error_msg = f"Unexpected error in process monitoring: {str(e)}"
                    self.after(0, lambda msg=error_msg: self.add_log(f"{msg}\n"))
                    break

            if not self.running or not self.process:
                return

            if work_dir is not None:
                err_log = os.path.join(work_dir, 'Err.log')
            else:
                err_log = 'Err.log'
            
            # Since stderr is merged to stdout, create empty error log
            error_message = ""
            with open(err_log, 'w') as f:
                f.write("stderr merged to stdout for real-time display\n")

            return_code = self.process.poll()
            if return_code != 0:
                self.after(0, lambda: messagebox.showwarning('Error', 
                    f'Process ended with error (return code: {return_code})\n\nError message:\n{error_message}'))
                self.go_button.config(text='Go', state='normal', relief='raised')
                self.running = False
                self.process = None
                return 0
            self.progress_bar['value'] = progress_now if progress_now is not None else 0

            self.add_log('\n'+'-'*42+'\n')

            elapsed_time = datetime.now() - begin_time
            formatted_time = timedelta2str(elapsed_time)
            begin_time_str = begin_time.strftime('%Y-%m-%d %H:%M:%S')
            end_time_str = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            info_str = f'Begin time: {begin_time_str}\nEnd time: {end_time_str}\nTime elapsed: {formatted_time}\n'
            self.add_log(info_str)
            self.add_log(self.break_line+self.finish_info)
            self.add_log("Please check Err.log if any error occurred.")
            
            if not error_message:
                messagebox.showinfo('Finish', f"{info_str}\n{self.finish_info}")

            self.go_button.config(text='Go', state='normal', relief='raised')
            self.running = False
            self.process = None

        # Start monitoring in a separate thread
        monitor_thread = threading.Thread(target=check_process)
        monitor_thread.daemon = True
        monitor_thread.start()

    @staticmethod
    def parse_child_process_log(log):
        if log == '':
            return log, None
        if log.find("%") != -1 and log.find("|") != -1:
            # tqdm log
            progress = int(log.split('%')[0].split(" ")[-1])
            log = log.split(":")[0]+f": {progress}%"
            return log, progress
        try:
            nums = log.split("(")[1].split(")")[0].split('/')
            now_section = 1
            all_section = 1
            if log.find('[') != -1 and log.find(']') != -1:
                section = log.split("[")[1].split("]")[0].split('/')
                now_section = int(section[0].strip())
                all_section = int(section[1].strip())
            progress = int(int(nums[0].strip()) / int(nums[1].strip()) * 100 / all_section
                           + 100 * (now_section - 1) / all_section)

            return log, progress
        except:
            return log, None


class SettingPage(ScrollableFrame):
    @dataclass
    class path_setting:
        name: str
        widget: PathWidget
        value: str = ''

    def __init__(self, root, config_filename='configs.toml'):
        super().__init__(root)

        path_list = ['Alchemd_Path',
                     'CAR_Path',
                     'AlchemdConvTools_Path',
                     'AlchemdToolkit_Path',
                     'GenambRBFE_Path',
                     'Watvina_Path',
                     'Wcc_Path']

        path_list_demo = [
            'AlchemdToolkit_Path'
        ]


        self.path_dict = {i: SettingPage.path_setting(name=i,
                                                      widget=PathWidget(self.scrollable_frame,
                                                                        text=i.replace('_', ' '))
                                               ) for i in path_list_demo}
        grid_row_count = 0
        path_setting_num = len(self.path_dict)

        for path_setting in self.path_dict.values():
            path_setting.widget.grid(row=grid_row_count, column=0, sticky='ew')
            grid_row_count += 1

        save_small_frame = tk.Frame(self.scrollable_frame)
        self.save_button = tk.Button(save_small_frame, text='Save', font=(AVAIL_FONT, 12), command=self.save_config)
        self.save_button.grid(row=0, column=0, sticky='wn', padx=10, pady=10)
        self.checking_status_label = tk.Label(save_small_frame, text=' ')
        self.checking_status_label.grid(row=0, column=1, sticky='w', pady=10)
        save_small_frame.grid(row=path_setting_num, sticky='wen', padx=10)

        info = textwrap.dedent("""\
        Alchemd (Alchemd GUI)

            _    _      _                        _ 
           / \  | | ___| |__   ___ _ __ ___   __| |
          / _ \ | |/ __| '_ \ / _ \ '_ ` _ \ / _` |
         / ___ \| | (__| | | |  __/ | | | | | (_| |
        /_/   \_\_|\___|_| |_|\___|_| |_|_|_|\__,_|


        Still waters run deep.
        
        Authors: Runduo Liu, Yilin Zhong, Yufen Yao,
                 Wanyi Huang, Hai-Bin Luo, Zhe Li \
        """

        # Author: <NAME>
        # Reference: <WEBSITE>
        # Other information: <Other>\
        # """
                               )
        print(info, flush=True)
        self.info_widget = tk.Label(self.scrollable_frame, text=info, font=("JetBrains Mono", 12), justify='left',
                                    anchor='w', relief='groove', bd=3)
        self.info_widget.grid(row=path_setting_num+1, column=0, sticky='wes', padx=5, pady=5, ipadx=0, ipady=0)

        self.grid_columnconfigure(0, weight=1)

        self.config_file = os.path.join(get_application_path(), config_filename)
        self.toml_dict = {}
        self.main_config_dict = None
        self.config_filename = config_filename

        if os.path.exists(self.config_file) and os.path.isfile(self.config_file):
            with open(self.config_file, 'r') as f:
                self.toml_dict = toml.load(self.config_file)
                self.main_config_dict = self.toml_dict['Main']
            self.update_variables()
        else:
            self.main_config_dict = {
                'Alchemd_Path': 'path to Alchemd',
                'CAR_Path': 'path to CAR',
                'AlchemdConvTools_Path': 'path to AlchemdConvTools',
                'AlchemdToolkit_Path': ALCHEMD_TOOLKIT_DIR,
                'GenambRBFE_Path': 'path to GenambRBFE',
                'Watvina_Path': os.path.join(ALCHEMD_TOOLKIT_DIR, 'dependencies', 'watvina'),
                'Wcc_Path': 'path to Wcc',
            }
            self.toml_dict['Main'] = self.main_config_dict
            self.update_variables()
            self.save_config()

    def update_variables(self):
        if self.main_config_dict is None:
            return 0
        for name, path_setting in self.path_dict.items():
            if name in self.main_config_dict.keys():
                path_setting.value = self.main_config_dict[name]
                path_setting.widget.update_text(path_setting.value)

    def save_config(self):
        alchemdTK_path_config_file = os.path.join(self.path_dict['AlchemdToolkit_Path'].widget.get_text(),
                                                  'core', self.config_filename)

        self.main_config_dict = {name: path_setting.widget.get_text() for name, path_setting in self.path_dict.items()}
        self.setup_paths()
        if not os.path.exists(alchemdTK_path_config_file):
            # Create file in a cross-platform way
            open(alchemdTK_path_config_file, 'a').close()
        # print(self.main_config_dict, flush=True)
        self.checking_status_label.config(text="Checking dependency availability...")
        self.update()
        failed_dependencies = self.test_dependencies()

        if len(failed_dependencies) > 0:
            _nums = len(failed_dependencies)
            msg = 'The following dependencies are invalid:\n'
            msg += "\n".join(failed_dependencies)
            print_info = f"Invalid dependencies: {failed_dependencies}" if _nums > 1\
                else f"{failed_dependencies[0]} is invalid."
            self.checking_status_label.config(text=print_info)
            self.update()
            messagebox.showwarning(title=f"Warning!", message=msg)
            return False
        else:
            self.checking_status_label.config(text="Dependencies are all valid and the config is saved.")
            self.update()
            messagebox.showinfo(title="Saved.", message="Dependencies are all valid.")

        self.toml_dict['Main'] = self.main_config_dict
        self.update_variables()
        with open(self.config_file, 'w') as f:
            toml.dump(self.toml_dict, f)
        with open(alchemdTK_path_config_file, 'w') as f:
            toml.dump(self.toml_dict, f)

    def setup_paths(self, from_tk_path=True):
        if from_tk_path:
            tk_path = self.main_config_dict['AlchemdToolkit_Path']
            dependencies_path = os.path.join(tk_path, 'dependencies')
            self.main_config_dict['Alchemd_Path'] = os.path.join(dependencies_path, 'Alchemd')
            self.main_config_dict['CAR_Path'] = os.path.join(dependencies_path, 'CAR-FEP')
            self.main_config_dict['GenambRBFE_Path'] = os.path.join(dependencies_path, 'genambRBFE')
            self.main_config_dict['AlchemdConvTools_Path'] = os.path.join(dependencies_path, 'AlchemConvTools')
            self.main_config_dict['Wcc_Path'] = os.path.join(dependencies_path, 'Weighted_cc-main')
            self.main_config_dict['Watvina_Path'] = os.path.join(dependencies_path, 'watvina')


    def test_dependencies(self) -> list:
        """
        check how many dependencies are available
        :return: dependencies not available
        """
        def _test_program_fail(path, program, is_bin=False, is_sh=False) -> bool:
            if is_bin:
                exec_cmd = f"{os.path.join(path, program)} -h >/dev/null 2>&1"
            elif is_sh:
                exec_cmd = f"/bin/bash {os.path.join(path, program)} -h >/dev/null 2>&1"
            else:
                exec_cmd = f"python {os.path.join(path, program)} -h >/dev/null 2>&1"
            return_code = os.system(exec_cmd)
            return_code >>= 8

            if return_code != 0:
                return True
            return False

        failed_dependencies = []
        if not os.path.exists(os.path.join(self.main_config_dict['GenambRBFE_Path'], 'genambRBFE')):
            failed_dependencies.append('GenambRBFE')
        if _test_program_fail(self.main_config_dict['Watvina_Path'], 'watvina', is_bin=True):
            failed_dependencies.append('Watvina')
        if _test_program_fail(self.main_config_dict['Wcc_Path'], 'wcc_main.py'):
            failed_dependencies.append('Wcc')
        # if _test_program_fail(self.main_config_dict['AlchemdToolkit_Path'], 'core/prepare_file.py'):
        #     failed_dependencies.append("AlchemdToolkit")
        if _test_program_fail(self.main_config_dict['CAR_Path'], 'segmented_converge_control.py'):
            failed_dependencies.append("CAR")
        # if _test_program_fail(self.main_config_dict['Alchemd_Path'], 'openmm-FEP-run.py'):
        #     failed_dependencies.append('Alchemd')
        if _test_program_fail(self.main_config_dict['AlchemdConvTools_Path'], 'one_end_fe_aly.py'):
            failed_dependencies.append('AlchemConvTools')

        return failed_dependencies


class App:
    def __init__(self, root):
        outer_frame = tk.Frame(root)
        outer_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True, padx=5, pady=5)
        navi_bar = NaviBar(outer_frame)
        info_panel = InfoPanel(outer_frame, navi_bar.get_active_page_info)
        outer_frame.grid_columnconfigure(0, weight=1)
        outer_frame.grid_rowconfigure(0, weight=1)
        navi_bar.grid(row=0, column=0, sticky='nsew', padx=5, pady=5)
        info_panel.grid(row=1, column=0, sticky='sew', padx=5, pady=5)


class NaviBar(tk.Frame):
    def __init__(self, root):

        super().__init__(root)

        self.default_bg = self.cget("bg")
        self.now_page = ""
        self.job_toml_dict = {}

        self.content_frame = tk.Frame(self, relief='groove', borderwidth=2)

        self.pages = {"Prepare": PreparePage(self.content_frame),
                      "Run": RunPage(self.content_frame),
                      "Analyze": AnalyzePage(self.content_frame),
                      "Setting": SettingPage(self.content_frame)}

        self.content_frame.grid(row=1, column=0, sticky='nsew', pady=0, ipadx=0, ipady=0, columnspan=len(self.pages))

        self.button_list = {}
        self.page_list = {}
        self.main_config_dict = self.pages["Setting"].main_config_dict

        self.grid_rowconfigure(1, weight=1)

        for key, value in self.pages.items():
            self.add_button(key, value)

        self.switch_page(list(self.button_list.keys())[0])

    def switch_page(self, name):

        self.main_config_dict = self.pages["Setting"].main_config_dict

        button_active_list = list(self.button_list.keys())
        button_active_list.remove(name)
        page_hide_list = list(self.page_list.keys())
        page_hide_list.remove(name)

        if self.now_page != '':
            last_activate_page = self.page_list[self.now_page]
            if isinstance(last_activate_page, ABC_ConfigPage):
                last_activate_page.sync_args('w')
                self.job_toml_dict = last_activate_page.toml_dict
        self.now_page = name

        self.button_list[name].config(state=DISABLED, relief=FLAT, bg='#daeaf2')
        for button_name in button_active_list:
            self.button_list[button_name].config(state=ACTIVE, relief='raise', bg=self.default_bg)

        active_page = self.page_list[name]
        active_page.pack(fill='both', side='top', expand=True)
        active_page.switch_update()
        if isinstance(active_page, ABC_ConfigPage):
            active_page.sync_args('r', from_dict=self.job_toml_dict)
        if hasattr(active_page, 'script_base'):
            setattr(active_page, 'script_base', self.main_config_dict['AlchemdToolkit_Path'])
        for page_name in page_hide_list:
            self.page_list[page_name].pack_forget()

    def add_button(self, name, page_frame: tk.Frame):
        tmp_button = tk.Button(self, text=name, command=lambda: self.switch_page(name), font=(AVAIL_FONT, 11))
        self.columnconfigure(len(self.button_list), weight=1)
        tmp_button.grid(row=0, column=len(self.button_list), sticky="ew")
        self.button_list[name] = tmp_button
        self.page_list[name] = page_frame

    def get_active_page_info(self) -> SubmitInfo:
        now_page = self.page_list[self.now_page]

        if isinstance(now_page, ABC_ConfigPage):
            return now_page.report_submit_info()

        return SubmitInfo(state=False, page_type='info', extra_info='This is setting page.', work_dir=None)


if __name__ == '__main__':
    window = tk.Tk()
    window.title("Alchemd")
    window.geometry("600x970")
    window.minsize(width=600, height=400)

    App(window)
    window.mainloop()
