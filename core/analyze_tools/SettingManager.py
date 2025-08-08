import copy
import os.path
import toml
from typing import Any, Tuple
autofill_token = '<AUTOFILL>'

class LineToolKit:
    def __init__(self):
        pass

    @staticmethod
    def value_parser(value: str) -> Any:
        if value.lower() == 'true':
            return True
        if value.lower() == 'false':
            return False
        if value.lower() == 'none':
            return None
        if value.isdigit():
            return int(value)
        if value.replace('.', '').isdigit():
            return float(value)

        return value.strip('"').strip("'")

    @staticmethod
    def extract_notes(line: str):
        line = line.strip()
        note = None
        content = line
        if line.find('#') != -1:
            content = line[:line.find('#')]
            note = line[line.find('#'):]
        return content, note

    @staticmethod
    def is_blank_line(line: str):
        if len(line.strip()) == 0:
            return True
        return False

    @staticmethod
    def extract_values(line: str) -> Tuple[str, str]:
        line = line.strip()
        key = line.split('=')[0].strip()
        value = line.split('=')[1].strip()
        return key, value

    @staticmethod
    def is_title_line(line: str):
        line = line.strip()
        if line.startswith('[') and line.endswith(']'):
            return True
        return False

    @staticmethod
    def path_pattern_parse(pattern: str):

        _pair_idx = pattern.lower().find('<pair>')
        _edge_idx = pattern.lower().find('<edge>')
        if _pair_idx == -1 or _edge_idx == -1:
            raise ValueError('Invalid pattern. Should contain both <pair> and <edge>')

        # print(pattern)
        # print(_pair_idx, _edge_idx)

        path_to_pair = pattern[:_pair_idx - 1]
        path_pair_to_edge = pattern[_pair_idx + 7:_edge_idx - 1]
        path_after_edge = pattern[_edge_idx + 7:]
        # print(path_to_pair, path_pair_to_edge)
        return path_to_pair, path_pair_to_edge, path_after_edge

    @staticmethod
    def path_gen_from_pattern(pair, edge, pattern: str) -> str:
        """
        :param pair: pair name
        :param edge: edge name
        :param pattern: path pattern. Example: /path/to/{pair}/to/{edge}
        :return: path fill with pair and edge.
        """
        return_path = pattern.replace('<pair>', pair)
        return_path = return_path.replace('<edge>', edge)
        return return_path


class SectionSettings:
    autofill_token = autofill_token

    def __init__(self):
        self.title = None
        self.__props = {}
        self.__props_note = {}

    @classmethod
    def from_block(cls, block: str):
        tmp_settings = cls()
        for line in block.split('\n'):
            if LineToolKit.is_blank_line(line):
                continue

            line = line.strip()

            if line.startswith('#'):
                # A note occupies whole line.
                continue

            if LineToolKit.is_title_line(line):
                # Title line
                tmp_settings.title = line[1:-1]
                continue

            content, note = LineToolKit.extract_notes(line)
            key, value = LineToolKit.extract_values(content)
            tmp_settings.set_prop(key, LineToolKit.value_parser(value), note)

        return tmp_settings

    def get_prop(self, prop_name: str) -> Any:
        if self.has_prop(prop_name):
            return self.__props[prop_name]
        return None

    def set_prop(self, prop_name: str, value: Any, note: str = None):
        self.__props[prop_name] = value
        if note is not None:
            self.__props_note[prop_name] = note

    def __all_props_filled(self) -> bool:
        for value in self.__props.values():
            if value == autofill_token:
                print("Haven't filled all properties.")
                return False
        return True

    def gen_file(self, filename: str):
        if self.__all_props_filled():
            with open(filename, 'a') as f:
                f.write(f'[{self.title}]\n')
                for key, value in self.__props.items():
                    f.write(f'{key} = {value}\n')
                f.write('\n\n')

    def output_to_block(self):
        if self.__all_props_filled():
            return_block = f'[{self.title}]\n'
            for key, value in self.__props.items():
                return_block += f'{key} = {value}\n'
            return return_block

    def has_prop(self, prop_name: str) -> bool:
        return prop_name in self.__props.keys()

    def apply_with_dict(self, props: dict):
        for key, value in props.items():
            self.__props[key] = LineToolKit.value_parser(f'{value}')

    def __str__(self):
        return_str = f'[{self.title}]\n'
        for key, value in self.__props.items():
            return_str += f'{key}: {value}\n'
        return return_str


class AllSettings:
    # TODO: use toml to replace this.

    def __init__(self, sections=None):
        self.__settings = {}
        if sections is not None:
            for section in sections:
                setattr(self, section, None)

    def add_section(self, title: str, settings: SectionSettings):
        if isinstance(settings, SectionSettings):
            self.__settings[title] = settings
            setattr(self, title, settings)
        else:
            raise TypeError(f"Settings of {title} is not a SectionSettings")

    def gen_file(self, file, section_setting_dict: dict = None):
        if section_setting_dict is not None:
            tmp_all_settings = copy.deepcopy(self)
            tmp_all_settings.apply_settings_with_dict(section_setting_dict)
            for section_title, section_settings in tmp_all_settings.__settings.items():
                section_settings.gen_file(file)
        else:
            for section_title, section_settings in self.__settings.items():
                section_settings.gen_file(file)

    @classmethod
    def from_file(cls, file):
        tmp_all_settings = AllSettings([])

        with open(file, 'r') as f:
            setting_block = ''
            for line in f:
                if len(line.strip()) == 0:
                    continue
                if LineToolKit.is_title_line(line):
                    if setting_block != '':
                        title = setting_block.split('\n')[0][1:-1]
                        tmp_all_settings.add_section(title, SectionSettings.from_block(setting_block))
                        setting_block = line
                        continue
                setting_block += line

            if setting_block != '':
                title = setting_block.split('\n')[0][1:-1]
                tmp_all_settings.add_section(title, SectionSettings.from_block(setting_block))

        return tmp_all_settings

    def apply_settings_with_dict(self, section_setting_dict: dict):
        """
        :param section_setting_dict: A dict. key is section title and value is section settings dict.
        Like {'alchemical':{'lambdas_json': 'lambdas.json', }, }
        :return:
        """
        for section_title, section_settings in section_setting_dict.items():
            tmp_section = getattr(self, section_title)
            tmp_section.apply_with_dict(section_settings)

    def has_prop(self, section_name, property_name):
        if section_name in self.__settings.keys():
            return self.__settings[section_name].has_prop(property_name)
        return False


def get_default_path(key):
    alchemd_core_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    default_config_file = os.path.join(alchemd_core_path, 'configs.toml')
    if os.path.isfile(default_config_file) and os.path.exists(default_config_file):
        with open(default_config_file, 'r') as f:
            config = toml.load(f)
        if config is None:
            raise Exception('Config file is empty')
        main_config = config['Main']
        if key in main_config.keys():
            return main_config[key]
    raise FileExistsError(f'Config file is not exist. Check configs.toml in {os.path.dirname(os.path.abspath(__file__))}')


if __name__ == '__main__':
    print(get_default_path('Alchemd_Path'))