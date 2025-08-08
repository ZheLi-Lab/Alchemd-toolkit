import time
import datetime
from typing import List, Union, Tuple


class DateParser(object):
    def __init__(self):
        pass

    @staticmethod
    def str2datetime(date_str, offset_hour: int = 0) -> datetime.datetime:
        if len(date_str) == 19:  # 2024-08-01 09:18.14  total 19 chars
            return_time = datetime.datetime.strptime(date_str, "%Y-%m-%d %H:%M.%S")
            return_time = return_time + datetime.timedelta(hours=offset_hour)
            return return_time
        elif len(date_str) == 31:  # Fri, 23 Aug 2024 02:52:21 -0400 total 31 chars
            return_time = datetime.datetime.strptime(date_str[:25], "%a, %d %b %Y %H:%M:%S")
            offset_hour = int(date_str[-5:-2])
            return_time = return_time - datetime.timedelta(hours=offset_hour)
            return return_time
        else:
            raise Exception(f"len of date_str is not included: {len(date_str)} {date_str}")

    @staticmethod
    def _timedelta_convert(timedelta: datetime.timedelta) -> Tuple[int, int, int, int]:
        days = timedelta.days
        seconds = timedelta.total_seconds()
        minutes = int(seconds // 60 % 60)
        hours = int(seconds // 3600 % 24)
        seconds = int(seconds % 60)
        return days, hours, minutes, seconds
    
    @staticmethod
    def __element_parser(ele: Union[str,datetime.datetime]) -> datetime.datetime:
        if isinstance(ele, str):
            return_datetime = DateParser.str2datetime(ele)
        elif isinstance(ele,datetime.datetime):
            return_datetime = ele
        else:
            raise Exception(f"{type(time1)} is not avail.")
        return return_datetime
        
    @staticmethod
    def get_time_gap_str(time1, time2) -> str:
        t1 = DateParser.__element_parser(time1)
        t2 = DateParser.__element_parser(time2)
        gap = DateParser._timedelta_convert(t2 - t1)
        return_time_str = f'{gap[2]}m'
        if gap[1] > 0:
            return_time_str = f'{gap[1]}h' + return_time_str
        if gap[0] > 0:
            return_time_str = f'{gap[0]}d' + return_time_str
        return return_time_str

    @staticmethod
    def time_parser(time_str):
        d = 0
        h = 0
        if time_str.find('d') != -1:
            d = int(time_str[:time_str.find('d')])
            time_str = time_str[time_str.find('d') + 1:]
        if time_str.find('h') != -1:
            h = int(time_str[:time_str.find('h')]) + d * 24
            time_str = time_str[time_str.find('h') + 1:]
        m = int(time_str[:time_str.find('m')])
        return h, m

    @staticmethod
    def convert2min(time_str):
        if time_str is None:
            return None
        sh, sm = DateParser.time_parser(time_str)
        return sh * 60 + sm

    @staticmethod
    def convert_from_min(minute: int) -> str:
        h = minute // 60
        d = h // 24
        m = minute % 60
        h = h % 24
        return_text = f'{m}m'
        if h > 0:
            return_text = f'{h}h' + return_text
        if d > 0:
            return_text = f'{d}d' + return_text
        return return_text
