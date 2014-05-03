
# download_temp_data.py

import datetime
from doit.tools import timeout 

data_sets = ['Tmean', 'Sunshine']

def task_get_temp_data():
    """Downloads the raw temperature data from the Met Office"""

    return {
        'actions': ['wget -O %(targets)s http://www.metoffice.gov.uk/climate/uk/datasets/Tmean/ranked/UK.txt'],
        'targets': ['UK_Tmean_data.txt'],
    }

def task_reformat_data():
    """Reformats all raw files for easier analysis"""

    for data_type in data_sets:
        yield {
            'actions': ['python reformat_weather_data.py %(dependencies)s > %(targets)s'],
            'file_dep': ['UK_{}_data.txt'.format(data_type)],
            'targets': ['UK_{}_data.reformatted.txt'.format(data_type)],
            'name': 'UK_{}_data.txt'.format(data_type),
        }