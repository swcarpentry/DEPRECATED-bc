
# monthly_raw_data_update.py

import datetime
from doit.tools import timeout 

data_sets = ['Tmean', 'Sunshine']

def get_data_file_parameters(data_type):
    """Takes a string describing the type of climate data, returns url and file name for that data"""

    base_url = 'http://www.metoffice.gov.uk/climate/uk/datasets/{0}/ranked/UK.txt'
    data_url = base_url.format(data_type)
    data_target = 'UK_{0}_data.txt'.format(data_type)
    return data_url, data_target

def task_download_data():
    """Downloads all raw data files from the Met Office website"""

    for data_type in data_sets:
        data_url, data_target = get_data_file_parameters(data_type)
        yield {
            'actions': ['wget -O %(targets)s {0}'.format(data_url)],
            'targets': [ data_target ],
            'name' : data_type,
            'uptodate': [timeout(datetime.timedelta(weeks=4))],

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