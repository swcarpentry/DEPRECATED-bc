import os

def get_data_file_parameters(data_type):
    base_url = 'http://www.metoffice.gov.uk/climate/uk/datasets/{0}/ranked/UK.txt'
    data_url = base_url.format(data_type)
    data_target = 'UK_{0}_data.txt'.format(data_type)
    return data_url, data_target

def get_unstacked_target(data_file):
    base, ext = os.path.splitext(data_file)
    return '.'.join([base, 'unstacked', ext])

def task_data():
    for data_type in ['Tmean', 'Sunshine']:
        data_url, data_target = get_data_file_parameters(data_type)
        yield {
            'actions': ['wget -O %(targets)s {0}'.format(data_url)],
            'targets': [ data_target ],
            'name': 'download {0} data'.format(data_type),
        }

def task_unstack_sunshine_data():
    for data_type in ['Tmean', 'Sunshine']:
        data_url, data_dep = get_data_file_parameters(data_type)
        data_target = get_unstacked_target(data_dep)
        yield {
            'actions': ['python unstack_weather_data.py %(dependencies)s > %(targets)s'],
            'file_dep': [data_dep],
            'targets': [data_target],
            'name': 'unstack {0}'.format(data_dep)
        }


