def task_get_temp_data():
    data_url = 'http://www.metoffice.gov.uk/climate/uk/datasets/Tmean/ranked/UK.txt'
    data_target = 'UK_mean_temp.txt'
    return {
        'actions': ['wget -O %(targets)s {0}'.format(data_url)],
        'targets': [ data_target ],
    }

def task_get_sunshine_data():
    data_url = 'http://www.metoffice.gov.uk/climate/uk/datasets/Sunshine/ranked/UK.txt'
    data_target = 'UK_Sunshine.txt'
    return {
        'actions': ['wget -O %(targets)s {0}'.format(data_url)],
        'targets': [ data_target ],
    }

def task_unstack_sunshine_data():
    return {
        'actions': ['python unstack_weather_data.py %(dependencies)s > %(targets)s'],
        'file_dep': ['UK_Sunshine_data.txt'],
        'targets': ['UK_Sunshine_data.unstacked.txt'],
    }


