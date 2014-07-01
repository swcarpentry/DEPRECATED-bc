
# sub_tasks_no_name.py

data_sets = ['Tmean', 'Sunshine']

def task_reformat_data():
    """Reformats all raw files for easier analysis"""

    for data_type in data_sets:
        yield {
            'actions': ['python reformat_weather_data.py %(dependencies)s > %(targets)s'],
            'file_dep': ['UK_{}_data.txt'.format(data_type)],
            'targets': ['UK_{}_data.reformatted.txt'.format(data_type)],
        }
    