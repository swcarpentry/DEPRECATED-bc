
# two_tasks.py

def task_reformat_temperature_data():
    """Reformats the raw temperature data file for easier analysis"""
        
    return {
        'file_dep': ['UK_Tmean_data.txt'],
        'targets': ['UK_Tmean_data.reformatted.txt'],
        'actions': ['python reformat_weather_data.py UK_Tmean_data.txt > UK_Tmean_data.reformatted.txt'],
}

def task_reformat_sunshine_data():
    """Reformats the raw sunshine data file for easier analysis"""

    return {
        'file_dep': ['UK_Sunshine_data.txt'],
        'targets': ['UK_Sunshine_data.reformatted.txt'],
        'actions': ['python reformat_weather_data.py UK_Sunshine_data.txt > UK_Sunshine_data.reformatted.txt'],
    }