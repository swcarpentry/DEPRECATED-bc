
# one_task.py

def task_reformat_temperature_data():
    """Reformats the raw temperature data file for easier analysis"""
    
    return {
        'file_dep': ['UK_Tmean_data.txt'],
        'targets': ['UK_Tmean_data.reformatted.txt'],
        'actions': ['python reformat_weather_data.py UK_Tmean_data.txt > UK_Tmean_data.reformatted.txt'],
    }
