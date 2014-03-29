def task_get_temp_data():
    return {
        'actions': ['wget -O %(targets)s http://www.metoffice.gov.uk/climate/uk/datasets/Tmean/ranked/UK.txt'],
        'targets': ['UK_mean_temp.txt'],
    }
