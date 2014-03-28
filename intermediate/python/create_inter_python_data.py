"""Create the data for the Software Carpentry Intermediate Python lectures"""

import numpy as np
import pandas as pd

np.random.seed(26)
datasets = {'A1': [0, 0.5, 0.7, 10],
            'A2': [0, 0.5, 0.7, 50],
            'A3': [0, 0.5, 0.3, 50],
            'B1': [3, 0.7, 0.2, 50],
            'B2': [3, 0.7, 0.7, 50]}

def make_data(intercept, tempslope, rainfallslope, numyears):
    years = np.arange(2010 - numyears, 2011)
    temps = np.random.uniform(70, 90, len(years))
    rainfalls = np.random.uniform(100, 300, len(years))
    noise = 2 * np.random.randn(len(years))
    mosquitos = intercept + tempslope * temps + rainfallslope * rainfalls + noise
    return zip(years, temps, rainfalls, mosquitos)

def export_data(data, filename):
    df = pd.DataFrame(data, columns=['year', 'temperature', 'rainfall','mosquitos'])
    df.to_csv(filename, index=False, float_format='%.0f')

for site in datasets:
    data = make_data(*datasets[site])
    if site == 'A1':
        #create a shorter dataset for first example
        data = data[-10:]
    export_data(data, '%s_mosquito_data.csv' % site)
