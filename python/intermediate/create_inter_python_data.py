"""Create the data for the Software Carpentry Intermediate Python lectures"""

import numpy as np
import pandas as pd

np.random.seed(26)
years = np.arange(1960, 2011)
temps = np.random.uniform(70, 90, len(years))
rainfalls = np.random.uniform(100, 300, len(years))
noise = 2 * np.random.randn(len(years))
mosquitos = 0.5 * temps + 0.7 * rainfalls + noise

data = zip(years, temps, rainfalls, mosquitos)
df = pd.DataFrame(data, columns=['year', 'temperature', 'rainfall','mosquitos'])
df.to_csv('mosquito_data_A2.csv', index=False, float_format='%.0f')
df_short = df[-10:]
df_short.to_csv('mosquito_data_A1.csv', index=False, float_format='%.0f')
