import plot_temperature
import numpy as np

def test_output_read_csv_file():
    '''
    Test that read_csv_file gives the expected output
    '''
    year, temperature, rainfall, mosquitos = plot_temperature.read_csv_file('mosquito_data_A1.csv')
    year_key = np.arange(2001, 2011, 1)
    assert (year == year_key).all(), 'year doesn\'t match key'
    
def test_convert_fahrenheit_to_celsius():
    temp_c = plot_temperature.convert_fahrenheit_to_celsius(32)
    assert temp_c == 32, '32F != 0C'