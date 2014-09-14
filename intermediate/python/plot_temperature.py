import numpy as np
from matplotlib import pyplot

def read_csv_file(filename):
    '''
    This code will read in a CSV file of year, temp, rainfall, mosquitos
    '''
    assert type(filename) is str, 'filename must be a string'
    year, temperature, rainfall, mosquitos = np.genfromtxt(filename, unpack = True, skiprows = 1, delimiter = ',')
    return year, temperature, rainfall, mosquitos
    
def convert_fahrenheit_to_celsius(temp_in_f):
    '''
    This code will convert an array of temperatures from fahrenheit to celsius
    '''
    assert (type(temp_in_f) is float) or (type(temp_in_f) is int) or (isinstance(temp_in_f, np.ndarray)), 'temp must be float or int'
    temp_in_c = (temp_in_f - 32) * 5 /9.0
    return temp_in_c
    
def plot_data(x, y):
    '''
    This code will plot arrays of x, y with symbol 
    '''
    pyplot.plot(x, y, 'o')
    pyplot.savefig('mosquitos_data.pdf')
    pyplot.close()
    
if __name__ == "__main__":
	year, temperature, rainfall, mosquitos = read_csv_file('mosquito_data_A1.csv')
	temp_in_c = convert_fahrenheit_to_celsius(temperature)
	plot_data(year, mosquitos)
	print 'all done'