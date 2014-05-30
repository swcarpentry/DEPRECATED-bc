import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser(description='Reformats a met-office weather data file. Input data has one row per year and one column per month. Output data has a date column and a value column.')
parser.add_argument('data_file',metavar='DATA_FILE', help='Data file containing met office weather stats.')

def get_month_data(month_data):
    """Takes a two column DataFrame and returns a one column DataFrame where the index is a pandas period""" 
    
    # Each pair of columns contains data for a specific month. What month is this?
    month = month_data.columns[0]
        
    # Given a year, return a Period object representing this month in that year
    def get_period(year):
        return pd.Period('{0} {1}'.format(month,year))
    
    # Change the index of the dataframe to be the monthly period
    month_data.index = map(get_period, full_data.iloc[:,2])
    
    # Rename the columns
    month_data.columns = ['value', 'year']
    month_data.index.name = 'month'
        
    # Remove the year column, we don't need it anymore.
    return month_data.drop('year', axis=1)

def unstack_data(full_data):
    """
    Takes a dataframe with monthly columns and yearly rows. Returns a single column
    dataframe where each row is a specific month of a specific year.
    """

    # Loop over columns of the DataFrame in groups of 2
    # and feed them to get_month_data
    monthly_data = [ get_month_data(full_data.iloc[:, i:i+2]) for i in range(1,25,2)]

    # Add all the data for the individual months together
    unstacked_data = pd.concat(monthly_data)

    return unstacked_data.sort()

if __name__ == '__main__':

    args = parser.parse_args()

    full_data = pd.read_csv(args.data_file, delim_whitespace=True,
            skiprows=7)

    unstacked_data = unstack_data(full_data)

    # Write the new dataframe to stdout
    try:
        unstacked_data.to_csv(sys.stdout)

    # Don't fall over if we pipe the output to head
    except IOError:
        pass

    
