#!/usr/bin/env python

'''
Module/script to calculate mean number of sightings of a given animal in a 
given sightings csv file.
'''

import sys
import matplotlib.mlab as ml
import numpy as np


def get_sightings(filename, focusanimal):

    # Load table
    tab = ml.csv2rec(filename)

    # Standardize capitalization of focusanimal
    focusanimal = focusanimal.capitalize()

    # Find number of records and total count of animals seen
    isfocus = (tab['animal'] == focusanimal)
    totalrecs = np.sum(isfocus)

    if totalrecs == 0:
        meancount = 0
    else:
        meancount = np.mean(tab['count'][isfocus])

    # Return num of records and animals seen
    return totalrecs, meancount


def get_sightings_loop(filename, focusanimal):

    # Load table
    tab = ml.csv2rec(filename)

    # Standardize capitalization of focusanimal
    focusanimal = focusanimal.capitalize()

    # Loop through all records, countings recs and animals
    totalrecs = 0.
    totalcount = 0.
    for rec in tab:
        if rec['animal'] == focusanimal:
            totalrecs += 1
            totalcount += rec['count']

    if totalrecs==0:
        meancount = 0
    else:
        meancount = totalcount/totalrecs

    # Return num of records and animals seen
    return totalrecs, meancount

if __name__ == '__main__':
    #print sys.argv
    filename = sys.argv[1]
    focusanimal = sys.argv[2]
    print get_sightings(filename, focusanimal)
