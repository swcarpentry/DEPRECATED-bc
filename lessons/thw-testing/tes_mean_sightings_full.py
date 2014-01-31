#!/usr/bin/python

'''
Unit tests for mean_sightings module. Uses sightings_tab_sm.csv as test data.
'''

from mean_sightings import get_sightings

# Always use this file for testing
filename = 'sightings_tab_sm.csv'

# Start the unit tests
def test_owl_is_correct():
    owlrec, owlmean = get_sightings(filename, 'Owl')
    assert owlrec == 2, 'Number of records for owl is wrong'
    assert owlmean == 17, 'Mean sightings for owl is wrong'

def test_muskox_is_correct():
    oxrec, oxmean = get_sightings(filename, 'Muskox')
    assert oxrec == 2, 'Number of records for muskox is wrong'
    assert oxmean == 25.5, 'Mean sightings for muskox is wrong'

def test_animal_not_present():
    animrec, animmean = get_sightings(filename, 'NotPresent')
    print animrec, animmean
    assert animrec == 0, 'Animal missing should return zero records'
    assert animmean == 0, 'Animal missing should return zero mean'

def test_arg_capitalization_wrong():
    owlrec, owlmean = get_sightings(filename, 'oWl')
    assert owlrec == 2, 'Number of records for owl is wrong'
    assert owlmean == 17, 'Mean sightings for owl is wrong'
