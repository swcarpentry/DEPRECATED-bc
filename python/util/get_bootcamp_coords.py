#!/usr/bin/env python

'''Extract the latitude and longitude coordinates stored for each bootcamp
in the SWC site repository for the metadata in the index.html yaml file and
write the output to a csv file.
'''

import os
import argparse
import glob
import itertools

import yaml
try:
    from yaml import CLoader as YLoader
except ImportError:
    # fall back on Python implementation
    from yaml import Loader as YLoader


def get_lat_lon(path):
    '''Extract the lat and lon metadata from the yaml file 
    specified by `path`'''
    with open(path, 'rt') as f:
        data = yaml.load_all(f, Loader=YLoader)

        # Only grab the first yaml document that contains the 
        # bootcamp metadata. The second is the html
        meta_data = list(itertools.islice(data, 1))[0]

        latlng = meta_data['latlng']
        lat, lon = map(float, latlng.split(','))

    return lat, lon



parser = argparse.ArgumentParser()
parser.add_argument('-d', dest='base_dir', help='Location of SWC site bootcamp directory')
parser.add_argument('-o', dest='out_dir', help='Directory to write output')

args = parser.parse_args()

bc_files = glob.glob(os.path.join(args.base_dir, 'bootcamps/*', 'index.html'))
if not bc_files:
    raise Exception('No bootcamp directories detected. Please check input directory')


output = os.path.join(args.out_dir, 'swc_bc_coords.csv')
with open(output, 'w') as f:
    f.write('# Latitude, Longitude\n')
    for bcf in bc_files:
        lat, lon = get_lat_lon(bcf)
        f.write('{:10.6f},{:10.6f}\n'.format(lat, lon))
