import os
import argparse
import glob
import itertools
import sys
import numpy as np

import yaml
from yaml import CLoader as YLoader

parser = argparse.ArgumentParser()
parser.add_argument('-d', dest='base_dir', help='Location of SWC site bootcamp directory')
parser.add_argument('-o', dest='out_dir', help='Directory to write output')

args = parser.parse_args()

bc_files = glob.glob(os.path.join(args.base_dir, 'bootcamps/*', 'index.html'))
if not bc_files:
    print 'Warning: No bootcamp directories detected. Please check input directory'
    sys.exit()

coord_data = []

for bcf in bc_files:
    with open(bcf, 'rt') as f:
        data = yaml.load_all(f, Loader=YLoader)

        # Only grab the first yaml document that contains the bootcamp metadata. The second is the html
        meta_data = list(itertools.islice(data, 1))[0]

        latlng = meta_data['latlng']
        coord_data.append(map(float, latlng.split(',')))

# Save coordinate data to csv file
header = 'Latitude, Longitude'
np.savetxt(os.path.join(args.out_dir, 'swc_bc_coords.csv'), coord_data, delimiter=',', fmt='%10.6f', header=header)
