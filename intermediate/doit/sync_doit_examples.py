""" 
This script is intended to keep the example doit scripts
in doit_examples/ in sync with the contents of the
iPython notebooks used for teaching. It iterates over
all iPython notebooks in the current directory and looks
for cells that contain the doit magic. If the first
comment line contains a filename, it writes the contents
of that cell to the relevant file in doit_examples/
"""

import simplejson
import os
import glob

# Iterate over notebooks in this directory
for nbpath in glob.glob('0?-*.ipynb'):

    # Open notebook and load as json
    with open(nbpath, 'r') as nbtxt:
        nbdata = simplejson.load(nbtxt)
        
    # Iterate over cells
    for cell in nbdata['worksheets'][0]['cells']:
        
        # If a code cell, check if the first line starts with %%doit
        if cell['cell_type'] == 'code':
            lines = cell['input']
            if lines and lines[0][:6] == '%%doit':
                
                # If it does, find the first comment line and check that it looks like a filename
                for line in lines:
                    if line[0] == '#':
                    
                        if line[-4:-1] == '.py':

                            # Extract the filename
                            fname = line[1:].strip()
                            fpath = os.path.join('doit_examples', fname)

                            print 'Found an example. Writing to {0}'.format(fpath)
                            
                            # Write the contents of the cell to the filename
                            with open(fpath, 'w') as example_file:
                                example_file.writelines(lines[1:])
                            
                        break
