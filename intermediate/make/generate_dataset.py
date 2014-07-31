import numpy as np
import pandas as pd
import string

"""
Description
------------

Generates mock datasets for make lesson

Run `python generate_dataset.py` and six datasets will be generated.

"""


NUM_SPP = 3  # Number of species in dataset
SAMPS = np.random.random_integers(5, 20, NUM_SPP)  # Number of samples
BASE_MEAN = 10
species = list(string.lowercase)[:NUM_SPP]

for i in xrange(NUM_SPP):

    spp_col = np.repeat(species[i], SAMPS[i])
    length = np.random.normal(BASE_MEAN + 2*i, 1, size=SAMPS[i])
    weight = np.random.normal(2 * length + 5, 1, size=SAMPS[i])

    tdata1 = pd.DataFrame(zip(length, spp_col), columns=["length", "species"])
    tdata1.to_csv("data-1-%i.dat" % (i + 1), index=False)

    tdata2 = pd.DataFrame(zip(weight, spp_col), columns=["weight", "species"])
    tdata2.to_csv("data-2-%i.dat" % (i + 1), index=False)







