import numpy as np
import pandas as pd
# Ensure matplotlib doesn't try to open windows
# From http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import time

""" Makes a bar plot based on summary data """


# Command line arguments
figure_name = sys.argv[1]
summary_data = sys.argv[2]

# Read in data
tdata = pd.read_csv(summary_data)

# Set up values for bar plot
species = tdata[tdata.columns[0]]
x_pos = np.arange(len(species))
value = tdata['mean']
SE = tdata['std'] / np.sqrt(tdata['count'])

# Make and save plot
plt.bar(x_pos, value, yerr=SE, align='center', alpha=0.4)
plt.xticks(x_pos, species)
plt.ylabel('Value')
plt.xlabel('Species')
plt.title('Plot from %s on %s at %s' % (summary_data,
                                        time.strftime("%m/%d/%Y"),
                                        time.strftime("%H:%M:%S")))
plt.savefig(figure_name)


