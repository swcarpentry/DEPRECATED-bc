import sys
import pandas as pd

"""
Description
------------
Simple stats.py script computes and saves as output the data summary for
some number of data sets.

"""

summary_name = sys.argv[1]
data_files = sys.argv[2:]

# Try so pattern_rule.mk works
try:
    data_files.remove("stats.py")
except ValueError:
    pass

spp_list = []
summary_list = []

# Summarize each data_file
for td in data_files:

    tdata = pd.read_csv(td)
    spp_list.append(tdata['species'].unique()[0])
    summary_list.append(tdata.describe().T)

summary = pd.concat(summary_list)
summary.index = spp_list
summary.index.name = "species"

summary.to_csv(summary_name)

