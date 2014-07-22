import sys
import numpy as np
from matplotlib import pyplot as plt

usage_string = """
Results.py plots the real overlap between two bed files
versus the overlap when one of the results files is randomly
shuffled.

Usage: python random.py shuffled_results.txt real_results.txt
"""

def main():

    if not len(sys.argv) == 3:
        sys.exit(usage_string)

    script = sys.argv[0]
    shuffled_results_file = sys.argv[1]
    real_results_file = sys.argv[2]

    process(shuffled_results_file, real_results_file)


def process(shuffled_results_file, real_results_file):
    shuffled_data = np.loadtxt(shuffled_results_file)
    real_data = np.loadtxt(real_results_file)

    print 'Real data:'
    print real_data
    print
    print 'Shuffled data:'
    print shuffled_data

    plt.hist(shuffled_data, color='black', label='Shuffled Data')
    plt.axvline(real_data, color='red', ls='--', lw=2, label='Real Data')
    plt.xlabel('Nucleotides overlap')
    plt.ylabel('Number of shuffled results')
    plt.legend()
    plt.show()

main()
