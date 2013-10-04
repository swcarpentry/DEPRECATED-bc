import sys
import numpy as np
from optparse import OptionParser

def main():
    script, flags, filenames = handle_args(sys.argv)
    if filenames == []:
        process('stdin', None)
    elif flags.merge:
        process('all', filenames)
    else:
        for f in filenames:
            process(f, [f])

def handle_args(args):
    script, rest = args[0], args[1:]
    parser = OptionParser()
    parser.add_option('-m', '--merge', dest='merge', help='Merge data from all files',
                      default=False, action='store_true')
    options, args = parser.parse_args(args=rest)
    return script, options, args

def process(title, filenames):
    if filenames is None:
        data = np.loadtxt(sys.stdin, delimiter=',')
        display(title, data, 1)
    else:
        results = [np.loadtxt(f, delimiter=',') for f in filenames]
        assert all([x.shape == results[0].shape for x in results]), 'File sizes differ'
        for r in results[1:]:
            results[0] += r
        display(title, results[0], len(results))

def display(title, data, number):
    print title
    scaling = float(number * data.shape[1])
    densities = data.sum(1) / scaling
    for d in densities:
        print d

# Run the program.
main()
