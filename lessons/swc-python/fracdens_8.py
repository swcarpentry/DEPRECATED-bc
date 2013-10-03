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
        width, filled = count(sys.stdin)
        display(title, filled, width)
    else:
        results = [count(f) for f in filenames]
        width, filled = results[0]
        assert all([r[0] == width for r in results]), 'File widths unequal'
        for a in results[1:]:
            filled += a[1]
        display(title, filled, len(filenames) * width)

def count(source):
    if type(source) is str:
        reader = open(source, 'r')
    else:
        reader = source
    cells = np.array([list(x.strip()) for x in reader.readlines()], np.int32)
    if reader is not source:
        reader.close()
    return cells.shape[1], cells.sum(1)

def display(title, counts, scaling):
    print title
    for c in counts:
        print float(c) / scaling

# Run the program.
main()
