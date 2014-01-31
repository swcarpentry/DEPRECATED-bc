import sys
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
        number = 1
        width, filled = count(sys.stdin)
    else:
        number = len(filenames)
        with open(filenames[0], 'r') as source:
            width, filled = count(source)
        for f in filenames[1:]:
            new_width, new_filled = count(source)
            assert new_width == width, 'File widths are not the same'
            filled = combine(filled, new_filled)
    display(title, filled, number * width)

def count(source):
    result = []
    for line in source:
        line = line.strip()
        width = len(line)
        n = line.count('1')
        result.append(n)
    return width, result

def combine(left, right):
    assert len(left) == len(right), 'Data set lengths have unequal lengths'
    result = []
    for i in range(len(left)):
        result.append( left[i] + right[i] )
    return result

def display(title, counts, scaling):
    print title
    for c in counts:
        print float(c) / scaling

# Run the program.
main()
