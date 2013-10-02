import sys

def main():
    script, flags, filenames = handle_args(sys.argv)
    for f in filenames:
        process(f)

def handle_args(args):
    script, rest = args[0], args[1:]
    parser = OptionParser()
    parser.add_option('-m', '--merge', dest='merge', help='Merge data from all files',
                      default=False, action='store_true')
    options, args = parser.parse_args(args=rest)
    return script, options, args

def process(filename):
    print filename

# Run the program.
main()
