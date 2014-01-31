import sys
from optparse import OptionParser

def main():
    script, flags, filenames = handle_args(sys.argv)
    if flags.merge:
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
    print title
    print 'files:',
    for f in filenames:
        print f, # eventually replace this with real code
    print # make sure there's a newline at the end

# Run the program.
main()
