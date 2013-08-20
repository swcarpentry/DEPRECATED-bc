#!/usr/bin/env python
'''This is a skeleton python program.  It doesn't do anything, but it parses arguments and 
calls a subroutine.'''

def do_something(arg):
    sys.stderr.write("I'm doing something with argument %s!\n" % arg)
    return(42)

import sys, os
from optparse import OptionParser

if __name__ == '__main__':
    usage  = "usage: skeleton.py [<options>] [<arguments>]"
    parser = OptionParser(usage)
    parser.add_option("-i", "--input",  dest="infile", default=None, help="Input filename")
    parser.add_option("-o", "--output", dest="outfile", default=None, help="Output filename")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=True, help="Verbose")
    (opts, args) = parser.parse_args()
#    if len(args) == 0 :
#        sys.exit("No arguments supplied!")
#    argument = args[0]
    inputfile = opts.infile

#    if opts.verbose: 
#        sys.stdout.write("Argument: %s\n" % argument ) 

    do_something(inputfile)

    if opts.verbose: 
        sys.stdout.write("Done. \n")
