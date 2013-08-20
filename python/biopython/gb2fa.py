#!/usr/bin/env python
'''Converts genbank files to sequence-only nucleic acid fasta'''

from Bio import SeqIO
import sys

outputformat = "fasta"
if len(sys.argv) != 3 :
    print "wrong number of args"
    print "Usage: gb2fa.py <input.gbk> <output.fasta>"
    sys.exit()
print "Converting %s to %s " % (sys.argv[1], sys.argv[2])
generator = SeqIO.parse(sys.argv[1], "genbank")
outfile = open(sys.argv[2], "w")
for record in generator:
    outfile.write(record.format(outputformat))
