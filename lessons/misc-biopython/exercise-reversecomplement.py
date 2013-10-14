#!/usr/bin/env python
'''exercis-se_allsixframes.py   Write the code to translate sequences 
in all six frames.  ''' 

from Bio import SeqIO

def reversecomplement(record):
    ''' This function takes a SeqRecord object and returns its 
    reverse complement'''
#   your code goes here
    reversecomplementsequence = "N" * len(record.seq)
    return reversecomplementsequence 


#   Open a fastq file, goes through it record-by-record, and output
#   the sequence id, the sequence, and the translations 
generator = SeqIO.parse("data/test-sequences.fasta", "fasta")
for seqrecord in generator:
    reversesequence = reversecomplement(seqrecord)
    print ">%s\nORIG: %s" % (seqrecord.id, seqrecord.seq)
    print "REVC: %s" % reversesequence
