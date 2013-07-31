#!/usr/bin/env python
'''exercise-b-trim.py   Write the code to trim sequences of low-quality bases. 
You can expect to turn the data in data/tiny.fastq into the filtered data 
exactly like data/tiny.trimmed.fastq  '''

from Bio import SeqIO
import sys

def btrimmer(seqrecord):
    ''' This function takes a SeqRecord object containing fastq data and returns a Seq 
    object with low-quality bases  (bases with quality scores of 2 and below) 
    removed from the end of the read'''
#   your code goes here
    choppedsequence = seqrecord    #  This is a placeholder, it does not trim!
    return choppedsequence

#   This part opens a fastq file, goes through it record-by-record, calls btrimmer 
#   and writes fastq-formatted reuslts to standard out.
generator = SeqIO.parse("data/tiny.fastq", "fastq")
for fastqsequence in generator:
    choppedfastqsequence = btrimmer(fastqsequence)
    sys.stdout.write(choppedfastqsequence.format("fastq"))
