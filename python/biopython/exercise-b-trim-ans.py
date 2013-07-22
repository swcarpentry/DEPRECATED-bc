#!/usr/bin/env python
'''exercise-b-trim.py   Write the code to trim sequences of low-quality bases. 
You can expect to turn the data in data/tiny.fastq into the filtered data 
exactly like data/tiny.trimmed.fastq  '''

from Bio import SeqIO
import sys

def btrimmer(seqrecord):
    ''' This function takes a Seq object containing fastq data and returns a Seq 
    object with low-quality bases  (bases with quality scores of 2 and below) 
    removed from the end of the read'''
    i = len(seqrecord)-1
    while seqrecord.letter_annotations["phred_quality"][i] <= 2:
         i = i - 1
    choppedsequence = seq[0:i+1]    #  This does NOT do what you want
    return choppedsequence

#   This part opens a fastq file, goes through it record-by-record, calls btrimmer 
#   and writes fastq-formatted reuslts to standard out.
generator = SeqIO.parse("data/tiny.fastq", "fastq")
for fastqsequence in generator:
     choppedfastqsequence = btrimmer(fastqsequence)
     sys.stdout.write(choppedfastqsequence.format("fastq"))
