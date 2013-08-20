#!/usr/bin/env python
'''tiny-fastq-parser.py
opens an example FASTQ file and dumps the data fields'''

from Bio import SeqIO
generator = SeqIO.parse("data/tiny.fastq", "fastq")
for sequence in generator:
    print sequence.id
    print sequence.seq
    print sequence.letter_annotations["phred_quality"]

