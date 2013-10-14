#!/usr/bin/env python

from Bio import SeqIO
generator = SeqIO.parse("data/tiny.fastq", "fastq")
for sequence in generator:
     print sequence.id
     print sequence.seq
     print sequence.letter_annotations["phred_quality"]

