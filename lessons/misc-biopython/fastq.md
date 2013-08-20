---
layout: lesson
root: ../..
title: FASTQ Sequence Parsing
---
FASTQ is the de-facto standard data format for the output of modern high-throughput sequencing machines.
In addition to sequence ID and header, this format includes a quality symbol for each base.

One way to parse fastq is using exactly the same ```SeqIO.parse()``` method, just with ```fastq``` instead of ```fasta``` as the format parameter.
The sequence is in the ```seq``` attribute and the quality scores (as ints) is in the letter_annotations["phred_quality"] attribute.  

Another approach is to use ```FastqGeneralIterator```.  Unlike ```SeqIO.parse()``` it takes file handles (not file names) and has no format parameter (it only works for fastq).  It returns tuples with the sequence description, the sequence string, and the quality string without additional methods to interpret and format the results.   The following code snippet opens the file tiny.fastq and writes truncated versions of the data to standard out.  (Note that if any of the input sequences are less than 30 base pairs in length, this code breaks.)

```python
from Bio.SeqIO.QualityIO import FastqGeneralIterator
in_handle = open("tiny.fastq")
iterator = FastqGeneralIterator(in_handle)
for triplet in iterator:
     (description, sequence, quality) = triplet
     print "@%s\n%s\n+\n%s" % ( description, sequence[0:30], quality[0:30] )
```

This approach is faster, but doesn't put the FASTQ-format-interpretation methods at your disposal, so you have to handle the ASCII-to-quality decoding yourself.
