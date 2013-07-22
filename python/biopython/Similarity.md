### Similarity searching ###
Similarity searching is perhaps the fundamental operation of computational biology; comparisons between known sequences and novel sequences are the bread-and-butter of sequence interpretation.  
You can run BLAST on your laptop, on your department's server, or via the NCBI web interface.  
BLAST against large databases is an expensive operation, so if your computational plan requires running BLAST a million times, you probably need to re-think your plan.  
For small numbers of sequences and for high-value sequences (contigs, genomes)  BLAST is extremely popular.

```python
from Bio.Blast import NCBIWWW
mysterysequence = "GCACTTGTCTCCTGTTTACTCCCCTGAGCTTGAGGGGTTAACATGAAGGTCATCGATAGCAGGATAATAATACAGTA"
blastresults = NCBIWWW.qblast("blastn", "nr", sequence=mysterysequence ) 
print type(blastresults)
```

```python
from Bio.Blast import NCBIXML
result_handle = open("my_blast.xml")
blast_record = NCBIXML.read(result_handle)
```

```python
from Bio.Blast.Applications import NcbiblastxCommandline
help(NcbiblastxCommandline)
```

