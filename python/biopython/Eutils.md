
##What to expect##
We have to **get the data**, **get the data out of its container**, and **do something with the data**.  

###Get the data--NCBI's EUTILS###
One large data source is NCBI.  NCBI provides an interface to allow automated download of various data products using HTTP GET requests.

The documentation for this interface, called EFETCH, is here:
http://www.ncbi.nlm.nih.gov/books/NBK25499/

Before using EUTILS, we might want to know what kinds of things it can do:
#####Search engine for PUBMED: #####
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%22Life+with+6000+genes%22&retmax=100
#####Search engine for SRA#####
 http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=SRX015714
#####Search engine for Genbank#####
 http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=GFAJ1

Ok, you probably get the picture.  These queries return us lists of IDs that don't mean anything to us as humans, but that we can iterate over and retrieve automatically.

Once we've got lists of identifiers, we can retrieve data:
##### Pubmed abstracts:#####
This abstract http://www.ncbi.nlm.nih.gov/pubmed/8849441
can be retrieved in a machine-friendly format with this query:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=8849441&rettype=xml
##### SRA metadata bundles#####
This page of metadata about a dataset from Jeff Gordon's twin study  http://www.ncbi.nlm.nih.gov/sra/SRX015714 
can be retrieved from here:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=16175
##### Genbank records #####
Sequence data deposited by authors (with annotations, when provided by authors or the archive itself) can be retrieved in FASTA or genbank formats using the EUTILS suite.  The following URL should retrieve the human mitochondrial reference sequence (NC_012920) in genbank format:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=251831106&rettype=gb

##### Genome sequences #####
The following query will retrieve the genome of Candidatus Hodgkinia cicadicola (REFseq accession NC_012960.1) in fasta format:
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_012960.1&rettype=fasta

A table describing the supported formats is here:
http://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.chapter4_table1/?report=objectonly
Note that these queries return data in several different formats, some return XML-formatted data structures while
others return files in gb and fasta formats.  Of course, how you handle the data after retrieving
it is going to depend on the data format.

Biopython has subroutines that take EFETCH's options and parameters and return python objects.
This saves us from having to write code that directly talks to NCBI's EFETCH API, freeing us to spend our time elsewhere.  
We just need to find out how to use these subroutines.

Here is a example using the `Entrez.efetch`  biopython method to retrieve fasta sequences:
```python
#!/usr/bin/env python
import os,sys
from Bio import Entrez
from Bio import SeqIO

def downloadstuff(accessionno):
    filename = "%s.fasta" % accessionno        
    print "Trying efectch on %s, writing to %s" % ( accessionno, filename )
    if not os.path.isfile(filename):  
        net_handle = Entrez.efetch(db="nucleotide",id=accessionno,rettype="fasta", retmode="text")
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read() )
        out_handle.close()
        net_handle.close()
    else:
        print "skipping, %s already exists!" % filename

Entrez.email = "swc@example.com        #  Always tell NCBI who you are.  
Entrez.tool = "SoftwareCarpentryBootcamp"

accession = sys.argv[1]     # take the first program argument

downloadstuff(accession)   
```

