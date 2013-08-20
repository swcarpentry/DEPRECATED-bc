#!/usr/bin/env python
'''This script sends a efetch request to NCBI, requesting a genbank-formatted data file
and creates a .gbk file if successful
retrieve.py NC_000913  
should create NC_000913.gbk containing the annotated E. coli K12 reference genome '''

import os, sys
from Bio import Entrez

def downloadgbk(accessionno):
    filename = "%s.gbk" % accessionno       
    print "Trying efectch on %s, writing to %s" % ( accessionno, filename )
    if not os.path.isfile(filename):  
        net_handle = Entrez.efetch(db="nucleotide", id=accessionno, rettype="gb", retmode="text") 
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read() ) 
        out_handle.close()
        net_handle.close()
    else:
        print "skipping, %s already exists!" % filename

Entrez.email = "trimble@anl.gov"
Entrez.tool = "SoftwareCarpentryBootcamp"

if len(sys.argv) != 2:
    sys.exit("Usage: retrieve.py <accession number>")
accession = sys.argv[1]     # take the first program argument

downloadgbk(accession)
