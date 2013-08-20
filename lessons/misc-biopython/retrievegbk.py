#!/usr/bin/env python
'''This script sends a efetch request to NCBI, requesting a genbank-formatted data file
and creates a .gbk file if successful
retrieve.py NC_000913  
should create NC_000913.gbk containing the annotated E. coli K12 reference genome '''

import os, sys

def downloadgbk(accessionno):
    from Bio import Entrez
    Entrez.email = "swc@example.com    # Tell NCBI who you are!
    Entrez.tool = "SoftwareCarpentryBootcamp"

    filename = "%s.gbk" % accessionno       
    print "Trying efectch on %s, writing to %s" % ( accessionno, filename )
    if not os.path.isfile(filename):  
        net_handle = Entrez.efetch(db="nucleotide", id=accessionno, rettype="gb", retmode="text") 
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read()) 
        net_handle.close()
        out_handle.close()
    else:
        print "skipping, %s already exists!" % filename

def main():
    if len(sys.argv) != 2:    # check that exactly one argument was suppplied 
        sys.exit("Usage: retrieve.py <accession number>")
    accession = sys.argv[1]     # assign the first argument to accession
    downloadgbk(accession)      # call the subroutine downloadgbk

if __name__ == '__main__':
    main()
