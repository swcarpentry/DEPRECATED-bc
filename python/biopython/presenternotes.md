Notes on new Biopython module,  SWC boot camp at Tufts June 3-4, 2013, W. Trimble+A. Ahmadia

The initial biopython module included several demonstration examples of "when are we going to use this," several follow-along code snippets, and exercises giving partially-functional code with unimplemented subroutines.

The demonstration examples were	
* download biological data form a HTTP request / REST API and plot it
* produce a composite of 240,000 image files (from which biological big data derive) (http://www.mcs.anl.gov/~trimble/nodi/cell-lg.mp4) (flashy pictures)
The exercises/ examples were
* retrieval of data from NCBI's EUTILS REST API (retrieve data programmatically from a working subroutine)
* iterating over lines of a file and retrieving data programmatically
* looping over records in a fasta-formatted file 
* doing something to the fasta data records -- reverse complimentation method

We didn't get to 
* doing something to the fasta data records -- translation method
* FASTQ parsing / btrimming example and 
* Genbank parsing / format conversion

Backup demonstration examples (not used / not developed yet)
* programmatic concatenation of data files according to a table associating multiple data files with each sample
* "the demultiplexing task" -- sorting data records according string matching of sample label fields with values like ATCACG and GCCAAT

Audience members suggested potentially interesting applications included
* sequence alignment example (BWA bindings)
* multiple sequence alignment example (MUSCLE bindings)
* expensive similarity-search output parsing (BLAST bindings / parsing)

The first exercise (a moderately useful subroutine that is just a binding to a NIH data-delivery REST interface) was well recieved. 

We (Will Trimble + Aron Ahmaidia) went after the second exercise (looping over a data file and running a bio subroutine on data from each line) in a blank ipython notebook.  It seemed to work despite a few hiccups.

The reverse-complimentation method didn't go over a smoothly, in part because the task was in a sense silly ("So you could have done this in one line, calling reverse_complement()?"), and partly it seemed because there wasn't adequate introduction to the syntax.

The students looked tired after a morning of git and an afternoon of biopython, ipython notebook, and biopython debugging.  For this reason, flashy examples (and code that does Neat Stuff) might be good to stack at the end of the session, by which time some of the students are burned out.

