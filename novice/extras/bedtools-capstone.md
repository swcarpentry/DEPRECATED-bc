% bedtools Tutorial
% Aaron Quinlan
% November 22, 2013


Synopsis
========

Our goal is to work through examples that demonstrate how to 
explore, process and manipulate genomic interval files (e.g., BED, VCF, BAM) with the `bedtools` software package.

Some of our analysis will be based upon the Maurano et al exploration of DnaseI hypersensitivity sites in hundreds of primary tissue types.

    Maurano et al. Systematic Localization of Common Disease-Associated Variation in Regulatory DNA. Science. 2012. Vol. 337 no. 6099 pp. 1190-1195.

    www.sciencemag.org/content/337/6099/1190.short


This tutorial is merely meant as an introduction to whet your appetite. There are many, many more tools and options than presented here. We therefore encourage you to read the bedtools [documentation](http://bedtools.readthedocs.org/en/latest/).




\


Setup
=====
From the Terminal, create a new directory on your Desktop called "bedtools-demo".

    cd ~/Desktop
    mkdir bedtools-demo

Navigate into that directory.

    cd bedtools-demo

Download the sample BED files I have provided.

    curl -O http://quinlanlab.cs.virginia.edu/cshl2013/maurano.dnaseI.tgz
    curl -O http://quinlanlab.cs.virginia.edu/cshl2013/cpg.bed
    curl -O http://quinlanlab.cs.virginia.edu/cshl2013/exons.bed
    curl -O http://quinlanlab.cs.virginia.edu/cshl2013/gwas.bed
    curl -O http://quinlanlab.cs.virginia.edu/cshl2013/genome.txt

Now, we need to extract all of the 20 Dnase I hypersensitivity BED files from the "tarball" named
`maurano.dnaseI.tgz`.

    tar -zxvf maurano.dnaseI.tgz
    rm maurano.dnaseI.tgz

Let's take a look at what files we now have.

    ls -1


\


What are these files?
=========================
Your directory should now contain 23 BED files and 1 genome file. Twenty of these files (those starting with "f" for "fetal tissue") reflect Dnase I hypersensitivity sites measured in twenty different fetal tissue samples from the brain, heart, intestine, kidney, lung, muscle, skin, and stomach.

In addition: `cpg.bed` represents CpG islands in the human genome; `exons.bed` represents RefSeq exons from human genes; and `gwas.bed` represents human disease-associated SNPs that were identified in genome-wide association studies (GWAS).

The latter 3 files were extracted from the UCSC Genome Browser's [Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables?command=start).


\


The bedtools help
==================
To bring up the help, just type

    bedtools

As you can see, there are multiple "subcommands" and for bedtools to
work you must tell it which subcommand you want to use. Examples:

    bedtools intersect
    bedtools merge
    bedtools subtract


bedtools "intersect"
====================

The `intersect` command is the workhorse of the `bedtools` suite. It compares two BED/VCF/GFF files (or a BAM file and one of the aforementioned files) and identifies all the regions in the gemome where the features in the two files overlap (that is, share at least one base pair in common).

![](http://bedtools.readthedocs.org/en/latest/_images/intersect-glyph.png)

Default behavior
----------------
By default, `intersect` reports the intervals that represent overlaps between your two files.  To demonstrate, let's identify all of the CpG islands that overlap exons.

    bedtools intersect -a cpg.bed -b exons.bed | head -5
    chr1    29320   29370   CpG:_116
    chr1    135124  135563  CpG:_30
    chr1    327790  328229  CpG:_29
    chr1    327790  328229  CpG:_29
    chr1    327790  328229  CpG:_29

Reporting the original feature in each file.
--------------------------------------------
The `-wa` (write A) and `-wb` (write B) options allow one to see the original records from the A and B files that overlapped.  As such, instead of not only showing you *where* the intersections occurred, it shows you *what* intersected.

    bedtools intersect -a cpg.bed -b exons.bed -wa -wb \
    | head -5
    chr1    28735   29810   CpG:_116    chr1    29320   29370   NR  _024540_exon_10_0_chr1_29321_r    0   -
    chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_ 039983_exon_0_0_chr1_134773_r    0   -
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028322_exon_2_0_chr1_324439_f    0   +
    chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_ 028327_exon_3_0_chr1_327036_f    0   +
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028325_exon_2_0_chr1_324439_f    0   +

How many base pairs of overlap were there?
------------------------------------------
The `-wo` (write overlap) option allows one to also report the *number* of base pairs of overlap between the features that overlap between each of the files.

    bedtools intersect -a cpg.bed -b exons.bed -wo \
    | head -10
    chr1    28735   29810   CpG:_116    chr1    29320   29370   NR  _024540_exon_10_0_chr1_29321_r    0   -   50
    chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_ 039983_exon_0_0_chr1_134773_r    0   -   439
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028322_exon_2_0_chr1_324439_f    0   +   439
    chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_ 028327_exon_3_0_chr1_327036_f    0   +   439
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028325_exon_2_0_chr1_324439_f    0   +   439
    chr1    713984  714547  CpG:_60 chr1    713663  714068  NR_ 033908_exon_6_0_chr1_713664_r    0   -   84
    chr1    762416  763445  CpG:_115    chr1    762970  763155  NR  _015368_exon_0_0_chr1_762971_f    0   +   185
    chr1    762416  763445  CpG:_115    chr1    763177  763229  NR  _047525_exon_0_0_chr1_763178_f    0   +   52
    chr1    762416  763445  CpG:_115    chr1    762970  763155  NR  _047524_exon_0_0_chr1_762971_f    0   +   185
    chr1    762416  763445  CpG:_115    chr1    762970  763155  NR_047523_exon_0_0_chr1_762971_f    0   +   185

Counting the number of overlapping features.
--------------------------------------------
We can also count, for each feature in the "A" file, the number of overlapping features in the "B" file. This is handled with the `-c` option.

    bedtools intersect -a cpg.bed -b exons.bed -c \
    | head
    chr1    28735   29810   CpG:_116    1
    chr1    135124  135563  CpG:_30 1
    chr1    327790  328229  CpG:_29 3
    chr1    437151  438164  CpG:_84 0
    chr1    449273  450544  CpG:_99 0
    chr1    533219  534114  CpG:_94 0
    chr1    544738  546649  CpG:_171    0
    chr1    713984  714547  CpG:_60 1
    chr1    762416  763445  CpG:_115    10
    chr1    788863  789211  CpG:_28 9

\


Find features that DO NOT overlap
--------------------------------------------
Often we want to identify those features in our A file that **do not** overlap features in the B file. The `-v` option is your friend in this case.

    bedtools intersect -a cpg.bed -b exons.bed -v \
    | head
    chr1    437151  438164  CpG:_84
    chr1    449273  450544  CpG:_99
    chr1    533219  534114  CpG:_94
    chr1    544738  546649  CpG:_171
    chr1    801975  802338  CpG:_24
    chr1    805198  805628  CpG:_50
    chr1    839694  840619  CpG:_83
    chr1    844299  845883  CpG:_153
    chr1    912869  913153  CpG:_28
    chr1    919726  919927  CpG:_15


Require a minimal fraction of overlap.
--------------------------------------------
Recall that the default is to report overlaps between features in A and B so long as *at least one basepair* of overlap exists. However, the `-f` option allows you to specify what fraction of each feature in A should be overlapped by a feature in B before it is reported.

Let's be more strict and require 50% of overlap.

    bedtools intersect -a cpg.bed -b exons.bed \
    -wo -f 0.50 \
    | head
    chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_ 039983_exon_0_0_chr1_134773_r    0   -   439
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028322_exon_2_0_chr1_324439_f    0   +   439
    chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_ 028327_exon_3_0_chr1_327036_f    0   +   439
    chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_ 028325_exon_2_0_chr1_324439_f    0   +   439
    chr1    788863  789211  CpG:_28 chr1    788770  794826  NR_ 047525_exon_4_0_chr1_788771_f    0   +   348
    chr1    788863  789211  CpG:_28 chr1    788770  794826  NR_ 047524_exon_3_0_chr1_788771_f    0   +   348
    chr1    788863  789211  CpG:_28 chr1    788770  794826  NR_ 047523_exon_3_0_chr1_788771_f    0   +   348
    chr1    788863  789211  CpG:_28 chr1    788858  794826  NR_ 047522_exon_5_0_chr1_788859_f    0   +   348
    chr1    788863  789211  CpG:_28 chr1    788770  794826  NR_ 047521_exon_4_0_chr1_788771_f    0   +   348
    chr1    788863  789211  CpG:_28 chr1    788858  794826  NR_ 047520_exon_6_0_chr1_788859_f    0   +   348


\


bedtools "merge"
====================
Many datasets of genomic features have many individual features that overlap one another (e.g. aligments from a ChiP seq experiment). It is often useful to just cobine the overlapping into a single, contiguous interval. The bedtools `merge` command will do this for you.

![](http://bedtools.readthedocs.org/en/latest/_images/merge-glyph.png)


Input must be sorted
--------------------
The merge tool requires that the input file is sorted by chromosome, then by start position. This allows the merging algorithm to work very quickly without requiring any RAM.

If you run the merge command on a file that is not sorted, you will get an error. For example:

    bedtools merge -i exons.bed
    chr1    66999824    67000051
    chr1    67091529    67091593
    chr1    67098752    67098777
    chr1    67101626    67101698
    chr1    67105459    67105516
    chr1    67108492    67108547
    chr1    67109226    67109402
    chr1    67126195    67126207
    chr1    67133212    67133224
    chr1    67136677    67136702
    chr1    67137626    67137678
    chr1    67138963    67139049
    chr1    67142686    67142779
    chr1    67145360    67145435
    chr1    67147551    67148052
    chr1    67154830    67154958
    chr1    67155872    67155999
    chr1    67161116    67161176
    chr1    67184976    67185088
    chr1    67194946    67195102
    chr1    67199430    67199563
    chr1    67205017    67205220
    chr1    67206340    67206405
    chr1    67206954    67207119
    ERROR: input file: (exons.bed) is not sorted by chrom then start.
       The start coordinate at line 26 is less than the start at line 25

To correct this, you need to sort your BED using the UNIX `sort` utility.

    sort -k1,1 -k2,2n exons.bed > exons.sort.bed

Let's try again.
    
    bedtools merge -i exons.sort.bed | head -10

The result is a new set of intervals representing the merged set of intervals in the input. That is, if a base pair in the genome is covered by 10 features, it will now only be represented once in the output file.



