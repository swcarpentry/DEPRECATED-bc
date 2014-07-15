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


Count the number of overlapping intervals.
------------------------------------------
A more sophisticated approach would be to not only merge overlapping intervals, but also report the *number* of intervals that were integrated into the new, merged interval. One does this with the `-n` option.

    bedtools merge -i exons.sort.bed -n | head -20
    chr1    11873   12227   1
    chr1    12612   12721   1
    chr1    13220   14829   2
    chr1    14969   15038   1
    chr1    15795   15947   1
    chr1    16606   16765   1
    chr1    16857   17055   1
    chr1    17232   17368   1
    chr1    17605   17742   1
    chr1    17914   18061   1
    chr1    18267   18366   1
    chr1    24737   24891   1
    chr1    29320   29370   1
    chr1    34610   35174   2
    chr1    35276   35481   2
    chr1    35720   36081   2
    chr1    69090   70008   1
    chr1    134772  139696  1
    chr1    139789  139847  1
    chr1    140074  140566  1

Merging features that are close to one another.
-----------------------------------------------
With the `-d` (distance) option, one can also merge intervals that do not overlap, yet are close to one another. For example, to merge features that are no more than 1000bp apart, one would run:

    bedtools merge -i exons.sort.bed -d 1000 -n | head -20
    chr1    11873   18366   12
    chr1    24737   24891   1
    chr1    29320   29370   1
    chr1    34610   36081   6
    chr1    69090   70008   1
    chr1    134772  140566  3
    chr1    323891  328581  10
    chr1    367658  368597  3
    chr1    621095  622034  3
    chr1    661138  665731  3
    chr1    700244  700627  1
    chr1    701708  701767  1
    chr1    703927  705092  2
    chr1    708355  708487  1
    chr1    709550  709660  1
    chr1    713663  714068  1
    chr1    752750  755214  2
    chr1    761585  763229  10
    chr1    764382  764484  9
    chr1    776579  778984  1


\


bedtools "complement"
=====================
We often want to know which intervals of the genome are **NOT** "covered" by intervals in a given feature file. For example, if you have a set of ChIP-seq peaks, you may also want to know which regions of the genome are not bound by the factor you assayed. The `complement` addresses this task.

![](http://bedtools.readthedocs.org/en/latest/_images/complement-glyph.png)

As an example, let's find all of the non-exonic (i.e., intronic or intergenic) regions of the genome.  Note, to do this you need a ["genome"](http://bedtools.readthedocs.org/en/latest/content/general-usage.html#genome-file-format) file, which tells `bedtools` the length of each chromosome in your file.  *Consider why the tool would need this information...*

    bedtools complement -i exons.bed -g genome.txt \
    > non-exonic.bed
    head non-exonic.bed
    chr1    0   11873
    chr1    12227   12612
    chr1    12721   13220
    chr1    14829   14969
    chr1    15038   15795
    chr1    15947   16606
    chr1    16765   16857
    chr1    17055   17232
    chr1    17368   17605
    chr1    17742   17914


\


bedtools "genomecov"
====================
For many analyses, one wants to measure the genome wide coverage of a feature file. For example, we often want to know what fraction of the genome is covered by 1 feature, 2 features, 3 features, etc. This is frequently crucial when assessing the "uniformity" of coverage from whole-genome sequencing. This is done with the versatile `genomecov` tool.

![](http://bedtools.readthedocs.org/en/latest/_images/genomecov-glyph.png)

As an example, let's produce a histogram of coverage of the exons throughout the genome. Like the `merge` tool, `genomecov` requires pre-sorted data. It also needs a genome file as above.

    bedtools genomecov -i exons.sort.bed -g genome.txt

This should run for 3 minutes or so. At the end of your output, you should see something like:

    genome  0   3062406951  3137161264  0.976171
    genome  1   44120515    3137161264  0.0140638
    genome  2   15076446    3137161264  0.00480576
    genome  3   7294047 3137161264  0.00232505
    genome  4   3650324 3137161264  0.00116358
    genome  5   1926397 3137161264  0.000614057
    genome  6   1182623 3137161264  0.000376972
    genome  7   574102  3137161264  0.000183
    genome  8   353352  3137161264  0.000112634
    genome  9   152653  3137161264  4.86596e-05
    genome  10  113362  3137161264  3.61352e-05
    genome  11  57361   3137161264  1.82844e-05
    genome  12  52000   3137161264  1.65755e-05
    genome  13  55368   3137161264  1.76491e-05
    genome  14  19218   3137161264  6.12592e-06
    genome  15  19369   3137161264  6.17405e-06
    genome  16  26651   3137161264  8.49526e-06
    genome  17  9942    3137161264  3.16911e-06
    genome  18  13442   3137161264  4.28477e-06
    genome  19  1030    3137161264  3.28322e-07
    genome  20  6329    3137161264  2.01743e-06
    ...


\


Producing BEDGRAPH output
--------------------------
Using the `-bg` option, one can also produce BEDGRAPH output which represents the "depth" fo feature coverage for each base pair in the genome:

    bedtools genomecov -i exons.sort.bed -g genome.txt -bg | head -20
    chr1    11873   12227   1
    chr1    12612   12721   1
    chr1    13220   14361   1
    chr1    14361   14409   2
    chr1    14409   14829   1
    chr1    14969   15038   1
    chr1    15795   15947   1
    chr1    16606   16765   1
    chr1    16857   17055   1
    chr1    17232   17368   1
    chr1    17605   17742   1
    chr1    17914   18061   1
    chr1    18267   18366   1
    chr1    24737   24891   1
    chr1    29320   29370   1
    chr1    34610   35174   2
    chr1    35276   35481   2
    chr1    35720   36081   2
    chr1    69090   70008   1
    chr1    134772  139696  1

    
\


Sophistication through chaining multiple bedtools
=================================================
Analytical power in `bedtools` comes from the ability to "chain" together multiple tools in order to construct rather sophisicated analyses with very little programming - you just need **genome arithmetic**!  Have a look at the examples [here](http://bedtools.readthedocs.org/en/latest/content/advanced-usage.html).


\


Principal component analysis
=============================

We will use the bedtools implementation of a Jaccard statistic to meaure the similarity of two 
datasets. Briefly, the Jaccard statistic measures the ratio of the number of *intersecting* base 
pairs to the *total* number of base pairs in the two sets.  As such, the score ranges from 0.0 to 1.
0; lower values reflect lower similarity, whereas higher values reflect higher similarity.

Let's walk through an example: we would expect the Dnase hypersensivity sites to be rather similar 
between two samples of the **same** fetal tissue type.  Let's test:

    bedtools jaccard \
        -a fHeart-DS16621.hotspot.twopass.fdr0.05.merge.bed \
        -b fHeart-DS15839.hotspot.twopass.fdr0.05.merge.bed
    intersection    union   jaccard
    81269248    160493950   0.50637

But what about the similarity of two **different** tissue types?

    bedtools jaccard \
        -a fHeart-DS16621.hotspot.twopass.fdr0.05.merge.bed \
        -b fSkin_fibro_bicep_R-DS19745.hg19.hotspot.twopass.fdr0.05.merge.bed
    intersection    union   jaccard
    28076951    164197278   0.170995

Hopefully this demonstrates how the Jaccard statistic can be used as a simple statistic to reduce 
the dimensionality of the comparison between two large (e.g., often containing thousands or 
millions of intervals) feature sets.


\


A Jaccard statistic for all 400 pairwise comparisons.
------------------------------------------------------


We are going to take this a bit further and use the Jaccard statistic to measure the similarity of 
all 20 tissue samples against all other 20 samples.  Once we have a 20x20 matrix of similarities, 
we can use dimensionality reduction techniques such as hierarchical clustering or principal 
component analysis to detect higher order similarities among **all** of the datasets.


We will use GNU parallel to compute a Jaccard statistic for the 400 (20*20) pairwise comparisons 
among the fetal tissue samples.

But first, we need to install [GNU parallel](http://www.gnu.org/software/parallel/).

    brew install parallel

Next, we need to install a tiny script I wrote for this analysis.

    curl -O http://quinlanlab.cs.virginia.edu/cshl2013/make-matrix.py


Now, we can use `parallel` to, you guessed it, compute the 400 pairwise Jaccard statistics in parallel using as many processors as you have available.

    parallel "bedtools jaccard -a {1} -b {2} \
             | awk 'NR>1' \
             | cut -f 3 \
             > {1}.{2}.jaccard" \
             ::: `ls *.merge.bed` ::: `ls *.merge.bed`

This command will create a single file containing the pairwise Jaccard measurements from all 400 tests.

    find . \
        | grep jaccard \
        | xargs grep "" \
        | sed -e s"/\.\///" \
        | perl -pi -e "s/.bed./.bed\t/" \
        | perl -pi -e "s/.jaccard:/\t/" \
        > pairwise.dnase.txt

A bit of cleanup to use more intelligible names for each of the samples.

    cat pairwise.dnase.txt \
    | sed -e 's/.hotspot.twopass.fdr0.05.merge.bed//g' \
    | sed -e 's/.hg19//g' \
    > pairwise.dnase.shortnames.txt
 
Now let's make a 20x20 matrix of the Jaccard statistic. This will allow the data to play nicely with R.

    awk 'NF==3' pairwise.dnase.shortnames.txt \
    | awk '$1 ~ /^f/ && $2 ~ /^f/' \
    | python make-matrix.py \
    > dnase.shortnames.distance.matrix
 
Let's also make a file of labels for each dataset so that we can label each dataset in our R plot.

    cut -f 1 dnase.shortnames.distance.matrix | cut -f 1 -d "-" | cut -f 1 -d "_" > labels.txt
 
Now start up R. (This assumes you have installed the `ggplot2` package).

    R

You should see something very similar to this:


    R version 2.15.1 (2012-06-22) -- "Roasted Marshmallows"
    Copyright (C) 2012 The R Foundation for Statistical Computing
    ISBN 3-900051-07-0
    Platform: x86_64-apple-darwin12.0.0 (64-bit)
    
    R is free software and comes with ABSOLUTELY NO WARRANTY.
    You are welcome to redistribute it under certain conditions.
    Type 'license()' or 'licence()' for distribution details.
    
      Natural language support but running in an English locale
    
    R is a collaborative project with many contributors.
    Type 'contributors()' for more information and
    'citation()' on how to cite R or R packages in publications.
    
    Type 'demo()' for some demos, 'help()' for on-line help, or
    'help.start()' for an HTML browser interface to help.
    Type 'q()' to quit R.
    
    >

No paste these commands into the R console:

    library(ggplot2)
    library(RColorBrewer)
    blues <- colorRampPalette(c('dark blue', 'light blue'))
    greens <- colorRampPalette(c('dark green', 'light green'))
    reds <- colorRampPalette(c('pink', 'dark red'))
     
    setwd("~/Desktop/bedtools-demo")
    x <- read.table('dnase.shortnames.distance.matrix')
    labels <- read.table('labels.txt')
    ngroups <- length(unique(labels))
    pca <- princomp(x)
    qplot(pca$scores[,1], pca$scores[,2], color=labels[,1],     geom="point", size=1) +
      scale_color_manual(values = c(blues(4), greens(5), reds(5))) 

You should see this:

![](http://quinlanlab.cs.virginia.edu/cshl2013/pca.png)

Et voila. Note that PCA was used in this case as an example of what PCA does for the CSHL Adv. Seq. course. Heatmaps are a more informative visualization in this case.



