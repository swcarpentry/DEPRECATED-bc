---
layout: lesson
root: ../../..
title: Capstone example Python and Bedtools
---

Synopsis
========

Our goal is to work through examples that demonstrate how to explore, process and manipulate bed files with the `bedtools` software package.

Bed files are commonly used to represent features in genomes, for example promotors, genes, centromeres etc. - essentially anything that can be identified with a particular genomic location. If you are working with published data from other groups (for example, transcription factor binding sites identified by ChIP-seq) you will likely find that the data is provided in bed format. One of the most commonly asked questions of bioinformaticians is whether two different sets of features overlap with each other, so today we'll be showing you how to answer that question.

This tutorial is meant as an introduction to whet your appetite. There are many, many more tools and options than presented here. We therefore encourage you to read the bedtools [documentation](http://bedtools.readthedocs.org/en/latest/).

During this lesson, you'll be working on some scripts and data which are in their own github repository. In our git lesson, we covered the situation where you have been given write access to a repository. Now we're going to cover what to do if you don't have write access to a repository - this is called forking.

Forking the github repository
=============================

The data and code that you need to work on are located at <http://github.com/rbeagrie/bedtools-example>. You are going to want to make changes to the code, but this repository is owned by someone else. Instead of creating a new project, you want to [fork](../../gloss.html#fork) it; i.e., clone it on GitHub. You can do this using the GitHub web interface:

<img src="img/git-fork-ui.png" alt="The Fork Button" />

If you visit the [bedtools-example repository](http://github.com/rbeagrie/bedtools-example) and click the "fork" button, you will have your own copy of the repository under your own github account - so there should now be a page that looks like https://github.com/YOUR_USERNAME/bedtools-example.

Setup
=====

From the Terminal, clone your fork of the repository into a new folder on your Desktop called "bedtools-example".

~~~
$ cd ~/Desktop
$ git clone git@github.com:YOUR_USERNAME_GOES_HERE/bedtools-example.git
~~~
{:class="in"}


Navigate into that directory.

~~~
$ cd bedtools-example
~~~
{:class="in"}

Let's take a look at what files we now have.

~~~
$ ls
~~~
{:class="in"}

~~~
data scripts
~~~
{:class="out"}

So we have two directories, one containing data and the other containing some scripts. We're going to start with the data folder, so lets look at what's in there.

~~~
$ ls data
~~~
{:class="in"}

~~~
cpg.bed exons.bed genome.txt gwas.bed
~~~
{:class="out"}


What are these files?
---------------------

The data directory should contain 3 BED files and 1 genome file. 

`data/cpg.bed` represents CpG islands in the human genome; `data/exons.bed` represents RefSeq exons from human genes; and `data/gwas.bed` represents human disease-associated SNPs that were identified in genome-wide association studies (GWAS). `data/genome.txt` file contains the lengths of all the human chromosomes.

These 4 files were extracted from the UCSC Genome Browser's [Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables?command=start).

Introducing BedTools
====================

All bed files contain three columns which list the chromosome, the starting coordinate and the end co-ordinate. They may also include some other information, such as a name for each region, but the first three columns are always the same. For example, have a look at `data/cpg.bed`:

~~~
head data/cpg.bed 
~~~
{:class="in"}
~~~
chr1    28735   29810   CpG:_116
chr1    135124  135563  CpG:_30
chr1    327790  328229  CpG:_29
chr1    437151  438164  CpG:_84
chr1    449273  450544  CpG:_99
chr1    533219  534114  CpG:_94
chr1    544738  546649  CpG:_171
chr1    713984  714547  CpG:_60
chr1    762416  763445  CpG:_115
chr1    788863  789211  CpG:_28
~~~
{:class="out"}

The fourth column here is the name of each CpG island - which in this case is just a numerical ID assigned by the folks at UCSC.

The bedtools help
-----------------

You should have installed bedtools as part of the setup for this bootcamp. In which case, you can bring up the help by typing:

~~~
$ bedtools
~~~
{:class="in"}

If that didn't work, you probably don't have bedtools installed properly. Have a look at <http://bedtools.readthedocs.org/en/latest/content/installation.html> for some help.

If it did work, you should see that there are multiple "subcommands". For bedtools to work you must tell it which subcommand you want to use. Examples:

~~~
$ bedtools intersect
$ bedtools merge
$ bedtools subtract
~~~
{:class="in"}

This should bring up the specific help text associated with each tool.

bedtools "intersect"
====================

The `intersect` command is the workhorse of the `bedtools` suite. It compares two bed files and identifies all the regions in the genome where the features in the two files overlap (that is, share at least one base pair in common).

![](http://bedtools.readthedocs.org/en/latest/_images/intersect-glyph.png)

Default behavior
----------------
By default, `intersect` reports the intervals that represent overlaps between your two files. This is represented by the top line of the image above (the one which reads "A intersect B"). As an example, you may know that CpG islands often mark transcription start sites. Imagine that we want to find all the exons in the genome which might contain a promotor. One way of doing this would be to identify all of the CpG islands that overlap exons. For example, let's list the first five overlaps between CpG islands and exons.

~~~
$ bedtools intersect -a data/cpg.bed -b data/exons.bed | head
~~~
{:class="in"}
~~~
chr1    29320   29370   CpG:_116
chr1    135124  135563  CpG:_30
chr1    327790  328229  CpG:_29
chr1    327790  328229  CpG:_29
chr1    327790  328229  CpG:_29
chr1    713984  714068  CpG:_60
chr1    762970  763155  CpG:_115
chr1    763177  763229  CpG:_115
chr1    762970  763155  CpG:_115
chr1    762970  763155  CpG:_115
~~~
{:class="out"}

Notice how the output lists the name of the CpG island which is overlapping, but not the name of the exon. If we switch the order of files, we get all the extra columns from the exons.bed file, including the name (fourth column). The fifth column is a score (not used in this context) and the sixth is the strand.

~~~
$ bedtools intersect -a data/exons.bed -b data/cpg.bed | head -5 
~~~
{:class="in"}
~~~
chr1  50489434  50489626  NM_032785_exon_13_0_chr1_50489435_r    0  -
chr1  16767166  16767348  NM_018090_exon_0_0_chr1_16767167_f     0  +
chr1  33546713  33546895  NM_052998_exon_0_0_chr1_33546714_f     0  +
chr1  33546988  33547109  NM_052998_exon_1_0_chr1_33546989_f     0  +
chr1  33547201  33547243  NM_052998_exon_2_0_chr1_33547202_f     0  +
chr1  16767166  16767270  NM_001145278_exon_0_0_chr1_16767167_f  0  +
chr1  16767166  16767348  NM_001145277_exon_0_0_chr1_16767167_f  0  +
chr1  8384389   8384719   NM_001080397_exon_0_0_chr1_8384390_f   0  +
chr1  8385895   8386102   NM_001080397_exon_2_0_chr1_8385878_f   0  +
chr1  8390268   8390800   NM_001080397_exon_3_0_chr1_8390269_f   0  +
~~~
{:class="out"}


Reporting the original feature in each file.
--------------------------------------------
The `-wa` (write A) and `-wb` (write B) options allow one to see the original records from the A and B files that overlapped.  As such, instead of only showing you *where* the intersections occurred, it shows you *what* intersected. This is represented by the second line of the image above (the one which reads "A intersect B (-wa)"). A

~~~
$ bedtools intersect -a cpg.bed -b exons.bed -wa -wb | head -5
~~~
{:class="in"}
~~~
chr1    28735   29810   CpG:_116    chr1    29320   29370   NR_024540_exon_10_0_chr1_29321_r    0   -
chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_039983_exon_0_0_chr1_134773_r    0   -
chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_028322_exon_2_0_chr1_324439_f    0   +
chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_028327_exon_3_0_chr1_327036_f    0   +
chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_028325_exon_2_0_chr1_324439_f    0   +
~~~
{:class="out"}

The first three columns of the output show us where there are overlaps between exons and cpg islands. The fourth column contains the name of the CpG island which is overlapping. The rest of the output is the full line from the exons.bed file.

Counting the number of overlapping features.
--------------------------------------------
We can also count, for each feature in the "A" file, the number of overlapping features in the "B" file. This is handled with the `-c` option.

~~~
$ bedtools intersect -a data/cpg.bed -b data/exons.bed -c  | head
~~~
{:class="in"}
~~~
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
~~~
{:class="out"}

Again, we see that the order of the bed files is important. We get an output line for every feature in `data/cpg.bed`, even if there is no overlap with `data/exons.bed`. If we reverse the order of the two files again, we get a different output:

~~~
$ bedtools intersect -a data/exons.bed -b data/cpg.bed -c  | head
~~~
{:class="in"}
~~~
chr1    66999824    67000051    NM_032291_exon_0_0_chr1_66999825_f  0   +   0
chr1    67091529    67091593    NM_032291_exon_1_0_chr1_67091530_f  0   +   0
chr1    67098752    67098777    NM_032291_exon_2_0_chr1_67098753_f  0   +   0
chr1    67101626    67101698    NM_032291_exon_3_0_chr1_67101627_f  0   +   0
chr1    67105459    67105516    NM_032291_exon_4_0_chr1_67105460_f  0   +   0
chr1    67108492    67108547    NM_032291_exon_5_0_chr1_67108493_f  0   +   0
chr1    67109226    67109402    NM_032291_exon_6_0_chr1_67109227_f  0   +   0
chr1    67126195    67126207    NM_032291_exon_7_0_chr1_67126196_f  0   +   0
chr1    67133212    67133224    NM_032291_exon_8_0_chr1_67133213_f  0   +   0
chr1    67136677    67136702    NM_032291_exon_9_0_chr1_67136678_f  0   +   0
~~~
{:class="out"}

Now we get one line for every exon. In this case, the first ten exons actually don't overlap with any CpG islands.

Require a minimal fraction of overlap.
--------------------------------------------
Recall that the default is to report overlaps between features in A and B so long as *at least one basepair* of overlap exists. However, the `-f` option allows you to specify what fraction of each feature in A should be overlapped by a feature in B before it is reported.

Let's be more strict and require 50% of overlap.

~~~
$ bedtools intersect -a data/cpg.bed -b data/exons.bed -c -f 0.50 | head
~~~
{:class="in"}
~~~
chr1    28735   29810   CpG:_116    0
chr1    135124  135563  CpG:_30 1
chr1    327790  328229  CpG:_29 3
chr1    437151  438164  CpG:_84 0
chr1    449273  450544  CpG:_99 0
chr1    533219  534114  CpG:_94 0
chr1    544738  546649  CpG:_171    0
chr1    713984  714547  CpG:_60 0
chr1    762416  763445  CpG:_115    0
chr1    788863  789211  CpG:_28 7
~~~
{:class="out"}

If you compare this to the result you see without the -f flag you will see that the very first CpG island is now reported as having 0 overlaps. This means it must have overlapped an exon by less than 50%.

Find features that DO NOT overlap
--------------------------------------------
Often we want to identify those features in our A file that **do not** overlap features in the B file. The `-v` option is your friend in this case. This is represented in the last line of the image above.

~~~
$ bedtools intersect -a cpg.bed -b exons.bed -v  | head
~~~
{:class="in"}
~~~
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
~~~
{:class="out"}

How many base pairs of overlap were there?
------------------------------------------
The `-wo` (write overlap) option allows one to also report the *number* of base pairs of overlap between the features that overlap between each of the files.

~~~
$ bedtools intersect -a cpg.bed -b exons.bed -wo  | head -10
~~~
{:class="in"}
~~~
chr1    28735   29810   CpG:_116    chr1    29320   29370   NR_024540_exon_10_0_chr1_29321_r    0   -   50
chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_039983_exon_0_0_chr1_134773_r    0   -   439
chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_028322_exon_2_0_chr1_324439_f    0   +   439
chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_028327_exon_3_0_chr1_327036_f    0   +   439
chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_028325_exon_2_0_chr1_324439_f    0   +   439
chr1    713984  714547  CpG:_60 chr1    713663  714068  NR_033908_exon_6_0_chr1_713664_r    0   -   84
chr1    762416  763445  CpG:_115    chr1    762970  763155  NR_015368_exon_0_0_chr1_762971_f    0   +   185
chr1    762416  763445  CpG:_115    chr1    763177  763229  NR_047525_exon_0_0_chr1_763178_f    0   +   52
chr1    762416  763445  CpG:_115    chr1    762970  763155  NR_047524_exon_0_0_chr1_762971_f    0   +   185
chr1    762416  763445  CpG:_115    chr1    762970  763155  NR_047523_exon_0_0_chr1_762971_f    0   +   185
~~~
{:class="out"}

Do exons and CpG islands overlap significantly?
===============================================

OK, now we've explored the capabilities of BedTools a little, let's use it to answer a biological question. Let's try to answer if CpG islands and exons overlap significantly - that is, do they overlap more than we would expect by random chance?

Going back to base pairs overlap
--------------------------------

In order to start answering this question, we need to get the base pairs of overlap between all of the elements. From the section above, we know we can get this using the -wo flag:

~~~
$ bedtools intersect -a cpg.bed -b exons.bed -wo | head -10
~~~
{:class="in"}
~~~
chr1  28735   29810   CpG:_116    chr1    29320   29370   NR_024540_exon_10_0_chr1_29321_r    0   -   50
chr1    135124  135563  CpG:_30 chr1    134772  139696  NR_039983_exon_0_0_chr1_134773_r    0   -   439
chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_028322_exon_2_0_chr1_324439_f    0   +   439
chr1    327790  328229  CpG:_29 chr1    327035  328581  NR_028327_exon_3_0_chr1_327036_f    0   +   439
chr1    327790  328229  CpG:_29 chr1    324438  328581  NR_028325_exon_2_0_chr1_324439_f    0   +   439
chr1    713984  714547  CpG:_60 chr1    713663  714068  NR_033908_exon_6_0_chr1_713664_r    0   -   84
chr1    762416  763445  CpG:_115    chr1    762970  763155  NR_015368_exon_0_0_chr1_762971_f    0   +   185
chr1    762416  763445  CpG:_115    chr1    763177  763229  NR_047525_exon_0_0_chr1_763178_f    0   +   52
chr1    762416  763445  CpG:_115    chr1    762970  763155  NR_047524_exon_0_0_chr1_762971_f    0   +   185
chr1    762416  763445  CpG:_115    chr1    762970  763155  NR_047523_exon_0_0_chr1_762971_f    0   +   185
~~~
{:class="out"}

Can you spot anything odd about this output? Lines 3-5 indicate that one of our CpG islands is overlapping with three alternative versions of the same exon. For the purposes of our analysis though, this should really count as a single overlap. What we need to do is merge overlapping exons in exons.bed together.

bedtools "merge"
====================

Many datasets of genomic features have many individual features that overlap one another (e.g. alignments from a ChiP seq experiment). It is often useful to just combine the overlapping into a single, contiguous interval. The bedtools `merge` command will do this for you.

![](http://bedtools.readthedocs.org/en/latest/_images/merge-glyph.png)


Input must be sorted
--------------------

The merge tool requires that the input file is sorted by chromosome, then by start position. This allows the merging algorithm to work very quickly without requiring any RAM.

If you run the merge command on a file that is not sorted, you will get an error. For example:

~~~
bedtools merge -i exons.bed
~~~
{:class="in"}
~~~
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
~~~
{:class="in"}
~~~
ERROR: input file: (exons.bed) is not sorted by chrom then start.
   The start coordinate at line 26 is less than the start at line 25
~~~
{:class="err"}

To correct this, you need to sort your BED. We could use the UNIX `sort` utility that we already covered, which would look like this:

~~~
$ sort -k1,1 -k2,2n exons.bed | head
~~~
{:class="in"}
~~~
chr1    11873   12227   NR_046018_exon_0_0_chr1_11874_f 0   +
chr1    12612   12721   NR_046018_exon_1_0_chr1_12613_f 0   +
chr1    13220   14409   NR_046018_exon_2_0_chr1_13221_f 0   +
chr1    14361   14829   NR_024540_exon_0_0_chr1_14362_r 0   -
chr1    14969   15038   NR_024540_exon_1_0_chr1_14970_r 0   -
~~~
{:class="out"}

But if you can't remember the right arguments to use for sort, the command `bedtools sort` will also do what we want:

~~~
$ bedtools sort -i exons.bed | head
~~~
{:class="in"}
~~~
chr1    11873   12227   NR_046018_exon_0_0_chr1_11874_f 0   +
chr1    12612   12721   NR_046018_exon_1_0_chr1_12613_f 0   +
chr1    13220   14409   NR_046018_exon_2_0_chr1_13221_f 0   +
chr1    14361   14829   NR_024540_exon_0_0_chr1_14362_r 0   -
chr1    14969   15038   NR_024540_exon_1_0_chr1_14970_r 0   -
~~~
{:class="out"}

Now instead of saving the sorted results to a separate file, we can actually pipe the result into bedtools merge:
    
~~~
$ bedtools sort -i exons.bed | bedtools merge | head
~~~
{:class="in"}
~~~
chr1    11873   12227
chr1    12612   12721
chr1    13220   14829
chr1    14969   15038
chr1    15795   15947
~~~
{:class="out"}

The result is a new set of intervals representing the merged set of intervals in the input. That is, if a base pair in the genome is covered by 10 features, it will now only be represented once in the output file. Notice how the merged file has less columns than the original file - we've lost the name of the exon and the strand.

Let's save the output of bedtools merge to a file and recalculate the overlap:

~~~
$ bedtools sort -i exons.bed | bedtools merge > exons.merged.bed
$ bedtools intersect -a cpg.bed -b exons.merged.bed -wo | head
~~~
{:class="in"}
~~~
chr1    28735   29810   CpG:_116    chr1    29320   29370   50
chr1    135124  135563  CpG:_30 chr1    134772  139696  439
chr1    327790  328229  CpG:_29 chr1    324438  328581  439
chr1    713984  714547  CpG:_60 chr1    713663  714068  84
chr1    762416  763445  CpG:_115    chr1    761585  762902  486
chr1    762416  763445  CpG:_115    chr1    762970  763155  185
chr1    762416  763445  CpG:_115    chr1    763177  763229  52
chr1    788863  789211  CpG:_28 chr1    788770  794826  348
chr1    854765  854973  CpG:_16 chr1    854714  854817  52
chr1    858970  861632  CpG:_257    chr1    861120  861180  60
~~~
{:class="out"}

This gives us the coordinates of the overlap, the name of the CpG island, the coordinates of the CpG island and the overlap in base pairs. Notice how the repetition in line 3-5 is now gone.

#Challenge: stranded merge

All we actually need from the output of `bedtools intersect` is the overlap column. We can use the unix `cut` command to extract just the 8th field (or column) using the -f flag.

~~~
$ bedtools intersect -a cpg.bed -b exons.merged.bed -wo | cut -f 8 | head
~~~
{:class="in"}
~~~
50
439
439
84
486
185
52
348
52
60
~~~
{:class="out"}

Now we need a program to read in this column of numbers and give us the sum. Hopefully this should ring some bells - The `readings.py` program we wrote earlier will do this with a few small changes.

~~~
$ bedtools intersect -a cpg.bed -b exons.merged.bed -wo | cut -f 8 | python readings.py --sum
~~~
{:class="in"}
~~~
7235722.0
~~~
{:class="out"}

Let's save this result to a file so that we remember it later:

~~~
$ bedtools intersect -a cpg.bed -b exons.merged.bed -wo | cut -f 8 | python readings.py --sum > cpg_result.txt
~~~
{:class="in"}


So we know what the overlap is between CpG islands and exons. How do we know whether this is significant? How do we know what we would expect by random chance?

BedTools includes a command called shuffle which will randomize one of our input files. Let's try randomizing the CpG islands:

~~~
bedtools shuffle -i cpg.bed -g genome.txt  | head 
~~~
{:class="in"}
~~~
chr7    154365006   154366081   CpG:_116
chr2    104893603   104894042   CpG:_30
chr11   11786369    11786808    CpG:_29
chr1    116263073   116264086   CpG:_84
chr10   10674197    10675468    CpG:_99
chrX    7680868 7681763 CpG:_94
chr6    34093178    34095089    CpG:_171
chr1    88263374    88263937    CpG:_60
chr1    150416635   150417664   CpG:_115
chr6    39650957    39651305    CpG:_28
~~~
{:class="out"}

Notice how the CpG islands have been shuffled to random chromosomes. This means that each chromosome will have roughly the same number of CpG islands, regardless of length. A better random model would be to keep the distribution of CpG islands over chromosomes and just shuffle the positions of each CpG on each chromosome. We can do this using the -chrom flag.

~~~
$ bedtools shuffle -chrom -i cpg.bed -g genome.txt | head 
~~~
{:class="in"}
~~~
chr1    2367858 2368933 CpG:_116
chr1    225323760   225324199   CpG:_30
chr1    37368198    37368637    CpG:_29
chr1    13872987    13874000    CpG:_84
chr1    91280571    91281842    CpG:_99
chr1    74656135    74657030    CpG:_94
chr1    64079046    64080957    CpG:_171
chr1    54034535    54035098    CpG:_60
chr1    237633659   237634688   CpG:_115
chr1    51954093    51954441    CpG:_28
~~~
{:class="out"}

Now let's save that result to a file, and count the overlap for the shuffled features.

~~~
$ bedtools shuffle -chrom -i cpg.bed -g genome.txt > cpg.shuffled.bed
$ bedtools intersect -a cpg.shuffled.bed -b exons.merged.bed -wo | cut -f 8 | python readings.py --sum
~~~
{:class="in"}
~~~
681980.0
~~~
{:class="out"}

Great, so our shuffled file overlaps much less with exons than real CpG islands do. How do we know if the difference is significant? One approach is to make lots of different randomly shuffled files, and count the overlap for each one. If we do this 100 times, and we never see an overlap that is as big as the one we see for the original `data/cpg.bed` file, then we can say that the overlap between exons and CpG islands is significant with a P-value (roughly) < 0.01.

First, prepare a directory to hold our shuffled results:

~~~
$ mkdir shuffled_cpg
~~~
{:class="in"}

In the `scripts` folder, you should see a file called `scripts/do_shuffle.sh`. This is a bash script that contains a loop from 1 to 100. It runs bedtools shuffle 100 times, and saves the results to a folder. Have a look at the contents of the file:

~~~
$ cat scripts/do_shuffle.sh
~~~
{:class="in"}
~~~
#! /bin/bash


for i in {1..100}
do
    bedtools shuffle -chrom -i data/cpg.bed -g data/genome.txt > shuffled_cpg/cpg.shuffled${i}.bed
done
~~~
{:class="out"}

Run the script, and then check to make sure you have some shiny new files in the `shuffled_cpg` folder:

~~~
$ scripts/do_suffle.sh
$ ls shuffled_cpg
~~~
{:class="in"}
~~~
cpg.shuffled100.bed  cpg.shuffled21.bed  cpg.shuffled33.bed  cpg.shuffled45.bed  cpg.shuffled57.bed  cpg.shuffled69.bed  cpg.shuffled80.bed  cpg.shuffled92.bed
cpg.shuffled10.bed   cpg.shuffled22.bed  cpg.shuffled34.bed  cpg.shuffled46.bed  cpg.shuffled58.bed  cpg.shuffled6.bed   cpg.shuffled81.bed  cpg.shuffled93.bed
cpg.shuffled11.bed   cpg.shuffled23.bed  cpg.shuffled35.bed  cpg.shuffled47.bed  cpg.shuffled59.bed  cpg.shuffled70.bed  cpg.shuffled82.bed  cpg.shuffled94.bed
cpg.shuffled12.bed   cpg.shuffled24.bed  cpg.shuffled36.bed  cpg.shuffled48.bed  cpg.shuffled5.bed   cpg.shuffled71.bed  cpg.shuffled83.bed  cpg.shuffled95.bed
cpg.shuffled13.bed   cpg.shuffled25.bed  cpg.shuffled37.bed  cpg.shuffled49.bed  cpg.shuffled60.bed  cpg.shuffled72.bed  cpg.shuffled84.bed  cpg.shuffled96.bed
cpg.shuffled14.bed   cpg.shuffled26.bed  cpg.shuffled38.bed  cpg.shuffled4.bed   cpg.shuffled61.bed  cpg.shuffled73.bed  cpg.shuffled85.bed  cpg.shuffled97.bed
cpg.shuffled15.bed   cpg.shuffled27.bed  cpg.shuffled39.bed  cpg.shuffled50.bed  cpg.shuffled62.bed  cpg.shuffled74.bed  cpg.shuffled86.bed  cpg.shuffled98.bed
cpg.shuffled16.bed   cpg.shuffled28.bed  cpg.shuffled3.bed   cpg.shuffled51.bed  cpg.shuffled63.bed  cpg.shuffled75.bed  cpg.shuffled87.bed  cpg.shuffled99.bed
cpg.shuffled17.bed   cpg.shuffled29.bed  cpg.shuffled40.bed  cpg.shuffled52.bed  cpg.shuffled64.bed  cpg.shuffled76.bed  cpg.shuffled88.bed  cpg.shuffled9.bed
cpg.shuffled18.bed   cpg.shuffled2.bed   cpg.shuffled41.bed  cpg.shuffled53.bed  cpg.shuffled65.bed  cpg.shuffled77.bed  cpg.shuffled89.bed
cpg.shuffled19.bed   cpg.shuffled30.bed  cpg.shuffled42.bed  cpg.shuffled54.bed  cpg.shuffled66.bed  cpg.shuffled78.bed  cpg.shuffled8.bed
cpg.shuffled1.bed    cpg.shuffled31.bed  cpg.shuffled43.bed  cpg.shuffled55.bed  cpg.shuffled67.bed  cpg.shuffled79.bed  cpg.shuffled90.bed
cpg.shuffled20.bed   cpg.shuffled32.bed  cpg.shuffled44.bed  cpg.shuffled56.bed  cpg.shuffled68.bed  cpg.shuffled7.bed   cpg.shuffled91.bed
~~~
{:class="out"}

Now we need to calculate the overlap between `data/exons.bed` and each of our shuffled bed files. The command for the first file is:

~~~
$ bedtools intersect -a shuffled_cpg/cpg.shuffled1.bed -b data/exons.merged.bed -wo | cut -f 8 | python scripts/readings.py --sum
~~~
{:class="in"}
~~~
630845.0
~~~
{:class="out"}

See if you can adapt this command and turn it into a for loop which runs over all of our shuffled results.

~~~
$ for filename in shuffled_cpg/*
> do
>     bedtools intersect -a $filename -b data/exons.merged.bed -wo | cut -f 8 | python scripts/readings.py --sum
> done | head
~~~
{:class="in"}
~~~
651838.0
655178.0
633355.0
636853.0
668339.0
651838.0
655178.0
633355.0
636853.0
668339.0
~~~
{:class="out"}

Save that to a file:

~~~
$ for filename in shuffled_cpg/*; 
> do
>     bedtools intersect -a $filename -b data/exons.merged.bed -wo | cut -f 8 | python scripts/readings.py --sum
> done > shuffled_results.txt
~~~
{:class="in"}


Plotting our results with python
--------------------------------

Also in your `scripts` directory, there should be a file called results.py. Open it and have a look.

~~~
$ cat scripts/results.py
~~~
{:class="in"}
~~~
import sys
import numpy as np


def main():
    script = sys.argv[0]
    shuffled_results_file = sys.argv[1]
    real_results_file = sys.argv[2]

    process(shuffled_results_file, real_results_file)


def process(shuffled_results_file, real_results_file):
    shuffled_data = np.loadtxt(shuffled_results_file)
    real_data = np.loadtxt(real_results_file)

    print 'Real data:'
    print real_data
    print
    print 'Shuffled data:'
    print shuffled_data

main()
~~~
{:class="out"}

Which prints:

~~~
python results.py shuffled_results.txt cpg_result.txt
~~~
{:class="in"}
~~~
Real data:
7235722.0

Shuffled data:
[ 640473.  638611.  644321.  642422.  634641.  660745.  645816.  635427.
  595615.  644916.  662412.  645489.  666141.  674597.  637099.  621515.
  647847.  653944.  641051.  685423.  658610.  686554.  618233.  655932.
  670053.  641412.  617865.  651189.  639806.  658123.  639381.  644652.
  667240.  689363.  627791.  625137.  635577.  643151.  616453.  633041.
  629223.  645209.  629201.  639179.  649602.  638849.  667827.  637550.
  652560.  647235.  710669.  626332.  689819.  646094.  631575.  633863.
  657661.  642538.  648691.  660292.  649780.  643179.  615128.  653863.
  610901.  613489.  624788.  680903.  617416.  654761.  663484.  672619.
  615038.  630960.  622431.  634951.  659257.  649931.  633901.  612363.
  639325.  683887.  656753.  690935.  661310.  692022.  633441.  644006.
  652187.  643428.  643711.  682621.  607918.  674996.  674909.  625800.
  642090.  662203.  667386.  660448.]
~~~
{:class="out"}

What happens if we call the script with no arguments?

~~~
$ python scripts/results.py
~~~
{:class="in"}
~~~
Traceback (most recent call last):
  File "results.py", line 23, in <module>
    main()
  File "results.py", line 7, in main
    shuffled_results_file = sys.argv[1]
IndexError: list index out of range
~~~
{:class="err"}

Not very helpful. Let's add a usage description. Modify the python file so that it looks like this:

~~~
import sys
import numpy as np

usage_string = """
Results.py plots the real overlap between two bed files
versus the overlap when one of the results files is randomly
shuffled.

Usage: python results.py shuffled_results.txt real_results.txt
"""

def main():

    if not len(sys.argv) == 3:
        sys.exit(usage_string)
~~~

Let's try plotting a histogram of the shuffled data, first you'll need to add `from matplotlib import pyplot as plt` at the top of the file, then add this at the bottom of the process function.

~~~
    print 'Shuffled data:'
    print shuffled_data

    plt.hist(shuffled_data, color='black')
    plt.savefig('cpg_overlap_exons.png')
~~~

Have a look at the image you just made. At the moment we're only plotting the shuffled results, so let's add the real overlap as a vertical red line. Update your results.py file so that it looks like this: 

~~~
import sys 
import numpy as np
from matplotlib import pyplot as plt 

usage_string = """ 
Results.py plots the real overlap between two bed files
versus the overlap when one of the results files is randomly
shuffled.

Usage: python results.py shuffled_results.txt real_results.txt
"""

def main():

    if not len(sys.argv) == 3:
        sys.exit(usage_string)

    script = sys.argv[0]
    shuffled_results_file = sys.argv[1]
    real_results_file = sys.argv[2]

    process(shuffled_results_file, real_results_file)


def process(shuffled_results_file, real_results_file):
    shuffled_data = np.loadtxt(shuffled_results_file)
    real_data = np.loadtxt(real_results_file)

    print 'Real data:'
    print real_data
    print
    print 'Shuffled data:'
    print shuffled_data

    plt.hist(shuffled_data, color='black')
    plt.axvline(real_data, color='red', ls='--', lw=2)
    plt.savefig('cpg_overlap_exons.png')

main()
~~~

Finally, let's add some axis labels and a legend:

~~~
import sys
import numpy as np
from matplotlib import pyplot as plt

usage_string = """
Results.py plots the real overlap between two bed files
versus the overlap when one of the results files is randomly
shuffled.

Usage: python results.py shuffled_results.txt real_results.txt
"""

def main():

    if not len(sys.argv) == 3:
        sys.exit(usage_string)

    script = sys.argv[0]
    shuffled_results_file = sys.argv[1]
    real_results_file = sys.argv[2]

    process(shuffled_results_file, real_results_file)


def process(shuffled_results_file, real_results_file):
    shuffled_data = np.loadtxt(shuffled_results_file)
    real_data = np.loadtxt(real_results_file)

    print 'Real data:'
    print real_data
    print
    print 'Shuffled data:'
    print shuffled_data

    plt.hist(shuffled_data, color='black', label='Shuffled Data')
    plt.axvline(real_data, color='red', ls='--', lw=2, label='Real Data')
    plt.xlabel('Nucleotides overlap')
    plt.ylabel('Number of shuffled results')
    plt.legend()
    plt.savefig('cpg_overlap_exons.png')

main()
~~~


