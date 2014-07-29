# RNAseq Analysis Example

This is an introduction to RNAseq analysis for use at Software Carpentry bootcamps that have covered novice R. It involves reading in some count data from an RNAseq experiment, exploring the data using base R functions and then analysis with the package DESeq2.

# Install required CRAN packages

First, install some packages that you'll use.


```r
install.packages("gplots")
install.packages("ggplot2")
install.packages("calibrate")
```

# Introduction and data import

The analysis of an RNAseq experiment begins with sequencing reads. These then need to be aligned to a reference genome or transcriptome. There are many different alignment tools available, but the process of alignment is both computationally intensive and time-consuming, so we won't cover it today. Once reads are aligned, the number of reads mapped to each gene can be counted. Again, there are several ways of doing this. The best way to find out about the tools that are available and suitable for your research is to look for recent review papers that compare the different tools.

The data for this tutorial comes from a PLOS ONE paper, [Genome-Wide Transcriptional Profiling of Skin and Dorsal Root Ganglia after Ultraviolet-B-Induced Inflammation](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0093338)[1], and the raw data can be downloaded from the [Gene Expression Omnibus database (GEO)](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54413). 

This data has already been downloaded and aligned to the human genome. The command line tool [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) was used to count reads mapped to human genes from the Ensembl annotation (available for download [here](http://www.ensembl.org/info/data/ftp/index.html)). 

The output from this tool is provided in the `counts.txt` file. Have a look at this file in the shell, using `head`.

Import the data into R as a `data.frame` and examine it again. You can set the arguments of `read.table` to import the first row as a header giving the column names, and the first column as row names. 


```r
# Filename with output from featureCounts
countFile <- "data/counts.txt"
# Read in the data
countData <- read.table(countFile, header=TRUE, row.names=1)
head(countData)
```

```
##                                  Chr
## DDX11L1                      1;1;1;1
## WASH7P     1;1;1;1;1;1;1;1;1;1;1;1;1
## MIR1302-10                     1;1;1
## FAM138A                        1;1;1
## OR4G4P                           1;1
## OR4G11P                            1
##                                                                                    Start
## DDX11L1                                                          11869;12595;12975;13221
## WASH7P     14363;14970;15796;16607;16854;17233;17498;17602;17915;18268;24734;29321;29534
## MIR1302-10                                                             29554;30267;30976
## FAM138A                                                                34554;35245;35721
## OR4G4P                                                                       52473;54830
## OR4G11P                                                                            62948
##                                                                                      End
## DDX11L1                                                          12227;12721;13052;14412
## WASH7P     14829;15038;15947;16765;17055;17368;17504;17742;18061;18379;24891;29370;29806
## MIR1302-10                                                             30039;30667;31109
## FAM138A                                                                35174;35481;36081
## OR4G4P                                                                       53312;54936
## OR4G11P                                                                            63887
##                               Strand Length
## DDX11L1                      +;+;+;+   1756
## WASH7P     -;-;-;-;-;-;-;-;-;-;-;-;-   2073
## MIR1302-10                     +;+;+   1021
## FAM138A                        -;-;-   1219
## OR4G4P                           +;+    947
## OR4G11P                            +    940
##            ctl1.fastq_tophat.accepted_hits.bam
## DDX11L1                                      0
## WASH7P                                       0
## MIR1302-10                                   0
## FAM138A                                      0
## OR4G4P                                       0
## OR4G11P                                      0
##            ctl2.fastq_tophat.accepted_hits.bam
## DDX11L1                                      0
## WASH7P                                       0
## MIR1302-10                                   0
## FAM138A                                      0
## OR4G4P                                       0
## OR4G11P                                      0
##            ctl3.fastq_tophat.accepted_hits.bam
## DDX11L1                                      0
## WASH7P                                       0
## MIR1302-10                                   0
## FAM138A                                      0
## OR4G4P                                       0
## OR4G11P                                      0
##            uvb1.fastq_tophat.accepted_hits.bam
## DDX11L1                                      0
## WASH7P                                       0
## MIR1302-10                                   0
## FAM138A                                      0
## OR4G4P                                       0
## OR4G11P                                      0
##            uvb2.fastq_tophat.accepted_hits.bam
## DDX11L1                                      0
## WASH7P                                       0
## MIR1302-10                                   0
## FAM138A                                      0
## OR4G4P                                       0
## OR4G11P                                      0
##            uvb3.fastq_tophat.accepted_hits.bam
## DDX11L1                                      0
## WASH7P                                       0
## MIR1302-10                                   0
## FAM138A                                      0
## OR4G4P                                       0
## OR4G11P                                      0
```

```r
colnames(countData)
```

```
##  [1] "Chr"                                
##  [2] "Start"                              
##  [3] "End"                                
##  [4] "Strand"                             
##  [5] "Length"                             
##  [6] "ctl1.fastq_tophat.accepted_hits.bam"
##  [7] "ctl2.fastq_tophat.accepted_hits.bam"
##  [8] "ctl3.fastq_tophat.accepted_hits.bam"
##  [9] "uvb1.fastq_tophat.accepted_hits.bam"
## [10] "uvb2.fastq_tophat.accepted_hits.bam"
## [11] "uvb3.fastq_tophat.accepted_hits.bam"
```

```r
class(countData)
```

```
## [1] "data.frame"
```

The data.frame contains information about genes (one gene per row) with the gene positions in the first five columns and then information about the number of reads aligning to the gene in each experimental sample. There are three replicates for control (column names starting with "ctl") and three for samples treated with ultraviolet-B light (starting "uvb"). We don't need the information on gene position for this analysis, just the counts for each gene and sample, so we can remove it from the data frame.


```r
# Remove first five columns (chr, start, end, strand, length)
countData <- countData[ ,-(1:5)]
head(countData)
```

```
##            ctl1.fastq_tophat.accepted_hits.bam
## DDX11L1                                      0
## WASH7P                                       0
## MIR1302-10                                   0
## FAM138A                                      0
## OR4G4P                                       0
## OR4G11P                                      0
##            ctl2.fastq_tophat.accepted_hits.bam
## DDX11L1                                      0
## WASH7P                                       0
## MIR1302-10                                   0
## FAM138A                                      0
## OR4G4P                                       0
## OR4G11P                                      0
##            ctl3.fastq_tophat.accepted_hits.bam
## DDX11L1                                      0
## WASH7P                                       0
## MIR1302-10                                   0
## FAM138A                                      0
## OR4G4P                                       0
## OR4G11P                                      0
##            uvb1.fastq_tophat.accepted_hits.bam
## DDX11L1                                      0
## WASH7P                                       0
## MIR1302-10                                   0
## FAM138A                                      0
## OR4G4P                                       0
## OR4G11P                                      0
##            uvb2.fastq_tophat.accepted_hits.bam
## DDX11L1                                      0
## WASH7P                                       0
## MIR1302-10                                   0
## FAM138A                                      0
## OR4G4P                                       0
## OR4G11P                                      0
##            uvb3.fastq_tophat.accepted_hits.bam
## DDX11L1                                      0
## WASH7P                                       0
## MIR1302-10                                   0
## FAM138A                                      0
## OR4G4P                                       0
## OR4G11P                                      0
```

```r
colnames(countData)
```

```
## [1] "ctl1.fastq_tophat.accepted_hits.bam"
## [2] "ctl2.fastq_tophat.accepted_hits.bam"
## [3] "ctl3.fastq_tophat.accepted_hits.bam"
## [4] "uvb1.fastq_tophat.accepted_hits.bam"
## [5] "uvb2.fastq_tophat.accepted_hits.bam"
## [6] "uvb3.fastq_tophat.accepted_hits.bam"
```

We can rename the columns to something shorter and a bit more readable.


```r
# Manually
c("ctl1", "ctl2", "ctl3", "uvb1", "uvb2", "uvb3")
# Using paste
?paste
paste("ctl", 1:3)
paste("ctl", 1:3, sep="")
?paste0
paste0("ctl", 1:3)
c(paste0("ctl", 1:3), paste0("uvb", 1:3))
```

An easier way to do this, especially for files with many columns, is to use the `gsub` command to strip out the extra information. This is also more robust to introduced errors, for example if the column order changes at some point in the future or you add additional replicates.


```r
# Using gsub -- robust
?gsub
gsub(pattern=".fastq_tophat.accepted_hits.bam", replacement="", x=colnames(countData))
```

```
## [1] "ctl1" "ctl2" "ctl3" "uvb1" "uvb2" "uvb3"
```

```r
colnames(countData) <- gsub(pattern=".fastq_tophat.accepted_hits.bam", replacement="", x=colnames(countData))
head(countData)
```

```
##            ctl1 ctl2 ctl3 uvb1 uvb2 uvb3
## DDX11L1       0    0    0    0    0    0
## WASH7P        0    0    0    0    0    0
## MIR1302-10    0    0    0    0    0    0
## FAM138A       0    0    0    0    0    0
## OR4G4P        0    0    0    0    0    0
## OR4G11P       0    0    0    0    0    0
```

## Exercise 1
Find the gene with the highest expression in any sample -- remember, each row is a gene. Extract the expression data for this gene for all samples. In which sample does it have the highest expression? 

What is the function of the gene? Can you suggest why this is the top expressed gene?

Hint 1: use the `apply` function from the introductory R lessons.

Hint 2: try `?which.max`.



# Data investigation using base R

We can investigate this data a bit more using some of the basic R functions before going on to use more sophisticated analysis tools.

First make a copy of the data, because we'll need it later. We will work on the copy. We will calculate the mean for each gene for each condition and plot them.


```r
countData2 <- countData #make a copy

# get Control columns
colnames(countData2)
```

```
## [1] "ctl1" "ctl2" "ctl3" "uvb1" "uvb2" "uvb3"
```

```r
?grep #grep searches for matches to a pattern 
grep("ctl", colnames(countData2))
```

```
## [1] 1 2 3
```

```r
ctlCols <- grep("ctl", colnames(countData2))
head(countData2[,ctlCols])
```

```
##            ctl1 ctl2 ctl3
## DDX11L1       0    0    0
## WASH7P        0    0    0
## MIR1302-10    0    0    0
## FAM138A       0    0    0
## OR4G4P        0    0    0
## OR4G11P       0    0    0
```

```r
head(apply(countData2[, ctlCols], 1, mean))
```

```
##    DDX11L1     WASH7P MIR1302-10    FAM138A     OR4G4P    OR4G11P 
##          0          0          0          0          0          0
```

```r
# here we'll use rowMeans instead, it's a convenient shortcut, and also faster!
countData2$ctlMean <- rowMeans(countData2[, ctlCols])

# same for uvb
uvbCols <- grep("uvb", colnames(countData2))
countData2$uvbMean <- rowMeans(countData2[, uvbCols])
```

Plot the mean expression of each gene in control against the UVB sample mean. Are there any outliers?


```r
plot(countData2$ctlMean, countData2$uvbMean)
```

![plot of chunk plot_means](./analysis_files/figure-html/plot_means.png) 


```r
library("ggplot2")
ggplot(countData2, aes(x=ctlMean, y=uvbMean)) + geom_point()
```

![plot of chunk ggplot_means](./analysis_files/figure-html/ggplot_means.png) 

## Exercise 2
How could you make this plot more informative and look more professional? 

Hint: try using a log scale. You can also changing colours, transparencies, sizes, or shapes of points. 

`?par` will give you information on lots of graphical parameters that can be set. Help for ggplot2 can be found [here](http://docs.ggplot2.org/current/).





There are many more options you can use to alter the appearance of these plots.

# Find candidate differentially expressed genes

We can find candidate differentially expressed genes by looking for genes with a large change between control and UVB samples. A common threshold used is log2 fold change more than 2 or less than -2. We will calculate log2 fold change for all the genes and colour the genes with log2 fold change of more than 2 or less than -2 on the plot.

First, check for genes with a mean expression of 0. Putting zeroes into the log2 fold change calculation will produce NAs, so we might want to remove these genes. Note: this is for mathematical reasons, although different software may produce different results when you try to do `log2(0)`

`TRUE` and `FALSE` can also be represented as 1 and 0. This is useful for getting the total number of observations for which a condition is true. 


```r
TRUE == 0
```

```
## [1] FALSE
```

```r
TRUE == 1
```

```
## [1] TRUE
```

```r
FALSE == 0
```

```
## [1] TRUE
```

This can be applied to testing whether genes have a mean expression of more than zero.


```r
head(countData2$ctlMean)
```

```
## [1] 0 0 0 0 0 0
```

```r
head(countData2$ctlMean > 0)
```

```
## [1] FALSE FALSE FALSE FALSE FALSE FALSE
```

```r
head(as.numeric(countData2$ctlMean > 0))
```

```
## [1] 0 0 0 0 0 0
```

When we call `sum(countData2$ctlMean > 0)`, we're really asking, "how many genes have a mean above 0 in the control group?"


```r
# discuss: why to remove zeroes (NAs produced)
sum(countData2$ctlMean > 0)
```

```
## [1] 474
```

```r
sum(countData2$uvbMean > 0)
```

```
## [1] 528
```

```r
nrow(countData2)
```

```
## [1] 56638
```

```r
countData2 <- subset(countData2, (countData2$ctlMean > 0 | countData2$uvbMean > 0))
# explain: | operator meaning OR in this context?
nrow(countData2)
```

```
## [1] 565
```


```r
# explain: what is fold change? why do we use log2 to quantify the fold change?
# discuss: Inf / -Inf may be produced in some cases. Concept of adding pseudocounts.
countData2$log2FC <- log2(countData2$uvbMean / countData2$ctlMean)
# again, reinforce that summing a logical vector gives you the number of 
# occurences of TRUE.
sum(countData2$log2FC > 2)
```

```
## [1] 170
```

```r
sum(countData2$log2FC < -2)
```

```
## [1] 46
```

Make a new column to store this information in.


```r
countData2$outlier <- FALSE
countData2$outlier[countData2$log2FC > 2] <- TRUE
countData2$outlier[countData2$log2FC < -2] <- TRUE
```


```r
plot(countData2$ctlMean, countData2$uvbMean, log="xy", pch=16)
```

```
## Warning: 91 x values <= 0 omitted from logarithmic plot
## Warning: 37 y values <= 0 omitted from logarithmic plot
```

```r
points(countData2$ctlMean[countData2$outlier==TRUE], countData2$uvbMean[countData2$outlier==TRUE], col="red", pch=16)
```

![plot of chunk plot_outliers](./analysis_files/figure-html/plot_outliers.png) 


```r
ggplot(countData2, aes(x=ctlMean, y=uvbMean, colour=outlier)) + geom_point() + scale_x_log10() + scale_y_log10() + theme_bw()
```

![plot of chunk ggplot_outliers](./analysis_files/figure-html/ggplot_outliers.png) 

What do you notice about the positions of the outliers on these plots? How would you interpret this?

# DESeq2 analysis

DESeq2 is an R package for analysis of RNAseq data. It is available from [Bioconductor](http://www.bioconductor.org/). Bioconductor is a project to provide tools for analysing high-throughput genomic data including RNA-seq, ChIP-seq and arrays. You can explore Bioconductor packages [here](http://www.bioconductor.org/packages/release/BiocViews.html#___Software). 


```r
# install and have a break to check everyone is up to date?
# explain bioconductor?
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```


```r
library("DESeq2")
```

```
## Loading required package: GenomicRanges
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Loading required package: IRanges
## Loading required package: GenomeInfoDb
## Loading required package: Rcpp
## Loading required package: RcppArmadillo
```

```
## Warning: package 'RcppArmadillo' was built under R version 3.1.1
```

```r
citation("DESeq2")
```

```
## 
##   Michael I Love, Wolfgang Huber and Simon Anders (2014):
##   Moderated estimation of fold change and dispersion for RNA-Seq
##   data with DESeq2. bioRxiv preprint
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {Moderated estimation of fold change and dispersion for RNA-Seq data with DESeq2},
##     author = {Michael I Love and Wolfgang Huber and Simon Anders},
##     year = {2014},
##     journal = {bioRxiv},
##     doi = {10.1101/002832},
##     url = {http://dx.doi.org/10.1101/002832},
##   }
```

It requires the count data to be in matrix form, and an additional dataframe describing the structure of the experiment.


```r
# countData is currently a data.frame, but DESeq2 expects its input to be in 
# matrix format, so we will convert our countData to a matrix.
class(countData)
```

```
## [1] "data.frame"
```

```r
countData <- as.matrix(countData)
class(countData)
```

```
## [1] "matrix"
```

```r
head(countData)
```

```
##            ctl1 ctl2 ctl3 uvb1 uvb2 uvb3
## DDX11L1       0    0    0    0    0    0
## WASH7P        0    0    0    0    0    0
## MIR1302-10    0    0    0    0    0    0
## FAM138A       0    0    0    0    0    0
## OR4G4P        0    0    0    0    0    0
## OR4G11P       0    0    0    0    0    0
```

```r
# construct colData dataframe
# three replicates of control and UVB.
colData <- data.frame(condition=c(rep("ctl", 3), rep("uvb",3)), row.names=colnames(countData))
```

DESeq works on a particular type of object called a DESeqDataSet.


```r
# introduce how DESeq2 works - type of object it works on etc
# instantiate the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
dds
```

```
## class: DESeqDataSet 
## dim: 56638 6 
## exptData(0):
## assays(1): counts
## rownames(56638): DDX11L1 WASH7P ... MT-TT MT-TP
## rowData metadata column names(0):
## colnames(6): ctl1 ctl2 ... uvb2 uvb3
## colData names(1): condition
```

Run the DESeq pipeline on this object. [Describe pipeline steps?]
Get results and have a look at them


```r
dds <- DESeq(dds)
```

```
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
```

```r
# Get differential expression results
res <- results(dds)
head(res)
```

```
## log2 fold change (MAP): condition uvb vs ctl 
## Wald test p-value: condition uvb vs ctl 
## DataFrame with 6 rows and 6 columns
##             baseMean log2FoldChange     lfcSE      stat    pvalue
##            <numeric>      <numeric> <numeric> <numeric> <numeric>
## DDX11L1            0             NA        NA        NA        NA
## WASH7P             0             NA        NA        NA        NA
## MIR1302-10         0             NA        NA        NA        NA
## FAM138A            0             NA        NA        NA        NA
## OR4G4P             0             NA        NA        NA        NA
## OR4G11P            0             NA        NA        NA        NA
##                 padj
##            <numeric>
## DDX11L1           NA
## WASH7P            NA
## MIR1302-10        NA
## FAM138A           NA
## OR4G4P            NA
## OR4G11P           NA
```

```r
table(res$padj<0.05)
```

```
## 
## FALSE  TRUE 
##   337    50
```

```r
# Order by adjusted p-value
res <- res[order(res$padj), ]
head(res)
```

```
## log2 fold change (MAP): condition uvb vs ctl 
## Wald test p-value: condition uvb vs ctl 
## DataFrame with 6 rows and 6 columns
##         baseMean log2FoldChange     lfcSE      stat    pvalue      padj
##        <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
## ANXA3     312.32          3.910    0.3475    11.251 2.287e-29 8.849e-27
## SFRP2    1711.27         -3.356    0.3407    -9.852 6.735e-23 1.303e-20
## MT2P1     167.36          4.067    0.4582     8.878 6.834e-19 8.816e-17
## IL8        31.70          5.976    0.9456     6.320 2.620e-10 2.535e-08
## PTPN13    358.34         -1.386    0.2362    -5.869 4.388e-09 3.396e-07
## CXCL1      42.47          3.691    0.6448     5.724 1.040e-08 6.705e-07
```

Combine DEseq results with the original counts data. Write significant results to a file.


```r
resData <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
head(resData)
```

```
##   Row.names baseMean log2FoldChange  lfcSE   stat    pvalue      padj
## 1     ANXA3   312.32          3.910 0.3475 11.251 2.287e-29 8.849e-27
## 2     SFRP2  1711.27         -3.356 0.3407 -9.852 6.735e-23 1.303e-20
## 3     MT2P1   167.36          4.067 0.4582  8.878 6.834e-19 8.816e-17
## 4       IL8    31.70          5.976 0.9456  6.320 2.620e-10 2.535e-08
## 5    PTPN13   358.34         -1.386 0.2362 -5.869 4.388e-09 3.396e-07
## 6     CXCL1    42.47          3.691 0.6448  5.724 1.040e-08 6.705e-07
##      ctl1     ctl2    ctl3   uvb1   uvb2   uvb3
## 1   26.33   35.131   48.83 378.95 783.44 601.27
## 2 1847.65 4865.598 2683.82 382.63 325.05 162.88
## 3   19.02    8.107   23.63 438.34 373.97 141.08
## 4    0.00    0.000    0.00  66.12  44.18  79.93
## 5  485.69  535.067  538.65 203.25 205.13 182.26
## 6   11.70    2.702    0.00  89.99  50.49  99.91
```

```r
names(resData)[1] <- "GeneID"
head(resData)
```

```
##   GeneID baseMean log2FoldChange  lfcSE   stat    pvalue      padj    ctl1
## 1  ANXA3   312.32          3.910 0.3475 11.251 2.287e-29 8.849e-27   26.33
## 2  SFRP2  1711.27         -3.356 0.3407 -9.852 6.735e-23 1.303e-20 1847.65
## 3  MT2P1   167.36          4.067 0.4582  8.878 6.834e-19 8.816e-17   19.02
## 4    IL8    31.70          5.976 0.9456  6.320 2.620e-10 2.535e-08    0.00
## 5 PTPN13   358.34         -1.386 0.2362 -5.869 4.388e-09 3.396e-07  485.69
## 6  CXCL1    42.47          3.691 0.6448  5.724 1.040e-08 6.705e-07   11.70
##       ctl2    ctl3   uvb1   uvb2   uvb3
## 1   35.131   48.83 378.95 783.44 601.27
## 2 4865.598 2683.82 382.63 325.05 162.88
## 3    8.107   23.63 438.34 373.97 141.08
## 4    0.000    0.00  66.12  44.18  79.93
## 5  535.067  538.65 203.25 205.13 182.26
## 6    2.702    0.00  89.99  50.49  99.91
```

```r
sig <- subset(resData, padj<0.05)
dir.create("results")
```

```
## Warning: 'results' already exists
```

```r
write.table(sig, file="results/sig.txt", sep="\t") # tab delim data
```

You can open this file in Excel or any text editor (try it now).

# Data Visualization

We can also do some exploratory plotting of the data.


```r
plotDispEsts(dds, main="Dispersion plot")
```

![plot of chunk plot_dispersion](./analysis_files/figure-html/plot_dispersion.png) 


```r
# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
plotPCA(rld)
```

![plot of chunk plot_heatmaps](./analysis_files/figure-html/plot_heatmaps1.png) 

```r
# Sample distance heatmap
head(assay(rld))
```

```
##            ctl1 ctl2 ctl3 uvb1 uvb2 uvb3
## DDX11L1       0    0    0    0    0    0
## WASH7P        0    0    0    0    0    0
## MIR1302-10    0    0    0    0    0    0
## FAM138A       0    0    0    0    0    0
## OR4G4P        0    0    0    0    0    0
## OR4G11P       0    0    0    0    0    0
```

```r
assay(rld)[1:5,1:5]
```

```
##            ctl1 ctl2 ctl3 uvb1 uvb2
## DDX11L1       0    0    0    0    0
## WASH7P        0    0    0    0    0
## MIR1302-10    0    0    0    0    0
## FAM138A       0    0    0    0    0
## OR4G4P        0    0    0    0    0
```

```r
t(assay(rld))[1:5,1:5]
```

```
##      DDX11L1 WASH7P MIR1302-10 FAM138A OR4G4P
## ctl1       0      0          0       0      0
## ctl2       0      0          0       0      0
## ctl3       0      0          0       0      0
## uvb1       0      0          0       0      0
## uvb2       0      0          0       0      0
```

```r
dist(t(assay(rld)))
```

```
##       ctl1  ctl2  ctl3  uvb1  uvb2
## ctl2 12.18                        
## ctl3 12.49 13.61                  
## uvb1 15.07 16.22 17.29            
## uvb2 20.01 19.73 20.43 15.05      
## uvb3 18.70 18.94 20.21 13.75 11.83
```

```r
as.matrix(dist(t(assay(rld))))
```

```
##       ctl1  ctl2  ctl3  uvb1  uvb2  uvb3
## ctl1  0.00 12.18 12.49 15.07 20.01 18.70
## ctl2 12.18  0.00 13.61 16.22 19.73 18.94
## ctl3 12.49 13.61  0.00 17.29 20.43 20.21
## uvb1 15.07 16.22 17.29  0.00 15.05 13.75
## uvb2 20.01 19.73 20.43 15.05  0.00 11.83
## uvb3 18.70 18.94 20.21 13.75 11.83  0.00
```

```r
sampleDists <- as.matrix(dist(t(assay(rld))))
heatmap(sampleDists)
```

![plot of chunk plot_heatmaps](./analysis_files/figure-html/plot_heatmaps2.png) 

```r
# better heatmap with gplots
library("gplots")
```

```
## Warning: package 'gplots' was built under R version 3.1.1
```

```
## KernSmooth 2.23 loaded
## Copyright M. P. Wand 1997-2009
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:IRanges':
## 
##     space
## 
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
heatmap.2(sampleDists)
```

![plot of chunk plot_heatmaps](./analysis_files/figure-html/plot_heatmaps3.png) 

```r
heatmap.2(sampleDists, col=colorpanel(64, "steelblue", "white"), key=FALSE, trace="none")
```

![plot of chunk plot_heatmaps](./analysis_files/figure-html/plot_heatmaps4.png) 

```r
heatmap.2(sampleDists, col=colorpanel(64, "black", "white"), key=FALSE, trace="none")
```

![plot of chunk plot_heatmaps](./analysis_files/figure-html/plot_heatmaps5.png) 

```r
heatmap.2(sampleDists, col=colorpanel(64, "red", "black", "green"), key=FALSE, trace="none")
```

![plot of chunk plot_heatmaps](./analysis_files/figure-html/plot_heatmaps6.png) 

```r
heatmap.2(sampleDists, col=colorpanel(64, "red", "white", "blue"), key=FALSE, trace="none")
```

![plot of chunk plot_heatmaps](./analysis_files/figure-html/plot_heatmaps7.png) 


```r
# Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")
```

![plot of chunk plot_pval_hist](./analysis_files/figure-html/plot_pval_hist.png) 



```r
# These are the plots that are most recognisable from papers
# MA Plot
par(pch=16)
with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x"))
```

```
## Warning: 56073 x values <= 0 omitted from logarithmic plot
```

```r
with(subset(res, padj<.05), points(baseMean, log2FoldChange, col="red", pch=16))
library("calibrate")
```

```
## Loading required package: MASS
```

```r
?textxy
res$Gene <- rownames(res)
with(subset(res, padj<.05), textxy(baseMean, log2FoldChange, labs=Gene, cex=1, col=2))
```

![plot of chunk MA_plot](./analysis_files/figure-html/MA_plot.png) 


```r
# Volcano plot
# Set point character
par(pch=16)
with(res, plot(log2FoldChange, -log10(pvalue), main="Volcano plot"))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), col="green"))
# Add legend
legend("topleft", legend=c("FDR<0.05", "|LFC|>1", "both"), pch=16, col=c("red","orange","green"))
# Label points
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=1))
```

![plot of chunk volcano_plot](./analysis_files/figure-html/volcano_plot.png) 

# References

1. Dawes JM, Antunes-Martins A, Perkins JR, Paterson KJ, Sisignano M, et al. (2014) Genome-Wide Transcriptional Profiling of Skin and Dorsal Root Ganglia after Ultraviolet-B-Induced Inflammation. PLoS ONE 9(4): e93338. doi: 10.1371/journal.pone.0093338 
