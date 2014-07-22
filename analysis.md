# RNAseq Analysis Example

This is an introduction to RNAseq analysis for use at Software Carpentry bootcamps that have covered novice R. It involves reading in some count data from an RNAseq experiment, exploring the data using base R functions and then analysis with the package DESeq2.

## Install required CRAN packages

First, install some packages that you'll use.


```r
install.packages("gplots")
install.packages("ggplot2")
install.packages("calibrate")
```

Import the data as a `data.frame` and examine it. The data is stored in a text file with the first line being a header giving the column names, and the row names in the first column. 


```r
## Filename with output from featureCounts
countfile <- "data/counts.txt"
## Read in the data
countdata <- read.table(countfile, header=TRUE, row.names=1)
head(countdata)
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
colnames(countdata)
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
class(countdata)
```

```
## [1] "data.frame"
```
It contains information about genes (one gene per row) with the gene positions in the first five columns and then information about the number of reads aligning to the gene in each experimental sample. We don't need the information on gene position, so we can remove it from the data frame.


```r
# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,-(1:5)]
head(countdata)
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
colnames(countdata)
```

```
## [1] "ctl1.fastq_tophat.accepted_hits.bam"
## [2] "ctl2.fastq_tophat.accepted_hits.bam"
## [3] "ctl3.fastq_tophat.accepted_hits.bam"
## [4] "uvb1.fastq_tophat.accepted_hits.bam"
## [5] "uvb2.fastq_tophat.accepted_hits.bam"
## [6] "uvb3.fastq_tophat.accepted_hits.bam"
```

We can rename the columns to something a bit more readable.

```r
## Manually
c("ctl1", "ctl2", "ctl3", "uvb1", "uvb2", "uvb3")
## Using paste
?paste
paste("ctl", 1:3)
paste("ctl", 1:3, sep="")
?paste0
paste0("ctl", 1:3)
c(paste0("ctl", 1:4), paste0("uvb", 1:5))
```

Using `gsub` is a more reproducible way to do this.

```r
## Using gsub -- reproducible
?gsub
gsub(pattern=".fastq_tophat.accepted_hits.bam", replacement="", x=colnames(countdata))
```

```
## [1] "ctl1" "ctl2" "ctl3" "uvb1" "uvb2" "uvb3"
```

```r
colnames(countdata) <- gsub(pattern=".fastq_tophat.accepted_hits.bam", replacement="", x=colnames(countdata))
head(countdata)
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

### Exercise 1
Find the gene with the highest expression in any sample. Extract  the expression data for this gene for all samples. In which sample does it have the highest expression? What is the function of the gene? Can you suggest why this is the top expressed gene?


```r
max(apply(countdata, 1, max)) #max expression is 7013
```

```
## [1] 7013
```

```r
which.max(apply(countdata, 1, max)) #gene is EEF1A1P9
```

```
## EEF1A1P9 
##    13514
```

```r
countdata[13514, ] #get other sample data - max is in uvb1
```

```
##          ctl1 ctl2 ctl3 uvb1 uvb2 uvb3
## EEF1A1P9 3570 3788 4345 7013 4217 3630
```

```r
#this is a pseudogene - maybe an artefact of only aligning reads to a single chromosome?
```

## Data investigation using base R

We can investigate this data a bit more using some of the basic R functions before going on to use more sophisticated analysis tools.

We will calculate the mean for each gene for each condition. First make a copy of the data, because we'll need it later. We will work on the copy.


```r
countdata2 <- countdata

#get Control columns
colnames(countdata2)
```

```
## [1] "ctl1" "ctl2" "ctl3" "uvb1" "uvb2" "uvb3"
```

```r
grep("ctl", colnames(countdata2))
```

```
## [1] 1 2 3
```

```r
ctlCols <- grep("ctl", colnames(countdata2))
head(countdata2[,ctlCols])
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
countdata2$ctlMean <- apply(countdata2[, ctlCols], 1, mean)

#same for uvb
uvbCols <- grep("uvb", colnames(countdata2))
countdata2$uvbMean <- apply(countdata2[, uvbCols], 1, mean)
```

Plot the mean expression of each gene in control against the UVB sample mean. Look for outliers.


```r
plot(countdata2$ctlMean, countdata2$uvbMean)
```

![plot of chunk plot_means](./analysis_files/figure-html/plot_means.png) 


```r
library(ggplot2)
ggplot(countdata2, aes(x=ctlMean, y=uvbMean)) + geom_point()
```

![plot of chunk ggplot_means](./analysis_files/figure-html/ggplot_means.png) 

### Exercise 2
How could you make this plot more informative and look more professional? Hint: try using a log scale. Try changing colours, transparencies, sizes, or shapes of points. 

`help(par)` will give you information on lots of graphical parameters that can be set. Help for ggplot2 can be found [here](http://docs.ggplot2.org/current/).


```r
plot(countdata2$ctlMean, countdata2$uvbMean, log="xy")
```

```
## Warning: 56164 x values <= 0 omitted from logarithmic plot
## Warning: 56110 y values <= 0 omitted from logarithmic plot
```

![plot of chunk exercise2_1](./analysis_files/figure-html/exercise2_1.png) 


```r
ggplot(countdata2, aes(x=ctlMean, y=uvbMean)) + geom_point() + scale_x_log10() + scale_y_log10() + theme_bw()
```

![plot of chunk exercise2_2](./analysis_files/figure-html/exercise2_2.png) 
There are lots more options you can use to alter the appearance of these plots.

##Find candidate differentially expressed genes

We can find candidate differentially expressed genes by looking for genes with a large change between control and UVB samples. A common threshold used is log2 fold change more than 2 fold. We will calculate log2 fold change for all the genes and colour the genes with log2 fold change more than 2 fold on the plot.


```r
sum(countdata2$ctlMean > 0)
```

```
## [1] 474
```

```r
sum(countdata2$uvbMean > 0)
```

```
## [1] 528
```

```r
nrow(countdata2)
```

```
## [1] 56638
```

```r
countdata2 <- subset(countdata2, (countdata2$ctlMean > 0 | countdata2$uvbMean > 0))
#explain | operator meaning OR in this context?
nrow(countdata2)
```

```
## [1] 565
```


```r
countdata2$log2FC <- log2(countdata2$uvbMean / countdata2$ctlMean)
sum(countdata2$log2FC > 2)
```

```
## [1] 170
```

```r
sum(countdata2$log2FC < -2)
```

```
## [1] 46
```
Make a new column to store this information in.


```r
countdata2$outlier <- FALSE
countdata2$outlier[countdata2$log2FC > 2] <- TRUE
countdata2$outlier[countdata2$log2FC < -2] <- TRUE
```


```r
plot(countdata2$ctlMean, countdata2$uvbMean, log="xy")
```

```
## Warning: 91 x values <= 0 omitted from logarithmic plot
## Warning: 37 y values <= 0 omitted from logarithmic plot
```

```r
points(countdata2$ctlMean[countdata2$outlier==TRUE], countdata2$uvbMean[countdata2$outlier==TRUE], col="red")
```

![plot of chunk plot_outliers](./analysis_files/figure-html/plot_outliers.png) 


```r
ggplot(countdata2, aes(x=ctlMean, y=uvbMean, colour=outlier)) + geom_point() + scale_x_log10() + scale_y_log10() 
```

![plot of chunk ggplot_outliers](./analysis_files/figure-html/ggplot_outliers.png) 

## DESeq2 analysis

DESeq2 is an R package for analysis of RNAseq data. It is available from [Bioconductor](http://www.bioconductor.org/). Bioconductor is a project to provide tools for analysing high-throughput genomic data including RNA-seq, ChIP-seq and arrays. You can explore Bioconductor packages [here](http://www.bioconductor.org/packages/release/BiocViews.html#___Software). 


```r
#install and have a break to check everyone is up to date?
#explain bioconductor?
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```


```r
library(DESeq2)
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
# Convert to matrix
class(countdata)
```

```
## [1] "data.frame"
```

```r
countdata <- as.matrix(countdata)
class(countdata)
```

```
## [1] "matrix"
```

```r
head(countdata)
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
# construct coldata dataframe
#three replicates of control and UVB.
coldata <- data.frame(condition=c(rep("ctl", 3), rep("uvb",3)), row.names=colnames(countdata))
```

DESeq works on a particular type of object called a DESeqDataSet.


```r
#introduce how DESeq2 works - type of object it works on etc
# instantiate the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
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
## Order by adjusted p-value
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
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
head(resdata)
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
names(resdata)[1] <- "GeneID"
head(resdata)
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
sig <- subset(resdata, padj<0.05)
write.table(sig, file="results/sig.txt", sep="\t") #tab delim data
```

## Data Visualization

We can also do some exploratory plotting of the data.


```r
plotDispEsts(dds, main="Dispersion plot")
```

![plot of chunk plot_dispersion](./analysis_files/figure-html/plot_dispersion.png) 


```r
## Regularized log transformation for clustering/heatmaps, etc
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
## better heatmap with gplots
library(gplots)
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
## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")
```

![plot of chunk plot_pval_hist](./analysis_files/figure-html/plot_pval_hist.png) 



```r
#These are the plots that are most recognisable from papers
# MA Plot
par(pch=16)
with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x"))
```

```
## Warning: 56073 x values <= 0 omitted from logarithmic plot
```

```r
with(subset(res, padj<.05), points(baseMean, log2FoldChange, col="red", pch=16))
library(calibrate)
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
## Set point character
par(pch=16)
with(res, plot(log2FoldChange, -log10(pvalue), main="Volcano plot"))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), col="green"))
## Add legend
legend("topleft", legend=c("FDR<0.05", "|LFC|>1", "both"), pch=16, col=c("red","orange","green"))
## Label points
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=1))
```

![plot of chunk volcano_plot](./analysis_files/figure-html/volcano_plot.png) 
