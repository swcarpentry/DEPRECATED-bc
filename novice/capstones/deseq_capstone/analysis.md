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
countFile <- "data/counts.txt"
## Read in the data
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
## DDX11L1                                    337
## WASH7P                                     319
## MIR1302-10                                 313
## FAM138A                                    321
## OR4G4P                                     322
## OR4G11P                                    343
##            ctl2.fastq_tophat.accepted_hits.bam
## DDX11L1                                    322
## WASH7P                                     377
## MIR1302-10                                 348
## FAM138A                                    333
## OR4G4P                                     311
## OR4G11P                                    339
##            ctl3.fastq_tophat.accepted_hits.bam
## DDX11L1                                    325
## WASH7P                                     299
## MIR1302-10                                 328
## FAM138A                                    340
## OR4G4P                                     364
## OR4G11P                                    335
##            uvb1.fastq_tophat.accepted_hits.bam
## DDX11L1                                    131
## WASH7P                                     213
## MIR1302-10                                 250
## FAM138A                                    253
## OR4G4P                                     378
## OR4G11P                                    418
##            uvb2.fastq_tophat.accepted_hits.bam
## DDX11L1                                    158
## WASH7P                                     194
## MIR1302-10                                 247
## FAM138A                                    234
## OR4G4P                                     315
## OR4G11P                                    450
##            uvb3.fastq_tophat.accepted_hits.bam
## DDX11L1                                    136
## WASH7P                                     203
## MIR1302-10                                 263
## FAM138A                                    276
## OR4G4P                                     405
## OR4G11P                                    437
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
It contains information about genes (one gene per row) with the gene positions in the first five columns and then information about the number of reads aligning to the gene in each experimental sample. We don't need the information on gene position, so we can remove it from the data frame.


```r
# Remove first five columns (chr, start, end, strand, length)
countData <- countData[ ,-(1:5)]
head(countData)
```

```
##            ctl1.fastq_tophat.accepted_hits.bam
## DDX11L1                                    337
## WASH7P                                     319
## MIR1302-10                                 313
## FAM138A                                    321
## OR4G4P                                     322
## OR4G11P                                    343
##            ctl2.fastq_tophat.accepted_hits.bam
## DDX11L1                                    322
## WASH7P                                     377
## MIR1302-10                                 348
## FAM138A                                    333
## OR4G4P                                     311
## OR4G11P                                    339
##            ctl3.fastq_tophat.accepted_hits.bam
## DDX11L1                                    325
## WASH7P                                     299
## MIR1302-10                                 328
## FAM138A                                    340
## OR4G4P                                     364
## OR4G11P                                    335
##            uvb1.fastq_tophat.accepted_hits.bam
## DDX11L1                                    131
## WASH7P                                     213
## MIR1302-10                                 250
## FAM138A                                    253
## OR4G4P                                     378
## OR4G11P                                    418
##            uvb2.fastq_tophat.accepted_hits.bam
## DDX11L1                                    158
## WASH7P                                     194
## MIR1302-10                                 247
## FAM138A                                    234
## OR4G4P                                     315
## OR4G11P                                    450
##            uvb3.fastq_tophat.accepted_hits.bam
## DDX11L1                                    136
## WASH7P                                     203
## MIR1302-10                                 263
## FAM138A                                    276
## OR4G4P                                     405
## OR4G11P                                    437
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

An easier way to do this, especially for files with many columns, is to use the `gsub` command to strip out the extra information. This is also more robust to introduced errors, for example if the column order changes at some point in the future.

```r
## Using gsub -- reproducible
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
## DDX11L1     337  322  325  131  158  136
## WASH7P      319  377  299  213  194  203
## MIR1302-10  313  348  328  250  247  263
## FAM138A     321  333  340  253  234  276
## OR4G4P      322  311  364  378  315  405
## OR4G11P     343  339  335  418  450  437
```

### Exercise 1
Find the gene with the highest expression in any sample. Extract  the expression data for this gene for all samples. In which sample does it have the highest expression? What is the function of the gene? Can you suggest why this is the top expressed gene?


```r
max(apply(countData, 1, max)) #max expression is 7013
```

```
## [1] 450
```

```r
topGene <- which.max(apply(countData, 1, max)) #gene is EEF1A1P9
countData[topGene, ] #get other sample data - max is in uvb1
```

```
##         ctl1 ctl2 ctl3 uvb1 uvb2 uvb3
## OR4G11P  343  339  335  418  450  437
```

```r
#this is a pseudogene - maybe an artefact of only aligning reads to a single chromosome?
```

## Data investigation using base R

We can investigate this data a bit more using some of the basic R functions before going on to use more sophisticated analysis tools.

We will calculate the mean for each gene for each condition. First make a copy of the data, because we'll need it later. We will work on the copy.


```r
countData2 <- countData

#get Control columns
colnames(countData2)
```

```
## [1] "ctl1" "ctl2" "ctl3" "uvb1" "uvb2" "uvb3"
```

```r
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
## DDX11L1     337  322  325
## WASH7P      319  377  299
## MIR1302-10  313  348  328
## FAM138A     321  333  340
## OR4G4P      322  311  364
## OR4G11P     343  339  335
```

```r
apply(countData2[, ctlCols], 1, mean) 
```

```
##    DDX11L1     WASH7P MIR1302-10    FAM138A     OR4G4P    OR4G11P 
##      328.0      331.7      329.7      331.3      332.3      339.0
```

```r
# here we'll use rowMeans instead, its a convenient shortcut, and also faster!
countData2$ctlMean <- rowMeans(countData2[, ctlCols])

#same for uvb
uvbCols <- grep("uvb", colnames(countData2))
countData2$uvbMean <- rowMeans(countData2[, uvbCols])
```

Plot the mean expression of each gene in control against the UVB sample mean. Look for outliers.


```r
plot(countData2$ctlMean, countData2$uvbMean)
```

![plot of chunk plot_means](./analysis_files/figure-html/plot_means.png) 


```r
library(ggplot2)
ggplot(countData2, aes(x=ctlMean, y=uvbMean)) + geom_point()
```

![plot of chunk ggplot_means](./analysis_files/figure-html/ggplot_means.png) 

### Exercise 2
How could you make this plot more informative and look more professional? Hint: try using a log scale. Try changing colours, transparencies, sizes, or shapes of points. 

`help(par)` will give you information on lots of graphical parameters that can be set. Help for ggplot2 can be found [here](http://docs.ggplot2.org/current/).


```r
plot(countData2$ctlMean, countData2$uvbMean, log="xy")
```

![plot of chunk exercise2_1](./analysis_files/figure-html/exercise2_1.png) 


```r
ggplot(countData2, aes(x=ctlMean, y=uvbMean)) + geom_point() + scale_x_log10() + scale_y_log10() + theme_bw()
```

![plot of chunk exercise2_2](./analysis_files/figure-html/exercise2_2.png) 
There are lots more options you can use to alter the appearance of these plots.

##Find candidate differentially expressed genes

We can find candidate differentially expressed genes by looking for genes with a large change between control and UVB samples. A common threshold used is log2 fold change more than 2 fold. We will calculate log2 fold change for all the genes and colour the genes with log2 fold change more than 2 fold on the plot.


```r
sum(countData2$ctlMean > 0)
```

```
## [1] 6
```

```r
sum(countData2$uvbMean > 0)
```

```
## [1] 6
```

```r
nrow(countData2)
```

```
## [1] 6
```

```r
countData2 <- subset(countData2, (countData2$ctlMean > 0 | countData2$uvbMean > 0))
#explain | operator meaning OR in this context?
nrow(countData2)
```

```
## [1] 6
```


```r
countData2$log2FC <- log2(countData2$uvbMean / countData2$ctlMean)
sum(countData2$log2FC > 2)
```

```
## [1] 0
```

```r
sum(countData2$log2FC < -2)
```

```
## [1] 0
```
Make a new column to store this information in.


```r
countData2$outlier <- FALSE
countData2$outlier[countData2$log2FC > 2] <- TRUE
countData2$outlier[countData2$log2FC < -2] <- TRUE
```


```r
plot(countData2$ctlMean, countData2$uvbMean, log="xy")
points(countData2$ctlMean[countData2$outlier==TRUE], countData2$uvbMean[countData2$outlier==TRUE], col="red")
```

![plot of chunk plot_outliers](./analysis_files/figure-html/plot_outliers.png) 


```r
ggplot(countData2, aes(x=ctlMean, y=uvbMean, colour=outlier)) + geom_point() + scale_x_log10() + scale_y_log10() 
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
## DDX11L1     337  322  325  131  158  136
## WASH7P      319  377  299  213  194  203
## MIR1302-10  313  348  328  250  247  263
## FAM138A     321  333  340  253  234  276
## OR4G4P      322  311  364  378  315  405
## OR4G11P     343  339  335  418  450  437
```

```r
# construct colData dataframe
#three replicates of control and UVB.
colData <- data.frame(condition=c(rep("ctl", 3), rep("uvb",3)), row.names=colnames(countData))
```

DESeq works on a particular type of object called a DESeqDataSet.


```r
#introduce how DESeq2 works - type of object it works on etc
# instantiate the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
```

```
## Warning: some variables in design formula are characters, converting to
## factors
```

```r
dds
```

```
## class: DESeqDataSet 
## dim: 6 6 
## exptData(0):
## assays(1): counts
## rownames(6): DDX11L1 WASH7P ... OR4G4P OR4G11P
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
```

```
## Warning: the parametric fit of dispersion estimates over the mean of counts
## failed, which occurs when the trend is not well captured by the
## function y = a/x + b. A local regression fit is automatically performed,
## and the analysis can continue. You can specify fitType='local' or 'mean'
## to avoid this message if re-running the same data.
## When using local regression fit, the user should examine plotDispEsts(dds)
## to make sure the fitted line is not sharply curving up or down based on
## the position of individual points.
## Warning: Estimated rdf < 1.0; not estimating variance
```

```
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
## DDX11L1        225.4      -0.800662   0.11188  -7.15620 8.295e-13
## WASH7P         261.6      -0.316571   0.10258  -3.08609 2.028e-03
## MIR1302-10     289.2      -0.001925   0.09442  -0.02039 9.837e-01
## FAM138A        290.4      -0.005237   0.09441  -0.05546 9.558e-01
## OR4G4P         354.2       0.498655   0.09741   5.11892 3.073e-07
## OR4G11P        397.6       0.725214   0.07570   9.58062 9.647e-22
##                 padj
##            <numeric>
## DDX11L1    2.488e-12
## WASH7P     3.042e-03
## MIR1302-10 9.837e-01
## FAM138A    9.837e-01
## OR4G4P     6.146e-07
## OR4G11P    5.788e-21
```

```r
table(res$padj<0.05)
```

```
## 
## FALSE  TRUE 
##     2     4
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
##             baseMean log2FoldChange     lfcSE      stat    pvalue
##            <numeric>      <numeric> <numeric> <numeric> <numeric>
## OR4G11P        397.6       0.725214   0.07570   9.58062 9.647e-22
## DDX11L1        225.4      -0.800662   0.11188  -7.15620 8.295e-13
## OR4G4P         354.2       0.498655   0.09741   5.11892 3.073e-07
## WASH7P         261.6      -0.316571   0.10258  -3.08609 2.028e-03
## MIR1302-10     289.2      -0.001925   0.09442  -0.02039 9.837e-01
## FAM138A        290.4      -0.005237   0.09441  -0.05546 9.558e-01
##                 padj
##            <numeric>
## OR4G11P    5.788e-21
## DDX11L1    2.488e-12
## OR4G4P     6.146e-07
## WASH7P     3.042e-03
## MIR1302-10 9.837e-01
## FAM138A    9.837e-01
```

Combine DEseq results with the original counts data. Write significant results to a file.


```r
resData <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
head(resData)
```

```
##    Row.names baseMean log2FoldChange   lfcSE     stat    pvalue      padj
## 1    OR4G11P    397.6       0.725214 0.07570  9.58062 9.647e-22 5.788e-21
## 2    DDX11L1    225.4      -0.800662 0.11188 -7.15620 8.295e-13 2.488e-12
## 3     OR4G4P    354.2       0.498655 0.09741  5.11892 3.073e-07 6.146e-07
## 4     WASH7P    261.6      -0.316571 0.10258 -3.08609 2.028e-03 3.042e-03
## 5 MIR1302-10    289.2      -0.001925 0.09442 -0.02039 9.837e-01 9.837e-01
## 6    FAM138A    290.4      -0.005237 0.09441 -0.05546 9.558e-01 9.837e-01
##    ctl1  ctl2  ctl3  uvb1  uvb2  uvb3
## 1 313.1 288.2 292.5 480.9 541.6 469.3
## 2 307.6 273.7 283.8 150.7 190.2 146.1
## 3 293.9 264.4 317.9 434.9 379.1 435.0
## 4 291.2 320.5 261.1 245.1 233.5 218.0
## 5 285.7 295.8 286.4 287.6 297.3 282.5
## 6 293.0 283.1 296.9 291.1 281.6 296.4
```

```r
names(resData)[1] <- "GeneID"
head(resData)
```

```
##       GeneID baseMean log2FoldChange   lfcSE     stat    pvalue      padj
## 1    OR4G11P    397.6       0.725214 0.07570  9.58062 9.647e-22 5.788e-21
## 2    DDX11L1    225.4      -0.800662 0.11188 -7.15620 8.295e-13 2.488e-12
## 3     OR4G4P    354.2       0.498655 0.09741  5.11892 3.073e-07 6.146e-07
## 4     WASH7P    261.6      -0.316571 0.10258 -3.08609 2.028e-03 3.042e-03
## 5 MIR1302-10    289.2      -0.001925 0.09442 -0.02039 9.837e-01 9.837e-01
## 6    FAM138A    290.4      -0.005237 0.09441 -0.05546 9.558e-01 9.837e-01
##    ctl1  ctl2  ctl3  uvb1  uvb2  uvb3
## 1 313.1 288.2 292.5 480.9 541.6 469.3
## 2 307.6 273.7 283.8 150.7 190.2 146.1
## 3 293.9 264.4 317.9 434.9 379.1 435.0
## 4 291.2 320.5 261.1 245.1 233.5 218.0
## 5 285.7 295.8 286.4 287.6 297.3 282.5
## 6 293.0 283.1 296.9 291.1 281.6 296.4
```

```r
sig <- subset(resData, padj<0.05)
dir.create("results")
```

```
## Warning: 'results' already exists
```

```r
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
##             ctl1  ctl2  ctl3  uvb1  uvb2  uvb3
## DDX11L1    7.947 7.888 7.906 7.636 7.724 7.624
## WASH7P     8.078 8.126 8.026 7.998 7.977 7.947
## MIR1302-10 8.170 8.187 8.171 8.173 8.189 8.165
## FAM138A    8.186 8.170 8.192 8.183 8.167 8.191
## OR4G4P     8.369 8.321 8.406 8.565 8.493 8.565
## OR4G11P    8.500 8.463 8.470 8.718 8.785 8.705
```

```r
assay(rld)[1:5,1:5]
```

```
##             ctl1  ctl2  ctl3  uvb1  uvb2
## DDX11L1    7.947 7.888 7.906 7.636 7.724
## WASH7P     8.078 8.126 8.026 7.998 7.977
## MIR1302-10 8.170 8.187 8.171 8.173 8.189
## FAM138A    8.186 8.170 8.192 8.183 8.167
## OR4G4P     8.369 8.321 8.406 8.565 8.493
```

```r
t(assay(rld))[1:5,1:5]
```

```
##      DDX11L1 WASH7P MIR1302-10 FAM138A OR4G4P
## ctl1   7.947  8.078      8.170   8.186  8.369
## ctl2   7.888  8.126      8.187   8.170  8.321
## ctl3   7.906  8.026      8.171   8.192  8.406
## uvb1   7.636  7.998      8.173   8.183  8.565
## uvb2   7.724  7.977      8.189   8.167  8.493
```

```r
dist(t(assay(rld)))
```

```
##         ctl1    ctl2    ctl3    uvb1    uvb2
## ctl2 0.10004                                
## ctl3 0.08184 0.13488                        
## uvb1 0.43456 0.45264 0.40115                
## uvb2 0.39578 0.42692 0.37861 0.13597        
## uvb3 0.44902 0.46984 0.40809 0.05509 0.15427
```

```r
as.matrix(dist(t(assay(rld))))
```

```
##         ctl1   ctl2    ctl3    uvb1   uvb2    uvb3
## ctl1 0.00000 0.1000 0.08184 0.43456 0.3958 0.44902
## ctl2 0.10004 0.0000 0.13488 0.45264 0.4269 0.46984
## ctl3 0.08184 0.1349 0.00000 0.40115 0.3786 0.40809
## uvb1 0.43456 0.4526 0.40115 0.00000 0.1360 0.05509
## uvb2 0.39578 0.4269 0.37861 0.13597 0.0000 0.15427
## uvb3 0.44902 0.4698 0.40809 0.05509 0.1543 0.00000
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
```

![plot of chunk volcano_plot](./analysis_files/figure-html/volcano_plot.png) 

```r
## Label points
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=1))
```
