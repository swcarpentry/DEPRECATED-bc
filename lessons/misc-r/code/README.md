The individual scripts we wrote on Day 1 of an R-flavored SWC Boot Camp, run in two instances: NESCent and Duke in May 2013.

Here the files are referred to w/o reflecting the creating of subdirectories, which happened on Day 2, and is the current state of the repo.

Instructor has done whatever is necessary to ensure that RStudio and R launches as "virginally" as possible. Suppress usage of custom `.Rprofile`. Ensure will launch with empty workspace and history. Ensure will launch with user's home directory as working directory.

block01: Students fired up RStudio. Basic exploration of the environment RStudio provides. Entered commands live in the R console. Discussed organizing a project, especially an R analytical project, and R's notions of workspace and working directory. Created an RStudio project to use for the remainder of the bootcamp. Sent commands from the History to the source editor. Eventually saved these as the stand-alone script `block01_toyExample.R`. Practiced grooming code and using RStudio's facilities for sending code from the source editor to the console. The R markdown file `block01_basicsWorkspaceWorkingDirProject.rmd` (in `code/`) creates the html file `block01_basicsWorkspaceWorkingDirProject.html`; that gives Jenny's "script" for block 1. 
  * inputs: none
  * code: `block01_preProject.R`, `block01_postProject.R`, `block01_toyExample.R`
  * outputs: `avgX.txt`, `niftyPlot.pdf`
  
block02: Basic care and feeding of the most common R objects. Special emphasis on `data.frames`. `read.table` and friends for import. Bit of figure-making with the `lattice` package. Using the `subset()` and `with()` functions and the `data=` and `subset=` arguments found in many functions to do computations _in situ_ with added bonus of readable code. How to access various bits of various R objects, i.e. indexing. Accurate transcript can be found in the script `block02_careFeedingData.R`, which is NOT meant to be run as a whole -- it's for interactive use.
  * inputs: `gapminderDataFiveYear.txt`
  * code: `block02_careFeedingData.R`
  * output: none
  
block03: Data aggregation = doing something repetitive for various logical bits of an R object. E.g. taking means of rows in a matrix, computing statistical summaries for variables in a `data.frame`, fitting a model to sub-`data.frames` induced by separating the Gapminder data out by country. Used the `apply` family of functions in base R and also introduced the add-on package `plyr`. Accurate transcript can be found in the script `block03_dataAggregation.R`, which is NOT meant to be run as a whole -- it's for interactive use.
  * inputs: `gapminderDataFiveYear.txt`
  * code: `block03_dataAggregation.R`
  * output: none

block04: Sort of a capstone "putting it all together" piece. Revisiting country specific linear models of life expectancy against year. Before writing those results to file for later use, reordering the continent factor rationally (based on rate of life expectancy gains) and dropping Oceania (too few countries). Purpose was to demonstrate typical hygiene for factor variables. Different ways to write rectangular data to file with various pros/cons: `write.table`, `dput`, `saveRDS`, which have natural relationships with `read.table`, `dget`, `readRDS`. Accurate transcript of live work can be found in the script `block04_puttingAllTogether.R`, which is NOT meant to be run as a whole -- it's for interactive use. We did package some of our work nicely as scripts that could be `knit` and/or `source`'d or put into a pipeline (see below).
  * `block04_puttingAllTogether.R`
    - inputs: `gapminderDataFiveYear.txt`
    - output: none
  * `01_countrySpecificInterceptSlope.R`
    - inputs: `gapminderDataFiveYear.txt`
    - output: `gCoef.txt`, `gCoef.rds`
  * `02_slopeComparisonAsiaVsAmericas.R` (we used this to demonstrate the super-lightweight dynamic report generation capability of RStudio: "File --> Compile notebook", also available as a button; later accomplished same from command line in a `Makefile`)
    - inputs: `gCoef.rds`
    - outputs (after compiling notebook): `02_slopeComparisonAsiaVsAmericas.html`
  * `03_slopeComparisonAsiaVsAmericas.R` (we used this to demonstrate how a stand-alone script could leave files behind for later use, such as a PDF and the results of a two-sample t-test; essentially equivalent to `02_slopeComparionsAsiaVsAmericas.R` but optimized for running in a hands-off way using old-school techniques, e.g. `sink()`)
      - inputs: `gCoef.rds`
      - outputs: `slopes_AsiaVsAmericas.pdf`, `02_slopeComparisonAsiaVsAmericas_fromSink.txt`
      
block99: Challenges issued for further work.

