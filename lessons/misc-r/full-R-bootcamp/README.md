# R bootcamp material



This `master` set is based on material from a [2-day bootcamp](https://github.com/swcarpentry/2013-10-09-canberra) ([@karthik](https://github.com/karthik)) ran in Canberra, Australia in October 2013. The material is meant to cover the full range of topics typically covered in a regular Python bootcamp. 

Shell and Git are intentionally missing from this set. Please consult other lesson folders for those topics.

| Topic | Materials |
| ------ | -------- |
| Basics | Introduction to R, data types, best practices, seeking help, Using the RStudio IDE |
| Functions | Basics of control flow, scoping rules and functions in R |
| Data Manipulation | A full introduction to the apply family, dealing with IO in R, and a full hands on example of cleaning messy data in R |
| Data Visualization | A complete introduction to `ggplot2` |
| Testing | Documentation with `roxygen2`, Unit testing with `testthat`|
| Reproducible Research | knitr, make |


## Package installation

Please install the following packages:

```coffee
install.packages("devtools")
# You'll need other non-R dependencies before you can install devtools. Please see the additional_software.md page for more instructions.
install.packages(c("reshape2", "plyr", "ggplot2", "knitr", "testthat", "assertthat", stringr", "pander"))
```

## Notes for instructor teaching with this material.

* Please pull from the `bc` repo and decide what topics within R you would like to cover. Then retain the necessary folders and delete the rest. It might also help to number section folders in the order in which you might cover the material.

* Next, knit all of the `Rmd` files in each folder to make sure they parse correctly. Please avoid doing this the night before you teach in case any material needs to be updated to keep up with package changes.

* Commit the rendered files to your repo (NOT `bc`) so you have material to show on the projector while you teach. You might also consider doing a sed replacement for **```r** to  **coffee** for better syntax highlighting. 

* If you decide to teach Make in the context of `knitr`, please ensure that Make is correctly installed on all computers (see additional software under R-basics) section. Also keep in mind that the material here is intentionally sparse. There will be material in other lesson folders with more details that you might want to pull in.

