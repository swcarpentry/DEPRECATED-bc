plyr
========================================================

plyr is an R Package for Split-Apply-Combine workflows.  It's functional programming model encourages writing reusable functions which can be called across varied datasets and frees the analyst/scientist from needing to manage for loop indecies.

# plyr Functions

plyr has functions for operating on lists, data frames and arrays.  Each function performs:

1. A splitting operation
2. An "apply" where a function is called on each split in turn.
3. A recombination of output data as a single data object.

# Functions continued
The functions are named based on which type of object they expect as input (array, list or data frame) and which type of data structure should be returned as output.

This gives us 9 core functions **ply.  There are an adittional three functions which will only perform the split and apply steps, and not any combine step.  They're named by their input data type and represent null output by a `_` (see table)

|INPUT\OUTPUT|Array|Data frame|List|Discarded|
|------------|:---:|:--------:|:--:|--------:|
|Array|aaply|adply|alply|a_ply|
|Data frame|daply|ddply|dlply|d_ply|
|list|laply|ldply|llply|l_ply|

# What can we do with these functions?


```r
library(plyr)
iris_mean <- ddply(iris, .(Species), function(df) {
    data.frame(mean.petal_length = mean(df$Petal.Length))
})
iris_mean
```

```
##      Species mean.petal_length
## 1     setosa             1.462
## 2 versicolor             4.260
## 3  virginica             5.552
```


## Debrief

### What we did

1. Told plyr we wanted to send in a data.frame of iris data
2. Told it to split (or subset) the iris data.frame by "Species".
3. Wrote an anonymous function to process `df` for each sub-data frame.

### What plyr did for us

1. Coaxed input to a `data.frame` if it wasn't already
2. Performed subsetting operation on unique values of `Species`
3. In-turn fed those three subsets to our user defined function.
4. Took the output and coaxed it back to a `data.frame` if it wasn't already.

The steps that plyr performed for us are typically very fragile pieces of code in many scientific analysis workflows. plyr is able to give us an *abstraction layer* where we don't have to think about the mechanics of accounting for data types and specific loop indecies.

In this way plyr is a building block for simplifying complex scientific analysis code.

# A larger data set

Using the ozone dataset that comes with plyr we can explore a more complicated and multi-dimensional example.  The ozone data set is a 3-dimensional array with the first two dimensions being latitude and longitude, while the third is a monthly time series over 6 years. 


```r
overview_loc <- adply(ozone, c(1, 2), summary)
head(overview_loc)
```

```
##     lat   long Min. 1st Qu. Median Mean 3rd Qu. Max.
## 1 -21.2 -113.8  242     258    268  268     280  288
## 2 -18.7 -113.8  240     258    265  266     274  284
## 3 -16.2 -113.8  242     256    263  263     270  278
## 4 -13.7 -113.8  242     256    260  260     266  274
## 5 -11.2 -113.8  242     254    260  259     264  274
## 6  -8.7 -113.8  242     252    258  258     262  276
```

For each of the 576 locations we can see some summary univariate statistics across the time dimension (based on margins of c(1,2) ).

If we wanted to see the univariate statistics for all of space, at each of the 72 timesteps, we could run with a margin of 3.


```r
summarize_arrayvals <- function(a) {
    summary(as.vector(a))
}
overview_temporal <- adply(ozone, 3, summarize_arrayvals)
head(overview_temporal)
```

```
##   time Min. 1st Qu. Median Mean 3rd Qu. Max.
## 1    1  242     248    252  259     260  312
## 2    2  238     248    252  258     258  334
## 3    3  240     250    254  261     260  338
## 4    4  242     248    252  259     258  334
## 5    5  242     250    256  263     266  348
## 6    6  246     256    260  269     278  330
```


Notice how the output contains the output of summary with a cbind of the dimension that was used for splitting (lat,lon) and (time) respectively.  plyr is doing much of the heavy lifting and allowing you to simply write a function that operates on the splits in a generalized way.


# Aditional capabilities

## Multicore 

All plyr functions also have a `.parallel` method which will parallelize the inner loops of plyr operations.  In order to use `.parallel` you must use a package like `doMC` which is appropriate for your operating system. 

## Progress bars

By specifying `.progress = 'text' ` to your plyr function parameters, plyr will display a progress bar for long-running plyr operations.

## print
The `*_ply` functions have a `.print` option which will automatically wrap your function in a `print()` statement so outputs get written to display devices such as PDFs, PNGs etc.


# More information
The author of plyr, Hadely Wickham, has an excellent and freely availble paper called ["The Split-Apply-Combine Strategy for Data Analysis"](http://www.jstatsoft.org/v40/i01/paper).  It is an excellent reference for further learning and examples.



