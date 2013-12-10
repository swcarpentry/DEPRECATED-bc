# Introduction to data.table

What is the `data.table` library for and why would you want to use it?  Doesn't base R come with data frames build in already? Turns out that there are some things that can be done MUCH faster and more easily with data.table. 

## Installation


```r
install.packages("data.table")
```



```r
library(data.table)
```


## Read in a file

Almost every analysis project in R will require the loading of data. `data.table` can read in files more efficiently than the standard functions in R, including `read.csv()`, `read.table()`, or `scan()`. This means that large files load fast.  Let's take a look at the function `fread()`. 


```r
# First we make a large data frame (10,000 rows long) and write it to a csv
# file, which we'll read in in the next step
df <- data.frame(a = 1:10^4, b = 1:10^4)
write.table(df, file = "test.csv", sep = ",", row.names = FALSE)

# instead of using read.csv(), we use fread(), a function provided in the
# data.frame library
fread("test.csv")
```

```
##            a     b
##     1:     1     1
##     2:     2     2
##     3:     3     3
##     4:     4     4
##     5:     5     5
##    ---            
##  9996:  9996  9996
##  9997:  9997  9997
##  9998:  9998  9998
##  9999:  9999  9999
## 10000: 10000 10000
```

```r

# We can track the amount of time a command takes to execute using the
# system.time() function
system.time(read.csv("test.csv"))
```

```
##    user  system elapsed 
##   0.013   0.000   0.014
```

```r
system.time(fread("test.csv"))
```

```
##    user  system elapsed 
##   0.001   0.000   0.001
```

```r

# much faster!

# and, just to prove they're identical once read in, we compare the two data
# sets to one another with identical()
identical(data.frame(fread("test.csv")), read.csv("test.csv"))
```

```
## [1] TRUE
```


## Combine data.frames

`data.table` can do more than just read in files though.  Another often-completed task is combining two data.frames.  Let compare the base R approach to the `data.table` version.


```r
# First we'll create two distinct data frames to join together
df1 <- data.frame(a = 1:4, b = letters[1:4])
df2 <- data.frame(a = 5:8, b = letters[5:8])

# put them into a list
l <- list(df1, df2)

# first, the base R version:
do.call(rbind, l)
```

```
##   a b
## 1 1 a
## 2 2 b
## 3 3 c
## 4 4 d
## 5 5 e
## 6 6 f
## 7 7 g
## 8 8 h
```

```r

# then the data.frame version, rbindlist()
rbindlist(l)
```

```
##    a b
## 1: 1 a
## 2: 2 b
## 3: 3 c
## 4: 4 d
## 5: 5 e
## 6: 6 f
## 7: 7 g
## 8: 8 h
```


This too is much faster with the data.table `rbindlist()` than with the base R `rbind()`.  To show the difference, we'll do the same procedure with data.frames 10 million rows long.


```r
df3 <- data.frame(a = 1:10^7, b = 1:10^7)
df4 <- data.frame(a = 1:10^7, b = 1:10^7)
l <- list(df3, df4)

system.time(do.call(rbind, l))
```

```
##    user  system elapsed 
##   2.132   0.451   2.583
```

```r
system.time(rbindlist(l))
```

```
##    user  system elapsed 
##   0.028   0.002   0.030
```

```r
# as you can see, rbindlist() is more than 80 times faster than base R's
# rbind() function
```


## Fast manipulation

Data table inherits from data.frame.  That means that `data.table` objects can also be used in functions that require data frames.  In addition, however, data.table objects can also be accessed with a completely different syntax, that is in many ways similar to SQL queries.  

The general syntax is (where DT is a data.table):    
`DT[where,select|update,group by][having][order by][ ]...[ ]`

This syntax enables all sorts of powerful manipulations and subsetting of data.tables, but one simple application of this that we describe here is that using this `data.table` syntax, you can perform fast manipulation of subsets of `data.table` objects.  In base R, the `tapply()` function is often used to calculate sums, means, or other statistics on the data in one column based on the categorical values in a second column, but in `data.table` this is built in to the syntax.  You can check out [this link](http://datatable.r-forge.r-project.org/datatable-faq.pdf) for more info.


```r
# first we make a data.table object and fill it with some random numbers
DT <- data.table(x = rep(letters, 3847), v = rnorm(100022))

# calculate the sum of the v column (leaving the first argument blank
# applies the sum() function over all rows)
DT[, sum(v)]
```

```
## [1] -123.8
```

```r

# calculate the sum of the values in the v column, grouped by the
# categorical values in the x column
DT[, sum(v), by = x]
```

```
##     x        V1
##  1: a    3.8061
##  2: b   46.1137
##  3: c  -35.6644
##  4: d  -18.8847
##  5: e -107.8458
##  6: f   75.0386
##  7: g  -50.6009
##  8: h  -87.8698
##  9: i  -50.3990
## 10: j  100.4895
## 11: k  -65.5616
## 12: l   35.7656
## 13: m  -93.0530
## 14: n    0.3912
## 15: o  -21.5754
## 16: p   63.9354
## 17: q   77.9314
## 18: r  -27.9210
## 19: s   -1.5271
## 20: t  -70.6666
## 21: u -155.4469
## 22: v   65.6844
## 23: w   67.6411
## 24: x  143.0250
## 25: y   37.0549
## 26: z  -53.6828
##     x        V1
```

```r

# and check out how fast it is!  More than 6 times faster than tapply().
# Because both functions are pretty fast, we replicate sthe function calls
# 10 times to show the difference.
system.time(replicate(10, tapply(DT$v, DT$x, sum)))
```

```
##    user  system elapsed 
##   0.233   0.027   0.260
```

```r
system.time(replicate(10, DT[, sum(v), by = x]))
```

```
##    user  system elapsed 
##   0.035   0.004   0.039
```

