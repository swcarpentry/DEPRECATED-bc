# Introduction to data.table

What is the `data.table` package for and why would you want to use it?  Doesn't base R come with data frames build in already? Turns out that there are some things that can be done MUCH faster and more easily with data.table. 

## Installation


```r
install.packages(c("data.table", "microbenchmark"))
```



```r
library(data.table)
library(microbenchmark)
```


## Read in a file

Almost every analysis project in R will require the loading of data. `data.table` can read in files more efficiently than the standard functions in R, including `read.csv()`, `read.table()`, or `scan()`. This means that large files load fast.  Let's take a look at the function `fread()`. 


```r
# First we make a data frame (10,000 rows long) and write it to a csv file,
# which we'll read in in the next step
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
# microbenchmark() function
(op <- microbenchmark(read.csv("test.csv"), fread("test.csv"), times = 100))
```

```
## Unit: microseconds
##                  expr     min    lq median    uq   max neval
##  read.csv("test.csv") 19461.3 20927  22422 24442 38694   100
##     fread("test.csv")   916.9  1055   1220  1355  2034   100
```

```r

# Visualize the difference in timing
plot(op)
```

![plot of chunk fread](figure/fread.png) 

```r

# much faster!

# and, just to prove they're identical once read in, we compare the two data
# sets to one another with identical()
identical(data.frame(fread("test.csv")), read.csv("test.csv"))
```

```
## [1] TRUE
```


If you did this with a larger dataset, say 10^6 rows long (ca. 13.8 MB), then the speed up becomes really meaningful (almost 60x faster):

```
(op <- microbenchmark(
  read.csv("test.csv"),
  fread("test.csv")))

Unit: milliseconds
                 expr        min         lq    median        uq       max neval
 read.csv("test.csv") 2844.94893 3282.67950 3440.3038 3708.6769 5758.1040   100
    fread("test.csv")   87.05953   94.19449  131.1998  169.6652  429.7068   100
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

# then the data.table version, rbindlist()
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
df3 <- data.frame(a = 1:10^6, b = 1:10^6)
df4 <- data.frame(a = 1:10^6, b = 1:10^6)
l <- list(df3, df4)

system.time(do.call(rbind, l))
```

```
##    user  system elapsed 
##   0.508   0.063   0.572
```

```r
system.time(rbindlist(l))
```

```
##    user  system elapsed 
##   0.007   0.000   0.008
```

```r

microbenchmark(do.call(rbind, l), rbindlist(l))
```

```
## Unit: milliseconds
##               expr     min     lq median    uq   max neval
##  do.call(rbind, l) 305.671 490.96 543.48 609.5 959.8   100
##       rbindlist(l)   6.157  12.16  14.34  18.4 351.0   100
```

```r
# as you can see, rbindlist() is many times faster than base R's rbind()
# function
```


## Fast manipulation

Data table inherits from data.frame.  That means that `data.table` objects can also be used in functions that require data frames.  In addition, however, data.table objects can also be accessed with a completely different syntax, that is in many ways similar to SQL queries.  

The general syntax is (where DT is a data.table):    
`DT[where,select|update,group by][having][order by][ ]...[ ]`

This syntax enables all sorts of powerful manipulations and subsetting of data.tables, but one simple application of this that we describe here is that using this `data.table` syntax, you can perform fast manipulation of subsets of `data.table` objects.  In base R, the `tapply()` function is often used to calculate sums, means, or other statistics on the data in one column based on the categorical values in a second column, but in `data.table` this is built in to the syntax.  You can check out [this link](http://datatable.r-forge.r-project.org/datatable-faq.pdf) for more info.

We'll need the `diamonds` dataset from the `ggplot2` package.


```r
install.packages("ggplot2")
```




```r
# first we make a data.table object and fill it with some random numbers
library(ggplot2)
data(diamonds)
head(diamonds)
```

```
##   carat       cut color clarity depth table price    x    y    z
## 1  0.23     Ideal     E     SI2  61.5    55   326 3.95 3.98 2.43
## 2  0.21   Premium     E     SI1  59.8    61   326 3.89 3.84 2.31
## 3  0.23      Good     E     VS1  56.9    65   327 4.05 4.07 2.31
## 4  0.29   Premium     I     VS2  62.4    58   334 4.20 4.23 2.63
## 5  0.31      Good     J     SI2  63.3    58   335 4.34 4.35 2.75
## 6  0.24 Very Good     J    VVS2  62.8    57   336 3.94 3.96 2.48
```

```r
diamonds_dt <- data.table(diamonds)
# DT <- data.table(x = rep(letters, 3847), v = rnorm(100022))

# calculate the sum of the Sepal.Length column (leaving the first argument
# blank applies the sum() function over all rows)
diamonds_dt[, sum(price)]
```

```
## [1] 212135217
```

```r

# calculate the sum of the values in the Sepal.Length column, grouped by the
# categorical values in the Species column
diamonds_dt[, sum(price), by = color]
```

```
##    color       V1
## 1:     E 30142944
## 2:     I 27608146
## 3:     J 14949281
## 4:     H 37257301
## 5:     F 35542866
## 6:     G 45158240
## 7:     D 21476439
```

```r

# More than 3 times faster than tapply(). Because both functions are pretty
# fast, we replicate sthe function calls 100 times to show the difference.
microbenchmark(tapply = tapply(diamonds$price, diamonds$color, sum), data.table = diamonds_dt[, 
    sum(price), by = color], times = 100)
```

```
## Unit: milliseconds
##        expr    min     lq median     uq   max neval
##      tapply 17.390 19.728 21.609 25.240 43.16   100
##  data.table  5.457  5.868  6.438  7.231 34.44   100
```



## Some common tasks in both Base R and data.table

Task | Base R  | data.table
------------- | ------------- | -------------
Creation | data.frame(...) | data.table(...)
Select column | x[,"column"]  | x[,column]
First 6 rows | head(x)  | head(x)
Grouping | tapply(x$col1,x$col2,sum) | x[,sum(col1),by=col2]
Combine two data.frame's | rbind(df1,df2) | rbindlist(list(df1,df2))
Merge data.frame's | merge(df1,df2) | DT1[DT2]
