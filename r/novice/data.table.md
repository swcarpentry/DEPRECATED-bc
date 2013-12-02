# Introduction to data.table

## Installation


```r
install.packages("data.table")
```



```r
library(data.table)
```


## Read in a file

`data.table` can read in files more efficiently thatn the standard functions in R, including `read.csv`, `read.table`, or `scan`. Let's take a look at the function `fread`. 


```r
df <- data.frame(a = 1:10^4, b = 1:10^4)
write.table(df, file = "test.csv", sep = ",", row.names = FALSE)
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

system.time(read.csv("test.csv"))
```

```
##    user  system elapsed 
##   0.033   0.001   0.038
```

```r
system.time(fread("test.csv"))
```

```
##    user  system elapsed 
##   0.002   0.000   0.002
```

```r
identical(data.frame(fread("test.csv")), read.csv("test.csv"))
```

```
## [1] TRUE
```


## Combine data.frames


```r
df1 <- data.frame(a = 1:4, b = letters[1:4])
df2 <- data.frame(a = 5:8, b = letters[5:8])
l <- list(df1, df2)
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

```r
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


This is much faster than `rbind`


```r
df3 <- data.frame(a = 1:10^7, b = 1:10^7)
df4 <- data.frame(a = 1:10^7, b = 1:10^7)
l <- list(df3, df4)
system.time(do.call(rbind, l))
```

```
##    user  system elapsed 
##   4.024   0.675   4.881
```

```r
system.time(rbindlist(l))
```

```
##    user  system elapsed 
##   0.062   0.002   0.065
```


## Fast manipulation


```r
DT <- data.table(x = rep(letters, 3847), v = rnorm(100022))
DT[, sum(v)]
```

```
## [1] -576.1
```

```r
DT[, sum(v), by = x]
```

```
##     x       V1
##  1: a   -6.244
##  2: b   17.962
##  3: c  -30.130
##  4: d -132.280
##  5: e   31.365
##  6: f  -48.264
##  7: g  -73.862
##  8: h  -87.046
##  9: i    3.829
## 10: j  -20.351
## 11: k   56.968
## 12: l    3.250
## 13: m  -47.952
## 14: n  -89.082
## 15: o   21.362
## 16: p -107.552
## 17: q    4.441
## 18: r  -44.574
## 19: s   37.492
## 20: t   41.824
## 21: u   14.918
## 22: v  -90.491
## 23: w  -57.714
## 24: x   73.694
## 25: y   19.280
## 26: z  -66.937
##     x       V1
```

```r

system.time(replicate(10, tapply(DT$v, DT$x, sum)))
```

```
##    user  system elapsed 
##   0.396   0.053   0.450
```

```r
system.time(replicate(10, DT[, sum(v), by = x]))
```

```
##    user  system elapsed 
##   0.069   0.008   0.076
```

