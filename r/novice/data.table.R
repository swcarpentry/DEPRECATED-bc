
## ----install, eval=FALSE-------------------------------------------------
## install.packages("data.table")


## ----load----------------------------------------------------------------
library(data.table)


## ----fread---------------------------------------------------------------
# First we make a large data frame (10,000 rows long) and write it to a csv file, which we'll read in in the next step
df <- data.frame(a = 1:10^4, b = 1:10^4)
write.table(df, file="test.csv", sep=",", row.names=FALSE)

# instead of using read.csv(), we use fread(), a function provided in the data.frame library
fread("test.csv")

# We can track the amount of time a command takes to execute using the system.time() function
system.time( read.csv("test.csv") )
system.time( fread("test.csv") )

# much faster!

# and, just to prove they're identical once read in, we compare the two data sets to one another with identical()
identical(data.frame(fread("test.csv")), read.csv("test.csv"))


## ----rbindlist-----------------------------------------------------------
# First we'll create two distinct data frames to join together 
df1 <- data.frame(a = 1:4, b = letters[1:4])
df2 <- data.frame(a = 5:8, b = letters[5:8])

# put them into a list 
l <- list(df1, df2)

# first, the base R version:
do.call(rbind, l)

# then the data.frame version, rbindlist()
rbindlist(l)



## ------------------------------------------------------------------------
df3 <- data.frame(a = 1:10^7, b = 1:10^7)
df4 <- data.frame(a = 1:10^7, b = 1:10^7)
l <- list(df3, df4)

system.time( do.call(rbind, l) )
system.time( rbindlist(l) )
# as you can see, rbindlist() is more than 80 times faster than base R's rbind() function


## ------------------------------------------------------------------------
# first we make a data.table object and fill it with some random numbers
DT <- data.table(x = rep(letters, 3847), v = rnorm(100022)) 

# calculate the sum of the v column (leaving the first argument blank applies the sum() function over all rows)
DT[,sum(v)]

# calculate the sum of the values in the v column, grouped by the categorical values in the x column
DT[,sum(v),by=x]

# and check out how fast it is!  More than 6 times faster than tapply(). Because both functions are pretty fast, we replicate sthe function calls 10 times to show the difference.
system.time( replicate(10, tapply(DT$v,DT$x,sum)) )
system.time( replicate(10, DT[,sum(v),by=x]) )


