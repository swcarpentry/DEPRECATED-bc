## this is an ugly script -- not meant to be "source'd"

## purpose is to collect representative, exhaustive code that might
## come up in this block

## importing data into a data.frame and then exploring it

## data import from URL
gdURL <- "http://www.stat.ubc.ca/~jenny/notOcto/STAT545A/examples/gapminder/data/gapminderDataFiveYear.txt"
gDat <- read.delim(gdURL)
gDat # not such a great idea ... too big!
str(gDat) # your main function for inspecting an object
## reading data from the web is a bit of a novelty
## more important: read from a locally stored plain text file

gDat <- read.delim("data/gapminderDataFiveYear.txt")
## read.table is the main workhorse for data import
## read.delim is merely a wrapper around that with
## certain arguments set to specific values
gDat <- read.delim("data/gapminderDataFiveYear.txt",
                   header = TRUE, sep = "\t")
str(gDat)
head(gDat, n = 10)
tail(gDat)
## tail often reveals more than head ... but even better is to look at some
## randomly selected rows!

## writing a function to peek into a data.frame
sample(x = nrow(gDat), size = 6)
sort(sample(x = nrow(gDat), size = 6))
gDat[sort(sample(x = nrow(gDat), size = 6)), ]
peek <- function(df) df[sort(sample(x = nrow(df), size = 6)), ]
peek
peek(gDat)
## for repeated use, use by other people, etc., one would want to upgrade this
## function! e.g. check if df is a data.frame

names(gDat)
dim(gDat)
nrow(gDat)
ncol(gDat)
head(rownames(gDat))
length(gDat)
str(gDat)
summary(gDat)
## dimnames(gDat) # not run; doesn't work well with so many rows

## student question:
## how would we make year categorical?
gDat$yearCat <- factor(gDat$year)
str(gDat)
summary(gDat)

## smell-testing dataset always includes making graphs!
## even though we aren't tackling graphics per se today,
## I can't resist!

## plotting packages: use lattice and/or ggplot2!
## JB knows lattice best, so that's what we'll use
library(lattice)
xyplot(lifeExp ~ year, gDat)
xyplot(lifeExp ~ gdpPercap, gDat)
xyplot(lifeExp ~ gdpPercap, gDat,
       subset = country == "Colombia")
xyplot(lifeExp ~ gdpPercap, gDat,
       subset = country == "Colombia",
       type = c("p", "r"))
xyplot(lifeExp ~ gdpPercap | continent, gDat,
       subset = year == 2007)
xyplot(lifeExp ~ gdpPercap, gDat,
       group = continent,
       subset = year == 2007, auto.key = TRUE)

## back to vetting a recently imported data.frame ...
str(gDat)
summary(gDat)
gDat$lifeExp # dollar sign is how you access 1 variable
summary(gDat$lifeExp)
densityplot(~ lifeExp, gDat)

## factors: the R objects you love to hate
## JB's advice: use them! learn how to handle safely
table(gDat$continent)
summary(gDat$continent)
levels(gDat$continent)
str(gDat)
## note that the actual *values* are integers, not the character codes,
## e.g. Africa or Europe, that are more user-visible
## never ever ever forget that factors are SPECIAL
nlevels(gDat$continent)
barchart(table(gDat$continent))
dotplot(table(gDat$continent), type = "h")

## if you want just some rows and/or just some variables, for inspection or to
## assign as a new object, use subset()
subset(gDat, subset = country == "Cambodia")
subset(gDat, subset = country %in% c("Japan", "Belgium"))
subset(gDat, subset = year == 1952)
subset(gDat, subset = country == "Uruguay",
       select = c(country, year, lifeExp))
subset(gDat, subset = country != "United States",
       select = c(country, lifeExp))

## notice the similarity of interface for subset(), xyplot(), lm(), ...
## this is not a coincidence!

## it is generally a good idea to take advantage of the "data" argument where
## you can specify a data.frame where the variables you want to work with are to
## be found

## functions with a "data" argument also often offer a "subset" argument

xyplot(lifeExp ~ year, gDat,
       subset = country == "Cambodia")
subset(gDat, subset = country == "Cambodia",
       select = c(country, year, lifeExp))
myFit <- lm(lifeExp ~ year, gDat,
            subset = country == "Cambodia")
summary(myFit)
(minYear <- min(gDat$year))
myFit <- lm(lifeExp ~ I(year - minYear), gDat,
            subset = country == "Cambodia")
summary(myFit)

## if a function does not have a "data" argument, with() can be handy
## if you are addicted to attach(), which is to be avoided!, with() is your
## methadone!
with(gDat, mean(lifeExp))
mean(gDat$lifeExp)

with(subset(gDat, country == "Cambodia"),
     cor(lifeExp, gdpPercap))
## much nicer than this:
foo <- gDat[gDat$country == "Cambodia", ]
cor(foo$lifeExp, foo$gdpPercap)

## sometimes a data.frame is just not what you need

## quick tour of other flavors of R objects and some important, sometimes
## surprising features of the language

## vectors are the garden variety R object
x <- 3 * 4
x
is.vector(x)
length(x)
x[1]
x[2]
x[2] <- 10
x
x[5] <- 3
x
x <- rnorm(25)
x

## R lurves to do vectorized computations
x^2

(y <- 1:3)
(z <- 3:7)
y + z
## recycling happens!
(y  <- 1:10)
(z <- 1:5)
y + z
## you are not always warned!

## atomic vectors MUST hold items of same "flavor"
## R will coerce to a least commong "flavor" if it must
x <- c("cabbage", pi, 0.3, TRUE)
x
class(x)

## lists are more general than vectors
x <- list("cabbage", pi, 0.3, TRUE)
x
class(x)

## how you can index vectors (and, with some modifications, lists, arrays, etc.)
x <- round(rnorm(8), 2)
set.seed(1) # how to make things *repeatably* random
x <- round(rnorm(8), 2)
x
names(x) <- letters[seq_along(x)]
x
x[4]
x[c(3, 6, 7)]
x[-1]
x < 0
x[x < 0]
which(x < 0)
x[c(TRUE, FALSE)]
x["b"]
x[c("b","f")]
