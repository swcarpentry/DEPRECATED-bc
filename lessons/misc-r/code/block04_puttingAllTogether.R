## this is an ugly script -- not meant to be "source'd"
gDat <- read.delim("data/gapminderDataFiveYear.txt")
str(gDat)

## function that returns estimated intercept and slope from linear regression of
## lifeExp on year
## anticipated input: Gapminder data for one country
## see block03 script for development of this
jFun <- function(z) {
  jCoef <- coef(lm(lifeExp ~ I(year - yearMin), z))
  names(jCoef) <- c("intercept", "slope")
  return(jCoef)
}

## get intercept and slope for each country
library(plyr)
yearMin <- min(gDat$year)
gCoef <- ddply(gDat, .(country, continent), jFun)
str(gCoef)
tail(gCoef)

## reordering the continent factor rationally (vs. alphabetically)
library(lattice)
bwplot(slope ~ continent, gCoef)
## alphabetical order is rarely justified

gCoef$continent <- reorder(gCoef$continent, gCoef$slope)
str(gCoef)
bwplot(slope ~ continent, gCoef)
levels(gCoef$continent)
## much more logical
## tables and plots sorted this way will give more insight!

## drop Oceania ... too few countries
gCoef <- subset(gCoef, continent != "Oceania")
levels(gCoef$continent)
## Oceania is still there as factor level?!?
subset(gCoef, continent == "Oceania")
## but no observation have continent Oceania!
gCoef <- droplevels(gCoef)
## gets rid of unused factor levels
levels(gCoef$continent)

## how to store data.frames for posterity

## plain text is always good!
## we slowly built up the arguments below
## you may not want all of these?
write.table(gCoef, "results/gCoef.txt", quote = FALSE, row.names = FALSE, sep = "\t")

## bad news about plain text storage:
## upon re-import via read.table, factors levels revert to alphabetical order
levels(gCoef$continent) # rational level order
rm(gCoef)
gCoef <- read.delim("results/gCoef.txt")
levels(gCoef$continent) # back to alphabetical!!

## let's re-reorder the continent factor levels
gCoef$continent <- reorder(gCoef$continent, gCoef$slope)
bwplot(slope ~ continent, gCoef)
levels(gCoef$continent)

## how to store a data.frame in a plain text file WITH factor level order being
## preserved?
dput(gCoef, "results/gCoef_DPUT.txt")
rm(gCoef)
gCoef <- dget("results/gCoef_DPUT.txt")
levels(gCoef$continent)
## success!
## but the file written by dput() is not easy to browse or open
## in Excel or anywhere outside of R, for that matter

## if we abandon plain text, there is an R-specific binary format
saveRDS(gCoef, "results/gCoef.rds")
rm(gCoef)
gCoef <- readRDS("results/gCoef.rds")
levels(gCoef$continent)
## success!

## bottom line:
## I usually save as plain text AND via (dput|saveRDS|save)
## plain text is for long-term, language-agnostic storage
## *.rds is for my short-/medium-term reuse with R
