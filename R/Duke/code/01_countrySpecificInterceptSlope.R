library(plyr) # ddply()
gDat <- read.delim("data/gapminderDataFiveYear.txt")

## function that returns estimated intercept and slope from linear regression of
## lifeExp on year
## anticipated input: Gapminder data for one country
yearMin <- min(gDat$year)
jFun <- function(z) {
  jCoef <- coef(lm(lifeExp ~ I(year - yearMin), z))
  names(jCoef) <- c("intercept", "slope")
  return(jCoef)
}

gCoef <- ddply(gDat, .(country, continent), jFun)

## reorder continent factor levels to reflect rate of life expectancy gains
## (slowest growth to largest)
gCoef$continent <-
  reorder(gCoef$continent, gCoef$slope)

## drop Oceania ... too few countries
gCoef <- droplevels(subset(gCoef, continent != "Oceania"))

## store in plain text
write.table(gCoef, "results/gCoef.txt", quote = FALSE,
            row.names = FALSE, sep = "\t")

## store in R-specific binary format
## will preserve factor level order
saveRDS(gCoef, "results/gCoef.rds")
