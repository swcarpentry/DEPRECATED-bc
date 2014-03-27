
## ----, echo=FALSE, message=FALSE-----------------------------------------
# Set default to hide all results and figures.
# Can change this manually in individual chunks.
# Must load knitr so opts_chunk is in search path.
library(knitr)
opts_chunk$set(results="hide", fig.show="hide", fig.keep="none")


## ----irisbase------------------------------------------------------------
# Load some data and look at the first few lines
data(iris)
head(iris)

# Make a basic scatter plot
plot(iris$Sepal.Length, iris$Petal.Length)


## ----installggplot2, eval=FALSE------------------------------------------
## # Only need to do this once
## install.packages("ggplot2")


## ----loadggplot2---------------------------------------------------------
library(ggplot2)


## ----diamondshead--------------------------------------------------------
data(diamonds)
head(diamonds)
str(diamonds)


## ----basediamond---------------------------------------------------------
plot(diamonds$carat, diamonds$price)


## ----qplot1--------------------------------------------------------------
qplot(carat, price, data=diamonds)


## ----clarcolor-----------------------------------------------------------
qplot(carat, price, data=diamonds, col=clarity)


## ----facetclar-----------------------------------------------------------
qplot(carat, price, data=diamonds, facets = ~clarity)


## ----facetcol------------------------------------------------------------
qplot(carat, price, data=diamonds, facets = ~color)


## ----facetclarcol--------------------------------------------------------
qplot(carat, price, data=diamonds, facets = clarity~color)


## ----facetclarcol_colcut-------------------------------------------------
qplot(carat, price, data=diamonds, facets=clarity~color, col=cut)


