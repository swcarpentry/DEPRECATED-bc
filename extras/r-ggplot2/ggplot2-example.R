
## ----, echo=FALSE, message=FALSE, eval=TRUE------------------------------
library(knitr)
# opts_chunk$set(results="hide", fig.show="hide", fig.keep="none")
opts_chunk$set(message=FALSE)


## ----irisbase------------------------------------------------------------
# Load some data and look at the first few lines
data(iris)
head(iris)

# Make a basic scatter plot
with(iris, plot(Sepal.Length, Petal.Length))


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
with(diamonds, plot(carat, price))


## ----qplot1--------------------------------------------------------------
qplot(carat, price, data = diamonds)


## ----clarcolor-----------------------------------------------------------
qplot(carat, price, data = diamonds, col = clarity)


## ----facetclar-----------------------------------------------------------
qplot(carat, price, data = diamonds, facets = ~ clarity)


## ----facetcol------------------------------------------------------------
qplot(carat, price, data = diamonds, facets = ~ color)


## ----facetclarcol, fig.height=10, fig.width=10---------------------------
qplot(carat, price, data = diamonds, facets = clarity ~ color)


## ----facetclarcol_colcut, fig.height=10, fig.width=10--------------------
qplot(carat, price, data = diamonds, facets = clarity ~ color, col = cut)


## ----ggdiamonds----------------------------------------------------------
# Using the qplot convenience function:
# qplot(carat, price, data = diamonds)

# Using ggplot:
ggplot(diamonds, aes(carat, price)) + geom_point()


## ----ggclarcol-----------------------------------------------------------
# Using the qplot convenience function:
# qplot(carat, price, data = diamonds, col = clarity) 
# Using ggplot:
ggplot(diamonds, aes(carat, price, col=clarity)) + geom_point()


## ----ggfacet, fig.height=10, fig.width=10--------------------------------
# Using the qplot convenience function:  
# qplot(carat, price, data = diamonds, facets = clarity ~ color)   
# Using ggplot:
ggplot(diamonds, aes(carat, price)) + geom_point() + facet_grid(clarity ~ color)


## ----ggalpha-------------------------------------------------------------
ggplot(diamonds, aes(carat, price)) + geom_point(alpha=1/5)


## ----gghexbin------------------------------------------------------------
ggplot(diamonds, aes(carat, price)) + geom_hex()


## ----density2d-----------------------------------------------------------
ggplot(diamonds, aes(carat, price)) + stat_density2d()


## ----stripplot-----------------------------------------------------------
ggplot(iris, aes(Species, Sepal.Length)) + geom_point()


## ----jitter--------------------------------------------------------------
ggplot(iris, aes(Species, Sepal.Length)) + geom_jitter()
ggplot(iris, aes(Species, Sepal.Length)) + geom_jitter(aes(shape=Species), size=3)
ggplot(iris, aes(Species, Sepal.Length)) + geom_jitter(aes(shape=Species, colour=Species), size=3)


## ----boxplot-------------------------------------------------------------
ggplot(iris, aes(Species, Sepal.Length)) + geom_boxplot()


## ----violinplot----------------------------------------------------------
ggplot(iris, aes(Species, Sepal.Length)) + geom_violin()


## ----violinpoints--------------------------------------------------------
ggplot(iris, aes(Species, Sepal.Length)) + geom_violin(alpha=.5) + geom_jitter(aes(shape=Species, colour=Species), size=3)


## ----irisloess-----------------------------------------------------------
ggplot(iris, aes(Sepal.Width, Sepal.Length, color=Species)) + geom_point() + geom_smooth()


## ----irislm--------------------------------------------------------------
ggplot(iris, aes(Sepal.Width, Sepal.Length, color=Species)) + geom_point() + geom_smooth(method="lm")


## ----irisfacet-----------------------------------------------------------
ggplot(iris, aes(Sepal.Width, Sepal.Length, color=Species)) + geom_point() + facet_wrap(~Species) + geom_smooth(method="lm")


## ------------------------------------------------------------------------
library(plyr)
str(baseball)
head(baseball)


## ----hrpoints------------------------------------------------------------
ggplot(baseball, aes(year, hr)) + geom_point()


## ----baseballmean--------------------------------------------------------
ggplot(baseball, aes(year, hr)) + stat_summary(fun.y=mean, geom="line")


## ----baseballmeansd------------------------------------------------------
mystats <- function(x) {
  data.frame(ymin=mean(x)-sd(x), ymax=mean(x)+sd(x))
}
ggplot(baseball, aes(year, hr)) + stat_summary(fun.y=mean, geom="line") + stat_summary(fun.data="mystats", geom="ribbon", alpha=.2, fill="steelblue")


## ----pricehisto----------------------------------------------------------
ggplot(diamonds, aes(price)) + geom_histogram(binwidth=200)


## ----ggfillhisto---------------------------------------------------------
ggplot(diamonds, aes(price, fill=clarity)) + geom_histogram(position="fill", binwidth=200)


## ----ggboxplots, fig.width=12--------------------------------------------
ggplot(diamonds, aes(cut, price)) + geom_boxplot(aes(fill=color)) + scale_y_log10()


## ----ggdepthdensity, warning=FALSE---------------------------------------
g <- ggplot(diamonds, aes(depth, fill=cut)) 
g <- g + geom_density(alpha=1/4) 
g <- g + xlim(55, 70)
g <- g + ggtitle("Table Depths by Cut Quality")
g <- g + xlab("Table Depth") + ylab("Density")
g


## ----ggsave, eval=FALSE--------------------------------------------------
## ggsave(filename="~/Desktop/table-depth-density.png", plot=g)


