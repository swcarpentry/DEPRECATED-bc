# ESA data visualization workshop
# This script contains all the code used to generate figures for the PDF. No need to peek in here 
# unless you're lost or having trouble completing the exercises.

# -------------------------------------

# Ignore these few lines, they are options to make sure code appears correctly on the slide.

## ----setup, include=FALSE------------------------------------------------
opts_chunk$set(cache=TRUE, message=FALSE)
# smaller font size for chunks
opts_chunk$set(size = 'footnotesize')
options(width = 60)
knit_hooks$set(inline = function(x) { 
  if (is.numeric(x)) return(knitr:::format_sci(x, 'latex')) 
  knitr:::hi_latex(x) 
}) 


## ----load, echo=FALSE, results='hide', warning=FALSE, message=FALSE------
require(ggplot2)
require(reshape2)
require(plyr)


## ----installation, tidy=FALSE, echo=TRUE, eval=FALSE---------------------
## install.packages("ggplot2", dependencies = TRUE)
## install.packages("plyr")
## install.packages("ggthemes")
## install.packages("reshape2")
## install.packages("gridExtra")
## install.packages("devtools")
## # Then a few packages to acquire data from the web to visualize
## install.packages("rfisheries")
## install.packages("rgbif")
## install.packages("taxize")
## # optional
## install_github("rWBclimate", "ropensci")


## ----some_data, tidy=FALSE, echo=TRUE------------------------------------
head(iris)


## ----plyrexample , fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold',  tidy=FALSE----
iris[1:2, ]
# Note the use of the '.' function to allow 'Species' to be used 
# without quoting
ddply(iris, .(Species), summarize, 
      mean.Sep.Wid = mean(Sepal.Width, na.rm = TRUE))


## ----reshapeexample , fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold',  tidy=FALSE----
iris[1:2, ]
df  <- melt(iris, id.vars = "Species")
df[1:2, ]


## ----reshapeexample2 , fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold',  tidy=FALSE----
df[1:2, ]
dcast(df, Species ~ variable, mean)


## ----data_summary, echo=TRUE,  tidy=FALSE--------------------------------
myplot <- ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width))
summary(myplot)


## ----first_plotb , eval=FALSE, tidy=FALSE--------------------------------
## ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width))
##  + geom_point()
## 
## myplot <- ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width))
## myplot + geom_point()


## ----first_plot , fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center',  tidy=FALSE----
ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width)) +
geom_point()


## ----first_plot_size , fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center',  tidy=FALSE----
ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width)) +
geom_point(size = 3)


## ----first_plot_color , fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center',  tidy=FALSE----
ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
geom_point(size = 3)


## ----first_plot_shape , fig.width=6, fig.height=3.5, out.width='.75\\linewidth', fig.show='hold', fig.align='center',  tidy=FALSE----
ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
geom_point(aes(shape = Species), size = 3)
# Why aes(shape = Species)?


## ----subset_data, tidy=FALSE, eval=FALSE---------------------------------
## # Make a small sample of the diamonds dataset
## d2 <- diamonds[sample(1:dim(diamonds)[1], 1000), ]


## ----ex1, echo = FALSE, out.width='.75\\linewidth', fig.height=4, fig.align='center'----
d2 <- diamonds[sample(1:dim(diamonds)[1], 1000), ]
ggplot(d2, aes(carat, price, color = color)) + geom_point() + theme_gray()


## ----boxplots1 , fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center', tidy=FALSE----
library(MASS)
ggplot(birthwt, aes(factor(race), bwt)) + geom_boxplot()


## ----boxplots2, echo=TRUE, tidy=FALSE------------------------------------
myplot <- ggplot(birthwt, aes(factor(race), bwt)) + geom_boxplot()
summary(myplot)


## ----facetgrid1, eval = TRUE, fig.width=4, fig.height=3, out.width='.75\\linewidth', fig.align='center', tidy=FALSE----
ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
geom_point() +
facet_grid(Species ~ .)


## ----facet_grid2, eval = TRUE, fig.width=4, fig.height=3, out.width='.75\\linewidth', fig.align='center', tidy= FALSE----
ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
geom_point() +
facet_grid(. ~ Species)


## ----facet_wrap, eval = TRUE, fig.width=4, fig.height=3, out.width='.75\\linewidth', fig.align='center', tidy = FALSE----
ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
geom_point() +
facet_wrap( ~ Species) # notice lack of .


## ----color_list2, eval = FALSE, fig.width=4, fig.height=6, out.width='.75\\linewidth', fig.align='center'----
## aes(color = variable) # mapping
## color = "black" # setting
## 
## # Or add it as a scale
## scale_fill_manual(values = c("color1", "color2"))


## ----color_list, eval = FALSE, fig.width=4, fig.height=5, out.width='.75\\linewidth'----
## library(RColorBrewer)
## display.brewer.all()


## ----barcolors, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center', tidy=FALSE----
df  <- melt(iris, id.vars = "Species")
ggplot(df, aes(Species, value, fill = variable)) +
geom_bar(stat = "identity", position = "dodge") +
scale_fill_brewer(palette = "Set1")


## ----facetgridcolors, eval = TRUE, fig.width=4, fig.height=3, out.width='.75\\linewidth', fig.align='center', tidy=FALSE----
ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
geom_point() +
facet_grid(Species ~ .) +
scale_color_manual(values = c("red", "green", "blue"))


## ----boxplots3, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center', tidy=FALSE----
library(MASS)
ggplot(birthwt, aes(factor(race), bwt)) +
geom_boxplot(width = .2) +
scale_y_continuous(labels = (paste0(1:4, " Kg")),
breaks = seq(1000, 4000, by = 1000))


## ----scale_list, eval=FALSE, tidy=FALSE----------------------------------
## scale_fill_discrete(); scale_colour_discrete()
## scale_fill_hue(); scale_color_hue()
## scale_fill_manual();  scale_color_manual()
## scale_fill_brewer(); scale_color_brewer()
## scale_linetype(); scale_shape_manual()


## ----histogr, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center',  tidy=FALSE----
h <- ggplot(faithful, aes(x = waiting))
h + geom_histogram(binwidth = 30, colour = "black")


## ----histogra, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center',  tidy=FALSE----
h <- ggplot(faithful, aes(x = waiting))
h + geom_histogram(binwidth = 8, fill = "steelblue",
colour = "black")


## ----line_setup, echo=FALSE----------------------------------------------
setwd('~/Github/esa_data_viz/Intro_lecture')


## ----linea , fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center', tidy=FALSE----
climate <- read.csv("../data/climate.csv", header = T)
ggplot(climate, aes(Year, Anomaly10y)) +
geom_line()


## ----lineb, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center',  tidy=FALSE----
ggplot(climate, aes(Year, Anomaly10y)) +
geom_ribbon(aes(ymin = Anomaly10y - Unc10y,
ymax = Anomaly10y + Unc10y),
fill = "blue", alpha = .1) +
geom_line(color = "steelblue")


## ----ex2, echo = FALSE, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.align='center', fig.show='hold'----
cplot <- ggplot(climate, aes(Year, Anomaly10y))
cplot <- cplot + geom_line(size = 0.7, color = "black")
cplot <- cplot + geom_line(aes(Year, Anomaly10y + Unc10y), linetype = "dashed", size = 0.7, color = "red")
cplot <- cplot + geom_line(aes(Year, Anomaly10y - Unc10y), linetype = "dashed", size = 0.7, color = "red")
cplot + theme_gray()


## ----barone, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center', tidy=FALSE----
ggplot(iris, aes(Species, Sepal.Length)) +
geom_bar(stat = "identity")


## ----bartwo, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center', tidy=FALSE----
df  <- melt(iris, id.vars = "Species")
ggplot(df, aes(Species, value, fill = variable)) +
geom_bar(stat = "identity")


## ----barthree, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center', tidy=FALSE----
ggplot(df, aes(Species, value, fill = variable)) +
geom_bar(stat = "identity", position = "dodge")


## ----barthree2, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center', tidy=FALSE----
ggplot(df, aes(Species, value, fill = variable)) + 
geom_bar(stat = "identity", position="dodge", color="black")


## ----ex3, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center', echo = FALSE----
ggplot(d2, aes(clarity, fill = cut)) +
geom_bar(position = "dodge",stat = "bin") + theme_gray()


## ----ex4, echo = FALSE, warning = FALSE, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.align='center', fig.show='hold'----
clim <- read.csv('../data/climate.csv', header = TRUE)
clim$sign <- ifelse(clim$Anomaly10y<0, FALSE, TRUE)
# or as simple as
# clim$sign <- clim$Anomaly10y < 0
ggplot(clim, aes(Year, Anomaly10y)) + geom_bar(stat = "identity", aes(fill = sign)) + theme_gray()


## ----densityone , eval=TRUE, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center',  tidy=FALSE----
ggplot(faithful, aes(waiting)) + geom_density()


## ----densityonefove , eval=TRUE, fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center',  tidy=FALSE----
ggplot(faithful, aes(waiting)) +
geom_density(fill = "blue", alpha = 0.1)


## ----densitytwo , , fig.width=6, fig.height=4, out.width='.75\\linewidth', fig.show='hold', fig.align='center', tidy=FALSE----
ggplot(faithful, aes(waiting)) +
geom_line(stat = "density")


## ----adding_stats,, eval = TRUE, fig.width=5, fig.height=4, out.width='.75\\linewidth', fig.align='center', tidy=FALSE----
ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
geom_point(aes(shape = Species), size = 3) +
geom_smooth(method = "lm")


## ----adding_stats2,, eval = TRUE, tidy=FALSE, fig.width=7, fig.height=4, out.width='.75\\linewidth', fig.align='center'----
ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
geom_point(aes(shape = Species), size = 3) +
geom_smooth(method = "lm") +
facet_grid(. ~ Species)


## ----theme_list, eval = FALSE, fig.width=4, fig.height=6, out.width='.75\\linewidth', fig.align='center'----
## + theme()
## # see ?theme() for more options


## ----facet_wrap_theme, eval = FALSE, fig.width=4, fig.height=6, out.width='.75\\linewidth', fig.align='center', tidy=FALSE----
## ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
## geom_point(size = 1.2, shape = 16) +
## facet_wrap( ~ Species) +
## theme(legend.key = element_rect(fill = NA),
## legend.position = "bottom",
## strip.background = element_rect(fill = NA),
## axis.title.y = element_text(angle = 0))


## ----facet_wrap_theme_execc, eval = TRUE, echo = FALSE, fig.width=4, fig.height=3, out.width='.75\\linewidth'----
ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
geom_point(size = 1.2, shape = 16) +
facet_wrap( ~ Species) +
theme(legend.key = element_rect(fill = NA),
legend.position = "bottom",
strip.background = element_rect(fill = NA),
axis.title.y = element_text(angle = 0))


## ----facet_wrap_theme_exec, eval = FALSE, echo = TRUE, fig.width=4, fig.height=3, out.width='.75\\linewidth', fig.align='center'----
## install.packages('ggthemes')
## library(ggthemes)
## # Then add one of these themes to your plot
##  + theme_stata()
##  + theme_excel()
##  + theme_wsj()
##  + theme_solarized()


## ----custom_plots, eval=FALSE, tidy=FALSE--------------------------------
## my_custom_plot <- function(df, title = "", ...) {
##     ggplot(df, ...) +
##     ggtitle(title) +
##     whatever_geoms() +
##     theme(...)
## }


## ----custom_plots2, eval=FALSE, tidy=FALSE-------------------------------
## plot1 <- my_custom_plot(dataset1, title = "Figure 1")


## ----pub0, eval = FALSE, out.width='.75\\linewidth'----------------------
## ggsave('~/path/to/figure/filename.png')


## ----pub1, eval = FALSE, out.width='.75\\linewidth'----------------------
## ggsave(plot1, file = "~/path/to/figure/filename.png")


## ----pub2, eval = FALSE, out.width='.75\\linewidth', tidy = FALSE--------
## ggsave(file = "/path/to/figure/filename.png", width = 6,
## height =4)


## ----pub3, eval = FALSE, out.width=".75\\linesewidth"--------------------
## ggsave(file = "/path/to/figure/filename.eps")
## ggsave(file = "/path/to/figure/filename.jpg")
## ggsave(file = "/path/to/figure/filename.pdf")


