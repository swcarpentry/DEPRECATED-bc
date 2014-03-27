# Data visualization with ggplot2




## Review

Remember an earlier lesson where we introduced basic plotting commands using built-in data and built-in plotting tools. For example, we made a few plots using Edgar Anderson's famous iris dataset, which measured petal and sepal length and width for several different species of flower.


```coffee
# Load some data and look at the first few lines
data(iris)
head(iris)
```

```
##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 1          5.1         3.5          1.4         0.2  setosa
## 2          4.9         3.0          1.4         0.2  setosa
## 3          4.7         3.2          1.3         0.2  setosa
## 4          4.6         3.1          1.5         0.2  setosa
## 5          5.0         3.6          1.4         0.2  setosa
## 6          5.4         3.9          1.7         0.4  setosa
```

```coffee

# Make a basic scatter plot
plot(iris$Sepal.Length, iris$Petal.Length)
```

![plot of chunk irisbase](figure/irisbase.png) 


When we call plot in this way, we are using built-in, or *base* graphics. R's base graphics are powerful and nearly infinitely customizable. 

## ggplot2 

### The Diamonds dataset

Now let's look at a bigger dataset. We're going to be using a data visualization package called ggplot2 for drawing the plots, and the ggplot2 package comes with some data we're going to use for this example. 

Recall how to install and load packages. Install the package if you haven't already:


```coffee
# Only need to do this once
install.packages("ggplot2")
```


Then load it:


```coffee
library(ggplot2)
```


Now let's load the diamonds dataset and take a look at the first few rows and it's structure with commands we learned previously. To learn more about this dataset you can also run `?diamonds`. 


```coffee
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

```coffee
str(diamonds)
```

```
## 'data.frame':	53940 obs. of  10 variables:
##  $ carat  : num  0.23 0.21 0.23 0.29 0.31 0.24 0.24 0.26 0.22 0.23 ...
##  $ cut    : Ord.factor w/ 5 levels "Fair"<"Good"<..: 5 4 2 4 2 3 3 3 1 3 ...
##  $ color  : Ord.factor w/ 7 levels "D"<"E"<"F"<"G"<..: 2 2 2 6 7 7 6 5 2 5 ...
##  $ clarity: Ord.factor w/ 8 levels "I1"<"SI2"<"SI1"<..: 2 3 5 4 2 6 7 3 4 5 ...
##  $ depth  : num  61.5 59.8 56.9 62.4 63.3 62.8 62.3 61.9 65.1 59.4 ...
##  $ table  : num  55 61 65 58 58 57 57 55 61 61 ...
##  $ price  : int  326 326 327 334 335 336 336 337 337 338 ...
##  $ x      : num  3.95 3.89 4.05 4.2 4.34 3.94 3.95 4.07 3.87 4 ...
##  $ y      : num  3.98 3.84 4.07 4.23 4.35 3.96 3.98 4.11 3.78 4.05 ...
##  $ z      : num  2.43 2.31 2.31 2.63 2.75 2.48 2.47 2.53 2.49 2.39 ...
```


From this we can see this dataset has prices of nearly 54,000 diamonds along with various features about the diamonds, such as the weight, the quality of the cut, the color, the clarity, and measurements of various dimensions. 

### Plotting with ggplot2

If we wanted to do some exploratory data analysis we might start by plotting the price versus the weight of the diamond.


```coffee
plot(diamonds$carat, diamonds$price)
```

![plot of chunk basediamond](figure/basediamond.png) 


As we would expect there is definitely a relationship between the size of the diamond and its cost, but how do the other variables (cut, color, clarity) affect the price? We could examine the interrelationships of all these variables using base R graphics, but it would become extremely cumbersome using base R graphics. 

ggplot2 is a widely used R package that extends R's visualization capabilities. It takes the hassle out of things like creating legends, mapping other variables to scales like color, or faceting plots into small multiples. We'll learn about what all these things mean shortly. To start with, let's produce the same plot as before, but this time using ggplot2's `qplot` function:


```coffee
qplot(carat, price, data = diamonds)
```

![plot of chunk qplot1](figure/qplot1.png) 


The syntax is very similar to R's base graphics where you specify what's on the x and y axes, then give it the name of the data frame you want to use. We see again the strong relationship between the size of the diamond and its price, but the relationship doesn't appear linear. How does the diamond's clarity affect the weight-price relationship? 

### Faceting and scaling

One option we could use is to color-code the points by their clarity. Here, we pass another `col=` argument with the variable we'd like to use for color-coding:


```coffee
qplot(carat, price, data = diamonds, col = clarity)
```

![plot of chunk clarcolor](figure/clarcolor.png) 


Examining the plot you can see that poor-clarity diamonds (included, small inclusions, etc) fetch a lower price per carat relative to more clear diamonds (very small inclusions, internally flawless, etc). We can see that ggplot2 color-codes the points using a sensible default color scheme, and automatically draws a legend on the side for us. This requires a good deal of extra error-prone coding using base graphics.

However, with 54,000 points on this plot, there is a good deal of overplotting that obscures how clarity affects the nature of the weight-price relationship. How else might we visualize this data? This is where a *series of small multiples* is helpful. The idea of *small multiples* was popularized by data visualization expert Edward Tufte. The idea is that you create a large grid of small plots, where each plot shows a particular *facet* of the data. Here, each plot in the grid might be price vs. carat for each particular clarity level. You explain to your audience the axes and how to interpret each plot only once, and the audience will immediately understand the rest of the plots.

This can be accomplished easily using ggplot2:


```coffee
qplot(carat, price, data = diamonds, facets = ~clarity)
```

![plot of chunk facetclar](figure/facetclar.png) 


Here, the `facets` argument expects a forumla object that's constructed with the `~` operator. Here, we've plotted the price vs. weight separately for each level of clarity. We can see what we suspected before. With dirty diamonds (included, and perhaps small inclusions), the weight-price relationship is linear or slightly quadratic. Large diamonds can be purchased rather cheaply. But for very clear diamonds (internally flawless), the relationship is quadratic or even exponential.

Let's examine the weight-price relationship for various color ratings:


```coffee
qplot(carat, price, data = diamonds, facets = ~color)
```

![plot of chunk facetcol](figure/facetcol.png) 


Here we see that for whiter diamonds (D, E, F) the price rises more quickly with increasing weight than for yellower diamonds (H, I, J).

We can further facet the plot across two different categorical variables using the same syntax:


```coffee
qplot(carat, price, data = diamonds, facets = clarity ~ color)
```

![plot of chunk facetclarcol](figure/facetclarcol.png) 


Here we see that the price per carat rises very steeply for very white, very clear diamonds, while the relationship is nearly linear for yellower, more flawed diamonds. We can see that a perfect white diamond averages around $15,000 while a yellow included diamond can be had for only around $2,000.

Finally, we can combind both color-coding and faceting in the same plot. Let's use the same faceting scheme as last time, but color the points by the quality of the diamond's cut.


```coffee
qplot(carat, price, data = diamonds, facets = clarity ~ color, col = cut)
```

![plot of chunk facetclarcol_colcut](figure/facetclarcol_colcut.png) 


This color-coding reveals that clearer, whiter diamonds *generally* have higher quality cuts, but the relationship doesn't appear strong, visually. Looking down the plot toward clearer diamonds you start to see more "Ideal" cuts than at the top, which are the more included diamonds. 

What we've done here in addition to faceting is map a feature of the data (here, the cut quality) onto a scale (here, color). This behavior will work differently depending on whether you're looking at categorical or continuous variables. We can also map features to other *scales* such as `size=`, `shape=`, `linetype=`, or even transparency using `alpha=`. All of these different scales can be combined with each other or with facets, and give you an extremely powerful and easy-to-use graphical toolbox for exploratory data analysis. 

### Exercise

Now, install the `xx` package, and load the `xx` dataset. Using the techniques we've learned here, use faceting and scaling options to explore how xx, xx, and xx affect the relationship between xx and xx.
