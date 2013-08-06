# John Blischak
# 2013-06-03
# jdblischak at gmail dot com
############################################################
rm(list = ls())
library(ggplot2)
data(mpg)

# Motivation: 
# Graphics are important not only to concisely convey the 
# results of an analysis, but also as diagnostic tools to 
# inform the statistical analyses being performed. ggplot2 
# is an R package that allows you to create high-quality,
# informative graphics in a flexible, intuitive fashion.

# Resources:
# Documentation - http://docs.ggplot2.org/current/
# Academic paper - http://vita.had.co.nz/papers/layered-grammar.pdf

ggplot() + 
  layer(data = mpg,
                 mapping = aes(x = cty, y = hwy, 
                               shape = factor(drv)),
                 geom = 'point', 
                 stat = 'identity',
                 position = 'identity') +
  facet_grid(. ~ year) + theme_bw() +
  labs(title = 'Comparing fuel economy in the city versus the highway\nfor cars with different drive trains',
       x = 'Fuel economy in the city (mpg)',
       y = 'Fuel economy on the highway (mpg)')

ggplot() + 
  layer(data = mpg,
        mapping = aes(x = cty, fill = factor(drv)),
        geom = 'bar', 
        stat = 'bin',
        position = 'stack') +
  facet_grid(year ~ .) +
  labs(title = 'Distribution of fuel economy for city driving for cars with different drive trains',
       x = 'Fuel economy in the city (mpg)',
       y = 'Number of cars')

# Data and Mappings
# First, we need to supply ggplot2 with data for the graph.
# The data must be in a data.frame. Furthermore, we need to
# indicate which columns of the data frame should be mapped
# to specific aesthetics. Aesthetics are the "stuff" of the 
# graph, e.g. the x- and y- coordinates, color, shape, and
# size of a point, etc. As an example, we will be using 
# a dataset from the EPA on the fuel economy of popular
# car models (?mpg for more details).
base1 <- ggplot(data = mpg, 
                mapping = aes(x = cty, y = hwy, 
                              color = factor(drv)))
base1
str(base1)

base2 <- ggplot(data = mpg, 
                mapping = aes(x = cty, fill = factor(drv)))
str(base2)

# Geometric objects
# The geom specifies the type of plot that will be drawn,
# e.g. geom_point creates a scatterplot and geom_bar a
# barplot.
base1 + geom_point()
base1 + geom_line()
base1 + geom_point() + geom_line()
base1 + geom_path()

base2 + geom_histogram()
base2 + geom_density(alpha = .5)

# Statistical transformations 
# A stat performs a statistical operation on the data and
# then graphs the statistical summary as the output. This is
# especially useful for adding fitted lines.
base1 + stat_smooth()
base1 + stat_smooth(method = 'lm')

base2 + stat_bin()
base2 + stat_bin() + geom_bar()

# Question: 
# Why does stat_bin create a histogram exactly like
# geom_histogram or stat_bin + geom_bar?

# Positions
# Position arguments control the interaction of geometric
# objects on the graph, e.g. stacked vs. side-by-side
# barplots.
base1 + geom_point(position = 'jitter')

base2 + geom_histogram(position = 'stack')
base2 + geom_histogram(position = 'dodge')

# Faceting 
# Faceting splits the data into subsets and creates a
# separate graph for each subset of the data.
base1 + geom_point() + facet_grid(. ~ year)

base2 + geom_histogram() + facet_grid(drv ~ .)

# Challenge:
# Facet the histogram of city fuel economy with the rows
# as the number of cylinders and the columns as the year.
base2 + geom_histogram() + facet_grid(cyl ~ year)

# Scales
# Scale arguments control how the data is converted into
# aesthetics, e.g. change the limits of the x and y axes or
# change the spectrum of a color aesthetic.
base1 + geom_point() + scale_y_continuous(limits=c(20, 30))

base2 + geom_histogram(position = 'dodge') + 
  scale_x_continuous(limits=c(10, 20))

# Challenge:
# Using one of the scale functions, convert the axes of
# base1 (highway vs. city mpg) to log10 scale.

base1 + geom_point() + scale_x_log10() + scale_y_log10()

# Themes
# Theme arguments control the appearance of graph itself, 
# e.g. changing the font of the axis labels or changing the
# color of the graph's background.
base1 + geom_point() + 
  theme(panel.background = element_rect('yellow'))

base2 + geom_histogram() + theme(legend.position = 'left')

# Challenge:
# For the scatter plot, change the major and minor grid line
# to red.
base1 + geom_point() + 
  theme(panel.grid.major = element_line(color = 'red'),
        panel.grid.minor = element_line(color = 'red'))

# Convenience functions: 
# I have purposely been verbose in order to explain the
# grammar of graphics behind ggplot2. However, once you are
# comfortable with the basics you will want to start using
# the available convenience functions, e.g. qplot, 
# last_plot(), xlim, and ylim.

# Other useful packages:
# ggplot2 was started by Hadley Wickham. He has also started
# packages for high-level data manipulation called plyr and
# reshape2. These complement ggplot2 by facilitating the
# manipulation of your data into a data.frame.


# Final Exercise: 
# Plot city mpg (cty) versus engine
# displacement (displ). Map the number of cylinders (cyl) to
# the size of the points. Limit the y-axis to 15-25 mpg.
# Move the legend to the bottom of the plot. Use the labs
# function to add a title and better describe the x and y
# axis labels.

ggplot(data = mpg, 
       mapping = aes(x = displ, y = cty, size = cyl)) + 
  geom_point() +
  theme(legend.position = 'bottom') + ylim(15, 25) +
  labs(title = 'City fuel economy versus engine size',
       x = 'Engine size (liters)',
       y = 'City fuel economy (mpg)')

# Time permitting:
# If time permits, experiment with some of your own data to 
# see the power of ggplot2 (use read.table to pull data into
# R). Or you can use one of the datasets provided by the 
# ggplot2 package.



  