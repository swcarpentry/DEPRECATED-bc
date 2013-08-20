## use File --> Compile Notebook on me
## or click notebook button in top right of editor pane

## or get same result from the command line:
## Rscript -e "knitr::stitch_rhtml('code/02_slopeComparisonAsiaVsAmericas.R', output = 'prose/02_slopeComparisonAsiaVsAmericas.html')"

## either will require knitr package to work
## to install: install.packages("knitr")

library(lattice)

## the working directory will differ for live execution vs. command line execution vs.
## using Compile Notebook
## you will need to toggle the comment status of these lines to get
## what you want

## for Compile Notebook:
## str(gCoef <- readRDS("../results/gCoef.rds"))

## for live execution or command line execution, when working directory set
## correctly to top-level directory for this project
str(gCoef <- readRDS("results/gCoef.rds"))

hDat <-
  droplevels(subset(gCoef,
                    continent %in% c("Asia", "Americas")))
str(hDat)

## Figure appears fine via RStudio's Compile Notebook feature
dotplot(slope ~ continent, hDat)
## JENNY TO DO:
## This figure is NOT appearing in the .html file when we knitr::stitch() from
## the command line via "Rscript -e" or from a Makefile
## Figure that out and fix it!
## Suspect is related to figure/02-slopeComparisonAsiaVsAmericas-Rhtmlauto-report.png
## which IS being created. Suspect I need to tell knitr::stitch() to make a self-contained html file or
## make sure this png goes to the right place.
## this may eventually be helpful:
## opts_chunk$set(fig.path="prose/figs/")

t.test(slope ~ continent, hDat)
