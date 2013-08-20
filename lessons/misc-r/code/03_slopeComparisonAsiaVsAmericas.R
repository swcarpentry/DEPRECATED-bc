## basically a copy of 02_slopeComparisonAsiaVsAmericas.R
## that does not anticipate/require the "Compile Notebook" treatment
library(lattice)

str(gCoef <- readRDS("results/gCoef.rds"))
hDat <-
  droplevels(subset(gCoef,
                    continent %in% c("Asia", "Americas")))
str(hDat)

pdf("figs/slopes_AsiaVsAmericas.pdf")
dotplot(slope ~ continent, hDat)
dev.off()

sink("results/slopes_AsiaVsAmericas.txt")
t.test(slope ~ continent, hDat)
sink()