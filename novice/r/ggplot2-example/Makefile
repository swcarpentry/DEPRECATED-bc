OUTS := $(patsubst %.Rmd,%.md,$(wildcard *.Rmd)) $(patsubst %.Rmd,%.R,$(wildcard *.Rmd))

all: $(OUTS)

clean:
	rm -rf $(OUTS) figure/

%.md: %.Rmd
	Rscript -e 'library(knitr); opts_chunk$$set(results="hide", fig.show="hide"); knit("$*.Rmd")'
	gsed -i 's/```r/```coffee/g' $*.md

%.R: %.Rmd
	Rscript -e "library(knitr); purl('$*.Rmd')"
