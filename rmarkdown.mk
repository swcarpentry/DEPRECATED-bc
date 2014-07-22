#======================================================================
# Create HTML versions of R Markdown files.
#
# This is stored in a separate file so that rebuilding the site won't
# inadvertently trigger an attempt to turn R Markdown files into HTML
# for people who don't want to install R, knitr, and rmarkdown.  Run
# this using 'make -f rmarkdown.mk' (which runs this directly) or
# 'make rmarkdown' (which triggers a run from the main Makefile).
# ======================================================================

#----------------------------------------------------------------------
# Specify the default target before any other targets are defined so
# that we're sure which one Make will choose.
#----------------------------------------------------------------------

all : rmarkdown

#----------------------------------------------------------------------
# Rules.
#----------------------------------------------------------------------

# R Markdown files.  Add patterns here to convert files stored in
# other locations.
RMARKDOWN_SRC = \
	$(wildcard novice/r/??-*.Rmd) \

# Files converted to HTML.
RMARKDOWN_TX = $(patsubst %.Rmd,%.html,$(RMARKDOWN_SRC))

# Convert a .Rmd to .html.
%.html : %.Rmd
	Rscript -e "library(rmarkdown); \
                    render('$<', \
                    output_format = html_document(theme = NULL))"
	# Hack to add back in the YAML header for processing with
	# jekyll.
	head -n 4 $< > /tmp/yaml
	cat /tmp/yaml $@ > /tmp/lesson_w_yaml
	mv /tmp/lesson_w_yaml $@
	rm /tmp/yaml

# Target.
rmarkdown : $(RMARKDOWN_TX)
