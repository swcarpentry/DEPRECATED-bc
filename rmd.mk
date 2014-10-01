#======================================================================
# Create Markdown versions of R Markdown files.
#
# This is stored in a separate file so that rebuilding the site won't
# inadvertently trigger an attempt to turn R Markdown files into
# Markdown for people who don't want to install R and knitr. Run this
# using 'make -f rmd.mk' (which runs this directly) or 'make rmd'
# (which triggers a run from the main Makefile).
# ======================================================================

#----------------------------------------------------------------------
# Specify the default target before any other targets are defined so
# that we're sure which one Make will choose.
#----------------------------------------------------------------------

all : rmd

#----------------------------------------------------------------------
# Rules.
#----------------------------------------------------------------------

# Chunk options for knitr
CHUNK_OPTS = novice/r/chunk_options.R

# R Markdown files.  Add patterns here to convert files stored in
# other locations.
RMD_SRC = \
	$(wildcard novice/r/??-*.Rmd)

# Files converted to Markdown.
RMD_TX = $(patsubst %.Rmd,%.md,$(RMD_SRC))

# Convert a .Rmd to .md.
%.md: %.Rmd $(CHUNK_OPTS)
	cd $$(dirname $<) && \
        Rscript -e "knitr::knit('$$(basename $<)', \
                                output = '$$(basename $@)')"

# Target.
rmd : $(RMD_TX)
