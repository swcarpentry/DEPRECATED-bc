#-----------------------------------------------------------
# Re-make lecture materials.
#
# We use Jekyll to compile HTML and Markdown files into their final
# form, and nbconvert to translate IPython Notebooks to theirs.  The
# problem is that Jekyll always erases and re-creates the output
# directory, so any compiled notebooks we put there are lost when we
# run it.  The solution is to cache compiled notebooks in a temporary
# directory, and copy those compiled pages into the output directory
# as needed.
#-----------------------------------------------------------

# Directories.
OUT = _site
TMP = tmp
LINK_OUT = /tmp/bc-links

# Source and destination Markdown/HTML pages.
MARKDOWN_SRC = \
	LICENSE.md \
	NEW_MATERIAL.md \
	bib.md \
	gloss.md \
	rules.md \
	$(wildcard bash/novice/*.md) \
	$(wildcard git/novice/*.md) \
	$(wildcard python/novice/*.md) \
	$(wildcard sql/novice/*.md)
MARKDOWN_DST = \
	$(patsubst %.md,$(OUT)/%.html,$(MARKDOWN_SRC))

# Source, cached, and destination Notebook files/HTML pages.
NOTEBOOK_SRC = \
	$(wildcard bash/novice/*.ipynb) \
	$(wildcard git/novice/*.ipynb) \
	$(wildcard python/novice/*.ipynb) \
	$(wildcard sql/novice/*.ipynb)
NOTEBOOK_TMP = \
	$(patsubst %.ipynb,$(TMP)/%.html,$(NOTEBOOK_SRC))
NOTEBOOK_DST = \
	$(patsubst %.ipynb,$(OUT)/%.html,$(NOTEBOOK_SRC))

# Mark cached versions of compiled notebooks as SECONDARY so that GNU
# Make won't delete them after rebuilding.
.SECONDARY : $(NOTEBOOK_TMP)

#-----------------------------------------------------------

# Default action: show available commands (marked with double '#').
all : commands

## commands : show all commands
commands :
	@grep -E '^##' Makefile | sed -e 's/## //g'

## check    : build site.
#  We know we're done when the compiled IPython Notebook files are
#  in the output directory.
check : $(NOTEBOOK_DST)

# Cannot create final versions of compiled notebook files until Jekyll
# has re-created the output directory.
$(NOTEBOOK_DST) : $(OUT)

# Copy cached versions of compiled notebook files into output directory.
$(OUT)/%.html : $(TMP)/%.html
	cp $< $@

# Build HTML versions of Markdown source files using Jekyll.  This always
# erases and re-creates the output directory.
$(OUT) : $(MARKDOWN_SRC)
	jekyll -t build -d $(OUT)

# Build HTML versions of IPython Notebooks.  This is slow, so we cache
# the results in a temporary directory.
$(TMP)/%.html : %.ipynb
	@mkdir -p $$(dirname $@)
	ipython nbconvert --output="$(subst .html,,$@)" "$<"

## links    : check links
# Depends on linklint, an HTML link-checking module from
# http://www.linklint.org/, which has been put in bin/linklint.
# Look in output directory's 'error.txt' file for results.
links :
	@bin/linklint -doc $(LINK_OUT) -textonly -root $(OUT) /@

## fixme    : find places where fixes are needed.
fixme :
	grep -n FIXME $$(find -f bash git python sql -type f -print | grep -v .ipynb_checkpoints)

## clean    : clean up
clean :
	rm -rf $(OUT) $(TMP) $$(find . -name '*~' -print) $$(find . -name '*.pyc' -print)

## show     : show variables
show :
	@echo "MARKDOWN_SRC" $(MARKDOWN_SRC)
	@echo "MARKDOWN_DST" $(MARKDOWN_DST)
	@echo "NOTEBOOK_SRC" $(NOTEBOOK_SRC)
	@echo "NOTEBOOK_DST" $(NOTEBOOK_DST)
