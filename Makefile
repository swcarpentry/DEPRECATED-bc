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

# Source Markdown pages.
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

# Source, cached, and destination Notebook files/HTML pages.
NOTEBOOK_SRC = \
	$(wildcard bash/novice/*.ipynb) \
	$(wildcard git/novice/*.ipynb) \
	$(wildcard python/novice/*.ipynb) \
	$(wildcard sql/novice/*.ipynb)

NOTEBOOK_MD = \
	$(patsubst %.ipynb,%.md,$(NOTEBOOK_SRC))

HTML_DST = \
	$(patsubst %.md,$(OUT)/%.html,$(MARKDOWN_SRC)) \
	$(patsubst %.md,$(OUT)/%.html,$(NOTEBOOK_MD))

#-----------------------------------------------------------

# Default action: show available commands (marked with double '#').
all : commands

## commands : show all commands
commands :
	@grep -E '^##' Makefile | sed -e 's/## //g'

## check    : build site.
#  We know we're done when the compiled IPython Notebook files are
#  in the output directory.
check : $(OUT)/index.html

# Build HTML versions of Markdown source files using Jekyll.
$(OUT)/index.html : $(MARKDOWN_SRC) $(NOTEBOOK_MD)
	jekyll -t build -d $(OUT)
	mv $(OUT)/NEW_MATERIAL.html $(OUT)/index.html

# Build Markdown versions of IPython Notebooks.
%.md : %.ipynb
	ipython nbconvert --template=./swc.tpl --to=markdown --output="$(subst .md,,$@)" "$<"

## fixme    : find places where fixes are needed.
fixme :
	@grep -n FIXME $$(find -f bash git python sql -type f -print | grep -v .ipynb_checkpoints)

## gloss    : check glossary
gloss :
	@bin/gloss.py ./gloss.md $(MARKDOWN_DST) $(NOTEBOOK_DST)

## images   : create a temporary page to display images
images :
	@bin/make-image-page.py $(MARKDOWN_SRC) $(NOTEBOOK_SRC) > image-page.html
	@echo "Open ./image-page.html to view images"

## links    : check links
# Depends on linklint, an HTML link-checking module from
# http://www.linklint.org/, which has been put in bin/linklint.
# Look in output directory's 'error.txt' file for results.
links :
	@bin/linklint -doc $(LINK_OUT) -textonly -root $(OUT) /@

## clean    : clean up
clean :
	rm -rf $(OUT) $(TMP) $$(find . -name '*~' -print) $$(find . -name '*.pyc' -print)
	rm -f $(HTML_DST)

## show     : show variables
show :
	@echo "OUT" $(OUT)
	@echo "TMP" $(TMP)
	@echo "LINK_OUT" $(LINK_OUT)
	@echo "MARKDOWN_SRC" $(MARKDOWN_SRC)
	@echo "NOTEBOOK_SRC" $(NOTEBOOK_SRC)
	@echo "NOTEBOOK_MD" $(NOTEBOOK_MD)
	@echo "HTML_DST" $(HTML_DST)
