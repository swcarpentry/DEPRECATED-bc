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

# Auxiliary directories.
OUT = _site
TMP = tmp
LINK_OUT = /tmp/bc-links

# Source directories (add more here as lessons materialize).
NOVICE_DIR = bash git python sql
INTERMEDIATE_DIR = bash hash make

# Source and destination Markdown/HTML pages.
MARKDOWN_SRC = \
	LICENSE.md \
	NEW_MATERIAL.md \
	bib.md \
	gloss.md \
	rules.md \
	$(wildcard $(patsubst %,%/novice/*.md,$(NOVICE_DIR))) \
	$(wildcard $(patsubst %,%/intermediate/*.md,$(INTERMEDIATE_DIR)))
MARKDOWN_DST = \
	$(patsubst %.md,$(OUT)/%.html,$(MARKDOWN_SRC))

# Source, cached, and destination Notebook files/HTML pages.
NOTEBOOK_SRC = \
	$(wildcard $(patsubst %,%/novice/*.ipynb,$(NOVICE_DIR))) \
	$(wildcard $(patsubst %,%/intermediate/*.ipynb,$(INTERMEDIATE_DIR)))
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

# Build HTML versions of Markdown source files using Jekyll.  This
# always erases and re-creates the output directory, and replaces the
# template file for bootcamp websites with the NEW_MATERIAL.html index
# to the lessons.
$(OUT) : $(MARKDOWN_SRC)
	jekyll -t build -d $(OUT)
	mv $(OUT)/NEW_MATERIAL.html $(OUT)/index.html
	rm -rf $$(find _site -name '*_files' -type d -print)

# Build HTML versions of IPython Notebooks.  This is slow, so we cache
# the results in a temporary directory.
$(TMP)/%.html : %.ipynb
	@mkdir -p $$(dirname $@)
	ipython nbconvert --output="$(subst .html,,$@)" "$<"

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
	rm -rf $(OUT) $(TMP) $$(find . -name '*~' -print) $$(find . -name '*.pyc' -print) $$(find . -name '*_files' -type d -print)

## show     : show variables
show :
	@echo "MARKDOWN_SRC" $(MARKDOWN_SRC)
	@echo "MARKDOWN_DST" $(MARKDOWN_DST)
	@echo "NOTEBOOK_SRC" $(NOTEBOOK_SRC)
	@echo "NOTEBOOK_DST" $(NOTEBOOK_DST)
