#-----------------------------------------------------------
# Re-make lecture materials.
#-----------------------------------------------------------

OUT = _site

MARKDOWN_SRC = \
	LICENSE.md \
	gloss.md \
	rules.md \
	$(wildcard bash/novice/*.md) \
	$(wildcard git/novice/*.md) \
	$(wildcard python/novice/*.md) \
	$(wildcard sql/novice/*.md)
MARKDOWN_DST = \
	$(patsubst %.md,$(OUT)/%.html,$(MARKDOWN_SRC))

NOTEBOOK_SRC = \
	$(wildcard bash/novice/*.ipynb) \
	$(wildcard git/novice/*.ipynb) \
	$(wildcard python/novice/*.ipynb) \
	$(wildcard sql/novice/*.ipynb)
NOTEBOOK_DST = \
	$(patsubst %.ipynb,$(OUT)/%.html,$(NOTEBOOK_SRC))

#-----------------------------------------------------------

all : commands

## commands : show all commands
commands :
	@grep -E '^##' Makefile | sed -e 's/## //g'

## check    : build site.
check : $(NOTEBOOK_DST)
$(NOTEBOOK_DST) : $(OUT)

# Build HTML versions of Markdown source files.
$(OUT) : $(MARKDOWN_SRC)
	jekyll -t build -d $(OUT)

# Build HTML versions of IPython Notebooks (slow).
$(OUT)/%.html : %.ipynb
	@mkdir -p $$(dirname $@)
	ipython nbconvert --output="$(subst .html,,$@)" "$<"

## links    : check links
# Depends on linklint, an HTML link-checking module from
# http://www.linklint.org/, which has been put in bin/linklint.
links :
	@bin/linklint -doc /tmp/bc-links -textonly -root _site /@

## clean    : clean up
clean :
	rm -rf $(OUT) $$(find . -name '*~' -print) $$(find . -name '*.pyc' -print)

## show     : show variables
show :
	@echo "MARKDOWN_SRC" $(MARKDOWN_SRC)
	@echo "MARKDOWN_DST" $(MARKDOWN_DST)
	@echo "NOTEBOOK_SRC" $(NOTEBOOK_SRC)
	@echo "NOTEBOOK_DST" $(NOTEBOOK_DST)
