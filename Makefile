#-----------------------------------------------------------
# Re-make lecture materials.
#-----------------------------------------------------------

OUT = _site
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

## check    : build notebooks *after* Jekyll runs, because it wipes the output directory.
check : $(NOTEBOOK_DST)

$(NOTEBOOK_DST) : $(OUT)/README.md

$(OUT)/README.md :
	jekyll -t build -d $(OUT)

$(OUT)/%.html : %.ipynb
	@mkdir -p $$(dirname $@)
	ipython nbconvert --output="$(subst .html,,$@)" "$<"

## links    : check links
#  Depends on linklint, an HTML link-checking module from http://www.linklint.org/,
#  which has been put in bin/linklint.
links :
	@bin/linklint -doc /tmp/bc-links -textonly -root _site /@

## clean    : clean up
clean :
	rm -rf $(OUT) $$(find . -name '*~' -print) $$(find . -name '*.pyc' -print)

## show     : show variables
show :
	@echo "NOTEBOOK_SRC" $(NOTEBOOK_SRC)
	@echo "NOTEBOOK_DST" $(NOTEBOOK_DST)
