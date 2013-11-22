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

## check    : build locally into $(OUT) directory for checking
check : $(NOTEBOOK_DST)
	jekyll -t build -d $(OUT)

## clean    : clean up
clean :
	rm -rf $(OUT) $$(find . -name '*~' -print) $$(find . -name '*.pyc' -print)

## show     : show variables
show :
	@echo "NOTEBOOK_SRC" $(NOTEBOOK_SRC)
	@echo "NOTEBOOK_DST" $(NOTEBOOK_DST)

#-----------------------------------------------------------

# rule to make HTML pages from notebook files.
$(OUT)/%.html : %.ipynb
	@mkdir -p $$(dirname $@)
	ipython nbconvert --stdout $< > $@
