#======================================================================
# Create Markdown versions of IPython Notebooks.
#
# This is stored in a separate file so that rebuilding the site won't
# inadvertently trigger an attempt to turn IPython Notebok files into
# Markdown for people who don't want to install IPython.  Run this
# using 'make -f ipynb.mk' (which runs this directly) or 'make ipynb'
# (which triggers a run from the main Makefile).
#======================================================================

#----------------------------------------------------------------------
# Specify the default target before any other targets are defined so
# that we're sure which one Make will choose.
#----------------------------------------------------------------------

all : ipynb

#----------------------------------------------------------------------
# Rules.
#----------------------------------------------------------------------

# Templates.
IPYNB_TPL = _templates/ipynb.tpl

# IPython Notebooks.  Add patterns here to convert notebooks stored in
# other locations.
IPYNB_SRC = \
	$(wildcard novice/python/??-*.ipynb) \
	$(wildcard novice/sql/??-*.ipynb) \
	$(wildcard intermediate/doit/??-*.ipynb) \
	$(wildcard intermediate/python/??-*.ipynb) \
	$(wildcard intermediate/webdata/??-*.ipynb)

# Notebooks converted to Markdown.
IPYNB_TX = $(patsubst %.ipynb,%.md,$(IPYNB_SRC))

# Convert a .ipynb to .md.
%.md : %.ipynb $(IPYNB_TPL)
	ipython nbconvert --template=$(IPYNB_TPL) --to=markdown --output="$(subst .md,,$@)" "$<"

# Target.
ipynb : $(IPYNB_TX)
