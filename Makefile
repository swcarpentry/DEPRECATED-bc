#================================================================================
# Re-make lecture materials.
#================================================================================

# Directories.
WEBSITE = site
INSTALL = $(HOME)/sites/software-carpentry.org/v5
LINKS = /tmp/bc-links
COLLECTED = collected
PATCHED = patched

# Templates.
IPYNB_TPL = _templates/ipynb.tpl
BOOK_TPL = _templates/book.tpl
LESSON_TPL = _templates/lesson.tpl

# Files.
BOOK = $(WEBSITE)/book.html
INDEX = $(WEBSITE)/index.html

#--------------------------------------------------------------------------------
# Specify the default command.
#--------------------------------------------------------------------------------

# Default action: show available commands (marked with double '#').
all : commands

#--------------------------------------------------------------------------------
# Collect Markdown versions of IPython Notebooks.
#--------------------------------------------------------------------------------

# IPython Notebook source (split by directory for later interpolation).
IPYNB_SRC_PYTHON = $(sort $(wildcard novice/python/??-*.ipynb))
IPYNB_SRC_SQL = $(sort $(wildcard novice/sql/??-*.ipynb))

# Notebooks converted to Markdown.
IPYNB_COLLECTED_PYTHON = $(patsubst %.ipynb,$(COLLECTED)/%.md,$(IPYNB_SRC_PYTHON))
IPYNB_COLLECTED_SQL = $(patsubst %.ipynb,$(COLLECTED)/%.md,$(IPYNB_SRC_SQL))
IPYNB_COLLECTED = $(IPYNB_COLLECTED_PYTHON) $(IPYNB_COLLECTED_SQL)

# Convert from .ipynb to .md.
$(COLLECTED)/%.md : %.ipynb $(IPYNB_TPL)
	ipython nbconvert --template=$(IPYNB_TPL) --to=markdown --output="$(subst .md,,$@)" "$<"

#--------------------------------------------------------------------------------
# Markdown source (not including files generated from .ipynb).
#--------------------------------------------------------------------------------

# Raw Markdown.
MD_SRC = \
	contents.md \
	intro.md \
	team.md \
	novice/shell/index.md $(sort $(wildcard novice/shell/??-*.md)) \
	novice/git/index.md $(sort $(wildcard novice/git/??-*.md)) \
	novice/python/index.md \
	novice/sql/index.md \
	novice/extras/index.md $(sort $(wildcard novice/extras/??-*.md)) \
	novice/teaching/index.md  $(sort $(wildcard novice/teaching/??-*.md)) \
	novice/ref/index.md  $(sort $(wildcard novice/ref/??-*.md)) \
	bib.md \
	gloss.md \
	rules.md \
	LICENSE.md

# Collected Markdown files (just a copy).
MD_COLLECTED = $(patsubst %,$(COLLECTED)/%,$(MD_SRC))

# Copy over Markdown files.
$(COLLECTED)/%.md : %.md
	@mkdir -p $$(dirname $@)
	cp $< $@

#--------------------------------------------------------------------------------
# Patch image paths in files (directory by directory).
#--------------------------------------------------------------------------------

# Use string substitution to interpolate Markdown files generated from notebooks
# in the right order.
ALL_COLLECTED = $(subst $(COLLECTED)/novice/python/index.md,$(COLLECTED)/novice/python/index.md $(IPYNB_COLLECTED_PYTHON),\
                $(subst $(COLLECTED)/novice/sql/index.md,$(COLLECTED)/novice/sql/index.md $(IPYNB_COLLECTED_SQL),\
                $(MD_COLLECTED)))

# All patched files.
ALL_PATCHED = $(patsubst $(COLLECTED)/%,$(PATCHED)/%,$(ALL_COLLECTED))

$(PATCHED)/novice/shell/%.md : $(COLLECTED)/novice/shell/%.md
	@mkdir -p $$(dirname $@)
	sed -e 's!<img src="img!<img src="novice/shell/img!g' $< > $@
$(PATCHED)/novice/git/%.md : $(COLLECTED)/novice/git/%.md
	@mkdir -p $$(dirname $@)
	sed -e 's!<img src="img!<img src="novice/git/img!g' $< > $@
$(PATCHED)/novice/python/%.md : $(COLLECTED)/novice/python/%.md
	@mkdir -p $$(dirname $@)
	sed -e 's!<img src="img!<img src="novice/python/img!g' $< > $@
$(PATCHED)/novice/sql/%.md : $(COLLECTED)/novice/sql/%.md
	@mkdir -p $$(dirname $@)
	sed -e 's!<img src="img!<img src="novice/sql/img!g' $< > $@
$(PATCHED)/novice/extras/%.md : $(COLLECTED)/novice/extras/%.md
	@mkdir -p $$(dirname $@)
	sed -e 's!<img src="img!<img src="novice/extras/img!g' $< > $@

# Everything else can just be copied over.
$(PATCHED)/%.md : %.md
	@mkdir -p $$(dirname $@)
	cp $< $@

#--------------------------------------------------------------------------------
# Build the web site from the patched source files.
#--------------------------------------------------------------------------------

# All patched files in the website.
ALL_WEBSITE = $(patsubst $(PATCHED)/%.md,$(WEBSITE)/%.html,$(ALL_PATCHED))

$(BOOK) : $(INDEX) $(BOOK_TPL) bin/make-book.py
	@mkdir -p $$(dirname $@)
	python bin/make-book.py $(BOOK_TMP) \
	| pandoc --email-obfuscation=none --template=$(BOOK_TPL) -t html -o - \
	| sed -e 's!../../gloss.html#!#g:!g' \
	| sed -e 's!../gloss.html#!#g:!g' \
	> $@

$(INDEX) : $(ALL_WEBSITE)
	cp $(WEBSITE)/contents.html $@

$(WEBSITE)/%.html : $(PATCHED)/%.md
	@mkdir -p $$(dirname $@)
	pandoc --email-obfuscation=none --template=$(LESSON_TPL) -t html -o $@ $<

#--------------------------------------------------------------------------------
# Extra files used in web site (CSS, images, etc.).
#--------------------------------------------------------------------------------

# CSS files.
CSS_SRC = $(wildcard css/*.css) $(wildcard css/*/*.css)
CSS_OUT = $(patsubst %,$(WEBSITE)/%,$(CSS_SRC))

# Image files.
IMG_SRC = \
	$(wildcard novice/*/img/*.png) $(wildcard novice/*/img/*.svg) \
	$(wildcard img/*.png) $(wildcard img/slides/*.png)
IMG_OUT = $(patsubst %,$(WEBSITE)/%,$(IMG_SRC))

# Rules (repeated because there's no easy way to abstract across suffixes).
$(WEBSITE)/%.css : %.css
	@mkdir -p $$(dirname $@)
	cp $< $@
$(WEBSITE)/%.png : %.png
	@mkdir -p $$(dirname $@)
	cp $< $@
$(WEBSITE)/%.svg : %.svg
	@mkdir -p $$(dirname $@)
	cp $< $@

#--------------------------------------------------------------------------------
# All invokable targets.
#--------------------------------------------------------------------------------

## site     : build the whole site.
site : $(BOOK) $(CSS_OUT) $(IMG_OUT)

## commands : show all commands.
commands :
	@grep -E '^##' Makefile | sed -e 's/## //g'

## contribs : list contributors (uses .mailmap file).
contribs :
	git log --pretty=format:%aN | sort | uniq

## fixme    : find places where fixes are needed.
fixme :
	@grep -i -n FIXME $$(find -f shell git python sql -type f -print | grep -v .ipynb_checkpoints)

## gloss    : check glossary.
gloss :
	@bin/gloss.py ./gloss.md $(MARKDOWN_DST) $(NOTEBOOK_DST)

## valid    : check validity of HTML book.
# Depends on xmllint to check validity of generated pages.
# Also depends on linklint, an HTML link-checking module from
# http://www.linklint.org/, which has been put in bin/linklint.
# Look in output directory's 'error.txt' file for results.
valid : tmp-book.html
	xmllint --noout tmp-book.html 2>&1 | python bin/unwarn.py
	@bin/linklint -doc $(LINKS) -textonly -root $(WEBSITE) /@

## sterile  : _really_ clean up.
sterile : clean
	rm -rf $(COLLECTED) $(PATCHED)

## clean    : clean up all generated files.
clean : tidy
	rm -rf $(WEBSITE)

## tidy     : clean up intermediate files only.
tidy :
	rm -rf \
	$$(find . -name '*~' -print) \
	$$(find . -name '*.pyc' -print)

## show     : show variables for debugging
show :
	@echo OUT $(WEBSITE)
	@echo INSTALL $(INSTALL)
	@echo LINKS $(LINKS)
	@echo COLLECTED $(COLLECTED)
	@echo PATCHED $(PATCHED)
	@echo IPYNB_TPL $(IPYNB_TPL)
	@echo BOOK_TPL $(BOOK_TPL)
	@echo BOOK $(BOOK)
	@echo IPYNB_SRC_PYTHON $(IPYNB_SRC_PYTHON)
	@echo IPYNB_SRC_SQL $(IPYNB_SRC_SQL)
	@echo IPYNB_COLLECTED_PYTHON $(IPYNB_COLLECTED_PYTHON)
	@echo IPYNB_COLLECTED_SQL $(IPYNB_COLLECTED_SQL)
	@echo IPYNB_COLLECTED $(IPYNB_COLLECTED)
	@echo MD_SRC $(MD_SRC)
	@echo MD_COLLECTED $(MD_COLLECTED)
	@echo ALL_COLLECTED $(ALL_COLLECTED)
	@echo ALL_PATCHED $(ALL_PATCHED)

#--------------------------------------------------------------------------------
# Support.
#--------------------------------------------------------------------------------

# Stop GNU Make from deleting these temporary files.
.SECONDARY : $(ALL_COLLECTED) $(ALL_PATCHED)
