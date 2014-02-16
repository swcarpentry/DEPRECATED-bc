#-----------------------------------------------------------
# Re-make lecture materials.
#-----------------------------------------------------------

# Directories.
OUT = _site
LINK_OUT = /tmp/bc-links
BOOK = _book
INSTALL_DIR = $(HOME)/sites/software-carpentry.org/v5

# Source Markdown pages.
MARKDOWN_SRC = \
	LICENSE.md \
	contents.md \
	bib.md \
	gloss.md \
	rules.md \
	setup.md \
	team.md \
	$(wildcard shell/novice/*.md) \
	$(wildcard git/novice/*.md) \
	$(wildcard python/novice/*.md) \
	$(wildcard sql/novice/*.md) \
	$(wildcard teaching/novice/*.md)

NOTEBOOK_SRC = \
	$(wildcard python/novice/??-*.ipynb) \
	$(wildcard sql/novice/??-*.ipynb)

# Slides.
SLIDES_SRC = \
	$(wildcard slides/*.html)
SLIDES_DST = \
	$(patsubst %.html,$(OUT)/%.html,$(SLIDES_SRC))

NOTEBOOK_MD = \
	$(patsubst %.ipynb,%.md,$(NOTEBOOK_SRC))

HTML_DST = \
	$(patsubst %.md,$(OUT)/%.html,$(MARKDOWN_SRC)) \
	$(patsubst %.md,$(OUT)/%.html,$(NOTEBOOK_MD))

BOOK_SRC = \
	intro.md \
	shell/novice/index.md $(wildcard shell/novice/*-*.md) \
	git/novice/index.md $(wildcard git/novice/*-*.md) \
	python/novice/index.md $(patsubst %.ipynb,%.md,$(wildcard python/novice/??-*.ipynb)) \
	sql/novice/index.md $(patsubst %.ipynb,%.md,$(wildcard sql/novice/??-*.ipynb)) \
	extras/novice/index.md $(wildcard extras/novice/*-*.md) \
	teaching/novice/index.md  $(wildcard teaching/novice/*-*.md) \
	ref/novice/index.md  $(wildcard ref/novice/*-*.md) \
	bib.md \
	tmp/gloss.md \
	rules.md \
	team.md \
	LICENSE.md

BOOK_TMP = \
	$(patsubst %,tmp/%,$(BOOK_SRC))

BOOK_DST = $(OUT)/book.html

# Mark cached versions of compiled notebooks as SECONDARY so that GNU
# Make won't delete them after rebuilding.
.SECONDARY : $(NOTEBOOK_MD)

#-----------------------------------------------------------

# Default action: show available commands (marked with double '#').
all : commands

## quick    : build just the bootcamp home page.
quick :
	jekyll -t build -d $(OUT)

## install  : install.
install : $(OUT)/index.html
	rm -rf $(INSTALL_DIR)
	cp -r _site $(INSTALL_DIR)
	mv $(INSTALL_DIR)/contents.html $(INSTALL_DIR)/index.html

## site     : build site.
site : $(BOOK_DST)

$(BOOK_DST) : $(OUT)/index.html $(BOOK_TMP) _templates/book.tpl tmp/gloss.md bin/make-book.py
	python bin/make-book.py $(BOOK_TMP) \
	| pandoc --email-obfuscation=none --template=_templates/book.tpl -t html -o - \
	| sed -e 's!../../gloss.html#!#g:!g' \
	| sed -e 's!../gloss.html#!#g:!g' \
	> $@

# Build HTML versions of Markdown source files using Jekyll.
$(OUT)/index.html : $(MARKDOWN_SRC) $(NOTEBOOK_MD)
	jekyll -t build -d $(OUT)
	sed -i -e 's!img src="python/novice/!img src="!g' $(OUT)/python/novice/??-*.html

# Build Markdown versions of IPython Notebooks.
%.md : %.ipynb _templates/ipynb.tpl
	ipython nbconvert --template=_templates/ipynb.tpl --to=markdown --output="$(subst .md,,$@)" "$<"

# Patch targets and links in the glossary for inclusion in the book.
tmp/gloss.md : gloss.md
	@mkdir -p $$(dirname $@)
	sed -e 's!](#!](#g:!g' -e 's!<a name="!<a name="#g:!g' $< > $@

# Patch image paths in the sections.
tmp/shell/novice/%.md : shell/novice/%.md
	@mkdir -p $$(dirname $@)
	sed -e 's!<img src="img!<img src="shell/novice/img!g' $< > $@

tmp/git/novice/%.md : git/novice/%.md
	@mkdir -p $$(dirname $@)
	sed -e 's!<img src="img!<img src="git/novice/img!g' $< > $@

tmp/python/novice/%.md : python/novice/%.md
	@mkdir -p $$(dirname $@)
	sed -e 's!<img src="img!<img src="python/novice/img!g' $< > $@

tmp/sql/novice/%.md : sql/novice/%.md
	@mkdir -p $$(dirname $@)
	sed -e 's!<img src="img!<img src="sql/novice/img!g' $< > $@

tmp/extras/novice/%.md : extras/novice/%.md
	@mkdir -p $$(dirname $@)
	sed -e 's!<img src="img!<img src="extras/novice/img!g' $< > $@

# All other Markdown files used in the book.
tmp/%.md : %.md
	@mkdir -p $$(dirname $@)
	cp $< $@

#-----------------------------------------------------------

## commands : show all commands
commands :
	@grep -E '^##' Makefile | sed -e 's/## //g'

## contribs : list contributors (uses .mailmap file)
contribs :
	git log --pretty=format:%aN | sort | uniq

## fixme    : find places where fixes are needed.
fixme :
	@grep -i -n FIXME $$(find -f shell git python sql -type f -print | grep -v .ipynb_checkpoints)

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

## valid    : check validity of HTML book.
valid : tmp-book.html
	xmllint --noout tmp-book.html 2>&1 | python bin/unwarn.py

## clean    : clean up
clean : tidy
	@rm -rf $(OUT) $(NOTEBOOK_MD)

## tidy     : clean up intermediate files only
tidy :
	@rm -rf \
	image-page.html \
	tmp \
	$$(find . -name '*~' -print) \
	$$(find . -name '*.pyc' -print) \
	$$(find . -name '??-*_files' -type d -print)

## show     : show variables
show :
	@echo "OUT" $(OUT)
	@echo "TMP" $(TMP)
	@echo "LINK_OUT" $(LINK_OUT)
	@echo "MARKDOWN_SRC" $(MARKDOWN_SRC)
	@echo "NOTEBOOK_SRC" $(NOTEBOOK_SRC)
	@echo "NOTEBOOK_MD" $(NOTEBOOK_MD)
	@echo "HTML_DST" $(HTML_DST)
