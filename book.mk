#======================================================================
# Create all-in-one version of notes as HTML or PDF.
#
# This is stored in a separate file so that rebuilding the site won't
# inadvertently trigger an attempt to create the book.  Run this using
# 'make -f book.mk html' or 'make -f book.mk pdf', or use 'make html'
# or 'make pdf' (which triggers a run from the main Makefile).
#======================================================================

include vars.mk

#----------------------------------------------------------------------
# Convert SVG diagrams into PNG images for use in LaTeX.
#----------------------------------------------------------------------

INKSCAPE = /Applications/Inkscape.app/Contents/Resources/bin/inkscape
SVG_TO_PNG = $(INKSCAPE) --export-png
DIAGRAM_SRC = $(wildcard novice/*/img/*.svg)
DIAGRAM_DST = $(patsubst %.svg,%.png,$(DIAGRAM_SRC))

diagrams : $(DIAGRAM_DST)

%.png : %.svg
	$(SVG_TO_PNG) $@ $<

#----------------------------------------------------------------------
# All-in-one book version.
#----------------------------------------------------------------------

BOOK_HTML = $(SITE)/book.html
BOOK_TEX = $(SITE)/book.tex
BOOK_PDF = $(SITE)/book.pdf

# pdf : build the all-in-one PDF version of the material.
pdf : $(BOOK_PDF)

# html : build the all-in-one HTML version of the material.
html : $(BOOK_HTML)

# Build the temporary input for the book by concatenating relevant
# sections of Markdown files and then forcing a rebuild of the site.
$(BOOK_MD) :
	python bin/make-book.py $(MOST_SRC) > $@
	touch *.md
	make site

$(BOOK_HTML) : $(BOOK_MD) $(DIAGRAM_DST)
	make site
	sed -i -e 's@\.\./\.\./gloss.html#@#g:@g' $@
	sed -i -e 's@\.svg@\.png@g' $@
	sed -i -e 's@<h4 id="challenges.*">@<h4>@g' $@
	sed -i -e 's@<h4 id="key-points.*">@<h4>@g' $@
	sed -i -e 's@<h4 id="objectives.*">@<h4>@g' $@
	sed -i -e 's@<h4 id="next-steps.*">@<h4>@g' $@
	sed -i -e 's@<pre class="in"><code>@<pre class="in">@g' $@
	sed -i -e 's@<pre class="out"><code>@<pre class="out">@g' $@
	sed -i -e 's@<pre class="err"><code>@<pre class="err">@g' $@

$(BOOK_TEX) : $(BOOK_HTML)
	pandoc -f html -t latex \
	  --standalone --table-of-contents --no-highlight --ascii \
	  --template _templates/tex.tpl \
	  -o $@ $<
	sed -i -e 's@\\paragraph@\\mbox\{\}\\paragraph@g' $@
	sed -i -e 's@π@\$$\\pi\$$@' $@

$(BOOK_PDF) : $(BOOK_TEX)
	cd $(SITE) && pdflatex book
