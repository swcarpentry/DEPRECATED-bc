#======================================================================
# Makefile for Software Carpentry EPUB version of lessons.
#
# We use the HTML version generate by Jekyll to get better results
# with pandoc.
#======================================================================

EPUB = _site/book.epub
MOBI = _site/book.mobi
AZW3 = _site/book.azw3

CSS = _epub/stylesheet.css

EPUB_SOURCE = _site/book.html

## commands : show all commands.
commands :
	@grep -E '^##' epub.mk | sed -e 's/##//g'

## epub     : build epub version of lessons
epub: $(EPUB)

$(EPUB) : $(EPUB_SOURCE) $(CSS)
	pandoc -f html-native_divs -t epub -o $@ \
	    --standalone \
	    --epub-stylesheet=$(CSS) \
	    --epub-metadata=_epub/metadata.xml \
	    --epub-chapter-level=2 \
	    $(EPUB_SOURCE)
	mkdir -p _epub_tmp
	unzip -uo $@ -d _epub_tmp
	python bin/swc_fix_epub.py _epub_tmp
	cd _epub_tmp && zip ../$@ *.xhtml
	rm -rf _epub_tmp

## mobi     : build mobi version of lessons
mobi: $(MOBI)

$(MOBI): $(EPUB)
	ebook-convert $< $@ \
	    --mobi-ignore-margins \
	    --no-inline-toc

## azw      : build azw version of lessons
azw: $(AZW3)

$(AZW3): $(EPUB)
	ebook-convert $< $@ \
	    --no-inline-toc

## all      : build all version of lessons (epub, mobi, azw)
all: $(EPUB) $(MOBI) $(AZW3)
