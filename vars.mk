#======================================================================
# vars.mk : variables used in multiple Makefiles.
#======================================================================

# Output directory for local build.
SITE = _site

# All-in-one Markdown version of material.
BOOK_MD = ./book.md

# Source Markdown files.  These are listed in the order in which they
# appear in the final book-format version of the notes.
MOST_SRC = \
	 intro.md \
	 team.md \
	 novice/shell/index.md $(sort $(wildcard novice/shell/??-*.md)) \
	 novice/git/index.md $(sort $(wildcard novice/git/??-*.md)) \
	 novice/python/index.md $(sort $(wildcard novice/python/??-*.md)) \
	 novice/sql/index.md $(sort $(wildcard novice/sql/??-*.md)) \
	 novice/extras/index.md $(sort $(wildcard novice/extras/??-*.md)) \
	 novice/teaching/index.md  $(sort $(wildcard novice/teaching/??-*.md)) \
	 novice/ref/index.md  $(sort $(wildcard novice/ref/??-*.md)) \
	 bib.md \
	 gloss.md \
	 rules.md \
	 LICENSE.md

# All source pages (including things not in the book).
ALL_SRC = \
	contents.md \
	$(wildcard intermediate/python/*.md) \
	$(wildcard intermediate/doit/*.md) \
	$(MOST_SRC)

# Other files that the site depends on.
EXTRAS = \
       $(wildcard css/*.css) \
       $(wildcard css/*/*.css)

# Principal target files
INDEX = $(SITE)/index.html
