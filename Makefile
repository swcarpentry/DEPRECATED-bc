all : commands

## commands : show all commands
commands :
	@grep -E '^##' Makefile | sed -e 's/## //g'

## check    : build locally into _site directory for checking
check :
	jekyll -t build -d _site

## clean    : clean up
clean :
	rm -rf _site $$(find . -name '*~' -print)

## links      : check links
# depends on linklint, an html-link-checking module from http://www.linklint.org/

links :
	@linklint -doc /tmp/linkdoc -root _site /@

