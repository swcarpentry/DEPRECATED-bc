LATEX = pdflatex

ALL_TEX = $(wildcard *.tex)

all : book.pdf

book.pdf : $(ALL_TEX)
	$(LATEX) book
	$(LATEX) book

clean :
	rm -f book.pdf *.aux *.log *.toc *~
