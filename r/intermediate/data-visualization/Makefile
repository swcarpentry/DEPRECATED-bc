slides : ggplot.pdf 
all : slides


# This  rule accepts PDF targets with corresponding latex 
# source, and makes into PDFs.
ggplot.pdf : ggplot.Rnw data/*.csv
	Rscript -e "library(knitr); knit('ggplot.Rnw', quiet = TRUE)"
	pdflatex ggplot.tex


# Remove all intermediate outputs
clean :
	rm *.log *.nav *.out *.snm *.toc *.vrb *.aux *.tex
	rm -r cache
	rm -r figure

# Remove all PDF outputs then build them again
rebuild : clean slides