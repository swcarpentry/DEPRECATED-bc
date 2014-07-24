# Makefile
# Creates all files necessary for Intermediate Make lesson

help :
	@echo 'Did you mean to run make with no arguments?'
	@echo 'To create all the necessary files for this lesson, run'
	@echo ' make lesson'

lesson : paper.pdf

paper.pdf : paper.tex figure-1.svg figure-2.svg
	cat $^ > $@

figure-%.svg : summary-%.dat
	python create_figure.py $@ $^

figure-%.svg : summary-%.dat
	python create_figure.py $@ $^

summary-%.dat : data-%-*.dat stats.py
	python stats.py $@ $^

data-%.dat : generate_dataset.py
	python $^

clean :
	rm -f *.dat *.svg paper.pdf

.SECONDARY :
