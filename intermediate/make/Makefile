# Makefile
# Creates all files necessary for Intermediate Make lesson

help :
	@echo 'Did you mean to run make with no arguments?'
	@echo 'To create all the necessary files for this lesson, run'
	@echo ' make lesson'

lesson : summary-1.dat summary-2.dat

figure-%.svg : summary-%.dat
	python create_figure.py $@ $^

summary-%.dat : data-%-*.dat stats.py
	python stats.py $@ $^

data-%.dat : generate_dataset.py
	python $^

clean :
	rm -f *.dat *.svg
