# all-patterns.mk

paper.pdf : paper.tex figure-1.svg figure-2.svg
	cat $^ > $@

figure-%.svg : summary-%.dat
	python create_figure.py $@ $^

summary-%.dat : data-%-*.dat
	python stats.py $@ $^

summary-1.dat : stats.py
summary-2.dat : stats.py
