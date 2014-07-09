# false-dependencies.mk

paper.pdf : paper.tex figure-1.svg figure-2.svg
	cat $^ > $@

figure-%.svg : summary-%.dat
	python create_figure.py $@ $^

summary-%.dat : data-%-*.dat
	python stats.py $@ $^

data-*-*.dat : stats.py
	touch $@

.SECONDARY :
