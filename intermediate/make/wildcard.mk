# wildcard.mk

summary-1.dat : data-1-*.dat
	python stats.py $@ $^
