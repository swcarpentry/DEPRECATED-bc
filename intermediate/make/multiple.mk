# multiple.mk

summary-1.dat : data-1-1.dat data-1-2.dat data-1-3.dat
	python stats.py summary-1.dat data-1-1.dat data-1-2.dat data-1-3.dat
