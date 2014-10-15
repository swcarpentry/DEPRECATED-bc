# phony.mk

all : figure-1.svg figure-2.svg

figure-1.svg : summary-1.dat
	python create_figure.py figure-1.svg summary-1.dat

figure-2.svg : summary-2.dat
	python create_figure.py figure-2.svg summary-2.dat
