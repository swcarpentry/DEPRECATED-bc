# hello.mk

figure-1.svg : summary-1.dat
	python create_figure.py figure-1.svg summary-1.dat
