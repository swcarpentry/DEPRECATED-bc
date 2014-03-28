# double.mk

figure-1.svg : summary-1.dat
	sgr -N -r summary-1.dat > figure-1.svg

figure-2.svg : summary-2.dat
	sgr -N -r summary-2.dat > figure-2.svg
