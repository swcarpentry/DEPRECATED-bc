# hello.mk

figure-1.svg : summary-1.dat
	sgr -N -r summary-1.dat > figure-1.svg 
