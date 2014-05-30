# all-patterns.mk

paper.pdf : paper.wdp figure-1.svg figure-2.svg
	wdp2pdf $<

figure-%.svg : summary-%.dat
	sgr -N -r $@ $^

summary-%.dat : data-%-*.dat
	stats.py $@ $^

summary-1.dat : stats.py
summary-2.dat : stats.py
