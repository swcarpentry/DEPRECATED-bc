# pattern-rule.mk

figure-%.svg : summary-%.dat
	sgr -N -r $@ $^

summary-1.dat : data-1-*.dat
	stats.py $@ $^

summary-2.dat : data-2-*.dat
	stats.py $@ $^

summary-1.dat : stats.py
summary-2.dat : stats.py
