# with-macro.mk

STYLE_DIR=c:/papers/

paper.pdf : paper.wdp figure-1.svg figure-2.svg
	wdp2pdf --style ${STYLE_DIR}/euphoric.wps $<

figure-%.svg : summary-%.dat
	sgr -N -r -s ${STYLE_DIR}/euphoric.fig $@ $^

summary-%.dat : data-%-*.dat
	stats.py $@ $^

data-*-*.dat : stats.py
	touch $@
