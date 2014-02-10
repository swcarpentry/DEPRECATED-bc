# with-include.mk

include config.mk

WDP2PDF_FLAGS=--style ${STYLE_DIR}/euphoric.wps
SGR_FLAGS=-N -r -s ${STYLE_DIR}/euphoric.fig

paper.pdf : paper.wdp figure-1.svg figure-2.svg
	wdp2pdf ${WDP2PDF_FLAGS} $<

figure-%.svg : summary-%.dat
	sgr ${SGR_FLAGS} $@ $^

summary-%.dat : data-%-*.dat
	stats.py $@ $^

data-*-*.dat : stats.py
	touch $@
