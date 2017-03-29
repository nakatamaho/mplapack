LATEX=platex
DVI2PDF=dvipdfmx
DVI2PS=dvips
TEXFILE=userman

.SUFFIXES:	.tex .dvi .ps .pdf
.tex.dvi:	
	${LATEX} $< $@
.dvi.ps:	
	${DVI2PS} $< > $@

.dvi.pdf:	
	${DVI2PDF} $<

all: ${TEXFILE}.pdf #${TEXFILE}.ps 

clean:
	rm -f ${TEXFILE}.dvi ${TEXFILE}.pdf ${TEXFILE}.log ${TEXFILE}.out ${TEXFILE}.aux texput.log *~
