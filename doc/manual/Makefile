MANUAL=manual
DATE=`date "+%m%d"`
.PHONY: all clean

all: $(MANUAL).pdf

clean:
	rm -f $(MANUAL).pdf $(MANUAL).dvi $(MANUAL).aux $(MANUAL).log $(MANUAL).rtf *bbl *aux *log *~

$(MANUAL).pdf: $(MANUAL).dvi
	dvipdfmx $(MANUAL).dvi -o $(MANUAL).pdf
	mkdir -p PDF
	cp $(MANUAL).pdf PDF/$(MANUAL)_$(DATE).pdf

$(MANUAL).dvi: $(MANUAL).tex
	latex $(MANUAL).tex
	latex $(MANUAL).tex
	bibtex $(MANUAL)
	latex $(MANUAL).tex
	latex $(MANUAL).tex
