

FILENAME=cxs-software

#my_paper

all:	$(FILENAME).dvi
	dvips $(FILENAME).dvi
	dvipdf $(FILENAME).dvi
	acroread $(FILENAME).pdf

$(FILENAME).pdf: $(FILENAME).dvi
	dvipdf $(FILENAME).dvi

$(FILENAME).dvi: $(FILENAME).tex #$(FILENAME).bib
	latex $(FILENAME).tex
	latex $(FILENAME).tex
#	latex $(FILENAME).tex
#	bibtex $(FILENAME)



clean:
	rm -f $(FILENAME).dvi
	rm -f *~ *.aux *.log *.toc *.bak *.out *.blg *.bbl

clobber: clean
	rm -f $(FILENAME).ps $(FILENAME).pdf
