
SUBDIRS = figures 


all: pysbr.pdf 

pysbr.pdf: pysbr.tex pysbr.bbl pysbr.aux $(SUBDIRS)
	pdflatex -shell-escape $< 
	grep Rerun pysbr.log && pdflatex $< 
	grep Rerun pysbr.log && pdflatex $< 
	touch pysbr.bbl
	touch $@

pysbr.bbl: pysbr.aux references.bib myplainnat.bst
	bibtex pysbr

pysbr.aux: pysbr.tex $(SUBDIRS)
	pdflatex -shell-escape $<

subdirs: $(SUBDIRS)


$(SUBDIRS):
	$(MAKE) -C $@

cleanall: clean
	rm -f *.pdf *.acn *.glo *.ist *.toc *.spl *.dvi *.ps *.fff

clean:
	rm -f *.aux *.log *.out *.blg *.bbl
