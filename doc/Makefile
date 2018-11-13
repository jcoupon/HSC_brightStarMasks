# Sample Makefile for tex files.
# Replace 'file' with the base name of the tex file.
# commands:
#    'file'      creates dvi and ps file
#    bib         like all, but with updates bibtex info
#    clean       removes unneede files
#    tar         creates a backup of the dir + subdirs

source = FNS2015
DD     = $(source)-`date +"%Y%m%d"`.tgz

$(source):
	latex $(source).tex 

bib:
	latex $(source).tex
	bibtex $(source)
	latex $(source).tex
	latex $(source).tex

ps:
	dvips $(source).dvi -o $(source).ps

pdf:
	dvips -Ppdf -G0  $(source).dvi
	epstopdf $(source).ps   

pdflatex: 
	pdflatex $(source).tex

.PHONY: clean tar

clean:
	rm -f *~ $(source).b?? *.bbl *.blg *.aux *.log $(source).out $(source).dvi $(source).pdf $(source).toc

tar:
	tar czf $(DD) * --exclude=*.tgz
