SUBDIRS = examples

man_MANS = astrochem.1 plabun.1 plroute.1 chmconvert.1

EXTRA_DIST = $(man_MANS), astrochem.tex

if HAVE_PDFLATEX
if HAVE_BIBTEX

all: pdf

SUFFIXES = .pdf .tex

pdf: astrochem.pdf

astrochem.pdf: astrochem.tex examples/abundances.pdf
	$(PDFLATEX) astrochem.tex
	$(BIBTEX) astrochem.aux
	$(PDFLATEX) astrochem.tex
	$(PDFLATEX) astrochem.tex

CLEANFILES = astrochem.aux astrochem.log astrochem.bbl astrochem.blg astrochem.pdf

endif
endif

