#  Makefile.am - Automake file for Astrochem
#
#  Copyright (c) 2006-2011 Sebastien Maret
# 
#  This file is part of Astrochem.
#
#  Astrochem is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  Astrochem is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with Astrochem.  If not, see <http://www.gnu.org/licenses/>.

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

