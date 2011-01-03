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

AM_CFLAGS = -DPKGDATADIR=\"$(pkgdatadir)\"

AM_CFLAGS += ${OPENMP_CFLAGS}

bin_PROGRAMS = astrochem chmconvert

astrochem_SOURCES = astrochem.c astrochem.h input.c network.c rates.c solve.c output.c

chmconvert_SOURCES = chmconvert.py

chmconvert: chmconvert.py
	sed -e 's,\@PYTHON\@,${PYTHON},g' $? > $@
	chmod a+x chmconvert

if HAVE_PYMOD_BIGGLES 
if HAVE_PYMOD_NUMPY

bin_PROGRAMS += plabun plroute

plabun_SOURCES = plabun.py

plroute_SOURCES = plroute.py

plabun:	plabun.py
	sed -e 's,\@PYTHON\@,${PYTHON},g' $? > $@
	chmod a+x plabun

plroute: plroute.py
	sed -e 's,\@PYTHON\@,${PYTHON},g' $? > $@
	chmod a+x plroute

endif
endif