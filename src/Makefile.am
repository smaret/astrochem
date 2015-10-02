#  Makefile.am - Automake file for Astrochem
#
#  Copyright (c) 2006-2015 Sebastien Maret
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

AM_CFLAGS += ${OPENMP_CFLAGS} -Wall


lib_LTLIBRARIES = libastrochem.la
libastrochem_la_SOURCES = network.c rates.c solve.c input.c
libastrochem_la_LDFLAGS = -version-info 0:0:0
include_HEADERS = libastrochem.h

bin_PROGRAMS = astrochem chmconvert converttolegacy

astrochem_SOURCES = astrochem.c 
astrochem_CFLAGS = ${HDF5_CFLAGS} ${AM_CFLAGS}
astrochem_CPPFLAGS = ${HDF5_CPPFLAGS} ${AM_CPPFLAGS}
astrochem_LDFLAGS = ${HDF5_LDFLAGS} ${AM_LDFLAGS}
astrochem_LDADD = -lastrochem ${HDF5_LIBS} ${AM_LDFLAGS}

chmconvert_SOURCES = chmconvert.py
chmconvert$(EXEEXT): chmconvert.py
	cp $? $@
	chmod a+x $@

converttolegacy_SOURCES = converttolegacy.py
converttolegacy$(EXEEXT): converttolegacy.py
	cp $? $@
	chmod a+x $@


CLEANFILES = chmconvert.py  converttolegacy.py

if HAVE_PYMOD_BIGGLES

bin_PROGRAMS += plabun plroute

plabun_SOURCES = plabun.py

plroute_SOURCES = plroute.py

plabun$(EXEEXT): plabun.py
	cp $? $@
	chmod a+x $@

plroute$(EXEEXT): plroute.py
	cp $? $@
	chmod a+x $@

CLEANFILES += plabun.py plroute.py

endif

install-exec-hook:
	sed -i' ' -e '/@REMOVE_INST@/d' -e 's,$${prefix},${prefix},g' $(DESTDIR)$(bindir)/chmconvert
	chmod a+x $(DESTDIR)$(bindir)/chmconvert
	sed -i' ' -e '/@REMOVE_INST@/d' -e 's,$${prefix},${prefix},g' $(DESTDIR)$(bindir)/converttolegacy
	chmod a+x $(DESTDIR)$(bindir)/converttolegacy
if HAVE_PYMOD_BIGGLES
	sed -i' ' -e '/@REMOVE_INST@/d' -e 's,$${prefix},${prefix},g' $(DESTDIR)$(bindir)/plabun
	chmod a+x $(DESTDIR)$(bindir)/plabun
	sed -i' ' -e '/@REMOVE_INST@/d' -e 's,$${prefix},${prefix},g' $(DESTDIR)$(bindir)/plroute
	chmod a+x $(DESTDIR)$(bindir)/plroute
endif
