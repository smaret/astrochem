#  Makefile.am - Automake file for Astrochem
#
#  Copyright (c) 2006-2014 Sebastien Maret
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
libastrochem_la_SOURCES = input.c network.c rates.c solve.c output.c
libastrochem_la_LDFLAGS = -version-info 0:0:0
include_HEADERS = libastrochem.h

bin_PROGRAMS = astrochem apiuser chmconvert

astrochem_SOURCES = astrochem.c
astrochem_LDADD = -lastrochem

apiuser_SOURCES = apiuser.c
apiuser_LDADD = -lastrochem


chmconvert_SOURCES = chmconvert.py

chmconvert$(EXEEXT): chmconvert.py
	cp $? $@
	chmod a+x $@

CLEANFILES = chmconvert.py

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
	mv $(DESTDIR)$(bindir)/chmconvert $(DESTDIR)$(bindir)/chmconvert.tmp
	sed -i 's/libastrochemtools/astrochem.libastrochemtools/g'  $(DESTDIR)$(bindir)/chmconvert.tmp
	sed -e 's,$${prefix},${prefix},g' $(DESTDIR)$(bindir)/chmconvert.tmp > $(DESTDIR)$(bindir)/chmconvert
	chmod a+x $(DESTDIR)$(bindir)/chmconvert
	rm -f $(DESTDIR)$(bindir)/chmconvert.tmp
if HAVE_PYMOD_BIGGLES
	mv $(DESTDIR)$(bindir)/plabun $(DESTDIR)$(bindir)/plabun.tmp
	sed -e 's,$${prefix},${prefix},g' $(DESTDIR)$(bindir)/plabun.tmp > $(DESTDIR)$(bindir)/plabun
	chmod a+x $(DESTDIR)$(bindir)/plabun
	rm -f $(DESTDIR)$(bindir)/plabun.tmp
	mv $(DESTDIR)$(bindir)/plroute $(DESTDIR)$(bindir)/plroute.tmp
	sed -e 's,$${prefix},${prefix},g' $(DESTDIR)$(bindir)/plroute.tmp > $(DESTDIR)$(bindir)/plroute
	chmod a+x $(DESTDIR)$(bindir)/plroute
	rm -f $(DESTDIR)$(bindir)/plroute.tmp
endif
