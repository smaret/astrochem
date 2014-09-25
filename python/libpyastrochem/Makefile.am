#  Makefile.am - Automake file for Astrochem
#
#  Copyright (c) 2006-2013 Sebastien Maret
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

if HAVE_CYTHON

.pyx.c:
	@CYTHON@ -o $@ $<

CYTHON_CFLAGS = -pthread -fPIC -fwrapv -Wall -fno-strict-aliasing 
CYTHON_LDFLAGS = -pthread -fPIC -fwrapv -Wall -fno-strict-aliasing -shared 

pkgpyexec_LTLIBRARIES = libpyastrochem.la
libpyastrochem_la_SOURCES = libpyastrochem.pyx
libpyastrochem_la_CFLAGS = $(PYTHON_CPPFLAGS) $(CYTHON_CFLAGS)
libpyastrochem_la_LIBADD = $(top_builddir)/src/libastrochem.la
libpyastrochem_la_LDFLAGS = $(CYTHON_LDFLAGS) $(PYTHON_LDFLAGS)  -version-info 0:0:0 -module

endif

CLEANFILES = *.c
