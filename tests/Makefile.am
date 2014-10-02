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

AM_CFLAGS = -Wall

TESTS = input_test source_test source_dyn_test network_test rate_test full_solve_test abundances_test solver_test python_lib_test python_tools_test

check_PROGRAMS = input_test source_test source_dyn_test network_test rate_test full_solve_test abundances_test solver_test python_lib_test python_tools_test

input_test_SOURCES = input_test.c
input_test_LDADD = $(top_builddir)/src/libastrochem.la

source_test_SOURCES = source_test.c
source_test_LDADD = $(top_builddir)/src/libastrochem.la

source_dyn_test_SOURCES = source_dyn_test.c
source_dyn_test_LDADD = $(top_builddir)/src/libastrochem.la

network_test_SOURCES = network_test.c
network_test_LDADD = $(top_builddir)/src/libastrochem.la

rate_test_SOURCES = rate_test.c
rate_test_LDADD = $(top_builddir)/src/libastrochem.la

full_solve_test_SOURCES = full_solve_test.c
full_solve_test_LDADD = $(top_builddir)/src/libastrochem.la

abundances_test_SOURCES = abundances_test.c
abundances_test_LDADD = $(top_builddir)/src/libastrochem.la

solver_test_SOURCES = solver_test.c
solver_test_LDADD = $(top_builddir)/src/libastrochem.la

python_lib_test$(EXEEXT): python_lib_test.py
	cp $? $@
	chmod a+x $@

python_tools_test$(EXEEXT): python_tools_test.py
	cp $? $@
	chmod a+x $@

benchmark: benchmark.py
	cp $? $@
	chmod a+x $@
	./benchmark

CLEANFILES = input.ini source.mdl source_dyn.mdl network.chm *.abun
