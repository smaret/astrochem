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

AM_CFLAGS = -Wall -I$(top_builddir)/src/

TESTS = input_test source_test source_dyn_test network_test rate_test abundances_test solver_test python_lib_test python_tools_test astrochem_test

check_PROGRAMS = input_test source_test source_dyn_test network_test rate_test abundances_test solver_test python_lib_test python_tools_test astrochem_test

input_test_SOURCES = input_test.c
input_test_LDADD = $(top_builddir)/src/input.o $(top_builddir)/src/libastrochem.la

source_test_SOURCES = source_test.c
source_test_LDADD = $(top_builddir)/src/libastrochem.la  $(top_builddir)/src/input.o

source_dyn_test_SOURCES = source_dyn_test.c
source_dyn_test_LDADD = $(top_builddir)/src/libastrochem.la $(top_builddir)/src/input.o

network_test_SOURCES = network_test.c
network_test_LDADD = $(top_builddir)/src/libastrochem.la

rate_test_SOURCES = rate_test.c
rate_test_LDADD = $(top_builddir)/src/libastrochem.la

abundances_test_SOURCES = abundances_test.c
abundances_test_LDADD = $(top_builddir)/src/libastrochem.la

solver_test_SOURCES = solver_test.c
solver_test_LDADD = $(top_builddir)/src/libastrochem.la

astrochem_test$(EXEEXT): astrochem_test.py
	cp $? $@
	chmod a+x $@

python_lib_test$(EXEEXT): python_lib_test.py
	cp $? $@
	chmod a+x $@

python_tools_test$(EXEEXT): python_tools_test.py
	cp $? $@
	chmod a+x $@

benchmark_bin: benchmark.py
	cp $? $@
	chmod a+x $@

benchmark: benchmark_bin
	./benchmark_bin

benchmark_openmp_bin: benchmark_openmp.py
	cp $? $@
	chmod a+x $@

benchmark_openmp: benchmark_openmp_bin
	./benchmark_openmp_bin

CLEANFILES = astrochem_test.py python_lib_test.py			\
	python_tools_test.py benchmark.py benchmark_bin			\
	benchmark_openmp.py benchmark_openmp_bin input.ini source.mdl	\
	source_dyn.mdl network.chm *.h5
