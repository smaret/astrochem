AM_CFLAGS = -DDATADIR=\"$(datadir)\"

AM_CFLAGS += ${OPENMP_CFLAGS}

bin_PROGRAMS = astrochem

astrochem_SOURCES = astrochem.c astrochem.h input.c network.c rates.c solve.c output.c

if HAVE_PYTHON

bin_PROGRAMS += chmconvert plabun

chmconvert_SOURCES = chmconvert.py

chmconvert: chmconvert.py
	sed -e 's,\@PYTHON\@,${PYTHON},g' $? > $@
	chmod u+x chmconvert

plabun_SOURCES = plabun.py

plabun:	plabun.py
	sed -e 's,\@PYTHON\@,${PYTHON},g' $? > $@
	chmod u+x plabun

endif