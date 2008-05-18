AM_CFLAGS = -DDATADIR=\"$(datadir)\"

AM_CFLAGS += ${OPENMP_CFLAGS}

bin_PROGRAMS = astrochem chmconvert

astrochem_SOURCES = astrochem.c astrochem.h input.c network.c rates.c solve.c output.c

chmconvert_SOURCES = chmconvert.py

chmconvert: chmconvert.py
	sed -e 's,\@PYTHON\@,${PYTHON},g' $? > $@
	chmod u+x chmconvert

if HAVE_PYMOD_BIGGLES

bin_PROGRAMS += plabun

plabun_SOURCES = plabun.py

plabun:	plabun.py
	sed -e 's,\@PYTHON\@,${PYTHON},g' $? > $@
	chmod u+x plabun

endif