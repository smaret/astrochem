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