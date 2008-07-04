CHMCONVERT = ../src/chmconvert

nobase_dist_pkgdata_DATA = osu2003.chm osu2003.osu osu2005.osu osu2005.chm \
	osu2007.osu osu2007.chm udfa2006.udfa udfa2006.chm \
	udfa2006_dipole.udfa udfa2006_dipole.chm

osu2003.chm: osu2003.osu
	${CHMCONVERT} -o $@ $?

osu2005.chm: osu2005.osu
	${CHMCONVERT} -o $@ $?

osu2007.chm: osu2007.osu
	${CHMCONVERT} -o $@ $?

udfa2006.chm: udfa2006.udfa
	${CHMCONVERT} -o $@ $?

udfa2006_dipole.chm: udfa2006_dipole.udfa
	${CHMCONVERT} -o $@ $?

CLEANFILES = osu2003.chm osu2005.chm osu2007.chm udfa2006.chm udfa2006_dipole.chm

