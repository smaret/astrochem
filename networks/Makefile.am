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

CHMCONVERT = ../src/chmconvert

nobase_dist_pkgdata_DATA = osu2008.osu osu2008.chm udfa2006.udfa	\
	udfa2006.chm udfa2006_dipole.udfa udfa2006_dipole.chm

osu2008.chm: osu2008.osu
	${CHMCONVERT} -o $@ $?

udfa2006.chm: udfa2006.udfa
	${CHMCONVERT} -o $@ $?

udfa2006_dipole.chm: udfa2006_dipole.udfa
	${CHMCONVERT} -o $@ $?

CLEANFILES = osu2008.chm udfa2006.chm udfa2006_dipole.chm
