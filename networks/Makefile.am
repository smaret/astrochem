#  Makefile.am - Automake file for Astrochem
#
#  Copyright (c) 2006-2025 Sebastien Maret
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

EXTRA_DIST = osu2008.osu osu2009.osu

nobase_dist_pkgdata_DATA = osu2008.chm osu2009.chm

osu2008.chm: osu2008.osu
	${CHMCONVERT} -o $@ $?

osu2009.chm: osu2009.osu
	${CHMCONVERT} -o $@ $?

CLEANFILES = osu2008.chm osu2009.chm

