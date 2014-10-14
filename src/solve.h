/*
   solve.h - private solve Function prototypes, various constant and data
   structures for Astrochem.

   Copyright (c) 2006-2013 Sebastien Maret

   This file is part of Astrochem.

   Astrochem is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   Astrochem is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Astrochem.  If not, see <http://www.gnu.org/licenses/>.
   */

/* Various definitions and constants */

#ifndef _SOLVE_H_
#define _SOLVE_H_

int get_abundance_idx (const res_t * results, int cell_idx, int ts_idx,
                       int abund_idx);

int get_route_idx (const res_t * results, int cell_idx, int ts_idx,
                   int abund_idx, int route_idx);


#endif // _SOLVE_H_
