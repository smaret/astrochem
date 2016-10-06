/*
   rates.h - private rates Function prototypes, various constant and data
   structures for Astrochem.

   Copyright (c) 2006-2016 Sebastien Maret

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

#ifndef _RATES_H_
#define _RATES_H_

double
rate (double alpha, double beta, double gamm, int reaction_type,
      int reaction_no __attribute__ ((unused)), double nh, double av,
      double tgas, double tdust, double chi, double cosmic,
      double grain_size, double grain_abundance, double ice_abundance);

#endif // _RATES_H_
