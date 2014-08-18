/* 
   rate_test.c - Test the rate() function
   
   Copyright (c) 2006-2014 Sebastien Maret
   
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../src/libastrochem.h"
#include "../src/rates.h"

int equaltol (double value1, double value2, double reltol);

int
equaltol (double value1, double value2, double reltol)
{
  /* Checks if two double are equal within a given relative
     tolerance */
  
  if (value1 == 0)
    return 1;
  if (fabs ((value1 - value2) / value1) <= reltol)
    return 0;
  else
    return 1;
}
  
int
main (void)
{
  double alpha = 0;
  double beta = 0;
  double gamm = 0;
  int reaction_type = 0;
  int reaction_no = 0;
  double nh  = 1e4;
  double av = 10;
  double tgas = 10;
  double tdust = 10;
  double chi = 1;
  double cosmic = 1.3e-17;
  double grain_size = 1e-5; /* cm */
  double grain_abundance = 1.32e-12;
  double ice_abundance = 0.;
  double k;

  /* e(-) attachment on grains:
     e(-) + grain -> grain(-) */
  
  reaction_type = -1;
  alpha = 6.90e-15;
  beta = .5;
  gamm = 0;
  k = rate (alpha, beta, gamm, reaction_type, reaction_no, nh,
	    av, tgas, tdust, chi, cosmic, grain_size, grain_abundance,
	    ice_abundance);
  if (equaltol (9.536397e-4, k, 1e-6) == 1)
    return EXIT_FAILURE;

  /* H2 formation on grains:
     H + H -> H2 */
  
  reaction_type = 0;
  alpha = 4.95e-17;
  beta = .5;
  gamm = 0;
  k = rate (alpha, beta, gamm, reaction_type, reaction_no, nh,
	    av, tgas, tdust, chi, cosmic, grain_size, grain_abundance,
	    ice_abundance);
  if (equaltol (9.037422e-14, k, 1e-6) == 1)
    return EXIT_FAILURE;

  /* Cosmic-ray ionisation reaction:
     H2 + cosmic-ray -> H2(+) + e(-) */

  reaction_type = 1;
  alpha = 9.3e-01;

  k = rate (alpha, beta, gamm, reaction_type, reaction_no, nh,
	    av, tgas, tdust, chi, cosmic, grain_size, grain_abundance,
	    ice_abundance);
  if (equaltol (1.209000e-17, k, 1e-6) == 1)
    return EXIT_FAILURE;

  /* Ion-neutral reaction (no barrier):
     H3(+) + H2O -> H3O(+) + H */

  reaction_type = 2;
  alpha = 4.5e-09;
  beta = -.5;
  gamm = 0.;

  k = rate (alpha, beta, gamm, reaction_type, reaction_no, nh,
	    av, tgas, tdust, chi, cosmic, grain_size, grain_abundance,
	    ice_abundance);
  if (equaltol (2.464752e-8, k, 1e-6) == 1)
    return EXIT_FAILURE;

  /* Neutral-neutral reaction (with large barrier):
     H2 + OH -> H2O + H */

  reaction_type = 7;
  alpha = 8.4e-13;
  beta = 0;
  gamm = 1.04e3;

  k = rate (alpha, beta, gamm, reaction_type, reaction_no, nh,
	    av, tgas, tdust, chi, cosmic, grain_size, grain_abundance,
	    ice_abundance);
  if (equaltol (5.723387e-58, k, 1e-6) == 1)
    return EXIT_FAILURE;

  /* Photo-ionisation
     CO + uv-photon -> C + O */

  reaction_type = 13;
  alpha = 3.1e-11;
  beta = 0;
  gamm = 2.5;

  k = rate (alpha, beta, gamm, reaction_type, reaction_no, nh,
	    av, tgas, tdust, chi, cosmic, grain_size, grain_abundance,
	    ice_abundance);
  if (equaltol (4.305263e-22, k, 1e-6) == 1)
    return EXIT_FAILURE;

  /* Depletion
     CO -> CO(ice) */

  reaction_type = 20;
  alpha = 1;
  beta = 28;
  gamm = 0;

  k = rate (alpha, beta, gamm, reaction_type, reaction_no, nh,
	    av, tgas, tdust, chi, cosmic, grain_size, grain_abundance,
	    ice_abundance);
  if (equaltol (3.593005e-14, k, 1e-6) == 1)
    return EXIT_FAILURE;

  /* Thermal desorption
     CO(ice) -> CO */

  reaction_type = 21;
  alpha = 0;
  beta = 28;
  gamm = 1180;

  k = rate (alpha, beta, gamm, reaction_type, reaction_no, nh,
	    av, tgas, tdust, chi, cosmic, grain_size, grain_abundance,
	    ice_abundance);
  if (equaltol (8.239140e-40, k, 1e-2) == 1)
    return EXIT_FAILURE;

  /* Cosmic-ray desorption (from Hasegawa & Herbst 1993)
     CO(ice) + cosmic-ray -> CO */

  reaction_type = 22;
  alpha = 0;
  beta = 28;
  gamm = 1180;

  k = rate (alpha, beta, gamm, reaction_type, reaction_no, nh,
	    av, tgas, tdust, chi, cosmic, grain_size, grain_abundance,
	    ice_abundance);
  if (equaltol (2.194592e-14, k, 1e-6) == 1)
    return EXIT_FAILURE;

  /* Cosmic-ray desorption (from Bring & Johnson 2004)
     CO(ice) + cosmic-ray -> CO */

  reaction_type = 22;
  alpha = 4.3e-16;
  beta = 0;
  gamm = 0;

  k = rate (alpha, beta, gamm, reaction_type, reaction_no, nh,
	    av, tgas, tdust, chi, cosmic, grain_size, grain_abundance,
	    ice_abundance);
  if (equaltol (4.3e-16, k, 1e-6) == 1)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}

