/*
   rates.c - Compute the reaction rates.

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrochem.h"

double
rate (double alpha, double beta, double gamm, int reaction_type,
      int reaction_no __attribute__ ((unused)), double nh, double av,
      double tgas, double tdust, double chi, double cosmic,
      double grain_size, double grain_abundance, double ice_abundance)
{
  double k;                     /* Reaction rate (cm^-3 s^-1) */

  /* Reactions types and rates. The nomencalture is similar to the one
     of the Ohio State University database for astrochemistry, but it
     is extended to include depletion and desorption on/from the grain
     surfaces. */

  switch (reaction_type)
    {
    case -1:
      /* Gas-grain interaction (excluding depletion and desorption),
         Electron-grain recombination. */
      k = alpha * pow (tgas / 300, beta) * GAS_DUST_NUMBER_RATIO;
      break;

    case 0:
      /* H2 formation on grains */
      k = alpha * pow (tgas / 300, beta) * nh;
      break;

    case 1:
      /* Cosmic-ray ionization (direct process). Cosmic-ray induced
         photoreactions (indirect process)   */
      k = alpha * cosmic;
      break;

    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 14:
      /* Ion-molecule reactions, charge exchange reactions (2), Negative
         ion - neutral species reactions (3), Radiative association (4),
         Associative ejection (5), Neutral + neutral -> ion + electron
         (6), Neutral-neutral chemical reactions (7), Neutral-neutral
         radiative association (8), Dissociative recombination (9),
         Radiative recombination (10), Positive ion - negative ion
         recombination (11), Electron attachment (12), Others (14) */
      k = alpha * pow (tgas / 300, beta) * exp (-gamm / tgas);
      break;

    case 13:
      /* Photo-ionization, Photo-dissociation */
      k = chi * alpha * exp (-gamm * av);
      break;

    case 20:
      /* Depletion on the grains */
      {
        double thermal_veloc = pow (8 * CONST_CGSM_BOLTZMANN * tgas
                                    / (M_PI * beta * CONST_CGSM_MASS_PROTON),
                                    0.5);
        k =
          M_PI * pow (grain_size,
                      2) * alpha * thermal_veloc * grain_abundance * nh;
        break;
      }

    case 21:
      /* Thermal desorption */
      {
        double v0 = pow (2 * GRAIN_SITES_PER_CM2 * gamm * CONST_CGSM_BOLTZMANN
                         / (M_PI * M_PI * beta * CONST_CGSM_MASS_PROTON),
                         0.5);
        k = v0 * exp (-gamm / tdust);
        break;
      }

    case 22:
      /* Cosmic ray desorption */
      {
        if (alpha == 0.0)
          {
            double v0 =
              pow (2 * GRAIN_SITES_PER_CM2 * gamm * CONST_CGSM_BOLTZMANN /
                   (M_PI * M_PI * beta * CONST_CGSM_MASS_PROTON),
                   0.5);
            k = v0 * FRACTION_TIME_GRAIN_70K * exp (-gamm / 70.);
            break;
          }
        else
          {
            k = alpha;
            break;
          }
      }

    case 23:
      /* Photo-desorption */
      {
	if (grain_abundance != 0)
	  {
	    double x = (ice_abundance * nh) / (GRAIN_SITES_PER_CM2 * M_PI
					       * pow (grain_size, 2)
					       * grain_abundance * nh);
	    double Ypd = alpha * (1 - exp (-x / gamm));
	    k = chi * DRAINE_STANDARD_ISRF_FUV * exp (-2 * av) * M_PI
	      * pow (grain_size, 2) * grain_abundance * nh * Ypd;
	  }
	else
	  {
	    k = 0;
	  }
	break;
      }

    default:
      fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__,
               "unknown reaction type.\n");
      exit (1);
    }
  return k;
}
