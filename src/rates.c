/*
  Compute the reaction rates.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrochem.h"

#define TEMPERATURE_COSMIC_RAY 70
#define GAS_DUST_NUMBER_RATIO 7.57e+11

double 
rate(double alpha, double beta, double gamm, int reaction_type,
     int reaction_no, double av, double tgas, double tdust,
     double chi, double pdyield, double cosmic) 
{  
  double k; /* Reaction rate (cm^-3 s^-1) */

  /* Reactions types and rates. The nomencalture is similar to the one
     of the Ohio State University database for astrochemistry, but it
     is extended to include depletion and desorption on/from the grain
     surfaces. */

  /* FixMe: Should we check the values of rates with an exponential? */    
  
  switch (reaction_type)
    {
    case 0:
      /* Gas-grain interaction (excluding depletion and desorption), 
	 Electron-grain recombination. */
      /*k = alpha * pow (tgas / 300, beta) * GAS_DUST_NUMBER_RATIO;*/
      k = alpha * pow (tgas / 300, beta);
      break;
      
    case 1:
      /* Cosmic-ray ionization (direct process). Cosmic-ray induced 
	 photoreactions (indirect process)   */
      k = alpha * cosmic;
      break;
      
    case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 9: case 10: 
    case 11: case 12: case 14:
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
      /* Depletion on the grains 
	 FixMe: Only the sticking coefficient should go in here. The
	 grain size should be an input parameter. */
      {
	double thermal_veloc = pow (8 * CONST_CGSM_BOLTZMANN * tgas
				    / (M_PI * beta * CONST_CGSM_MASS_PROTON),
				    0.5);
	k = alpha * gamm * thermal_veloc;
	break;
      }
      
    case 21:
      /* Thermal desorption */
      {
	double nu0 = pow (2 * alpha * CONST_CGSM_BOLTZMANN / 
			  (pow (M_PI, 2) * CONST_CGSM_MASS_PROTON * beta),
			  0.5);
	k = nu0 * exp (-gamm / tdust );
	break;
      }
      
    case 22:
      /* Cosmic ray desorption */
      {
	double nu0 = pow (2 * alpha * CONST_CGSM_BOLTZMANN / 
			  (pow (M_PI, 2) * CONST_CGSM_MASS_PROTON * beta),
			  0.5);
	k = nu0 * exp (-gamm / TEMPERATURE_COSMIC_RAY);
	break;
      }
      
    case 23:
      /* Photo-desorption */
      k = chi * exp(-2 * av) * alpha * pdyield;
      break;
      
    default:
      fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__, 
	       "unknown reaction type.\n");
      exit (1);
    }
  return k;  
}

  
  
