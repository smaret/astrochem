/*
   solve.c - Build the ODE system and the jacobian matrix, and solve
   the system with CVODE.

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cvode/cvode.h>
#ifdef USE_LAPACK
#include <cvode/cvode_lapack.h>
#else
#include <cvode/cvode_dense.h>
#endif

#include "astrochem.h"
#include "rates.h"

static int f (realtype t, N_Vector y, N_Vector ydot, void *params);

static int jacobian (long int N, realtype t, N_Vector y, N_Vector fy,
                     DlsMat J, void *params, N_Vector tmp1,
                     N_Vector tmp2, N_Vector tmp3);

/*
  Right-hand side function of the ODE system.
*/

static int
f (realtype t __attribute__ ((unused)), N_Vector y, N_Vector ydot,
   void *params)
{
  int i;

  /* Get the parameters passed with *params. */

  double *reac_rates = ((params_t *) params)->reac_rates;
  const react_t *reactions = ((params_t *) params)->reactions;
  int n_reactions = ((params_t *) params)->n_reactions;
  int n_species = ((params_t *) params)->n_species;
  double nh = ((params_t *) params)->nh;
  double av = ((params_t *) params)->av;
  double tgas = ((params_t *) params)->tgas;
  double tdust = ((params_t *) params)->tdust;
  double chi = ((params_t *) params)->chi;
  double cosmic = ((params_t *) params)->cosmic;
  double grain_size = ((params_t *) params)->grain_size;
  double grain_abundance = ((params_t *) params)->grain_abundance;

  /* Loop over the reactions and build the right hand ODE system
     equations. */

  for (i = 0; i < n_species; i++)
    {
      NV_Ith_S (ydot, i) = 0;
    }

  for (i = 0; i < n_reactions; i++)
    {
      double y_product;

      /* Compute the production/destruction rate of the reaction. */

      y_product = 1;

      if (reactions[i].reaction_type == 0)
        {
          /* The formation of H2 is a first order reaction, so the
             product of the reactants needs to be divided by the H
             abundance. */
          y_product *= NV_Ith_S (y, reactions[i].reactants[0]);
          y_product *= reac_rates[i];
        }
      else if (reactions[i].reaction_type == 23)
        {
          /* Photo-desorption is a zeroth order reaction. However the
             reaction rate depends on the ice thickness, so it must be
             recomputed at each time step. */

          reac_rates[i] = rate (reactions[i].alpha,
                                reactions[i].beta,
                                reactions[i].gamma,
                                reactions[i].reaction_type,
                                reactions[i].reaction_no,
                                nh, av, tgas, tdust,
                                chi, cosmic,
                                grain_size,
                                grain_abundance,
                                NV_Ith_S (y, reactions[i].reactants[0]) / nh);
          y_product = reac_rates[i];
        }
      else
        {
          /* For other reactions, the production/destruction rate is
             the product of the reactants multiplied by the reaction rate. */
          int r;
          for( r = 0; r < MAX_REACTANTS; r++ )
            {
              if (reactions[i].reactants[r] != -1)
                y_product *= NV_Ith_S (y, reactions[i].reactants[r]);
            }
          y_product *= reac_rates[i];
        }
      /* Add a production term for each product of the reaction. */
      int p;
      for( p = 0 ; p < MAX_PRODUCTS; p++ )
        {
          if (reactions[i].products[p] != -1)
            NV_Ith_S (ydot, reactions[i].products[p]) += y_product;
        }
      /* Add a destruction term for each reactants of the reaction. */
      int r;
      for( r = 0; r < MAX_REACTANTS; r++ )
        {
          if (reactions[i].reactants[r] != -1)
            NV_Ith_S (ydot, reactions[i].reactants[r]) -= y_product;
        }
    }

  return (0);
}

/*
   Jacobian matrix.
   */

static int
jacobian (long int N __attribute__ ((unused)),
          realtype t __attribute__ ((unused)), N_Vector y,
          N_Vector fy __attribute__ ((unused)),
          DlsMat J, void *params,
          N_Vector tmp1 __attribute__ ((unused)),
          N_Vector tmp2 __attribute__ ((unused)),
          N_Vector tmp3 __attribute__ ((unused)))
{
  int i, j;

  /* Get the parameters passed with *params. */

  double *reac_rates = ((params_t *) params)->reac_rates;
  const react_t *reactions = ((params_t *) params)->reactions;
  int n_reactions = ((params_t *) params)->n_reactions;
  int n_species = ((params_t *) params)->n_species;
  double nh = ((params_t *) params)->nh;
  double av = ((params_t *) params)->av;
  double chi = ((params_t *) params)->chi;
  double grain_size = ((params_t *) params)->grain_size;
  double grain_abundance = ((params_t *) params)->grain_abundance;

  /* Compute the jacobian matrix. */

  for (i = 0; i < n_species; i++)
    {
      for (j = 0; j < n_species; j++)
        {
          DENSE_ELEM (J, i, j) = 0;
        }
    }

  for (i = 0; i < n_reactions; i++)
    {

      /* Compute the Jacobian matrix elements corresponding to
         the reaction. */

      double y_product;

      if (reactions[i].reaction_type == 0)
        {
          /* H2 formation on grains */

          DENSE_ELEM (J, reactions[i].reactants[0], reactions[i].reactants[0]) -=
           2 * reac_rates[i];
          DENSE_ELEM (J, reactions[i].products[0], reactions[i].reactants[0]) +=
           reac_rates[i];
        }
      else if (reactions[i].reaction_type == 23)
        {
          /* Photo-desorption */

	  if (grain_abundance != 0)
	    {
	      double jac_elem;

	      jac_elem =
		(chi * DRAINE_STANDARD_ISRF_FUV * exp (-2 * av) * reactions[i].alpha) /
		(GRAIN_SITES_PER_CM2 * reactions[i].gamma) *
		exp (-NV_Ith_S (y, reactions[i].reactants[0]) /
		     (GRAIN_SITES_PER_CM2 * M_PI * pow (grain_size, 2) *
		      grain_abundance * nh * reactions[i].gamma));

	      DENSE_ELEM (J, reactions[i].reactants[0], reactions[i].reactants[0]) -=
		jac_elem;
	      DENSE_ELEM (J, reactions[i].products[0], reactions[i].reactants[0]) +=
		jac_elem;
	    }
        }
      else
        {
          /* Other reactions */
          int r, r2, p;
          int nreactants = MAX_REACTANTS;
          for( r = 0; r < MAX_REACTANTS; r++ )
            {
              if( reactions[i].reactants[r] == -1 )
                {
                  nreactants = r;
                  break;
                }
            }
          for( r = 0; r < nreactants; r++ )
            {
              y_product =  reac_rates[i];
              for( r2 = 0; r2 < nreactants; r2++ )
                {
                  if ( r != r2 )
                    {
                      y_product *=  NV_Ith_S (y, reactions[i].reactants[ r2 ]);
                    }
                }

              for( r2 = 0; r2 < nreactants; r2++ )
                {
                  DENSE_ELEM (J, reactions[i].reactants[r2], reactions[i].reactants[r]) -= y_product;
                }
              for( p = 0; p < MAX_PRODUCTS; p++ )
                {
                  if(  reactions[i].products[p] != -1 )
                    {
                      DENSE_ELEM (J, reactions[i].products[p], reactions[i].reactants[r]) += y_product;
                    }
                }
            }
        }
    }

  return (0);
}



/**
 * Initialize the solver
 */
int solver_init( const cell_t* cell, const net_t* network, const phys_t* phys,
                 const double* abundances , double density, double abs_err, double rel_err,
                 astrochem_mem_t* astrochem_mem )
{
  astrochem_mem->density = density;

  /* Allocate the work array and fill it with the initial
     concentrations. Ignore species that are not in the
     network. */

  if (( astrochem_mem->y = N_VNew_Serial (network->n_species)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__,
               "array allocation failed.\n");
      return -1;
    }
  int i;
  for (i = 0; i < network->n_species; i++)
    {
      NV_Ith_S ( astrochem_mem->y, i) = abundances[i] * density;
    }

  /* Allocate an array for the reaction rates and compute them. */

  if ((astrochem_mem->params.reac_rates =
       malloc (sizeof (double) * (unsigned int) network->n_reactions)) ==
      NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
               __FILE__, __LINE__);
      return -1;
    }
  for (i = 0; i < network->n_reactions; i++)
    {
      astrochem_mem->params.reac_rates[i] = rate (network->reactions[i].alpha,
                                                  network->reactions[i].beta,
                                                  network->reactions[i].gamma,
                                                  network->reactions[i].reaction_type,
                                                  network->reactions[i].reaction_no,
                                                  cell->nh, cell->av, cell->tgas,
                                                  cell->tdust, phys->chi,
                                                  phys->cosmic,
                                                  phys->grain_size,
                                                  phys->grain_abundance, 0.);
    }

  /* Fill out a structure containing the parameters of the function
     defining the ODE system and the jacobian. */

  astrochem_mem->params.reactions = network->reactions;
  astrochem_mem->params.n_reactions = network->n_reactions;
  astrochem_mem->params.n_species = network->n_species;
  astrochem_mem->params.nh = cell->nh;
  astrochem_mem->params.av = cell->av;
  astrochem_mem->params.tgas = cell->tgas;
  astrochem_mem->params.tdust = cell->tdust;
  astrochem_mem->params.chi = phys->chi;
  astrochem_mem->params.cosmic = phys->cosmic;
  astrochem_mem->params.grain_size = phys->grain_size;
  astrochem_mem->params.grain_abundance = phys->grain_abundance;

  /* Define the ODE system and solve it using the Backward
     Differential Formulae method (BDF) with a Newton Iteration. The
     absolute error is multiplied by the density, because we compute
     concentrations and not abundances. */

  if (( astrochem_mem->cvode_mem = CVodeCreate (CV_BDF, CV_NEWTON)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: solver memory allocation failed.\n",
               __FILE__, __LINE__);
      return -1;
    }

  if ((CVodeInit ( astrochem_mem->cvode_mem, f, 0.0,  astrochem_mem->y ) != CV_SUCCESS)
      ||
      (CVodeSStolerances
       ( astrochem_mem->cvode_mem, rel_err,
         abs_err * density ) != CV_SUCCESS)
#ifdef USE_LAPACK
      || ((CVLapackDense ( astrochem_mem->cvode_mem, network->n_species) != CV_SUCCESS))
#else
      || ((CVDense ( astrochem_mem->cvode_mem, network->n_species) != CV_SUCCESS))
#endif
      || ((CVDlsSetDenseJacFn (  astrochem_mem->cvode_mem, jacobian) != CV_SUCCESS))
      || (CVodeSetUserData ( astrochem_mem->cvode_mem, &astrochem_mem->params) != CV_SUCCESS)
      || (CVodeSetMaxNumSteps (astrochem_mem->cvode_mem, CVODE_MXSTEPS) != CV_SUCCESS))
    {
      fprintf (stderr, "astrochem: %s:%d: solver initialization failed.\n",
               __FILE__, __LINE__);
      return -1;
    }
  return 0;
}

/*
  Solve the ODE system
 */

int solve( astrochem_mem_t* astrochem_mem, const net_t* network, double* abundances, double time , const cell_t* new_cell, int verbose )
{

  /* Computing new reac rate and abundance if cell parameter have evolved
     since last call */
  if( new_cell != NULL )
    {
      int i,j;
      for (j = 0; j < network->n_species; j++)
        {
          NV_Ith_S (astrochem_mem->y, j) *= new_cell->nh / astrochem_mem->params.nh;
        }

      for (i = 0; i < network->n_reactions; i++)
        {
          astrochem_mem->params.reac_rates[i] = rate (network->reactions[i].alpha,
                                                      network->reactions[i].beta,
                                                      network->reactions[i].gamma,
                                                      network->reactions[i].reaction_type,
                                                      network->reactions[i].reaction_no,
                                                      new_cell->nh, new_cell->av, new_cell->tgas,
                                                      new_cell->tdust, astrochem_mem->params.chi,
                                                      astrochem_mem->params.cosmic,
                                                      astrochem_mem->params.grain_size,
                                                      astrochem_mem->params.grain_abundance, 0.);
        }

      /* Fill out a structure containing the parameters of the function
         defining the ODE system and the jacobian. */

      astrochem_mem->params.nh = new_cell->nh;
      astrochem_mem->params.av = new_cell->av;
      astrochem_mem->params.tgas = new_cell->tgas;
      astrochem_mem->params.tdust = new_cell->tdust;
    }

  realtype t = 0.0;
  CVode ( astrochem_mem->cvode_mem, (realtype) time, astrochem_mem->y, &t, CV_NORMAL);

  /* Print the cell number, time and time step after each call. */

  if (verbose >= 2)
    {
      realtype h;

      CVodeGetLastStep (astrochem_mem->cvode_mem, &h);
      fprintf (stdout, "t = %8.2e  delta_t = %8.2e\n",
               (double) t / CONST_MKSA_YEAR,
               (double) h / CONST_MKSA_YEAR);
    }
  int i;
  for( i = 0; i < network->n_species; i++ )
    {
      abundances[i] =  (double) NV_Ith_S ( astrochem_mem->y , i ) / astrochem_mem->density ;
    }
  return 0;
}

/**
 * Close the solver
 */
void solver_close( astrochem_mem_t* astrochem_mem )
{
  if( astrochem_mem->params.reac_rates != NULL )
    {
      free (astrochem_mem->params.reac_rates);
    }
  if( astrochem_mem->y != NULL )
    {
      N_VDestroy_Serial (astrochem_mem->y);
    }
  if( astrochem_mem->cvode_mem != NULL )
    {
      CVodeFree (&astrochem_mem->cvode_mem);
    }
}

/**
 *  Allocate the abundance vector
 */
int alloc_abundances ( const net_t* network, double** abundances )
{
  if (( *abundances =
        malloc (sizeof (double) * network->n_species )) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
               __FILE__, __LINE__);
      return -1;
    }
  int i;
  for( i = 0 ; i < network->n_species ; i++ )
    {
      (*abundances)[i]=0;
    }
  return 0;
}

/**
 *  Initialise specific abundances
 */
int set_initial_abundances( const char** species, int n_initialized_species, const double* initial_abundances,
                            const net_t* network, double* abundances )
{
  int i,j;
  for( i = 0 ; i < network->n_alloc_species ; i++ )
    {
      for( j = 0 ; j < n_initialized_species ; j++ )
        {
          if( strcmp( network->species[i].name , species[j] ) == 0 )
            {
              abundances[i] = initial_abundances[j];
            }
        }
    }
  return 0;
}

/**
 * Free the abundances
 */
void free_abundances ( double* abundances )
{
  if( abundances != NULL )
    {
      free ( abundances);
    }
}
