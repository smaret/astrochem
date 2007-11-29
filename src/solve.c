/* 
   Build the ODE system and the jacobian matrix, and solve the
   system with CVODE.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <setjmp.h>
#include <string.h>
#include <math.h>

#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
 
#include "astrochem.h"

struct par {
  double *reac_rates;
  struct react *reactions;
  int n_reactions;
  int n_species;
};

static jmp_buf env;

static int f (realtype t, N_Vector y, N_Vector ydot, void *f_data);

static int jacobian (long int N, DenseMat J, realtype t, N_Vector y,
		     N_Vector fy, void *jacobian_data, N_Vector tmp1,
		     N_Vector tmp2, N_Vector tmp3);

void interrupt_handler (int sig);

/*
  Vector defining the ODE system.
*/

static int 
f (realtype t __attribute__ ((unused)), N_Vector y, N_Vector ydot,
   void *f_data)
{
  int i;
  
  /* Get the parameters passed with *f_data. */

  double *reac_rates = ((struct par *)f_data)->reac_rates;
  struct react *reactions = ((struct par *)f_data)->reactions;
  int n_reactions = ((struct par *)f_data)->n_reactions;
  int n_species = ((struct par *)f_data)->n_species;

  /* Loop over the reactions and build the right hand ODE system
     equations. */
  
  for (i = 0; i < n_species; i++)
    {
      NV_Ith_S (ydot, i) = 0;
    }

  for (i = 0; i < n_reactions; i++)
    {
      double y_product;

      /* Compute the product of the reactants multiplied by the
	 rate of the reaction. */

      y_product = 1;

      if (reactions[i].reactant1 != -1)
	y_product *= NV_Ith_S (y, reactions[i].reactant1);
      if (reactions[i].reactant2 != -1)
	y_product *= NV_Ith_S (y, reactions[i].reactant2);
      if (reactions[i].reactant3 != -1)
	y_product *= NV_Ith_S (y, reactions[i].reactant3);
      y_product *= reac_rates[i];

      /* Add a production term for each product of the reaction. */

      if (reactions[i].product1 != -1)
	NV_Ith_S (ydot, reactions[i].product1) += y_product;
      if (reactions[i].product2 != -1)
	NV_Ith_S (ydot, reactions[i].product2) += y_product;
      if (reactions[i].product3 != -1)
	NV_Ith_S (ydot, reactions[i].product3) += y_product;
      if (reactions[i].product4 != -1)
	NV_Ith_S (ydot, reactions[i].product4) += y_product;
      
      /* Add a destruction term for each reactant of the reaction. */

      if (reactions[i].reactant1 != -1)
	NV_Ith_S (ydot, reactions[i].reactant1) -= y_product;
      if (reactions[i].reactant2 != -1)
	NV_Ith_S (ydot, reactions[i].reactant2) -= y_product;
      if (reactions[i].reactant3 != -1)
	NV_Ith_S (ydot, reactions[i].reactant3) -= y_product;
    }
  
  return (0);
}

/*
  Jacobian matrix.
*/
	  
static int
jacobian (long int N __attribute__ ((unused)), DenseMat J, 
	  realtype t __attribute__ ((unused)), N_Vector y,
	  N_Vector fy __attribute__ ((unused)), void *jacobian_data,
	  N_Vector tmp1 __attribute__ ((unused)), 
	  N_Vector tmp2 __attribute__ ((unused)), 
	  N_Vector tmp3 __attribute__ ((unused)))
{
  int i, j;
  
  /* Get the parameters passed with *jacobian_data. */

  double *reac_rates = ((struct par *)jacobian_data)->reac_rates;
  struct react *reactions = ((struct par *)jacobian_data)->reactions;
  int n_reactions = ((struct par *)jacobian_data)->n_reactions;
  int n_species = ((struct par *)jacobian_data)->n_species;
    
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
      
      if ((reactions[i].reactant1 != -1) &&
	  (reactions[i].reactant2 == -1) &&
	  (reactions[i].reactant3 == -1))
	{
	  DENSE_ELEM (J, reactions[i].reactant1, reactions[i].reactant1) -=
	    reac_rates[i];
	  DENSE_ELEM (J, reactions[i].product1, reactions[i].reactant1) +=
	    reac_rates[i];
	  if (reactions[i].product2 != -1)
	    DENSE_ELEM (J, reactions[i].product2, reactions[i].reactant1) +=
	      reac_rates[i];
	  if (reactions[i].product3 != -1)
	    DENSE_ELEM (J, reactions[i].product3, reactions[i].reactant1) +=
	      reac_rates[i];
	  if (reactions[i].product4 != -1)
	    DENSE_ELEM (J, reactions[i].product4, reactions[i].reactant1) +=
	      reac_rates[i];
	}

      if ((reactions[i].reactant1 != -1) &&
	  (reactions[i].reactant2 != -1) &&
	  (reactions[i].reactant3 == -1))
	{
	  y_product = reac_rates[i] * NV_Ith_S (y, reactions[i].reactant2);

	  DENSE_ELEM (J, reactions[i].reactant1, reactions[i].reactant1) -=
	    y_product;
	  DENSE_ELEM (J, reactions[i].reactant2, reactions[i].reactant1) -=
	    y_product;
	  DENSE_ELEM (J, reactions[i].product1, reactions[i].reactant1) +=
	    y_product;
	  if (reactions[i].product2 != -1)
	    DENSE_ELEM (J, reactions[i].product2, reactions[i].reactant1) +=
	      y_product;
	  if (reactions[i].product3 != -1)
	    DENSE_ELEM (J, reactions[i].product3, reactions[i].reactant1) +=
	      y_product;
	  if (reactions[i].product4 != -1)
	    DENSE_ELEM (J, reactions[i].product4, reactions[i].reactant1) +=
	      y_product;

	  y_product = reac_rates[i] * NV_Ith_S (y, reactions[i].reactant1);

	  DENSE_ELEM (J, reactions[i].reactant1, reactions[i].reactant2) -=
	    y_product;
	  DENSE_ELEM (J, reactions[i].reactant2, reactions[i].reactant2) -=
	    y_product;
	  DENSE_ELEM (J, reactions[i].product1, reactions[i].reactant2) +=
	    y_product;
	  if (reactions[i].product2 != -1)
	    DENSE_ELEM (J, reactions[i].product2, reactions[i].reactant2) +=
	      y_product;
	  if (reactions[i].product3 != -1)
	    DENSE_ELEM (J, reactions[i].product3, reactions[i].reactant2) +=
	      y_product;
	  if (reactions[i].product4 != -1)
	    DENSE_ELEM (J, reactions[i].product4, reactions[i].reactant2) +=
	      y_product;
	}

      if ((reactions[i].reactant1 != -1) &&
	  (reactions[i].reactant2 != -1) &&
	  (reactions[i].reactant3 != -1))
	{
	  y_product = reac_rates[i] * NV_Ith_S (y, reactions[i].reactant2)
	    * NV_Ith_S (y, reactions[i].reactant3);

	  DENSE_ELEM (J, reactions[i].reactant1, reactions[i].reactant1) -=
	    y_product;
	  DENSE_ELEM (J, reactions[i].reactant2, reactions[i].reactant1) -=
	    y_product;
	  DENSE_ELEM (J, reactions[i].reactant3, reactions[i].reactant1) -=
	    y_product;
	  DENSE_ELEM (J, reactions[i].product1, reactions[i].reactant1) +=
	    y_product;
	  if (reactions[i].product2 != -1)
	    DENSE_ELEM (J, reactions[i].product2, reactions[i].reactant1) +=
	      y_product;
	  if (reactions[i].product3 != -1)
	    DENSE_ELEM (J, reactions[i].product3, reactions[i].reactant1) +=
	      y_product;
	  if (reactions[i].product4 != -1)
	    DENSE_ELEM (J, reactions[i].product4, reactions[i].reactant1) +=
	      y_product;

	  y_product = reac_rates[i] * NV_Ith_S (y, reactions[i].reactant1)
	    * NV_Ith_S (y, reactions[i].reactant3);

	  DENSE_ELEM (J, reactions[i].reactant1, reactions[i].reactant2) -=
	    y_product;
	  DENSE_ELEM (J, reactions[i].reactant2, reactions[i].reactant2) -=
	    y_product;
	  DENSE_ELEM (J, reactions[i].reactant3, reactions[i].reactant2) -=
	    y_product;
	  DENSE_ELEM (J, reactions[i].product1, reactions[i].reactant2) +=
	    y_product;
	  if (reactions[i].product2 != -1)
	    DENSE_ELEM (J, reactions[i].product2, reactions[i].reactant2) +=
	      y_product;
	  if (reactions[i].product3 != -1)
	    DENSE_ELEM (J, reactions[i].product3, reactions[i].reactant2) +=
	      y_product;
	  if (reactions[i].product4 != -1)
	    DENSE_ELEM (J, reactions[i].product4, reactions[i].reactant2) +=
	      y_product;

	  y_product = reac_rates[i] * NV_Ith_S (y, reactions[i].reactant1)
	    * NV_Ith_S (y, reactions[i].reactant2);

	  DENSE_ELEM (J, reactions[i].reactant1, reactions[i].reactant3) -=
	    y_product;
	  DENSE_ELEM (J, reactions[i].reactant2, reactions[i].reactant3) -=
	    y_product;
	  DENSE_ELEM (J, reactions[i].reactant3, reactions[i].reactant3) -=
	    y_product;
	  DENSE_ELEM (J, reactions[i].product1, reactions[i].reactant3) +=
	    y_product;
	  if (reactions[i].product2 != -1)
	    DENSE_ELEM (J, reactions[i].product2, reactions[i].reactant3) +=
	      y_product;
	  if (reactions[i].product3 != -1)
	    DENSE_ELEM (J, reactions[i].product3, reactions[i].reactant3) +=
	      y_product;
	  if (reactions[i].product4 != -1)
	    DENSE_ELEM (J, reactions[i].product4, reactions[i].reactant3) +=
	      y_product;
	}
    }

  return (0);
}

/*
  Handle interruption signal while solving equations.
*/

void
interrupt_handler (int sig __attribute__ ((unused)))
{
  longjmp (env, 1);
}

/*
  Solve the ODE system.
*/
 
int
solve (double chi, double pdyield, double cosmic,
       double abs_err, double rel_err,
       struct abund initial_abundances[],
       int n_initial_abundances, char *output_species[],
       int n_output_species, double av, double nh, double tgas,
       double tdust, struct react reactions[],
       int n_reactions, char *species[], int n_species,
       int shell_index, double tim[], int time_steps,
       double abundances[MAX_SHELLS][MAX_TIME_STEPS][MAX_OUTPUT_ABUNDANCES],
       int trace_routes, 
       struct rout routes[MAX_SHELLS][MAX_TIME_STEPS][MAX_OUTPUT_ABUNDANCES][N_OUTPUT_ROUTES],
       int verbose)
{
   realtype t = 0.0;
   struct par params;                  /* Parameters for f() and jacobian() */
   N_Vector y;                         /* Work array for the solver */
   void *cvode_mem;                    /* Memory space for the solver */
   double *reac_rates;                 /* Reaction rates */

   /* Allocate the work array and fill it with the initial
      concentrations. Ignore species that are not in the
      network. */
   
   if ((y = N_VNew_Serial (n_species)) == NULL)
     {
       fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__, 
		"array allocation failed.\n");
       exit(1);
     }
   {
     int i;

     for (i = 0; i < n_species; i++)
       {
	 NV_Ith_S (y, i) = 0;
       }
     for (i = 0; i < n_initial_abundances; i++)
       {
	 int spec_index;

	 spec_index = specie_index ((initial_abundances[i].specie),
				    species, n_species);
	 if (spec_index != -2)
	   NV_Ith_S (y, spec_index) = initial_abundances[i].abundance * nh;
       }
   }

   /* Allocate an array for the reaction rates and compute them. */
   
   if ((reac_rates = malloc (sizeof (double) * n_reactions)) == NULL)
     {
       fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
		__FILE__, __LINE__); 
       exit(1);
     }
   {
     int i;

     for (i = 0; i < n_reactions; i++)
       {
	 reac_rates[i] = rate (reactions[i].alpha,
			       reactions[i].beta,
			       reactions[i].gamma,
			       reactions[i].reaction_type,
			       reactions[i].reaction_no,
			       av, tgas, tdust,
			       chi, pdyield, cosmic);
       }
   }
   
   /* Allocate and fill out a structure containing the parameters of
      the function defining the ODE system and the jacobian. */
   
   if ((params.reactions = malloc (sizeof (species))) == NULL)
     {
       fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
		__FILE__, __LINE__); 
       exit(1);
     }
   params.reac_rates = reac_rates;
   params.reactions = reactions;
   params.n_reactions = n_reactions;
   params.n_species = n_species;

   /* Define the ODE system and solve it using the Backward
      Differential Formulae method (BDF) with a Newton Iteration. The
      absolute error is multiplied by the density, because we compute
      concentrations and not abundances. */
      
   if ((cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON)) == NULL)
     {
       fprintf (stderr, "astrochem: %s:%d: solver memory allocation failed.\n",
		__FILE__, __LINE__); 
       exit(1);
     }
   
   abs_err = abs_err * nh;
   if (CVodeMalloc (cvode_mem, f, 0.0, y, CV_SS, rel_err, &abs_err) != CV_SUCCESS)
     {
       fprintf (stderr, "astrochem: %s:%d: solver memory initialization failed.\n",
		__FILE__, __LINE__); 
       exit(1);
     }

   CVDense (cvode_mem, n_species);
   CVDenseSetJacFn (cvode_mem, jacobian, &params);
   CVodeSetFdata(cvode_mem, &params);
   
   {
     int i, j;
     void (*orig_interrupt_handler)(int);
     
     /* Solve the system for each time step. Return if an interruption
	signal is encountered. */
     
     /*orig_interrupt_handler = signal (SIGINT, interrupt_handler);

     if (setjmp(env) == 0)
       {
	 ;
       }
     else
       {
	 N_VDestroy_Serial (y);
	 CVodeFree (&cvode_mem);

	 return (1);
	 } */
     
     for (i = 0; i < time_steps; i++)
       {
	 CVode (cvode_mem, (realtype) tim[i], y, &t, CV_NORMAL);
	 
	 /* Print the shell number, time and time step after each call. */
	 
	 if (verbose >= 2)
	   {
	     realtype h;

	     CVodeGetLastStep (cvode_mem, &h);
	     fprintf (stdout, "shell = %3i  t = %8.2e  delta_t = %8.2e\n",
		      shell_index, (double) t / CONST_MKSA_YEAR, 
		      (double) h / CONST_MKSA_YEAR);
	   }
	 
	 /* Fill the array of abundances with the output species
	    abundances. Ignore species that are not in the
	    network. Abundance that are lower than MIN_ABUNDANCES are
	    set to 0. */
	 
	 for (j = 0; j < n_output_species; j++)
	   {
	     int spec_index;
	     
	     spec_index = specie_index (output_species [j], species, n_species);
	     if (spec_index != -2)
	       {
		 abundances[shell_index][i][j] = (double) NV_Ith_S (y, spec_index) / nh;
		 if (abundances[shell_index][i][j] < MIN_ABUNDANCE)
		   abundances[shell_index][i][j] = 0.;
	       }
	   }

	 /* Compute the rate of each formation/destruction route for
	    each output specie. */

	 if (trace_routes)
	   {
	     for (j = 0; j < n_output_species; j++)
	       {
		 int spec_index;
		 int k;
		 double r[MAX_REACTIONS][N_OUTPUT_ROUTES];
		 
		 spec_index = specie_index (output_species [j], species, n_species);
		 if (spec_index != -2)
		   {
		     for (k= 0; k < n_reactions; k++)
		       {
			 
		       }
		   }
	       }
	   }
       }

     /* Restore original signal handler */
     
     /* signal (SIGINT, orig_interrupt_handler); */
   }

   N_VDestroy_Serial (y);
   CVodeFree (&cvode_mem);

   return (0);
 
 }

