/* 
   network_test.c - Test the solve() function
   
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../src/astrochem.h"

/* FixMe: Why do we need to define these here? (segfault otherwise) */
double abundances[MAX_SHELLS][MAX_TIME_STEPS][MAX_OUTPUT_ABUNDANCES];
struct rout routes[MAX_SHELLS][MAX_TIME_STEPS][MAX_OUTPUT_ABUNDANCES][N_OUTPUT_ROUTES];

int
main (void)
{
  FILE *f;

  struct inp input_params;
  char chem_file[MAX_LINE]; 
  char source_file[MAX_LINE];
  char suffix[MAX_LINE];
  double chi;
  double cosmic;
  double grain_size;
  double grain_abundance;
  double ti;
  double tf;
  double abs_err;
  double rel_err;
  struct abund initial_abundances[MAX_INITIAL_ABUNDANCES];
  int n_initial_abundances;
  char *output_species[MAX_OUTPUT_ABUNDANCES];
  int n_output_species;
  int time_steps; 
  int trace_routes;

  struct mdl source_mdl;
  int n_shells;
  int shell[MAX_SHELLS];
  double av[MAX_SHELLS];
  double nh[MAX_SHELLS];
  double tgas[MAX_SHELLS];
  double tdust[MAX_SHELLS];
  int shell_index;

  struct net network;
  struct react reactions[MAX_REACTIONS];
  char *species[MAX_SPECIES];
  int n_reactions;
  int n_species;

  double tim[MAX_TIME_STEPS];

  int verbose = 0;

  /* Create the input.ini, source.mdl and network_chm files */

  f = fopen ("input.ini", "w");
  fprintf (f, "# This input file was created by solve_test\n");
  fprintf (f, "[files]\n");
  fprintf (f, "source = source.mdl\n");
  fprintf (f, "chem = network.chm\n");
  fprintf (f, "# Solver parameters\n");
  fprintf (f, "[solver]\n");
  fprintf (f, "ti = 1e-6\n");
  fprintf (f, "tf = 1e7\n");
  fprintf (f, "abs_err = 1e-15\n");
  fprintf (f, "rel_err = 1e-6\n");
  fprintf (f, "# Initial abundances\n");
  fprintf (f, "[abundances]\n");
  fprintf (f, "X = 1.0\n");
  fprintf (f, "Y = 0.0\n");
  fprintf (f, "# Output\n");
  fprintf (f, "[output]\n");
  fprintf (f, "time_steps = 128\n");
  fprintf (f, "abundances = X,Y\n");
  fclose (f);

  f = fopen ("source.mdl", "w");
  fprintf (f, "# This source model file was created by solve_test\n");
  fprintf (f, "0   20.0    1e+04    10.0    10.0\n");
  fclose (f);

  f = fopen ("network.chm", "w");
  fprintf (f, "# This network file was created by solve_test\n");
  fprintf (f, "X -> Y    1e-9    0    0    2    1\n");
  fclose (f);

  /* Read them */

  read_input ("input.ini", &input_params,verbose);

  read_source ("source.mdl", &source_mdl, verbose);

  read_network ("network.chm", &network, verbose);

  /* Solve the ODE system */

  {
    int i;

    for (i = 0; i < time_steps; i++)
      {   
	if (i < MAX_TIME_STEPS)
	  tim[i] = pow (10., log10 (ti) + i * (log10 (tf) - log10(ti)) 
			/ (time_steps - 1));
	else
	  return EXIT_FAILURE;
      }
  }

  shell_index = 0.;
  solve (chi, cosmic, grain_size, grain_abundance, 
	 abs_err, rel_err, initial_abundances,
	 n_initial_abundances, output_species,
	 n_output_species, av[shell_index],
	 nh[shell_index], tgas[shell_index],
	 tdust[shell_index], reactions, n_reactions,
	 species, n_species, shell_index, tim,
	 time_steps, abundances, trace_routes,
	 routes, verbose);


  /* Check the abundances */

  {
    int i;
    double x_abundance;
    double y_abundance;
    double x_abs_err;
    double x_rel_err;
    double y_abs_err;
    double y_rel_err;

    for (i = 0; i < time_steps; i++)
      {
	x_abundance = 1.0 * exp (-1e-9 * tim[i]);
	y_abundance = 1.0 - x_abundance;

	x_abs_err = fabs(abundances[0][i][0] - x_abundance);
	y_abs_err = fabs(abundances[0][i][1] - y_abundance);
	x_rel_err = x_abs_err / x_abundance;
	y_rel_err = y_abs_err / y_abundance;

	/* Errors accumulate after each time step, so the actual error
	   on the abundance is somewhat larger than the solver
	   relative tolerance. */
	
	if ((x_abs_err > abs_err) && (x_rel_err > rel_err * 5e2))
	  {
	    fprintf (stderr, "solve_test: %s:%d: incorrect abundance at t=%12.6e: expected %12.6e, got %12.6e.\n",
		     __FILE__, __LINE__, tim[i], x_abundance, abundances[0][i][0]); 
	    return EXIT_FAILURE;
	    }

	if ((y_abs_err > abs_err) && (y_rel_err > rel_err * 5e2))
	  {
	    fprintf (stderr, "solve_test: %s:%d: incorrect abundance at t=%12.6e: expected %12.6e, got %12.6e.\n",
		     __FILE__, __LINE__, tim[i], y_abundance, abundances[0][i][1]); 
	    return EXIT_FAILURE;
	  }
      }
  }
    
  return EXIT_SUCCESS;
}

