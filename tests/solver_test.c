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
#include "../src/libastrochem.h"

int
main (void)
{
  FILE *f;
  int cell_index;
  net_t network;
  int verbose = 0;

  /* Create the input.ini, source.mdl and network_chm files */

  f = fopen ("network.chm", "w");
  fprintf (f, "# This network file was created by full_solve_test\n");
  fprintf (f, "X -> Y    1e-9    0    0    2    1\n");
  fclose (f);

  /* Read them */
  read_network ("network.chm", &network, verbose);

  phys_t phys;
  phys.cosmic = COSMIC_DEFAULT;
  phys.chi = CHI_DEFAULT;
  phys.grain_size = GRAIN_SIZE_DEFAULT;
  phys.grain_abundance = 0;

  double abs_err, rel_err;
  abs_err = 1e-15;
  rel_err = 1e-6;

  const char* species[]  = { "X", "Y" };
  const double initial_abundances[] = { 1.0, 0.0 };

  double *abundances;
  alloc_abundances( &network, &abundances ); // Allocate the abundances array; it contains all species.
  set_initial_abundances(species, 2, initial_abundances, &network, abundances); // Set initial abundances

  double density = 1000;
  double av = 20;
  double temperature = 10;

  cell_t cell;
  cell.nh = density;
  cell.av = av;
  cell.tgas = temperature;
  cell.tdust = temperature; // Assume tgas = tdust in this specific case

  astrochem_mem_t astrochem_mem;

  if( solver_init( &cell, &network, &phys, abundances , density, abs_err, rel_err, &astrochem_mem ) != 0 )
    {
      return EXIT_FAILURE;
    }
  int i;
  double time = 0;
  for( i = 0; i < 1000 ; i++ )
    {

      time+= 10;
      solve( &astrochem_mem, &network, abundances, time, NULL,verbose);


      double x_abundance;
      double y_abundance;
      double x_abs_err;
      double x_rel_err;
      double y_abs_err;
      double y_rel_err;

      x_abundance = 1.0 * exp (-1e-9 *  time);
      y_abundance = 1.0 - x_abundance;
      x_abs_err = fabs( abundances[0] - x_abundance);
      y_abs_err = fabs( abundances[1] - y_abundance);
      x_rel_err = x_abs_err / x_abundance;
      y_rel_err = y_abs_err / y_abundance;

      /* Errors accumulate after each time step, so the actual error
         on the abundance is somewhat larger than the solver
         relative tolerance. */
      if ((x_abs_err > abs_err) && (x_rel_err > rel_err * 5e2))
        {
          fprintf (stderr, "full_solve_test: %s:%d: incorrect abundance at t=%12.15e: expected %12.15e, got %12.15e.\n",
                   __FILE__, __LINE__, time, x_abundance, abundances[0]);
          solver_close( &astrochem_mem );
          free_abundances( abundances );
          free_network (&network);
          return EXIT_FAILURE;
        }

      if ((y_abs_err > abs_err) && (y_rel_err > rel_err * 5e2))
        {
          fprintf (stderr, "full_solve_test: %s:%d: incorrect abundance at t=%12.6e: expected %12.6e, got %12.6e.\n",
                   __FILE__, __LINE__, time, y_abundance, abundances[1]);
          solver_close( &astrochem_mem );
          free_abundances( abundances );
          free_network (&network);
          return EXIT_FAILURE;
        }
    }
  solver_close( &astrochem_mem );
  free_abundances( abundances );
  free_network (&network);

  return EXIT_SUCCESS;
}
