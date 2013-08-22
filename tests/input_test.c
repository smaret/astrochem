/* 
   input_test.c - Test the read_input() function
   
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
#include "../src/astrochem.h"

int
main (void)
{
  FILE *f;
  inp_t input_params;
  int verbose = 0;

  /* Create the input.ini file */

  f = fopen ("input.ini", "w");
  fprintf (f, "# This input file was created by input_test\n");
  fprintf (f, "[files]\n");
  fprintf (f, "source = source.mdl\n");
  fprintf (f, "chem = ../../networks/osu2008.chm\n");
  fprintf (f, "# Physical paramaters\n");
  fprintf (f, "[phys]\n");
  fprintf (f, "cosmic = 1.3e-17\n");
  fprintf (f, "# Solver parameters\n");
  fprintf (f, "[solver]\n");
  fprintf (f, "ti = 1e-6\n");
  fprintf (f, "tf = 1e9\n");
  fprintf (f, "rel_err = 1e-6\n");
  fprintf (f, "# Initial abundances\n");
  fprintf (f, "[abundances]\n");
  fprintf (f, "H2      = 0.5\n");
  fprintf (f, "He      = 0.14\n");
  fprintf (f, "N       = 2.14e-5\n");
  fprintf (f, "O       = 1.76e-4\n");
  fprintf (f, "C(+)    = 7.30e-5\n");
  fprintf (f, "S(+)    = 8.00e-8\n");
  fprintf (f, "Si(+)   = 8.00e-8\n");
  fprintf (f, "Fe(+)   = 3.00e-9\n");
  fprintf (f, "Na(+)   = 2.00e-9\n");
  fprintf (f, "Mg(+)   = 7.00e-9\n");
  fprintf (f, "P(+)    = 3.00e-9\n");
  fprintf (f, "Cl(+)   = 4.00e-9\n");
  fprintf (f, "e(-)    = 7.32e-5\n");
  fprintf (f, "grain   = 1.32e-12\n");
  fprintf (f, "# Output\n");
  fprintf (f, "[output]\n");
  fprintf (f, "time_steps = 64\n");
  fprintf (f, "abundances = CO,C(+),C,e(-),OH,H3O(+),H,H2,HCO(+),CO(+),C4H,HCO(+),CH(+),CH\n");
  fprintf (f, "trace_routes = 1\n");
  fclose (f);

  /* Read it */

  read_input ("input.ini", &input_params, verbose);
  
  /* Check that the values are correct */
  
  if ((strcmp (input_params.files.source_file, "source.mdl") == 0) &&
      (strcmp (input_params.files.chem_file, "../../networks/osu2008.chm") == 0) &&
      (input_params.phys.chi == CHI_DEFAULT) &&
      (input_params.phys.cosmic == 1.3e-17) && 
      (input_params.phys.grain_abundance == 1.32e-12) && 
      (input_params.solver.ti == 1e-6 * CONST_MKSA_YEAR) &&
      (input_params.solver.tf == 1e9 * CONST_MKSA_YEAR) && 
      (input_params.solver.abs_err == ABS_ERR_DEFAULT) && 
      (input_params.solver.rel_err == 1e-6) && 
      (input_params.output.time_steps == 64) &&
      (input_params.abundances.n_initial_abundances == 14) &&
      (strcmp (input_params.abundances.initial_abundances[0].specie, "H2") == 0) &&
      (input_params.abundances.initial_abundances[0].abundance == 0.5) &&
      (strcmp (input_params.abundances.initial_abundances[1].specie, "He") == 0) &&
      (input_params.abundances.initial_abundances[1].abundance == 0.14) &&
      (strcmp (input_params.abundances.initial_abundances[2].specie, "N") == 0) &&
      (input_params.abundances.initial_abundances[2].abundance == 2.14e-5) &&
      (strcmp (input_params.abundances.initial_abundances[3].specie, "O") == 0) &&
      (input_params.abundances.initial_abundances[3].abundance == 1.76e-4) &&
      (strcmp (input_params.abundances.initial_abundances[4].specie, "C(+)") == 0) &&
      (input_params.abundances.initial_abundances[4].abundance == 7.30e-5) &&
      (input_params.output.n_output_species == 14) &&
      (strcmp (input_params.output.output_species[0], "CO") == 0) &&
      (strcmp (input_params.output.output_species[1], "C(+)") == 0) &&
      (strcmp (input_params.output.output_species[2], "C") == 0) &&
      (strcmp (input_params.output.output_species[3], "e(-)") == 0) &&
      (strcmp (input_params.output.output_species[4], "OH") == 0))
    {
      free_input_struct (&input_params);
      return EXIT_SUCCESS;
    }
  else
    {
      free_input_struct (&input_params);
      return EXIT_FAILURE;
    }
}

