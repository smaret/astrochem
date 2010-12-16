/* 
   tests/input_test.c - Test the read_input() function
   
   Copyright (c) 2006-2010 Sebastien Maret
   
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
main ()
{
  FILE *f;
  char *input_ini;

  char chem_file[MAX_LINE]; 
  char source_file[MAX_LINE];
  char suffix[MAX_LINE];
  double chi;
  double cosmic;
  double grain_size;
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

  int verbose = 0;

  input_ini = "[files]\n"
    "source = source.mdl\n"
    "chem = ../../networks/osu2008.chm\n"
    "# Physical paramaters\n"
    "[phys]\n"
    "cosmic = 1.3e-17\n"
    "# Solver parameters\n"
    "[solver]\n"
    "ti = 1e-6\n"
    "tf = 1e9\n"
    "rel_err = 1e-6\n"
    "# Initial abundances\n"
    "[abundances]\n"
    "H2      = 0.5\n"
    "He      = 0.14\n"
    "N       = 2.14e-5\n"
    "O       = 1.76e-4\n"
    "C(+)    = 7.30e-5\n"
    "S(+)    = 8.00e-8\n"
    "Si(+)   = 8.00e-8\n"
    "Fe(+)   = 3.00e-9\n"
    "Na(+)   = 2.00e-9\n"
    "Mg(+)   = 7.00e-9\n"
    "P(+)    = 3.00e-9\n"
    "Cl(+)   = 4.00e-9\n"
    "e(-)    = 7.32e-5\n"
    "# Output\n"
    "[output]\n"
    "time_steps = 64\n"
    "abundances = CO,C(+),C,e(-),OH,H3O(+),H,H2,HCO(+),CO(+),C4H,HCO(+),CH(+),CH\n"
    "trace_routes = 1\n";

  /* Create the input.ini file */

  f = fopen ("input.ini", "w");
  fprintf (f, "# This input file was created by input_test\n");
  fprintf (f, "%s", input_ini);
  fclose (f);

  /* Read it */

  read_input ("input.ini", chem_file, source_file, &chi, &cosmic,
	      &grain_size, &ti, &tf, &abs_err, &rel_err,
	      initial_abundances, &n_initial_abundances,
	      output_species, &n_output_species, &time_steps,
	      &trace_routes, suffix, verbose);
  
  /* Check that the values are correct */
  
  if ((strcmp (source_file, "source.mdl") == 0) &&
      (strcmp (chem_file, "../../networks/osu2008.chm") == 0) &&
      (chi == CHI_DEFAULT) &&
      (cosmic == 1.3e-17) && 
      (ti == 1e-6 * CONST_MKSA_YEAR) &&
      (tf == 1e9 * CONST_MKSA_YEAR) && 
      (abs_err == ABS_ERR_DEFAULT) && 
      (rel_err == 1e-6) && 
      (time_steps == 64) &&
      (n_initial_abundances == 13) &&
      (strcmp (initial_abundances[0].specie, "H2") == 0) &&
      (initial_abundances[0].abundance == 0.5) &&
      (strcmp (initial_abundances[1].specie, "He") == 0) &&
      (initial_abundances[1].abundance == 0.14) &&
      (strcmp (initial_abundances[2].specie, "N") == 0) &&
      (initial_abundances[2].abundance == 2.14e-5) &&
      (strcmp (initial_abundances[3].specie, "O") == 0) &&
      (initial_abundances[3].abundance == 1.76e-4) &&
      (strcmp (initial_abundances[4].specie, "C(+)") == 0) &&
      (initial_abundances[4].abundance == 7.30e-5) &&
      (n_output_species == 14) &&
      (strcmp (output_species[0], "CO") == 0) &&
      (strcmp (output_species[1], "C(+)") == 0) &&
      (strcmp (output_species[2], "C") == 0) &&
      (strcmp (output_species[3], "e(-)") == 0) &&
      (strcmp (output_species[4], "OH") == 0))
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

