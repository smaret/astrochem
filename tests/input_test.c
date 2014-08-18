/* 
   input_test.c - Test the read_input() function
   
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/libastrochem.h"

int
main (void)
{
  FILE *f;
  inp_t input_params;
  net_t network;
  int verbose = 0;

  /* Create the network.chm file */

  f = fopen ("network.chm", "w");
  fprintf (f, "# This network file was created by network test\n");
  fprintf (f, "# The reaction were extracted form the OSU 2008 network\n");
  fprintf (f, "N            + H                           -> H2                                                           4.95e-17  5.00e-01  0.00e+00  0    1\n");
  fprintf (f, "C(+)         + grain(-)                    -> C            + grain                                         4.90e-17  5.00e-01  0.00e+00  0    3\n");
  fprintf (f, "H3(+)        + H3O(+)                    -> S(+)           + H            + OH                          1.00e-16  5.00e-01  0.00e+00  0   13\n");
  fprintf (f, "C            + cosmic-ray                  -> Cl(+)         + e(-)                                          1.02e+03  0.00e+00  0.00e+00  1   15\n");
  fprintf (f, "CH5N         + cosmic-ray                  -> HCN          + H2           + H            + C4H               1.41e+03  0.00e+00  0.00e+00  1  176\n");
  fprintf (f, "CO(+)         + Fe                          -> Fe(+)        + C                                             2.60e-09  0.00e+00  0.00e+00  2  218\n");
  fprintf (f, "He(+)        + HNC                         -> C(+)         + N            + H            + He              4.43e-09 -5.00e-01  0.00e+00  2  735\n");
  fprintf (f, "C(-)         + NO                          -> CN(-)        + O                                             1.00e-09  0.00e+00  0.00e+00  3 3151\n");
  fprintf (f, "C(+)         + H                           -> CH(+)                                                        1.70e-17  0.00e+00  0.00e+00  4 3162\n");
  fprintf (f, "C(-)         + C                           -> C2           + e(-)                                          5.00e-10  0.00e+00  0.00e+00  5 3243\n");
  fprintf (f, "O            + CH                          -> HCO(+)       + e(-)                                          2.00e-11  4.40e-01  0.00e+00  6 3289\n");
  fprintf (f, "C            + CH                          -> C2           + H                                             6.59e-11  0.00e+00  0.00e+00  7 3290\n");
  fprintf (f, "C            + C                           -> CO           + photon                                        1.00e-17  0.00e+00  0.00e+00  8 3672\n");
  fprintf (f, "C2(+)        + e(-)                        -> Si(+)            + C                                             8.84e-08 -5.00e-01  0.00e+00  9 3688\n");
  fprintf (f, "C(+)         + e(-)                        -> C            + photon                                        4.40e-12 -6.10e-01  0.00e+00 10 4227\n");
  fprintf (f, "C(+)         + C(-)                        -> Na(+)            + C                                             2.30e-07 -5.00e-01  0.00e+00 11 4243\n");
  fprintf (f, "P(+)            + e(-)                        -> C(-)                                                         3.00e-15  0.00e+00  0.00e+00 12 4279\n");
  fprintf (f, "Mg(+)            + uv-photon                   -> C(+)         + e(-)                                          2.16e-10  0.00e+00  2.61e+00 13 4283\n");
  fclose (f);


  /* Create the input.ini file */

  f = fopen ("input.ini", "w");
  fprintf (f, "# This input file was created by input_test\n");
  fprintf (f, "[files]\n");
  fprintf (f, "source = source.mdl\n");
  fprintf (f, "chem = network.chm\n");
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

 /* Read the input file */
  read_input_file_names ("input.ini", &input_params.files, verbose);
  
  /* Read the chemical network file */
  read_network (input_params.files.chem_file, &network, verbose);
   
  /* Read the input file */
  read_input ("input.ini", &input_params, &network , verbose);
 
  
  /* Check that the values are correct */
    
  if ((strcmp (input_params.files.source_file, "source.mdl") == 0) &&
      (strcmp (input_params.files.chem_file, "network.chm") == 0) &&
      (input_params.phys.chi == CHI_DEFAULT) &&
      (input_params.phys.cosmic == 1.3e-17) && 
      (input_params.phys.grain_abundance == 1.32e-12) && 
      (input_params.solver.ti == 1e-6 * CONST_MKSA_YEAR) &&
      (input_params.solver.tf == 1e9 * CONST_MKSA_YEAR) && 
      (input_params.solver.abs_err == ABS_ERR_DEFAULT) && 
      (input_params.solver.rel_err == 1e-6) && 
      (input_params.output.time_steps == 64) &&
      (input_params.abundances.n_initial_abundances == 14) &&
      (strcmp (network.species_names[input_params.abundances.initial_abundances[0].species_idx], "H2") == 0) &&
      (input_params.abundances.initial_abundances[0].abundance == 0.5) &&
      (strcmp (network.species_names[input_params.abundances.initial_abundances[1].species_idx], "He") == 0) &&
      (input_params.abundances.initial_abundances[1].abundance == 0.14) &&
      (strcmp (network.species_names[input_params.abundances.initial_abundances[2].species_idx], "N") == 0) &&
      (input_params.abundances.initial_abundances[2].abundance == 2.14e-5) &&
      (strcmp (network.species_names[input_params.abundances.initial_abundances[3].species_idx], "O") == 0) &&
      (input_params.abundances.initial_abundances[3].abundance == 1.76e-4) &&
      (strcmp (network.species_names[input_params.abundances.initial_abundances[4].species_idx], "C(+)") == 0) &&
      (input_params.abundances.initial_abundances[4].abundance == 7.30e-5) &&
      (input_params.output.n_output_species == 14) &&
      (strcmp (network.species_names[input_params.output.output_species_idx[0]], "CO") == 0) &&
      (strcmp (network.species_names[input_params.output.output_species_idx[1]], "C(+)") == 0) &&
      (strcmp (network.species_names[input_params.output.output_species_idx[2]], "C") == 0) &&
      (strcmp (network.species_names[input_params.output.output_species_idx[3]], "e(-)") == 0) &&
      (strcmp (network.species_names[input_params.output.output_species_idx[4]], "OH") == 0))
    {
      free_input (&input_params);
      free_network (&network);
      return EXIT_SUCCESS;
    }
  else
    {
      free_input (&input_params);
      free_network (&network);
      return EXIT_FAILURE;
    }
}

