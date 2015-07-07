/* 
   network_test.c - Test the read_network() function
   
   Copyright (c) 2006-2015 Sebastien Maret
   
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
#include "libastrochem.h"
#include "network.h"

int
main (void)
{
  FILE *f;
  char chem_file[] = "network.chm";
  net_t network;
  
  int verbose = 1;

  /* Create the network.chm file */

  f = fopen ("network.chm", "w");
  fprintf (f, "# This network file was created by network test\n");
  fprintf (f, "# The reaction were extracted form the OSU 2008 network\n");
  fprintf (f, "H            + H                           -> H2                                                           4.95e-17  5.00e-01  0.00e+00  0    1\n");
  fprintf (f, "C(+)         + grain(-)                    -> C            + grain                                         4.90e-17  5.00e-01  0.00e+00  0    3\n");
  fprintf (f, "H3(+)        + grain(-)                    -> H2           + H            + grain                          1.00e-16  5.00e-01  0.00e+00  0   13\n");
  fprintf (f, "C            + cosmic-ray                  -> C(+)         + e(-)                                          1.02e+03  0.00e+00  0.00e+00  1   15\n");
  fprintf (f, "CH5N         + cosmic-ray                  -> HCN          + H2           + H            + H               1.41e+03  0.00e+00  0.00e+00  1  176\n");
  fprintf (f, "C(+)         + Fe                          -> Fe(+)        + C                                             2.60e-09  0.00e+00  0.00e+00  2  218\n");
  fprintf (f, "He(+)        + HNC                         -> C(+)         + N            + H            + He              4.43e-09 -5.00e-01  0.00e+00  2  735\n");
  fprintf (f, "C(-)         + NO                          -> CN(-)        + O                                             1.00e-09  0.00e+00  0.00e+00  3 3151\n");
  fprintf (f, "C(+)         + H                           -> CH(+)                                                        1.70e-17  0.00e+00  0.00e+00  4 3162\n");
  fprintf (f, "C(-)         + C                           -> C2           + e(-)                                          5.00e-10  0.00e+00  0.00e+00  5 3243\n");
  fprintf (f, "O            + CH                          -> HCO(+)       + e(-)                                          2.00e-11  4.40e-01  0.00e+00  6 3289\n");
  fprintf (f, "C            + CH                          -> C2           + H                                             6.59e-11  0.00e+00  0.00e+00  7 3290\n");
  fprintf (f, "C            + C                           -> C2           + photon                                        1.00e-17  0.00e+00  0.00e+00  8 3672\n");
  fprintf (f, "C2(+)        + e(-)                        -> C            + C                                             8.84e-08 -5.00e-01  0.00e+00  9 3688\n");
  fprintf (f, "C(+)         + e(-)                        -> C            + photon                                        4.40e-12 -6.10e-01  0.00e+00 10 4227\n");
  fprintf (f, "C(+)         + C(-)                        -> C            + C                                             2.30e-07 -5.00e-01  0.00e+00 11 4243\n");
  fprintf (f, "C            + e(-)                        -> C(-)                                                         3.00e-15  0.00e+00  0.00e+00 12 4279\n");
  fprintf (f, "O            + -13-CH                      -> H-13-CO(+)   + e(-)                                          2.00e-11  4.40e-01  0.00e+00  6 3289\n");
  fprintf (f, "CO                                         -> CO(ice)                                                      1.00e+00  2.80e+01  0.00e+00 20 10044\n");
  fclose (f);

  /* Read it */

  if( read_network(chem_file, &network, verbose) != EXIT_SUCCESS )
    {
      return EXIT_FAILURE;
    }

  /* Check that the values are correct */

  if ((network.n_reactions == 19) &&
      (network.n_species == 29) &&

      /* Reaction #1 */
      (network.reactions[0].reactants[0] == find_species("H", &network)) &&
      (network.reactions[0].reactants[1] ==  find_species("H", &network)) &&
      (network.reactions[0].reactants[2] == -1) &&
      (network.reactions[0].products[0] ==  find_species("H2", &network)) &&
      (network.reactions[0].products[1] == -1) &&
      (network.reactions[0].products[2] == -1) &&
      (network.reactions[0].products[3] == -1) &&
      (network.reactions[0].alpha == 4.95e-17) &&
      (network.reactions[0].beta == .5) &&
      (network.reactions[0].gamma == 0) &&
      (network.reactions[0].reaction_type == 0) &&
      (network.reactions[0].reaction_no == 1) &&

      /* Reaction #176 */
      (network.reactions[4].reactants[0] == find_species("CH5N", &network)) &&
      (network.reactions[4].reactants[1] == -1) &&
      (network.reactions[4].reactants[2] == -1) &&
      (network.reactions[4].products[0] ==  find_species("HCN", &network)) &&
      (network.reactions[4].products[1] ==  find_species("H2", &network)) &&
      (network.reactions[4].products[2] ==  find_species("H", &network)) &&
      (network.reactions[4].products[3] ==  find_species("H", &network) )&&
      (network.reactions[4].alpha == 1.41e3) &&
      (network.reactions[4].beta == 0) &&
      (network.reactions[4].gamma == 0) &&
      (network.reactions[4].reaction_type == 1) &&
      (network.reactions[4].reaction_no == 176) &&

      /* Reaction #4227 */
      (network.reactions[14].reactants[0] == find_species("C(+)", &network) )&&
      (network.reactions[14].reactants[1] == find_species("e(-)", &network) )&&
      (network.reactions[14].reactants[2] == -1) &&
      (network.reactions[14].products[0] == find_species("C", &network) )&&
      (network.reactions[14].products[1] == -1) &&
      (network.reactions[14].products[2] == -1) &&
      (network.reactions[14].products[3] == -1) &&
      (network.reactions[14].alpha == 4.40e-12) &&
      (network.reactions[14].beta == -.61) &&
      (network.reactions[14].gamma == 0) &&
      (network.reactions[14].reaction_type == 10) &&
      (network.reactions[14].reaction_no == 4227)  &&

      /* Mass and charge of a few species */

      (abs(network.species[find_species("C(+)", &network)].mass / UMA - 12) <= 0.01) &&
      (network.species[find_species("C(+)", &network)].charge == 1.0) &&
      (abs(network.species[find_species("HCO(+)", &network)].mass / UMA - 29) <= 0.01) &&
      (network.species[find_species("HCO(+)", &network)].charge == 1.0) &&
      (abs(network.species[find_species("H-13-CO(+)", &network)].mass / UMA - 30) <= 0.01) &&
      (network.species[find_species("H-13-CO(+)", &network)].charge == 1.0) &&
      (abs(network.species[find_species("CH5N", &network)].mass / UMA - 31) <= 0.01) &&
      (network.species[find_species("CH5N", &network)].charge == 0.0) &&
      (network.species[find_species("e(-)", &network)].charge == -1))
    {
      free_network (&network);
      return EXIT_SUCCESS;
    }
  else
    {
      free_network (&network);
      return EXIT_FAILURE;
    }
}

