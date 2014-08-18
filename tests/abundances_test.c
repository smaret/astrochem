/*
   network.test.c - Test the read_network.) function

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
#include "../src/libastrochem.h"

    int
main (void)
{
    FILE *f;
    net_t network;

    int verbose = 0;

    /* Create the network.chm file */
    f = fopen ("network.chm", "w");
    fprintf (f, "# This network file was created by full_solve_test\n");
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
    fprintf (f, "C            + uv-photon                   -> C(+)         + e(-)                                          2.16e-10  0.00e+00  2.61e+00 13 4283\n");
    fclose (f);

    /* Read it */
    read_network ("network.chm", &network, verbose);

    const char* species[]  = {"CH", "HCO(+)", "e(-)"};
    const double initial_abundances[] = {1e-4, 1e-9, 1e-3};

    double *abundances = NULL;
    alloc_abundances( &network, &abundances ); // Allocate the abundances array; it contains all species.
    if( abundances == NULL )
    {
        return EXIT_FAILURE;
    }
    int i;
    for( i = 0; i < network.n_species; i++ )
    {
        if( abundances[i] != 0 )
        {
            return EXIT_FAILURE;
        }
    }

    set_initial_abundances(species, 3, initial_abundances, &network, abundances); // Set initial abundances
    if(  abundances[7] != initial_abundances[2] )
    {
        return EXIT_FAILURE;
    }
    if(  abundances[22] != initial_abundances[0] )
    {
        return EXIT_FAILURE;
    }
    if(  abundances[23] != initial_abundances[1] )
    {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
