/* 
   source_test.c - Test the read_source() function
   
   Copyright (c) 2006-2011 Sebastien Maret
   
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
  char *source_mdl;

  int n_shells;
  int shell[MAX_SHELLS];
  double av[MAX_SHELLS];
  double nh[MAX_SHELLS];
  double tgas[MAX_SHELLS];
  double tdust[MAX_SHELLS];

  int verbose = 0;

  source_mdl = "# shell number, Av [mag], n(H) [cm^-3], Tgas [K], Tdust [K]\n"
    "0	 0.1	1e+02	15.0	12.0\n"
    "1	 1.0	1e+03	11.0	10.0\n"
    "2	10.0	1e+04	 8.0	 7.0\n";

  /* Create the input.ini file */

  f = fopen ("source.mdl", "w");
  fprintf (f, "# This source model file was created by source_test\n");
  fprintf (f, "%s", source_mdl);
  fclose (f);

  /* Read it */

  read_source ("source.mdl", shell, &n_shells, av, nh,
	       tgas, tdust, verbose);
  
  /* Check that the values are correct */

  if ((n_shells == 3) &&
      (av[0] == 0.1) &&
      (nh[0] == 1e2) &&
      (tgas[0] == 15) &&
      (tdust[0] == 12) &&
      (av[1] == 1.0) &&
      (nh[1] == 1e3) &&
      (tgas[1] == 11) &&
      (tdust[1] == 10) &&
      (av[2] == 10.0) &&
      (nh[2] == 1e4) &&
      (tgas[2] == 8) &&
      (tdust[2] == 7))
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

