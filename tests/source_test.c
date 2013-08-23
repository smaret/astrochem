/* 
   source_test.c - Test the read_source() function
   
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

  mdl_t source_mdl;
  int verbose = 0;

  /* Create the input.ini file */

  f = fopen ("source.mdl", "w");
  fprintf (f, "# This source model file was created by source_test\n");
  fprintf (f, "# shell number, Av [mag], n(H) [cm^-3], Tgas [K], Tdust [K]\n");
  fprintf (f, "0	 0.1	1e+02	15.0	12.0\n");
  fprintf (f, "1	 1.0	1e+03	11.0	10.0\n");
  fprintf (f, "2	10.0	1e+04	 8.0	 7.0\n");
  fclose (f);

  /* Read it */

  read_source ("source.mdl", &source_mdl, verbose);
  
  /* Check that the values are correct */
  if ((source_mdl.n_shells == 3) &&
      (source_mdl.shell[0].av == 0.1) &&
      (source_mdl.shell[0].nh == 1e2) &&
      (source_mdl.shell[0].tgas == 15) &&
      (source_mdl.shell[0].tdust == 12) &&
      (source_mdl.shell[1].av == 1.0) &&
      (source_mdl.shell[1].nh == 1e3) &&
      (source_mdl.shell[1].tgas == 11) &&
      (source_mdl.shell[1].tdust == 10) &&
      (source_mdl.shell[2].av == 10.0) &&
      (source_mdl.shell[2].nh == 1e4) &&
      (source_mdl.shell[2].tgas == 8) &&
      (source_mdl.shell[2].tdust == 7))
    {
      free_mdl(&source_mdl);
      return EXIT_SUCCESS;
    }

  else
  {
    free_mdl(&source_mdl);
    return EXIT_FAILURE;
  }
}

