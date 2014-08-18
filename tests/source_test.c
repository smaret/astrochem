/* 
   source_test.c - Test the read_source() function
   
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

  mdl_t source_mdl;
  inp_t fake;
  fake.output.time_steps = 128;
  fake.solver.ti = 0.000001  * CONST_MKSA_YEAR;
  fake.solver.tf = 10000000  * CONST_MKSA_YEAR;
  int verbose = 0;

  /* Create the input.ini file */

  f = fopen ("source.mdl", "w");
  fprintf (f, "# This source model file was created by source_test\n");
  fprintf (f, "# cell number, Av [mag], n(H) [cm^-3], Tgas [K], Tdust [K]\n");
  fprintf (f, "0	 0.1	1e+02	15.0	12.0\n");
  fprintf (f, "1	 1.0	1e+03	11.0	10.0\n");
  fprintf (f, "2	10.0	1e+04	 8.0	 7.0\n");
  fclose (f);

  /* Read it */

  read_source ("source.mdl", &source_mdl, &fake, verbose);

  /* Check that the values are correct */
  if ((source_mdl.n_cells == 3) && 
      (source_mdl.cell[0].av[0] == 0.1) &&
      (source_mdl.cell[0].nh[0] == 1e2) && 
      (source_mdl.cell[0].tgas[0] == 15) &&
      (source_mdl.cell[0].tdust[0] == 12) &&
      (source_mdl.cell[1].av[0] == 1.0) &&
      (source_mdl.cell[1].nh[0] == 1e3) &&
      (source_mdl.cell[1].tgas[0] == 11) &&
      (source_mdl.cell[1].tdust[0] == 10) &&
      (source_mdl.cell[2].av[0] == 10.0) &&
      (source_mdl.cell[2].nh[0] == 1e4) &&
      (source_mdl.cell[2].tgas[0] == 8) &&
      (source_mdl.cell[2].tdust[0] == 7) &&
      (source_mdl.ts.time_steps[10] - 332.988055 < 0.0001) && 
      (source_mdl.ts.time_steps[23] - 7130.784562 < 0.0001) &&
      (source_mdl.ts.time_steps[47] - 2040939.960351 < 0.0001) )
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

