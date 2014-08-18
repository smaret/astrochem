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
  int verbose = 0;

  /* Create the input.ini file */

  f = fopen ("source_dyn.mdl", "w");
  fprintf (f, "# Source model file example for a time-dependant source structure\n"
      "[times]\n"
      "   0    1.00e-06\n"
      "   1    1.24e-06\n"
      "   2    1.55e-06\n"
      "   3    1.92e-06\n"
      "   4    2.39e-06\n"
      "   5    2.97e-06\n"
      "   6    3.69e-06\n"
      "   7    4.59e-06\n"
      "   8    5.70e-06\n"
      "   9    7.09e-06\n"
      "  10    8.81e-06\n"
      "[cells]\n"
      "# cell (= shell) number, time index, Av [mag], nH [cm^-3], Tgas [K], Tdust [K]\n"
      "   0      0    0.00  1.00e+04   10.00   10.00\n"
      "   0      1    0.16  1.11e+04   11.50   11.50\n"
      "   0      2    0.31  1.24e+04   12.99   12.99\n"
      "   0      3    0.47  1.39e+04   14.49   14.49\n"
      "   0      4    0.63  1.55e+04   15.98   15.98\n"
      "   0      5    0.79  1.72e+04   17.48   17.48\n"
      "   0      6    0.94  1.92e+04   18.98   18.98\n"
      "   0      7    1.10  2.14e+04   20.47   20.47\n"
      "   0      8    1.26  2.39e+04   21.97   21.97\n"
      "   0      9    1.42  2.66e+04   23.46   23.46\n"
      "   0     10    1.57  2.97e+04   24.96   24.96\n"
      "   1      0    0.00  1.00e+04   10.00   10.00\n"
      "   1      1    0.16  1.11e+04   11.50   11.50\n"
      "   1      2    0.31  1.24e+04   12.99   12.99\n"
      "   1      3    0.47  1.39e+04   14.49   14.49\n"
      "   1      4    0.63  1.55e+04   15.98   15.98\n"
      "   1      5    0.79  1.72e+04   17.48   17.48\n"
      "   1      6    0.94  1.92e+04   18.98   18.98\n"
      "   1      7    1.10  2.14e+04   20.47   20.47\n"
      "   1      8    1.26  2.39e+04   21.97   21.97\n"
      "   1      9    1.42  2.66e+04   23.46   23.46\n"
      "   1     10    1.57  2.97e+04   24.96   24.96\n"
      "   2      0    0.00  1.00e+04   10.00   10.00\n"
      "   2      1    0.16  1.11e+04   11.50   11.50\n"
      "   2      2    0.31  1.24e+04   12.99   12.99\n"
      "   2      3    0.47  1.39e+04   14.49   14.49\n"
      "   2      4    0.63  1.55e+04   15.98   15.98\n"
      "   2      5    0.79  1.72e+04   17.48   17.48\n"
      "   2      6    0.94  1.92e+04   18.98   18.98\n"
      "   2      7    1.10  2.14e+04   20.47   20.47\n"
      "   2      8    1.26  2.39e+04   21.97   21.97\n"
      "   2      9    1.42  2.66e+04   23.46   23.46\n"
      "   2     10    1.57  2.97e+04   24.96   24.96\n");
  fclose (f);

  /* Read it */

  read_source ("./source_dyn.mdl", &source_mdl, &fake, verbose);
  /* Check that the values are correct */
  if ((source_mdl.n_cells == 3) &&
      (source_mdl.ts.n_time_steps == 11) &&
      (source_mdl.cell[0].av[0] == 0.0) &&
      (source_mdl.cell[0].nh[0] == 1e4) && 
      (source_mdl.cell[0].tgas[0] == 10) &&
      (source_mdl.cell[0].tdust[0] == 10) &&
      (source_mdl.cell[1].av[5] == 0.79) &&
      (source_mdl.cell[1].nh[5] == 1.72e4) &&
      (source_mdl.cell[1].tgas[5] == 17.48) &&
      (source_mdl.cell[1].tdust[5] == 17.48) &&
      (source_mdl.cell[2].av[10] == 1.57) &&
      (source_mdl.cell[2].nh[10] == 2.97e4) &&
      (source_mdl.cell[2].tgas[10] == 24.96) &&
      (source_mdl.cell[2].tdust[10] == 24.96) &&
      (source_mdl.ts.time_steps[5] / CONST_MKSA_YEAR - 2.97e-6 < 0.01e-6 ))
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

