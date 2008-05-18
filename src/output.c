/*
  output.c - Write the computed abundances in output files.

  Copyright (c) 2006-2008 Sebastien Maret
   
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "astrochem.h"

#define MAX_CHAR_FILENAME 64

void
output (int n_shells, double tim[], int time_steps,
	char *output_species[], int n_output_species,
	double abundances[MAX_SHELLS][MAX_TIME_STEPS][MAX_OUTPUT_ABUNDANCES],
	char *species[], int n_species, char *suffix, int verbose)
{
  FILE *f;
  char filename[MAX_CHAR_FILENAME];
  int i, j, k;

  if (verbose == 1)
    fprintf (stdout, "Writing abundances in output files... ");
  
  for (i = 0; i < n_output_species; i++)
    {
      /* Ignore species that are not in the network. */
      
      if (specie_index (output_species [i], species, n_species) != - 2)
	{

	  /* Open the file or exit if an error occurs. */

	  strncpy (filename, output_species[i], MAX_CHAR_FILENAME);
	  if (strlen (suffix) != 0)
	    {
	      strncat (filename, "_", MAX_CHAR_FILENAME - strlen (filename) - 1);
	      strncat (filename, suffix, MAX_CHAR_FILENAME - strlen (filename) - 1);
	    }
	  strncat (filename, ".abun", MAX_CHAR_FILENAME - strlen (filename) -1
		   - strlen (suffix));
	  if ((f = fopen (filename, "w")) == NULL)
	    {
	      fprintf (stderr, "astrochem: error: can't open %s\n", filename);
	      exit (1);
	    }

#ifdef OUTPUT_TIME_COLUMN

	  /* Write the abundances as in a function of time in different
	     columns, with a line for each shell. This is a good way to
	     write data if we have a lot of shells. However, it is
	     difficult to read if the source has only one (or a few)
	     shells */

	  /* Write the header. */
	  
	  fprintf (f, "# %s abundance computed by astrochem\n", output_species[i]);
	  fprintf (f, "# shell number / time [yr]\n");
	  fprintf (f, "#\n");
	  
	  /* Write the time */
	  
	  fprintf (f, "   ");
	  for (j = 0; j < time_steps; j++)
	    fprintf (f, "  %8.2e", tim[j] / CONST_MKSA_YEAR);
	  fprintf (f, "\n");
	  
	  /* Write the abundance as a function of time for each shell. */
	  
	  for (k = 0; k < n_shells; k++)
	    {
	      fprintf (f, "%3d", k);
	      for (j = 0; j < time_steps; j++)
		fprintf (f, "  %8.2e", abundances[k][j][i]);
	      fprintf (f, "\n");
	    }
	  
#else
	  
	  /* Write abundances as a function of time in different lines,
	     with a column for each shell. This is better if we have only
	     one (or a few) shells. */

	  /* Write the header. */
	  
	  fprintf (f, "# %s abundance computed by astrochem\n", output_species[i]);
	  fprintf (f, "# time [yr] / shell number\n");
	  fprintf (f, "#\n");
	  
	  /* Write the shell number */
	  
	  fprintf (f, "        ");
	  for (k = 0; k < n_shells; k++)
	    fprintf (f, "  %8d", k);
	  fprintf (f, "\n");
	  
	  /* Write the abundance as a function of time for each shell. */
	  
	  for (j = 0; j < time_steps; j++)
	    {
	      fprintf (f, "%8.2e", tim[j] / CONST_MKSA_YEAR);
	      for (k = 0; k < n_shells; k++)
		fprintf (f, "  %8.2e", abundances[k][j][i]);
	      fprintf (f, "\n");
	    } 
	  
#endif
     
	  fclose (f);
	}
    }
  if (verbose == 1)
    fprintf (stdout, "done.\n");
}
