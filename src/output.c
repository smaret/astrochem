/*
   output.c - Write the abundances and formation/destruction routes in
   output files.

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "astrochem.h"

#define MAX_CHAR_FILENAME 64

void
output (int n_shells, const inp_t * input_params, const mdl_t * source_mdl,
	const net_t * network, const res_t * results, int verbose)
{
  FILE *f;
  char filename[MAX_CHAR_FILENAME];
  int i, j, k, l;

  /* Write abundances in output files */

  if (verbose == 1)
    fprintf (stdout, "Writing abundances in output files... ");

  for (i = 0; i < input_params->output.n_output_species; i++)
    {
      /* Open the file or exit if an error occurs. */

      strncpy (filename,
	       network->species_names[input_params->
				      output.output_species_idx[i]],
	       MAX_CHAR_FILENAME);
      if (strlen (input_params->output.suffix) != 0)
	{
	  strncat (filename, "_", MAX_CHAR_FILENAME - strlen (filename) - 1);
	  strncat (filename, input_params->output.suffix,
		   MAX_CHAR_FILENAME - strlen (filename) - 1);
	}
      strncat (filename, ".abun", MAX_CHAR_FILENAME - strlen (filename) - 1
	       - strlen (input_params->output.suffix));
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

      fprintf (f, "# %s abundance computed by astrochem\n",
	       input_params->output.output_species[i]);
      fprintf (f, "# shell number / time [yr]\n");
      fprintf (f, "#\n");

      /* Write the time */

      fprintf (f, "   ");
      for (j = 0; j < source_mdl->n_time_steps; j++)
	fprintf (f, "  %8.2e", source_mdl->time_steps[j] / CONST_MKSA_YEAR);
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

      fprintf (f, "# %s abundance computed by astrochem\n",
	       network->species_names[input_params->
				      output.output_species_idx[i]]);
      fprintf (f, "# time [yr] / shell number\n");
      fprintf (f, "#\n");

      /* Write the shell number */

      fprintf (f, "        ");
      for (k = 0; k < n_shells; k++)
	fprintf (f, "  %8d", k);
      fprintf (f, "\n");

      /* Write the abundance as a function of time for each shell. */

      for (j = 0; j < input_params->output.time_steps; j++)
	{
	  fprintf (f, "%8.2e",
		   source_mdl->ts.time_steps[j] / CONST_MKSA_YEAR);
	  for (k = 0; k < n_shells; k++)
	    fprintf (f, "  %8.2e",
		     results->abundances[get_abundance_idx
					 (results, k, j, i)]);
	  fprintf (f, "\n");
	}

#endif

      fclose (f);
    }

  if (verbose == 1)
    fprintf (stdout, "done.\n");

  /* Write formation/destruction routes in output files */

  if (input_params->output.trace_routes == 1)
    {
      if (verbose == 1)
	fprintf (stdout,
		 "Writing formation/destruction routes in output files... ");

      for (i = 0; i < input_params->output.n_output_species; i++)
	{
	  /* Open the file or exit if an error occurs. */
	  strncpy (filename,
		   network->species_names[input_params->
					  output.output_species_idx[i]],
		   MAX_CHAR_FILENAME);
	  if (strlen (input_params->output.suffix) != 0)
	    {
	      strncat (filename, "_",
		       MAX_CHAR_FILENAME - strlen (filename) - 1);
	      strncat (filename, input_params->output.suffix,
		       MAX_CHAR_FILENAME - strlen (filename) - 1);
	    }
	  strncat (filename, ".rout",
		   MAX_CHAR_FILENAME - strlen (filename) - 1 -
		   strlen (input_params->output.suffix));
	  if ((f = fopen (filename, "w")) == NULL)
	    {
	      fprintf (stderr, "astrochem: error: can't open %s\n", filename);
	      exit (1);
	    }

	  /* Write the main formation routes as a function of time
	     and for each shell. For each formation/destruction
	     route, we write the reaction number and the reaction
	     rate. */

	  /* Write the header */

	  fprintf (f,
		   "# Main %s formation/destruction routes computed by astrochem\n",
		   network->species_names[input_params->
					  output.output_species_idx[i]]);
	  fprintf (f,
		   "# shell number  time [yr]  reaction number 1  reaction rate 1 [cm-3/s]... \n");
	  fprintf (f, "#\n");

	  /* Write the formation/destruction routes */

	  for (k = 0; k < n_shells; k++)
	    {
	      for (j = 0; j < input_params->output.time_steps; j++)
		{
		  fprintf (f, "%4i", k);
		  fprintf (f, "  %9.2e",
			   source_mdl->ts.time_steps[j] / CONST_MKSA_YEAR);
		  for (l = 0; l < N_OUTPUT_ROUTES; l++)
		    {
		      fprintf (f, "  %4i",
			       results->routes[get_route_idx
					       (results, k, j, i,
						l)].formation.reaction_no);
		      fprintf (f, " %9.2e",
			       results->routes[get_route_idx
					       (results, k, j, i,
						l)].formation.rate);
		    }
		  fprintf (f, "\n");
		  fprintf (f, "%4i", k);
		  fprintf (f, "  %9.2e",
			   source_mdl->ts.time_steps[j] / CONST_MKSA_YEAR);
		  for (l = 0; l < N_OUTPUT_ROUTES; l++)
		    {
		      fprintf (f, "  %4i",
			       results->routes[get_route_idx
					       (results, k, j, i,
						l)].destruction.reaction_no);
		      fprintf (f, " %9.2e",
			       -results->routes[get_route_idx
						(results, k, j, i,
						 l)].destruction.rate);
		    }
		  fprintf (f, "\n");
		}
	    }
	}

      if (verbose == 1)
	fprintf (stdout, "done.\n");
    }
}
