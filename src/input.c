/* 
   input.c - Read the input files needed by Astrochem.

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
#include <errno.h>
#include "astrochem.h"

/*
  Read the input file containing the parameters needed by the code:
  name of the chemistry network and source file names, physical 
  parameters, solver parameters, and initial abundances.
*/ 
  
void
read_input (const char *input_file, struct inp *input_params, int verbose)
{
  FILE *f;
  char line[MAX_LINE];
  char keyword[MAX_LINE];
  char parameter[MAX_LINE];
  char value[MAX_LINE];
  int  line_number = 0;
  int  i = 0, j = 0;
  errno = 0;
  
  /* Open the input file or exit if we can't open it. */

  if (verbose == 1)
    fprintf (stdout, "Reading input from %s.\n", input_file);
  f = fopen (input_file, "r" );
  if (!f)
    {
      fprintf (stderr, "astrochem: error: Can't open %s: %s\n", input_file,
	       strerror (errno));
      exit (EXIT_FAILURE);
    }

  /* Set the default values for parameters in the input file, in case
     the user didn't specify them. */
  
  strcpy (input_params->files.source_file, "");
  strcpy (input_params->files.chem_file, "");
  strcpy (&input_params->output.suffix, "");
  input_params->phys.chi = CHI_DEFAULT;
  input_params->phys.cosmic = COSMIC_DEFAULT;
  input_params->phys.grain_size = GRAIN_SIZE_DEFAULT;
  input_params->phys.grain_abundance = 0;
  input_params->solver.ti = TI_DEFAULT;
  input_params->solver.tf = TF_DEFAULT;
  input_params->solver.abs_err = ABS_ERR_DEFAULT;
  input_params->solver.rel_err = REL_ERR_DEFAULT;
  input_params->output.time_steps = TIME_STEPS_DEFAULT;
  input_params->output.trace_routes = TRACE_ROUTES_DEFAULT;

  /* Loop over the lines, and look for keywords (between brackets) and
     parameters/values (separated by "="). */

  while (fgets (line, MAX_LINE, f) != NULL)
    {
      line_number++;
      if (line[0] == '#') continue; /* Skip comments */
      if (sscanf (line, "[ %512[a-zA-Z] ]", keyword ) == 1)
	;
      else if (sscanf (line, "%s = %s", parameter, value ) == 2) {
	
	/* Source and chemical network files */
	
	if (strcmp (keyword, "files") == 0)
	  {
	    if (strcmp (parameter, "source") == 0)
	      strncpy (input_params->files.source_file, value, MAX_LINE);
	    else if (strcmp (parameter, "chem") == 0)
	      strncpy (input_params->files.chem_file, value, MAX_LINE);
	    else
	      input_error (input_file, line_number); /* Unknown parameter */
	  }

	/* Physical parameters */

	else if (strcmp (keyword, "phys") == 0)
	  {
	    if (strcmp (parameter, "chi") == 0)
	      input_params->phys.chi = atof (value);
	    else if (strcmp (parameter, "cosmic") == 0)
	      input_params->phys.cosmic = atof (value);
	    else if (strcmp (parameter, "grain_size") == 0)
	      input_params->phys.grain_size = atof (value) * 1e-4;   /* microns -> cm */
	    else
	      input_error (input_file, line_number);
	  }

	/* Solver parameters. Times are converted from
	   year to seconds */

	else if (strcmp (keyword, "solver") == 0)
	  {
	    if (strcmp (parameter, "ti") == 0)
	      input_params->solver.ti = atof (value) * CONST_MKSA_YEAR;
	    else if (strcmp (parameter, "tf") == 0)
	      input_params->solver.tf = atof (value) * CONST_MKSA_YEAR;
	    else if (strcmp (parameter, "abs_err") == 0)
	      input_params->solver.abs_err = atof (value);
	    else if (strcmp (parameter, "rel_err") == 0)
	      input_params->solver.rel_err = atof (value);
	    else
	      input_error (input_file, line_number);
	  }

	/* Initial abundances */
	
	else if (strcmp (keyword, "abundances") == 0)
	  {
	    if (i < MAX_INITIAL_ABUNDANCES)
	      {
		strncpy (input_params->abundances.initial_abundances[i].specie, parameter,
			 sizeof (input_params->abundances.initial_abundances[i].specie));
		input_params->abundances.initial_abundances[i].abundance = atof (value);

		/* Compute the total grain density */
		
		if ((strncmp (input_params->abundances.initial_abundances[i].specie, "grain",
			      sizeof (input_params->abundances.initial_abundances[i].specie)) == 0) ||
		    (strncmp (input_params->abundances.initial_abundances[i].specie, "grain(-)",
			      sizeof (input_params->abundances.initial_abundances[i].specie)) == 0) ||
		    (strncmp (input_params->abundances.initial_abundances[i].specie, "grain(+)",
			      sizeof (input_params->abundances.initial_abundances[i].specie)) == 0))
		  input_params->phys.grain_abundance += input_params->abundances.initial_abundances[i].abundance;
		
		i++;
	      } 
	    else 
	      {
		fprintf (stderr, "astrochem: error: the number of species "
			 "in %s exceed %i.\n", input_file, MAX_INITIAL_ABUNDANCES);
		exit (1);
	      }
	  }

	/* Output */
	
	else if (strcmp (keyword, "output") == 0)
	  {
	    if (strcmp (parameter, "time_steps") == 0)
		input_params->output.time_steps = atoi (value);
	    else if (strcmp (parameter, "abundances") == 0)
	      
	      /* Loop over the species (separated by a comma), and
		 copy them in the ouput_species array */
	      {
		const char delimiter[] = ",";
		char *output_specie;

		/* Structure initialization */

		int temp_i;
		for (temp_i = 0; temp_i < MAX_OUTPUT_ABUNDANCES; temp_i++)
		  {
		    input_params->output.output_species[temp_i]=NULL;
		  }

		if (j >= MAX_OUTPUT_ABUNDANCES)
		  {
		    fprintf (stderr, "astrochem: error: the number of species in output exceeds %i.\n",
			     MAX_OUTPUT_ABUNDANCES);
		    exit(1);
		  }
		if ((input_params->output.output_species[j] = malloc (sizeof (char) * MAX_CHAR_SPECIES)) == NULL)
		  {
		    fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
			     __FILE__, __LINE__); 
		    exit(1);
		  }
		output_specie = strtok (value, delimiter);
		strncpy (input_params->output.output_species[j], output_specie,
			 sizeof (char) * MAX_CHAR_SPECIES); 
		j++;
		while ((output_specie = strtok (NULL, delimiter)) != NULL)
		  {
		    if (j >= MAX_OUTPUT_ABUNDANCES)
		      {
			fprintf (stderr, "astrochem: error: the number of species in output exceeds %i.\n", 
				 MAX_OUTPUT_ABUNDANCES);
			exit(1);
		      }
		    if ((input_params->output.output_species[j] = malloc (sizeof (char) * MAX_CHAR_SPECIES)) == NULL)
		      {
			fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
				 __FILE__, __LINE__); 
			exit(1);
		      }
		    strncpy (input_params->output.output_species[j], output_specie,
			     sizeof (char) * MAX_CHAR_SPECIES);
		    j++;
		  }
	      }
	    else if (strcmp (parameter, "trace_routes") == 0)
	      input_params->output.trace_routes = atoi (value);
	    else if (strcmp (parameter, "suffix") == 0)
	      strncpy (&input_params->output.suffix, value, MAX_LINE);
	    else
	      input_error (input_file, line_number);
	  }
	
	/* Unknown or unspecified keyword */
     
	else
	  input_error (input_file, line_number-1); /* Unknown or unspecified keyword */
      } 
      
      /* Error while reading a parameter/values */

      else 
	input_error (input_file, line_number);
    }

  /* Close the file */
  
  fclose (f);
  
  input_params->abundances.n_initial_abundances = i;
  input_params->output.n_output_species = j;
  
  /* Check that the source file name and the chemical file name were specified
     in the input file. Also check that other parameters have acceptable 
     values. */

  if (strcmp (input_params->files.chem_file, "") == 0)
    {
      fprintf (stderr, "astrochem: error: no chemical network file specified in %s.\n",
	       input_file);
      exit (1);
    }
  if (strcmp (input_params->files.source_file, "") == 0)
    {
      fprintf (stderr, "astrochem: error: no source model file specified in %s.\n",
	       input_file);
      exit (1);
    }
  if (input_params->solver.ti <= 0.)
    {
      fprintf (stderr, "astrochem: warning: incorrect ti value specified in %s."
	       "Assuming default value.\n", input_file);
      input_params->solver.ti = TI_DEFAULT;
    }
  if (input_params->solver.tf < input_params->solver.ti)
    {
      fprintf (stderr, "astrochem: warning: incorrect tf value specified in %s."
	       "Assuming default value.\n", input_file);
      input_params->solver.tf = TF_DEFAULT;
    }
}

/*
  Display an error message and exit of an error is encountered while
  reading the input file.
*/

void 
input_error (const char *input_file, int line_number)
{
  fprintf (stderr, "astrochem: error: incorrect input in %s line %i.\n", 
	   input_file, line_number);
  exit (1);
}

/*
  Read the file containing the source model.
*/

void 
read_source (const char *source_file, struct mdl *source_mdl,const int verbose)
{
  FILE *f;
  char line[MAX_LINE];
  int line_number = 0;
  int i = 0;
  int col1;
  double col2, col3, col4, col5;

  /* Open the input file or exit if we can't open it */

  if (verbose == 1)
    fprintf (stdout, "Reading source model from %s.\n", source_file);
  f = fopen (source_file, "r" );
  if (!f)
    {
      fprintf (stderr, "astrochem: error: can't open %s.\n", source_file);
      exit (1);
    }
    
  /* Loop over the lines, and read the shell number, visual
     extinction, density, gas and dust temperature. */
  
  while (fgets (line, MAX_LINE, f) != NULL)
    {
      line_number++;
      if (line[0] == '#') continue; /* Skip comments */
      if (sscanf (line, "%d %lf %lf %lf %lf", 
		  &col1, &col2, &col3, &col4, &col5) == 5)
	{
	  if (i <= MAX_SHELLS - 1)
	    {
	      source_mdl->shell[i].av = col2;
	      source_mdl->shell[i].nh = col3;
	      source_mdl->shell[i].tgas = col4;
	      source_mdl->shell[i].tdust = col5;
	      i++;
	    }
	  else
	    {
	      fprintf (stderr, "astrochem: error: the number of shells in %s exceed %i.\n", 
		      source_file, MAX_SHELLS);
	      exit (1);
	    }
	} 
      else
	{
	  input_error (source_file, line_number);
	}
    }

  /* Close the file */
  
  fclose (f);

  /* Check the model has at least one shell */

  source_mdl->n_shells = i;
  if ( source_mdl->n_shells == 0 )
    {
      fprintf (stderr, "astrochem: error: no valid lines found in %s.\n",
	       source_file);
      exit (1);
    }
}

/* 
   Check that the initial_abundance and output species structure do
   not contain any specie that is not in the network. If they do, we
   simply display a warning message, since those will be ignored in
   solve() and output(). 
*/

void
check_species (struct abund initial_abundances[],
	       int n_initial_abundances, char *output_species[],
	       int n_output_species, char *species[], int n_species)
{
  int i;
  
  /* Check initial abundances */
  
  for (i=0; i < n_initial_abundances; i++)
    {
      if (specie_index (initial_abundances[i].specie, (const char* const*) species, n_species) == -2)
	  fprintf (stderr, "astrochem: warning: %s initial abundance given, "
		   "but is not in the network.\n", initial_abundances[i].specie);
    }

  /* Check output species */
  
  for (i=0; i < n_output_species; i++)
    {
      if (specie_index (output_species [i], (const char *const *) species, n_species) == -2)
	fprintf (stderr, "astrochem: warning: %s abundance requested, "
		 "but is not in the network.\n", output_species [i]);
    }
}

/* 
   Free the input structure.
*/

void
free_input_struct (struct inp * input_params)
{
  int i;

  for (i=0; i < MAX_OUTPUT_ABUNDANCES; i++)
    {
      if (input_params->output.output_species[i] != NULL)
        {
	  free (input_params->output.output_species[i]);
        }
    }
}
