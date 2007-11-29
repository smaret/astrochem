/* 
   Read the input files needed by astrochem.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "astrochem.h"

/*
  Read the input file containing the parameters needed by the code:
  name of the chemistry network and source file names, physical 
  parameters, solver parameters, and initial abundances.
*/ 
  
void 
read_input (char *input_file, char *chem_file, char *source_file,
	    double *chi, double *pdyield, double *cosmic,
	    double *ti, double *tf, double *abs_err,
	    double *rel_err, struct abund initial_abundances[],
	    int *n_initial_abundances, char *output_species[],
	    int *n_output_species, int *time_steps, int *trace_routes,
	    char *suffix, int verbose)
{
  FILE *f;
  char line[MAX_LINE];
  char keyword[MAX_LINE];
  char parameter[MAX_LINE];
  char value[MAX_LINE];
  int  line_number = 0;
  int  i = 0, j = 0;
  
  /* Open the input file or exit if we can't open it. */

  if (verbose == 1)
    fprintf (stdout, "Reading input from %s.\n", input_file);
  f = fopen (input_file, "r" );
  if (!f)
    {
      fprintf (stderr, "astrochem: error: can't open %s.\n", input_file);
      exit (1);
    }

  /* Set the default values for parameters in the input file, in case
     the user didn't specify them. */
  
  strcpy(source_file, "");
  strcpy(chem_file, "");
  strcpy(suffix, "");
  *chi = CHI_DEFAULT;
  *pdyield = PDYIELD_DEFAULT;
  *cosmic = COSMIC_DEFAULT;
  *ti = TI_DEFAULT;
  *tf = TF_DEFAULT;
  *abs_err = ABS_ERR_DEFAULT;
  *rel_err = REL_ERR_DEFAULT;
  *time_steps = TIME_STEPS_DEFAULT;
  *trace_routes = TRACE_ROUTES_DEFAULT;

  /* Loop over the lines, and look for keywords (between brackets) and
     parameters/values (separated by "="). */

  while (fgets (line, MAX_LINE, f) != NULL)
    {
      line_number++;
      if (line[0] == '#') continue; /* Skip comments */
      if (sscanf( line, "[ %512[a-zA-Z] ]", keyword ) == 1)
	;
      else if (sscanf (line, "%s = %s", parameter, value ) == 2) {

	/* Source and chemical network files */
	
	if (strcmp (keyword, "files") == 0)
	  {
	    if (strcmp (parameter, "source") == 0)
	      strncpy (source_file, value, MAX_LINE);
	    else if (strcmp (parameter, "chem") == 0)
	      strncpy (chem_file, value, MAX_LINE);
	    else
	      input_error (input_file, line_number); /* Unknown parameter */
	  }

	/* Physical parameters */

	else if (strcmp (keyword, "phys") == 0)
	  {
	    if (strcmp (parameter, "chi") == 0)
	      *chi = atof (value);
	    else if (strcmp (parameter, "pdyield") == 0)
	      *pdyield = atof (value);
	    else if (strcmp (parameter, "cosmic") == 0)
	      *cosmic = atof (value);
	    else
	      input_error (input_file, line_number);
	  }

	/* Solver parameters. Times are converted from
	   year to seconds */

	else if (strcmp (keyword, "solver") == 0)
	  {
	    if (strcmp (parameter, "ti") == 0)
	      *ti = atof (value) * CONST_MKSA_YEAR;
	    else if (strcmp (parameter, "tf") == 0)
	      *tf = atof (value) * CONST_MKSA_YEAR;
	    else if (strcmp (parameter, "abs_err") == 0)
	      *abs_err = atof (value);
	    else if (strcmp (parameter, "rel_err") == 0)
	      *rel_err = atof (value);
	    else
	      input_error (input_file, line_number);
	  }

	/* Initial abundances */
	
	else if (strcmp (keyword, "abundances") == 0)
	  {
	    if (i < MAX_INITIAL_ABUNDANCES)
	      {
		strncpy (initial_abundances[i].specie, parameter,
			 sizeof (initial_abundances[i].specie));
		initial_abundances[i].abundance = atof (value);
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
		*time_steps = atoi (value);
	    else if (strcmp (parameter, "abundances") == 0)
	      
	      /* Loop over the species (separated by a comma), and
		 copy them in the ouput_species array */
	      {
		const char delimiter[] = ",";
		char *output_specie;
		
		if ((output_species[j] = malloc (sizeof (char) * MAX_CHAR_SPECIES)) == NULL)
		  {
		    fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
			     __FILE__, __LINE__); 
		    exit(1);
		  }
		output_specie = strtok (value, delimiter);
		strncpy (output_species[j], output_specie,
			 sizeof (char) * MAX_CHAR_SPECIES); 
		j++;
		while ((output_specie = strtok (NULL, delimiter)) != NULL)
		  {
		    if ((output_species[j] = malloc (sizeof (char) * MAX_CHAR_SPECIES)) == NULL)
		      {
			fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
				 __FILE__, __LINE__); 
			exit(1);
		      }
		    strncpy (output_species[j], output_specie,
			     sizeof (char) * MAX_CHAR_SPECIES);
		    j++;
		  }
	      }
	    else if (strcmp (parameter, "time_steps") == 0)
	      *trace_routes = atoi (value);
	    else if (strcmp (parameter, "suffix") == 0)
	      strncpy (suffix, value, MAX_LINE);
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
  
  *n_initial_abundances = i;
  *n_output_species = j;
  
  /* Check that the source file name and the chemical file name were specified
     in the input file. Also check that other parameters have acceptable 
     values. */

  if (strcmp (chem_file, "") == 0)
    {
      fprintf (stderr, "astrochem: error: no chemical network file specified in %s.\n",
	       input_file);
      exit (1);
    }
  if (strcmp (source_file, "") == 0)
    {
      fprintf (stderr, "astrochem: error: no source model file specified in %s.\n",
	       input_file);
      exit (1);
    }
  if (*ti <= 0.)
    {
      fprintf (stderr, "astrochem: warning: incorrect ti value specified in %s."
	       "Assuming default value.\n", input_file);
      *ti = TI_DEFAULT;
    }
  if (*tf < *ti)
    {
      fprintf (stderr, "astrochem: warning: incorrect tf value specified in %s."
	       "Assuming default value.\n", input_file);
      *tf = TF_DEFAULT;
    }
}

/*
  Display an error message and exit of an error is encountered while
  reading the input file.
*/

void 
input_error (char *input_file, int line_number)
{
  fprintf (stderr, "astrochem: error: incorrect input in %s line %i.\n", 
	  input_file, line_number);
  exit (1);
}
  
/*
  Read the file containing the source model.
*/

void 
read_source (char *source_file, int shell[], int *n_shells,
	     double av[], double nh[], double tgas[], double tdust[],
	     int verbose)
{
  FILE *f;
  char line[MAX_LINE];
  int line_number = 0;
  int i = 0;
  int col1;
  double col2, col3, col4, col5;

  /* Open the input file or exit if we can't open it */

  if (verbose > 0)
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
	      shell[i] = col1;
	      av[i] = col2;
	      nh[i] = col3;
	      tgas[i] = col4;
	      tdust[i] = col5;
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

  *n_shells = i;
  if ( *n_shells == 0 )
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
      if (specie_index (initial_abundances[i].specie, species, n_species) == -2)
	  fprintf (stderr, "astrochem: warning: %s initial abundance given, "
		   "but is not in the network.\n", initial_abundances[i].specie);
    }

  /* Check output species */
  
  for (i=0; i < n_output_species; i++)
    {
      if (specie_index (output_species [i], species, n_species) == -2)
	fprintf (stderr, "astrochem: warning: %s abundance requested, "
		 "but is not in the network.\n", output_species [i]);
    }
}
