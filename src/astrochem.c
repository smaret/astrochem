/* 
   Astrochem - compute the abundances of chemical species in the
   interstellar medium as as function of time.
   
   Copyright (c) 2006-2013 Sebastien Maret

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "astrochem.h"


void usage (void);
void version (void);

int
main (int argc, char *argv[])
{
  inp_t input_params;
  mdl_t source_mdl;
  net_t network;
  res_t results;
  int cell_index;
  
  int verbose = 1;
  char *input_file;

  /* Parse options and command line arguments. Diplay help 
     message if no (or more than one) argument is given. */

  {
    int opt;
    
    static struct option longopts[] = {
      {"help",    no_argument, NULL, 'h'},
      {"version", no_argument, NULL, 'V'},
      {"verbose", no_argument, NULL, 'v'},
      {"quiet",   no_argument, NULL, 'q'},
      {0, 0, 0, 0}
    };
    
    while ((opt = getopt_long(argc, argv, "hVvq", longopts, NULL)) != -1)
      {
	switch (opt)
	  {
	  case 'h':
	    usage ();
	    exit (0);
	    break;
	  case 'V':
	    version ();
	    exit (0);
	    break;
	  case 'v':
	    verbose = 2;
	    break;
	  case 'q':
	    verbose = 0;
	    break;
	  default:
	    usage ();
	    exit (1);
	  }
      };
    argc -= optind;
    argv += optind;
    if (argc != 1) 
      {
	usage ();
	exit (1);
      }
    input_file = argv[0];
  }

 /* Read the input file */
  read_input_file_names (input_file, &input_params.files, verbose);
  
  /* Read the chemical network file */
  read_network (input_params.files.chem_file, &network, verbose);
   
  /* Read the input file */
  read_input (input_file, &input_params, &network , verbose);
  
  /* Read the source model file */
  read_source (input_params.files.source_file, &source_mdl, verbose);

  /* Check that the initial_abundance and output_species structure do
     not contain any specie that is not in the network. */

  /*check_species (input_params.abundances.initial_abundances, input_params.abundances.n_initial_abundances,
		 input_params.output.output_species, input_params.output.n_output_species, network.species, 
		 network.n_species);*/

  /* Allocate results */
   alloc_results( &results, input_params.output.time_steps, source_mdl.n_cells, input_params.output.n_output_species);


  /* Build the vector of time */

  {
    int i;
    
    for (i = 0; i <  input_params.output.time_steps; i++)
      {
	if (i < MAX_TIME_STEPS)
	    results.tim[i] = pow (10., log10 ( input_params.solver.ti) + i 
				   * (log10 (input_params.solver.tf) - log10 (input_params.solver.ti)) 
				   / (input_params.output.time_steps - 1));
	else
	  {
	    fprintf (stderr, "astrochem: error: the number of time" 
		     "steps in %s exceed %i.\n", input_file, MAX_TIME_STEPS);
        free_input (&input_params);
        free_mdl (&source_mdl);
        free_network (&network);
        free_results (&results);
	    return(EXIT_FAILURE);
	  }
      }
  }

  /* Solve the ODE system for each cell. */

#ifdef HAVE_OPENMP  
#pragma omp parallel shared (abundances) private (cell_index)
#endif

  {

#ifdef HAVE_OPENMP
#pragma omp for schedule (dynamic, 1) nowait
#endif

    for (cell_index = 0; cell_index < source_mdl.n_cells; cell_index++)
      {
	if (verbose >= 1)
	  fprintf (stdout, "Computing abundances in cell %d...\n", cell_index);
	solve (cell_index, &input_params, &source_mdl.cell[cell_index], &network, &results, verbose);
	if (verbose >= 1)
	  fprintf (stdout, "Done with cell %d.\n", cell_index);
      }
  }

  /* Write the abundances in output files */

  output (source_mdl.n_cells, &input_params, &network, &results, verbose);
  free_input (&input_params);
  free_mdl (&source_mdl);
  free_network (&network);
  free_results (&results);
  return (EXIT_SUCCESS);
}

/*
  Display help message.
*/

void
usage (void)
{
  fprintf (stdout, "Usage: astrochem [option...] [file]\n\n");
  fprintf (stdout, "Options:\n");
  fprintf (stdout, "   -h, --help         Display this help\n");
  fprintf (stdout, "   -V, --version      Print program version\n");
  fprintf (stdout, "   -v, --verbose      Verbose mode\n");
  fprintf (stdout, "   -q, --quiet        Suppress all messages\n");
  fprintf (stdout, "\n");
  fprintf (stdout, "See the astrochem(1) manual page for more information.\n");
  fprintf (stdout, "Report bugs to <%s>.\n", PACKAGE_BUGREPORT);
}

/*
  Display version.
*/

void
version (void)
{
  fprintf (stdout, "This is astrochem, version %s\n", PACKAGE_VERSION);
#ifdef HAVE_OPENMP
  fprintf (stdout, "OpenMP support enabled, ");
#else
  fprintf (stdout, "OpenMP support disabled, ");
#endif
#ifdef USE_LAPACK
  fprintf (stdout, "LAPACK support enabled.\n");
#else
  fprintf (stdout, "LAPACK support disabled.\n");
#endif
  fprintf (stdout, "Copyright (c) 2006-2013 Sebastien Maret\n");
  fprintf (stdout, "\n");
  fprintf (stdout, "This is free software. You may redistribute copies of it under the terms\n");
  fprintf (stdout, "of the GNU General Public License. There is NO WARRANTY, to the extent\n");
  fprintf (stdout, "permitted by law.\n");
}
