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

char chem_file[MAX_LINE];
char source_file[MAX_LINE];
char suffix[MAX_LINE];

double chi;
double cosmic;
double grain_size;
double grain_abundance;
double ti;
double tf;
double abs_err;
double rel_err;
struct abund initial_abundances[MAX_INITIAL_ABUNDANCES];
int n_initial_abundances;
char *output_species[MAX_OUTPUT_ABUNDANCES];
int n_output_species;
int time_steps;
int trace_routes;

int n_shells;
int shell[MAX_SHELLS];
double av[MAX_SHELLS];
double nh[MAX_SHELLS];
double tgas[MAX_SHELLS];
double tdust[MAX_SHELLS];
int shell_index;

struct react reactions[MAX_REACTIONS];
char *species[MAX_SPECIES];
int n_reactions;
int n_species;

double abundances[MAX_SHELLS][MAX_TIME_STEPS][MAX_OUTPUT_ABUNDANCES];
struct rout routes[MAX_SHELLS][MAX_TIME_STEPS][MAX_OUTPUT_ABUNDANCES][N_OUTPUT_ROUTES];
double tim[MAX_TIME_STEPS];

void usage (void);
void version (void);

int
main (int argc, char *argv[])
{  
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

  read_input (input_file, chem_file, source_file, &chi, &cosmic,
	      &grain_size, &grain_abundance, &ti, &tf, &abs_err, &rel_err,
	      initial_abundances, &n_initial_abundances,
	      output_species, &n_output_species, &time_steps,
	      &trace_routes, suffix, verbose);
  
  /* Read the source model file */

  read_source (source_file, shell, &n_shells, av, nh,
	       tgas, tdust, verbose);

  /* Read the chemical network file */

  read_network (chem_file, reactions, &n_reactions, 
		species, &n_species, verbose);

  /* Check that the initial_abundance and output_species structure do
     not contain any specie that is not in the network. */

  check_species (initial_abundances, n_initial_abundances,
		 output_species, n_output_species, species, 
		 n_species);

  /* Build the vector of time */

  {
    int i;
    
    for (i = 0; i < time_steps; i++)
      {
	if (i < MAX_TIME_STEPS)
	    tim[i] = pow (10., log10 (ti) + i * (log10 (tf) - log10(ti)) 
			   / (time_steps - 1));
	else
	  {
	    fprintf (stderr, "astrochem: error: the number of time" 
		     "steps in %s exceed %i.\n", input_file, MAX_TIME_STEPS);
	    exit (1);
	  }
      }
  }

  /* Solve the ODE system for each shell. */

#ifdef HAVE_OPENMP  
#pragma omp parallel shared (abundances) private (shell_index)
#endif

  {

#ifdef HAVE_OPENMP
#pragma omp for schedule (dynamic, 1) nowait
#endif

    for (shell_index = 0; shell_index < n_shells; shell_index++)
      {
	if (verbose >= 1)
	  fprintf (stdout, "Computing abundances in shell %d...\n", shell_index);

	  solve (chi, cosmic, grain_size, grain_abundance,
	       abs_err, rel_err, initial_abundances,
	       n_initial_abundances, output_species,
	       n_output_species, av[shell_index],
	       nh[shell_index], tgas[shell_index],
	       tdust[shell_index], reactions, n_reactions,
	       species, n_species, shell_index, tim,
	       time_steps, abundances, trace_routes,
	       routes, verbose);

	if (verbose >= 1)
	  fprintf (stdout, "Done with shell %d.\n", shell_index);
      }
  }

  /* Write the abundances in output files */

  output (n_shells, tim, time_steps, output_species, n_output_species,
	  abundances, species, n_species, trace_routes, routes, suffix,
	  verbose);

  exit (0);
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
