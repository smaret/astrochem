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

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "libastrochem.h"


void usage (void);
void version (void);

int
main (int argc, char *argv[])
{
  int verbose = 1;

  /* Parse options and command line arguments. Diplay help
     message if no (or more than one) argument is given. */

    {
      int opt;

      static struct option longopts[] = {
          {"help", no_argument, NULL, 'h'},
          {"version", no_argument, NULL, 'V'},
          {"verbose", no_argument, NULL, 'v'},
          {"quiet", no_argument, NULL, 'q'},
          {0, 0, 0, 0}
      };

      while ((opt = getopt_long (argc, argv, "hVvq", longopts, NULL)) != -1)
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
      if (argc != 0)
        {
          usage ();
          return 1;
        }
    }

  char *chem_file = "../networks/osu2009.chm";
  net_t network;
  read_network(chem_file, &network, verbose );

  phys_t phys;
  phys.cosmic = 1e-17;
  phys.chi = 0;
  phys.grain_size = GRAIN_SIZE_DEFAULT;
  phys.grain_abundance = 0;

  double abs_err, rel_err;
  abs_err = ABS_ERR_DEFAULT;
  rel_err = REL_ERR_DEFAULT;

  const char* species[]  = {"CO", "HCO(+)", "e(-)"};
  const double initial_abundances[] = {1e-4, 1e-9, 1e-9};

  double *abundances;
  alloc_abundances( &network, &abundances ); // Allocate the abundances array; it contains all species.
  set_initial_abundances(species, 3, initial_abundances, &network, abundances); // Set initial abundances

  double density = 1000;
  double av = 20;
  double temperature = 10;

  cell_t cell;
  cell.nh = &density;
  cell.av = &av;
  cell.tgas = &temperature;
  cell.tdust = &temperature; // Assume tgas = tdust in this specific case

  astrochem_mem_t astrochem_mem;

  if( solver_init( &cell, &network, &phys, abundances , density, abs_err, rel_err, &astrochem_mem ) != 0 )
    {
      return EXIT_FAILURE;
    }
  int i;
  double time = 0;
  for( i = 0; i< 1000000 ; i++)
    {
      time += 1e-6; // advance time
      solve( &astrochem_mem, &network, abundances, time, verbose);

      /* Do something with the results of abundances computations */
    }
  solver_close( &astrochem_mem );
  free_abundances( abundances );
  free_network (&network);
  return (EXIT_SUCCESS);
}

/*
   Display help message.
   */

void
usage (void)
{
  fprintf (stdout, "Usage: apiuser [option...]\n\n");
  fprintf (stdout, "Options:\n");
  fprintf (stdout, "   -h, --help         Display this help\n");
  fprintf (stdout, "   -V, --version      Print program version\n");
  fprintf (stdout, "   -v, --verbose      Verbose mode\n");
  fprintf (stdout, "   -q, --quiet        Suppress all messages\n");
  fprintf (stdout, "\n");
  fprintf (stdout,
           "See the astrochem(1) manual page for more information.\n");
  fprintf (stdout, "Report bugs to <%s>.\n", PACKAGE_BUGREPORT);
}

/*
   Display version.
   */

void
version (void)
{
  fprintf (stdout, "This is astrochem, version %s\n", PACKAGE_VERSION);
#ifdef USE_LAPACK
  fprintf (stdout, "LAPACK support enabled.\n");
#else
  fprintf (stdout, "LAPACK support disabled.\n");
#endif
  fprintf (stdout, "Copyright (c) 2006-2013 Sebastien Maret\n");
  fprintf (stdout, "\n");
  fprintf (stdout,
           "This is free software. You may redistribute copies of it under the terms\n");
  fprintf (stdout,
           "of the GNU General Public License. There is NO WARRANTY, to the extent\n");
  fprintf (stdout, "permitted by law.\n");
}
