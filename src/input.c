/*
   input.c - Read the input files needed by Astrochem.

   Copyright (c) 2006-2016 Sebastien Maret

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
#include <math.h>
#include <float.h>

#include "astrochem.h"
#include "input.h"

#include "network.h"

typedef enum
{
  R_STATIC = 0,  /**< Static Reading Mode, only cells >*/
  R_DYNAMIC = 1, /**< Dynamic Cells Reading Mode >*/
  R_TIMES = 2    /**< Dynamic Time Steps Reading Mode >*/
} SOURCE_READ_MODE;

int alloc_mdl (mdl_t * source_mdl, int n_cells, int n_time_steps);

int get_nb_active_line_section (const char *file, const char *section);

/*
  Read the input file containing the parameters needed by the code:
  name of the chemistry network and source file names, physical
  parameters, solver parameters, and initial abundances.
*/

int
read_input (const char *input_file, inp_t * input_params,
            const net_t * network, int verbose)
{
  FILE *f;
  char line[MAX_LINE];
  char keyword[MAX_LINE];
  char parameter[MAX_LINE];
  char value[MAX_LINE];
  int line_number = 0;
  int i = 0, j = 0;
  errno = 0;

  /* Open the input file or exit if we can't open it. */
  if (verbose == 1)
    fprintf (stdout, "Reading input from %s.\n", input_file);
  f = fopen (input_file, "r");
  if (!f)
    {
      fprintf (stderr, "astrochem: error: Can't open %s: %s\n", input_file,
               strerror (errno));
      return EXIT_FAILURE;
    }

  int n_output_species = 0;
  bool all_species = false;
  int n_initial_abundances =
   get_nb_active_line_section (input_file, "abundances");
  if( n_initial_abundances == -1 )
    {
      return EXIT_FAILURE;
    }
  while (fgets (line, MAX_LINE, f) != NULL)
    {
      if (strncmp (line, "abundances = ", 13) == 0)
        {
          char *localStr = &line[12];
          while (localStr != NULL)
            {
              localStr++;

              //Checking for all species
              if ((strncmp( localStr,"ALL", 3 ) == 0 ) ||
		  (strncmp( localStr,"all", 3 ) == 0))
                {
                  all_species = true;
                  n_output_species = network->n_species;
                  break;
                }
              localStr = strchr (localStr, ',');
              n_output_species++;
            }
        }
    }
  //Reset stream to beginning of file
  if (fseek (f, 0, SEEK_SET) != 0)
    {
      fprintf (stderr, "astrochem: error seeking begining of input file "
               "%s .\n", input_file);
      return EXIT_FAILURE;
    }

  /* Allocate a data structure for the input, and set the default
     values for parameters in the input file, in case the user doesn't
     specify them. */
  if( alloc_input (input_params, n_initial_abundances, n_output_species) != EXIT_SUCCESS )
    {
      return EXIT_FAILURE;
    }
  strcpy (input_params->files.source_file, "");
  strcpy (input_params->files.chem_file, "");
  strcpy (input_params->output.suffix, "");
  input_params->phys.chi = CHI_DEFAULT;
  input_params->phys.cosmic = COSMIC_DEFAULT;
  input_params->phys.grain_size = GRAIN_SIZE_DEFAULT;
  input_params->phys.grain_abundance = 0;
  input_params->phys.grain_mass_density = GRAIN_MASS_DENSITY_DEFAULT;
  input_params->phys.grain_gas_mass_ratio = GRAIN_GAS_MASS_RATIO_DEFAULT;
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
      if (line[0] == '#')
        continue;               /* Skip comments */
      if (sscanf (line, "[ %512[a-zA-Z] ]", keyword) == 1)
        ;
      else if (sscanf (line, "%s = %s", parameter, value) == 2)
        {

          /* Source and chemical network files */

          if (strcmp (keyword, "files") == 0)
            {
              if (strcmp (parameter, "source") == 0)
                strncpy (input_params->files.source_file, value, MAX_LINE);
              else if (strcmp (parameter, "chem") == 0)
                strncpy (input_params->files.chem_file, value, MAX_LINE);
              else
                {
                  fprintf (stderr, "astrochem: error: incorrect input in %s line %i.\n",
                           input_file, line_number);
                  return EXIT_FAILURE;
                }
            }

          /* Physical parameters */

          else if (strcmp (keyword, "phys") == 0)
            {
              if (strcmp (parameter, "chi") == 0)
                input_params->phys.chi = atof (value);
              else if (strcmp (parameter, "cosmic") == 0)
                input_params->phys.cosmic = atof (value);
              else if (strcmp (parameter, "grain_size") == 0)
                input_params->phys.grain_size = atof (value) * 1e-4;    /* microns -> cm */
              else if (strcmp (parameter, "grain_gas_mass_ratio") == 0)
                input_params->phys.grain_gas_mass_ratio = atof (value);
              else if (strcmp (parameter, "grain_mass_density") == 0)
                input_params->phys.grain_mass_density = atof (value) * 1e-3; /* kg/m^3-> g/cm^3 */

              else
                {
                  fprintf (stderr, "astrochem: error: incorrect input in %s line %i.\n",
                           input_file, line_number);
                  return EXIT_FAILURE;
                }
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
                {
                  fprintf (stderr, "astrochem: error: incorrect input in %s line %i.\n",
                           input_file, line_number);
                  return EXIT_FAILURE;
                }
            }

          /* Initial abundances */

          else if (strcmp (keyword, "abundances") == 0)
            {
              if (i < input_params->abundances.n_initial_abundances)
                {
                  int species_idx = find_species (parameter, network);
                  if (species_idx < 0)
                    {
                      fprintf (stderr,
                               "astrochem: warning: %s initial abundance given, "
                               "but is not in the network.\n", parameter);
                      input_params->abundances.n_initial_abundances--;
                    }
                  else
                    {
                      int k;
                      bool duplicated = false;
                      for( k = 0; k<i; k++ )
                        {
                          if(  input_params->abundances.initial_abundances[k].species_idx == species_idx )
                            {
                              duplicated = true;
                              fprintf (stderr,"astrochem: warning: duplicated initial abundances, keeping only the first"
                                       "initial abundances of %s.\n", parameter);
                              break;
                            }
                        }
                      if( duplicated )
                        {
                          continue;
                        }
                      input_params->abundances.
                       initial_abundances[i].species_idx = species_idx;
                      input_params->abundances.
                        initial_abundances[i].abundance = atof (value);

                      /* Compute the total grain density */
                      int g, gm, gp;
                      g = find_species ("grain", network);
                      gm = find_species ("grain(-)", network);
                      gp = find_species ("grain(+)", network);
                      if (input_params->abundances.
                          initial_abundances[i].species_idx == g
                          || input_params->abundances.
                          initial_abundances[i].species_idx == gm
                          || input_params->abundances.
                          initial_abundances[i].species_idx == gp)
                        input_params->phys.grain_abundance +=
                          input_params->abundances.
                          initial_abundances[i].abundance;
                      i++;
                    }
                }
              else
                {
                  fprintf (stderr,
                           "astrochem: error: the number of species is incorrect, != %i, file "
                           "%s may be corrupt .\n",
                           input_params->abundances.n_initial_abundances,
                           input_file);
                  return EXIT_FAILURE;
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
                  if( all_species )
                    {
                      //In case of all species , copy all specie idxs
                      int idx;
                      for( idx = 0; idx < network->n_species; idx++ )
                        {
                          input_params->output.output_species_idx[idx] = idx;
                        }
                    }
                  else
                    {
                      const char delimiter[] = ",";
                      char *output_specie;
                      output_specie = strtok (value, delimiter);
                      while( output_specie != NULL )
                        {
                          if (j >= input_params->output.n_output_species)
                            {
                              fprintf (stderr,
                                       "astrochem: error: the number of species in output exceeds %i, incoherent input file.\n",
                                       input_params->output.n_output_species);
                              break;
                            }
                          int species_idx = find_species (output_specie, network);
                          if (species_idx < 0)
                            {
                              fprintf (stderr,
                                       "astrochem: warning: %s abundance requested, "
                                       "but is not in the network.\n", output_specie);
                              input_params->output.n_output_species--;
                              output_specie = strtok ( NULL, delimiter);
                              continue;
                            }
                          int k;
                          bool duplicated = false;
                          for( k=0; k<j; k++ )
                            {
                              if( input_params->output.output_species_idx[k] == species_idx )
                                {
                                  fprintf (stderr,
                                           "astrochem: error: duplicated output species in input file : %s.\n", output_specie );
                                  input_params->output.n_output_species--;
                                  duplicated = true;
                                  break;
                                }
                            }
                          if ( duplicated )
                            {
                              output_specie = strtok ( NULL, delimiter);
                              continue;
                            }
                          input_params->output.output_species_idx[j] = species_idx;
                          j++;
                          output_specie = strtok ( NULL, delimiter);
                        }
                      if( j < input_params->output.n_output_species )
                        {
                          fprintf (stderr,
                                   "astrochem: error: the number of species in output inferior to %i, incoherent input file.\n",
                                   input_params->output.n_output_species);
                          input_params->output.n_output_species = j-1;
                        }
                    }
                }
              else if (strcmp (parameter, "trace_routes") == 0)
                input_params->output.trace_routes = atoi (value);
              else if (strcmp (parameter, "suffix") == 0)
                strncpy (input_params->output.suffix, value, MAX_LINE);
              else
                {
                  fprintf (stderr, "astrochem: error: incorrect input in %s line %i.\n",
                           input_file, line_number);
                  return EXIT_FAILURE;
                }
            }

          /* Unknown or unspecified keyword */

          else
            {
              fprintf (stderr, "astrochem: error: incorrect input in %s line %i.\n",
                       input_file, line_number);
              return EXIT_FAILURE;
            }
        }

      /* Error while reading a parameter/values */

      else
        {
          fprintf (stderr, "astrochem: error: incorrect input in %s line %i.\n",
                   input_file, line_number);
          return EXIT_FAILURE;
        }
    }

  /* Close the file */

  fclose (f);

  // Compute grain mass and grain abundance
  if( input_params->phys.grain_abundance != 0 )
    {
      fprintf (stderr,
               "astrochem: warning: use of a input grain abundance is deprecated, yet still supported. "
               "Please use the grain_gas_mass_ratio and grain_mass_density parameters instead.\n");
    }
  else
    {
      double grain_mass = (4./3. * M_PI * pow( input_params->phys.grain_size, 3) * input_params->phys.grain_mass_density);
      input_params->phys.grain_abundance = input_params->phys.grain_gas_mass_ratio * (CONST_CGSM_MASS_PROTON) / grain_mass;

      int g, gm, gp;
      g = find_species ("grain", network);
      gm = find_species ("grain(-)", network);
      gp = find_species ("grain(+)", network);
      if( g != -1 )
        {
          network->species[g].mass = grain_mass;
        }
      if( gm != -1 )
        {
          network->species[gm].mass = grain_mass;
        }
      if( gp != -1 )
        {
          network->species[gm].mass = grain_mass;
        }
    }

  /* Check that the source file name and the chemical file name were specified
     in the input file. Also check that other parameters have acceptable
     values. */

  if (strcmp (input_params->files.chem_file, "") == 0)
    {
      fprintf (stderr,
               "astrochem: error: no chemical network file specified in %s.\n",
               input_file);
      return EXIT_FAILURE;
    }
  if (strcmp (input_params->files.source_file, "") == 0)
    {
      fprintf (stderr,
               "astrochem: error: no source model file specified in %s.\n",
               input_file);
      return EXIT_FAILURE;
    }
  if (input_params->solver.ti <= 0.)
    {
      fprintf (stderr,
               "astrochem: warning: incorrect ti value specified in %s."
               "Assuming default value.\n", input_file);
      input_params->solver.ti = TI_DEFAULT;
    }
  if (input_params->solver.tf < input_params->solver.ti)
    {
      fprintf (stderr,
               "astrochem: warning: incorrect tf value specified in %s."
               "Assuming default value.\n", input_file);
      input_params->solver.tf = TF_DEFAULT;
    }

  // Check input abundances are neutral
  double gas_charge = 0;
  for( i = 0; i < input_params->abundances.n_initial_abundances; i++ )
    {
      gas_charge+= input_params->abundances.initial_abundances[i].abundance *
       network->species[ input_params->abundances.initial_abundances[i].species_idx ].charge;
    }
  if( gas_charge > DBL_EPSILON || gas_charge < -DBL_EPSILON )
    {
      fprintf (stderr, "astrochem: warning: gas charge is not neutral, charge = %e\n" ,gas_charge);
    }
  return EXIT_SUCCESS;
}

/*
  Read the file containing the source model.
*/

int
read_source (const char *source_file, mdl_t * source_mdl,
             const inp_t * input_params, const int verbose)
{
  FILE *f;
  char line[MAX_LINE];
  int line_number = 0;
  int n_cell = 0;
  int allocated = 0;
  SOURCE_READ_MODE mode = R_STATIC;     //0 static, 1 dynamic, 2 time_step reading
  double av, nh, tgas, tdust;



  /* Open the input file or exit if we can't open it */
  if (verbose == 1)
    fprintf (stdout, "Reading source model from %s.\n", source_file);
  f = fopen (source_file, "r");
  if (!f)
    {
      fprintf (stderr, "astrochem: error: can't open %s.\n", source_file);
      return EXIT_FAILURE;
    }
  int n_line = get_nb_active_line (source_file);
  if( n_line == -1 )
    {
      return EXIT_FAILURE;
    }
  /* Check the model has at least one cell */
  if (n_line == 0)
    {
      fprintf (stderr, "astrochem: error: no valid lines found in %s.\n",
               source_file);
      return EXIT_FAILURE;
    }

  /* Loop over the lines, and read the cell number, visual
     extinction, density, gas and dust temperature. */
  while (fgets (line, MAX_LINE, f) != NULL)
    {
      line_number++;
      if (line[0] == '#')
        continue;               /* Skip comments */

      if (allocated == 0)       //First line without comment determine if it is a static or dynamic source.
        {
          if (strncmp (line, "[times]", 7) == 0)        // Dynamic source file begin with [times] section
            {
              // Get the number of time steps, the number of cells*time_steps, check correctness
              int nts = get_nb_active_line_section (source_file, "times");
              if( nts == -1 )
                {
                  return EXIT_FAILURE;
                }
              int n_cells_times =
               get_nb_active_line_section (source_file, "cells");
              if( n_cells_times == -1 )
                {
                  return EXIT_FAILURE;
                }
              if (n_cells_times % nts != 0)
                {
                  fprintf (stderr,
                           "astrochem: %s: %d: error: incorrect format in source file %s .\n",
                           __FILE__, __LINE__, source_file);
                  return EXIT_FAILURE;
                }
              //Allocate the model
              if( alloc_mdl (source_mdl, n_cells_times / nts, nts) != EXIT_SUCCESS )
                {
                  return EXIT_FAILURE;
                }
              //Set reading flag to read time steps
              mode = R_TIMES;
              //Set source flag to dynamic
              source_mdl->mode = DYNAMIC;
              allocated = 1;
              continue;
            }
          else                  //Static source
            {
              //Allocate static model
              if( alloc_mdl (source_mdl, n_line, input_params->output.time_steps) != EXIT_SUCCESS )
                {
                  return EXIT_FAILURE;
                }
              //Set source flag to static
              source_mdl->mode = STATIC;
              allocated = 1;
              /* Build the vector of time */
              int i;
              for (i = 0; i < input_params->output.time_steps; i++)
                {
                  source_mdl->ts.time_steps[i] =
                    pow (10.,
                         log10 (input_params->solver.ti) +
                         i * (log10 (input_params->solver.tf) -
                              log10 (input_params->solver.ti)) /
                         (input_params->output.time_steps - 1));
                }
            }
        }
      if (strncmp (line, "[cells]", 7) == 0)    //Time to read the cells
        {
          //Set reading flag to read the dynamic cells
          mode = R_DYNAMIC;
          continue;
        }
      //Dynamic mode, reading time steps
      if (mode == R_TIMES)
        {
          int tmp_ts;
          double ts_val;
          sscanf (line, "%d %lf", &tmp_ts, &ts_val);    //Format and value not checked
          if (tmp_ts < source_mdl->ts.n_time_steps)
            {
              source_mdl->ts.time_steps[tmp_ts] = ts_val * CONST_MKSA_YEAR; // seconds
            }
          else
            {
              fprintf (stderr,
                       "astrochem: error: the time_steps idx %i in %s exceed %i.\n",
                       tmp_ts, source_file, source_mdl->ts.n_time_steps);
              return EXIT_FAILURE;
            }
        }
      //Dynamic mode, reading cells
      else if (mode == R_DYNAMIC)
        {
          int tmp_cell, tmp_ts;
          //Format and value not checked
          sscanf (line, "%d %d %lf %lf %lf %lf", &tmp_cell, &tmp_ts, &av, &nh,
                  &tgas, &tdust);
          if (tmp_ts < source_mdl->ts.n_time_steps
              || tmp_cell < source_mdl->n_cells)
            {
              source_mdl->cell[tmp_cell].av[tmp_ts] = av;
              source_mdl->cell[tmp_cell].nh[tmp_ts] = nh;
              source_mdl->cell[tmp_cell].tgas[tmp_ts] = tgas;
              source_mdl->cell[tmp_cell].tdust[tmp_ts] = tdust;
            }
          else
            {
              fprintf (stderr,
                       "astrochem: error: the time_steps idx %i in %s exceed %i or the cell idx %i exceed %i.\n",
                       tmp_ts, source_file, source_mdl->ts.n_time_steps,
                       tmp_cell, source_mdl->n_cells);
              return EXIT_FAILURE;
            }
        }
      //Static mode, reading cells
      else
        {
          if (n_cell < source_mdl->n_cells)
            {
              int tmp;
              if (sscanf (line, "%d %lf %lf %lf %lf",
                          &tmp, &av, &nh, &tgas, &tdust) != 5)
                {
                  fprintf (stderr,
                           "astrochem: %s: %d: error: incorrect format in source file %s .\n",
                           __FILE__, __LINE__, source_file);
                  return EXIT_FAILURE;
                }
              source_mdl->cell[n_cell].av[0] = av;
              source_mdl->cell[n_cell].nh[0] = nh;
              source_mdl->cell[n_cell].tgas[0] = tgas;
              source_mdl->cell[n_cell].tdust[0] = tdust;
              n_cell++;
            }
          else
            {
              fprintf (stderr,
                       "astrochem: error: the number of cells in %s exceed %i.\n",
                       source_file, source_mdl->n_cells);
              return EXIT_FAILURE;
            }
        }
    }
  /* Close the file */
  fclose (f);
  return EXIT_SUCCESS;
}

/*
  Allocate the structure containing the input parameters
*/

int
alloc_input (inp_t * input_params, int n_initial_abundances,
             int n_output_abundances)
{
  input_params->abundances.n_initial_abundances = n_initial_abundances;
  input_params->output.n_output_species = n_output_abundances;
  if ((input_params->abundances.initial_abundances =
       malloc (sizeof (abund_t) * n_initial_abundances)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
               __FILE__, __LINE__);
      return EXIT_FAILURE;
    }
  if ((input_params->output.output_species_idx =
       malloc (sizeof (char *) * n_output_abundances)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
               __FILE__, __LINE__);
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

/*
   Free the input structure.
 */
void
free_input (inp_t * input_params)
{
  free (input_params->output.output_species_idx);
  free (input_params->abundances.initial_abundances);
}

/*
  Allocate the structure containing the source model
*/

int
alloc_mdl (mdl_t * source_mdl, int n_cells, int n_time_steps)
{
  source_mdl->n_cells = n_cells;
  source_mdl->ts.n_time_steps = n_time_steps;
  if ((source_mdl->ts.time_steps =
       malloc (sizeof (double) * n_time_steps)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
               __FILE__, __LINE__);
      return EXIT_FAILURE;
    }
  /* We want the cells in source strcutres have the following layout in memory :

     [ cells[i] -> [  [ nh[0] , nh[1], .. , nh[n_time_steps-1] ] , [av[0], av[1], ..,] , [ tgas[0],.. ], [ tdust[0] ] ] , cells[i+1] -> ... ]
     although cells[i] et cells[i+1] could be pointing at two very different place, we want it to be memory contiguous


     We use a memory alignement technique, so  all cells will be contiguous and in each cells,
     nh[], av[], tgas[] and tdust[] will be contiguous
     First allocate a big data block wich will contain all data of all cells, ie all nh,av,tgas,tdust , for ech time steps, for each cells.
     Then allocate a array of struct of pointer wich will contain all pointer to the data.
     Finally point each pointer to the right block of data.
   */
  double *data;
  if ((data = malloc (4 * n_cells * n_time_steps * sizeof (double))) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
               __FILE__, __LINE__);
      return EXIT_FAILURE;
    }
  if ((source_mdl->cell = malloc (sizeof (cell_table_t) * n_cells)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
               __FILE__, __LINE__);
      return EXIT_FAILURE;
    }
  int i;
  for (i = 0; i < n_cells; i++)
    {
      source_mdl->cell[i].nh = &(data[4 * i * n_time_steps]);
      source_mdl->cell[i].av = &(data[(4 * i + 1) * n_time_steps]);
      source_mdl->cell[i].tgas = &(data[(4 * i + 2) * n_time_steps]);
      source_mdl->cell[i].tdust = &(data[(4 * i + 3) * n_time_steps]);
    }
  return EXIT_SUCCESS;
}

/*
   Free the source model structure
*/

void
free_mdl (mdl_t * source_mdl)
{
  // This is to free the big block of data
  free (source_mdl->cell[0].nh);
  free (source_mdl->cell);
  free (source_mdl->ts.time_steps);
}

/*
  Get the number of non commented line in a section of an input file
*/

int
get_nb_active_line_section (const char *file, const char *section)
{
  FILE *f;
  char line[MAX_LINE];
  int line_number = 0;
  f = fopen (file, "r");
  int section_flag = 0;
  int section_len = strlen (section);
  char full_section[section_len + 2];
  full_section[0] = '[';
  strcpy (&full_section[1], section);
  full_section[section_len + 1] = ']';
  if (!f)
    {
      fprintf (stderr, "astrochem: error: can't open %s.\n", file);
      return -1;
    }
  while (fgets (line, MAX_LINE, f) != NULL)
    {
      if (line[0] != '#')
        {
          if (section_flag == 0)
            {
              if (strncmp (line, full_section, section_len + 2) == 0)
                {
                  section_flag = 1;
                }
            }
          else
            {
              if (line[0] == '[')
                {
                  break;
                }
              else
                {
                  line_number++;
                }
            }
        }
    }
  fclose (f);
  return line_number;
}

/*
  Read the source model and network filenames from an input file
*/

int
read_input_file_names (const char *input_file, files_t * files, int verbose)
{
  FILE *f;
  char line[MAX_LINE];
  char keyword[MAX_LINE];
  char parameter[MAX_LINE];
  char value[MAX_LINE];
  int line_number = 0;
  errno = 0;

  /* Open the input file or exit if we can't open it. */

  if (verbose == 1)
    fprintf (stdout, "Reading input from %s.\n", input_file);
  f = fopen (input_file, "r");
  if (!f)
    {
      fprintf (stderr, "astrochem: error: Can't open %s: %s\n", input_file,
               strerror (errno));
      return EXIT_FAILURE;
    }

  strcpy (files->source_file, "");
  strcpy (files->chem_file, "");

  /* Loop over the lines, and look for keywords (between brackets) and
     parameters/values (separated by "="). */

  while (fgets (line, MAX_LINE, f) != NULL)
    {
      line_number++;
      if (line[0] == '#')
        continue;               /* Skip comments */
      if (sscanf (line, "[ %512[a-zA-Z] ]", keyword) == 1)
        ;
      else if (sscanf (line, "%s = %s", parameter, value) == 2)
        {
          /* Source and chemical network files */
          if (strcmp (keyword, "files") == 0)
            {
              if (strcmp (parameter, "source") == 0)
                strncpy (files->source_file, value, MAX_LINE);
              else if (strcmp (parameter, "chem") == 0)
                strncpy (files->chem_file, value, MAX_LINE);
              else
                {
                  fprintf (stderr, "astrochem: error: incorrect input in %s line %i.\n",
                           input_file, line_number);
                  return EXIT_FAILURE;
                }
            }
        }
    }
  fclose (f);
  return EXIT_SUCCESS;
}

/*
  Get the number of non commented lines in a file 
*/

int
get_nb_active_line (const char *file)
{
  FILE *f;
  char line[MAX_LINE];
  int line_number = 0;
  f = fopen (file, "r");
  if (!f)
    {
      fprintf (stderr, "astrochem: error: can't open %s.\n", file);
      return -1;
    }
  while (fgets (line, MAX_LINE, f) != NULL)
    {
      if (line[0] != '#')
        {
          line_number++;
        }
    }
  fclose (f);
  return line_number;
}
