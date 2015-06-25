/*
   Astrochem - compute the abundances of chemical species in the
   interstellar medium as as function of time.

   Copyright (c) 2006-2014 Sebastien Maret

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

#include "astrochem.h"
#include "input.h"
#include <hdf5.h>
#include <math.h>

#ifdef HAVE_OPENMP
/*Global variable, OpenMP lock*/
omp_lock_t  lock;
#endif

#define ABUNDANCE_DATASET_RANK 3
#define ROUTE_DATASET_RANK 4
void usage (void);
void version (void);

int
main (int argc, char *argv[])
{
  inp_t input_params;
  mdl_t source_mdl;
  net_t network;
  int cell_index;

  int verbose = 1;
  char *input_file;

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
              return EXIT_SUCCESS;
              break;
            case 'V':
              version ();
              return EXIT_SUCCESS;
              break;
            case 'v':
              verbose = 2;
              break;
            case 'q':
              verbose = 0;
              break;
            default:
              usage ();
              return EXIT_FAILURE;
            }
        };
      argc -= optind;
      argv += optind;
      if (argc != 1)
        {
          usage ();
          return EXIT_FAILURE;
        }
      input_file = argv[0];
    }

  /* Read the input file */
  if( read_input_file_names (input_file, &input_params.files, verbose) != EXIT_SUCCESS )
    {
      return EXIT_FAILURE;
    }

  /* Read the chemical network file */
  if( read_network (input_params.files.chem_file, &network, verbose) != EXIT_SUCCESS )
    {
      return EXIT_FAILURE;
    }

  /* Read the input file */
  if( read_input (input_file, &input_params, &network, verbose) != EXIT_SUCCESS )
    {
      return EXIT_FAILURE;
    }

  /* Read the source model file */
  if( read_source (input_params.files.source_file, &source_mdl, &input_params,
               verbose) != EXIT_SUCCESS )
    {
      return EXIT_FAILURE;
    }

  // Hdf5 files, datatype and dataspace
  hid_t       fid, datatype, dataspace, dataset, tsDataset, tsDataspace,  speciesDataset, speciesDataspace, speciesType;
  datatype = H5Tcopy(H5T_NATIVE_DOUBLE);

  hsize_t     dimsf[ ROUTE_DATASET_RANK ]={  source_mdl.n_cells, source_mdl.ts.n_time_steps, input_params.output.n_output_species, N_OUTPUT_ROUTES };
  fid = H5Fcreate( "astrochem_output.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 
  dataspace = H5Screate_simple( ABUNDANCE_DATASET_RANK, dimsf, NULL);

  // Add Atributes
  hid_t simpleDataspace = H5Screate(H5S_SCALAR);
  hid_t attrType = H5Tcopy(H5T_C_S1);
  H5Tset_size ( attrType, MAX_CHAR_FILENAME );
  H5Tset_strpad(attrType,H5T_STR_NULLTERM);
  hid_t attrNetwork = H5Acreate( fid, "chem_file", attrType, simpleDataspace, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite( attrNetwork, attrType, realpath( input_params.files.chem_file, NULL ) );
  H5Aclose( attrNetwork );
  hid_t attrModel = H5Acreate( fid, "source_file", attrType, simpleDataspace, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite( attrModel, attrType, realpath( input_params.files.source_file, NULL ) );
  H5Aclose( attrModel );

  H5Tclose( attrType );
  H5Sclose( simpleDataspace );

  // Define chunk property
  hsize_t     chunk_dims[ ROUTE_DATASET_RANK ] = { 1, 1, input_params.output.n_output_species, N_OUTPUT_ROUTES };
  hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(prop_id, ABUNDANCE_DATASET_RANK , chunk_dims);

  // Create dataset
  dataset = H5Dcreate(fid, "Abundances", datatype, dataspace, H5P_DEFAULT, prop_id, H5P_DEFAULT);

  int i;
  hid_t dataspaceRoute, route_t_datatype, r_t_datatype, route_prop_id, routeGroup;
  hid_t routeDatasets[ input_params.output.n_output_species ];
  if (input_params.output.trace_routes)
    {

      // Create route dataspace
      dataspaceRoute = H5Screate_simple( ROUTE_DATASET_RANK, dimsf, NULL);

      // Create route datatype
      r_t_datatype = H5Tcreate (H5T_COMPOUND, sizeof(r_t));
      H5Tinsert( r_t_datatype, "reaction_number", HOFFSET(r_t, reaction_no ), H5T_NATIVE_INT);
      H5Tinsert( r_t_datatype, "reaction_rate", HOFFSET(r_t, rate), H5T_NATIVE_DOUBLE);

      route_t_datatype = H5Tcreate (H5T_COMPOUND, sizeof(rout_t));
      H5Tinsert( route_t_datatype, "formation_rate", HOFFSET(rout_t, formation ), r_t_datatype );
      H5Tinsert( route_t_datatype, "destruction_rate", HOFFSET(rout_t, destruction ), r_t_datatype );

      // Define route chunk property
      route_prop_id = H5Pcreate(H5P_DATASET_CREATE);
      H5Pset_chunk( route_prop_id, ROUTE_DATASET_RANK, chunk_dims);


      // Create each named route dataset
      routeGroup = H5Gcreate( fid, "Routes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      char routeName[6] = "route_";
      char tempName[ MAX_CHAR_SPECIES + sizeof( routeName ) ];
      for( i = 0; i < input_params.output.n_output_species ; i++ )
        {
          strcpy( tempName, routeName );
          strcat( tempName, network.species[input_params.output.output_species_idx[i]].name );
          routeDatasets[i] = H5Dcreate( routeGroup, tempName, route_t_datatype, dataspaceRoute, H5P_DEFAULT, route_prop_id, H5P_DEFAULT);
        }
    }
  // Timesteps and species
  hsize_t n_ts = source_mdl.ts.n_time_steps;
  hsize_t n_species =  input_params.output.n_output_species ;
  tsDataspace = H5Screate_simple( 1, &n_ts, NULL);
  speciesDataspace = H5Screate_simple( 1, &n_species, NULL);
  speciesType = H5Tcopy (H5T_C_S1);
  H5Tset_size (speciesType, MAX_CHAR_SPECIES );
  H5Tset_strpad(speciesType,H5T_STR_NULLTERM);

  // Create ts and species datasets
  tsDataset = H5Dcreate(fid, "TimeSteps", datatype, tsDataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  speciesDataset = H5Dcreate(fid, "Species", speciesType, speciesDataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  double* convTs = (double*) malloc( sizeof(double)* source_mdl.ts.n_time_steps );
  for( i=0; i< source_mdl.ts.n_time_steps; i++ )
    {
      convTs[i] = source_mdl.ts.time_steps[i] / CONST_MKSA_YEAR;
    }

  H5Dwrite( tsDataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, convTs );


  char speciesName [ input_params.output.n_output_species ][ MAX_CHAR_SPECIES ];
  for( i = 0; i < input_params.output.n_output_species ; i++ )
    {
      strcpy( speciesName[i], network.species[input_params.output.output_species_idx[i]].name );
    }

  H5Dwrite( speciesDataset, speciesType, H5S_ALL, H5S_ALL, H5P_DEFAULT, speciesName );

  free( convTs );
  H5Dclose( tsDataset );
  H5Dclose( speciesDataset );
  H5Tclose( speciesType );
  H5Sclose( tsDataspace );
  H5Sclose( speciesDataspace );


#ifdef HAVE_OPENMP
  /*Initialize lock*/
  omp_init_lock(&lock);


  /* Solve the ODE system for each cell. */
    {
#pragma omp parallel for schedule (dynamic, 1)
#endif
      for (cell_index = 0; cell_index < source_mdl.n_cells; cell_index++)
        {
          if (verbose >= 1)
            fprintf (stdout, "Computing abundances in cell %d...\n",
                     cell_index);
          if( full_solve ( fid, dataset, routeDatasets, dataspace, dataspaceRoute, datatype, route_t_datatype, cell_index, &input_params, source_mdl.mode,
                       &source_mdl.cell[cell_index], &network, &source_mdl.ts, verbose) != EXIT_SUCCESS )
            {
	      exit (EXIT_FAILURE);
            }
          if (verbose >= 1)
            fprintf (stdout, "Done with cell %d.\n", cell_index);
        }
#ifdef HAVE_OPENMP
    }


  /*Finished lock mechanism, destroy it*/
  omp_destroy_lock(&lock);
#endif

  /*
   * Close/release hdf5 resources.
   */
  if (input_params.output.trace_routes)
    {


      for( i = 0; i <  input_params.output.n_output_species ; i++ )
        {
          H5Dclose(routeDatasets[i] );
        }
      H5Sclose(dataspaceRoute);
      H5Gclose(routeGroup);
      H5Pclose(route_prop_id);
      H5Tclose(r_t_datatype);
      H5Tclose(route_t_datatype);
    }
  H5Dclose(dataset);
  H5Pclose(prop_id);
  H5Sclose(dataspace);
  H5Tclose(datatype);

  H5Fclose(fid);

  free_input (&input_params);
  free_mdl (&source_mdl);
  free_network (&network);
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
  fprintf (stdout, "Copyright (c) 2006-2014 Sebastien Maret\n");
  fprintf (stdout, "\n");
  fprintf (stdout,
           "This is free software. You may redistribute copies of it under the terms\n");
  fprintf (stdout,
           "of the GNU General Public License. There is NO WARRANTY, to the extent\n");
  fprintf (stdout, "permitted by law.\n");
}

int
full_solve (hid_t fid, hid_t dataset, hid_t* routeDatasets, hid_t dataspace, hid_t routeDataspace, hid_t datatype, hid_t routeDatatype, int cell_index, const inp_t * input_params, SOURCE_MODE mode,
            const cell_table_t * cell, const net_t * network, const time_steps_t * ts,
            int verbose)
{

  double *abundances = NULL;
  alloc_abundances( network, &abundances ); // Allocate the abundances array; it contains all species.

  rout_t* routes = NULL;
  if (( routes =
        malloc (sizeof (rout_t) *  input_params->output.n_output_species * N_OUTPUT_ROUTES)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: routes allocation failed.\n",
               __FILE__, __LINE__);
      return EXIT_SUCCESS;
    }

  double* output_abundances = NULL;
  if (( output_abundances =
        malloc (sizeof (double) * input_params->output.n_output_species )) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
               __FILE__, __LINE__);
      return EXIT_FAILURE;
    }

  // Create the memory dataspace, selecting all output abundances
  hsize_t size = input_params->output.n_output_species;
  hid_t memDataspace = H5Screate_simple(1, &size, NULL);

  // Create the file dataspace, and prepare selection of a chunk of the file
  hid_t fileDataspace = H5Scopy(dataspace);
  hsize_t     count[3]={  1, 1,  input_params->output.n_output_species };

  hsize_t routeSize[2] = { input_params->output.n_output_species, N_OUTPUT_ROUTES };
  hsize_t     routeCount[4]={  1, 1,  input_params->output.n_output_species, N_OUTPUT_ROUTES };
  hid_t routeFileDataspace, routeMemDataspace;
  if (input_params->output.trace_routes)
    {
      // Create the route memory dataspace, selecting all output routes
      routeMemDataspace = H5Screate_simple(2, routeSize, NULL);
      // Create the route file dataspace, and prepare selection of a chunk of the file
      routeFileDataspace = H5Scopy(routeDataspace);
    }

  // Initializing abundance
#if 0 //Ultra complicated code
  const species_name_t* species = malloc( input_params->abundances.n_initial_abundances * sizeof(*species));
  double *initial_abundances = malloc( input_params->abundances.n_initial_abundances * sizeof(double) );

  int i;
  for( i = 0; i <  input_params->abundances.n_initial_abundances ; i++ )
    {
      strcpy( network->species_names[input_params->abundances.initial_abundances[i].species_idx ] , species[i] );
      initial_abundances[i] = input_params->abundances.initial_abundances[i].abundance;
    }
  set_initial_abundances( species, 3, initial_abundances, &network, abundances); // Set initial abundances
#else // same thing , without using api
  int i;
  for( i = 0; i <  input_params->abundances.n_initial_abundances ; i++ )
    {
      abundances[ input_params->abundances.initial_abundances[i].species_idx ] = input_params->abundances.initial_abundances[i].abundance;
    }
#endif


  double min_nh;                 /* Minimum density */

  /* Compute the minimum density to set the absolute tolerance of the
     solver */
  min_nh = cell->nh[0];
  if (mode == DYNAMIC)
    {
      int i;

      for (i = 1; i < ts->n_time_steps; i++)
        {
          if (cell->nh[i] < min_nh)
            {
              min_nh = cell->nh[i];
            }
        }
    }

  astrochem_mem_t astrochem_mem;
  cell_t cell_unik;
  cell_unik.av = cell->av[0];
  cell_unik.nh = cell->nh[0];
  cell_unik.tgas = cell->tgas[0];
  cell_unik.tdust = cell->tdust[0];
  if( solver_init( &cell_unik, network, &input_params->phys, abundances, min_nh, input_params->solver.abs_err,  input_params->solver.rel_err, &astrochem_mem ) != EXIT_SUCCESS )
    {
      return EXIT_FAILURE;
    }
  else
    {
      int i, j;

      /* Solve the system for each time step. */
      for (i = 0; i < ts->n_time_steps; i++)
        {


          if (i!=0 && mode == DYNAMIC)
            {
              cell_unik.av = cell->av[i];
              cell_unik.nh = cell->nh[i];
              cell_unik.tgas = cell->tgas[i];
              cell_unik.tdust = cell->tdust[i];

              if( solve( &astrochem_mem, network, abundances,  ts->time_steps[i], &cell_unik, verbose ) != EXIT_SUCCESS )
                {
                  return EXIT_FAILURE;
                }
            }
          else
            {
              if( solve( &astrochem_mem, network, abundances,  ts->time_steps[i], NULL, verbose ) != EXIT_SUCCESS )
                {
                  return EXIT_FAILURE;
                }
            }


          /* Fill the array of abundances with the output species
             abundances. Ignore species that are not in the
             network. Abundance that are lower than MIN_ABUNDANCES are
             set to 0. */

          for (j = 0; j < input_params->output.n_output_species; j++)
            {
              if (mode == STATIC)
                {
                  output_abundances[j] =
                   (double) NV_Ith_S (astrochem_mem.y, input_params->output.output_species_idx[j]) / cell->nh[0];
                }
              else
                {
                  output_abundances[j] =
                   (double) NV_Ith_S (astrochem_mem.y, input_params->output.output_species_idx[j]) / cell->nh[i];
                }
              if (output_abundances[j] < MIN_ABUNDANCE)
                output_abundances[j] = 0.;

#ifdef HAVE_OPENMP
              omp_set_lock(&lock);
#endif
              // Select a chunk of the file
              hsize_t     start[3]={  cell_index, i, 0 };
              H5Sselect_hyperslab( fileDataspace, H5S_SELECT_SET, start, NULL, count , NULL );

              // Write the chunk
              H5Dwrite(dataset, datatype, memDataspace, fileDataspace, H5P_DEFAULT,
                       output_abundances );

#ifdef HAVE_OPENMP
              omp_unset_lock(&lock);
#endif

            }

          /* Compute the rate of each formation/destruction route for
             each output specie. */

          if (input_params->output.trace_routes)
            {
              for (j = 0; j < input_params->output.n_output_species; j++)
                {
                  int k;
                  int l;
                  for (l = 0; l < N_OUTPUT_ROUTES; l++)
                    {
                      routes[ j*N_OUTPUT_ROUTES + l ].formation.rate = 0;
                      routes[ j*N_OUTPUT_ROUTES + l ].destruction.rate = 0;
                    }
                  for (k = 0; k < network->n_reactions; k++)
                    {
                      /* If the species is a product of the
                         reaction then compute the formation
                         rate. If the rate is greater than the
                         smallest rate in the formation route
                         structure, we add the current reaction
                         number and rate to that structure. */

                      bool specie_in_products = false;
                      int p;
                      for( p = 0; p < MAX_PRODUCTS; p++ )
                        {
                          if( network->reactions[k].products[p] ==  input_params->output.output_species_idx[j])
                            {
                              specie_in_products = true;
                              break;
                            }
                        }
                      if( specie_in_products )
                        {
                          r_t formation_route;
                          double min_rate;
                          unsigned int min_rate_index;
                          if (network->reactions[k].reaction_type == 0)
                            {
                              formation_route.rate = astrochem_mem.params.reac_rates[k];
                              formation_route.rate *=
                               NV_Ith_S (astrochem_mem.y, network->reactions[k].reactants[0]);
                            }
                          else if (network->reactions[k].reaction_type == 23)
                            {
                              formation_route.rate = astrochem_mem.params.reac_rates[k];
                            }
                          else
                            {
                              formation_route.rate = astrochem_mem.params.reac_rates[k];
                              int r;
                              for( r = 0; r < MAX_REACTANTS; r++ )
                                {
                                  if( network->reactions[k].reactants[r] != -1 )
                                    {
                                      formation_route.rate *=
                                       NV_Ith_S (astrochem_mem.y, network->reactions[k].reactants[r]);
                                    }
                                }
                            }
                          formation_route.reaction_no =
                           network->reactions[k].reaction_no;
                          min_rate = routes[ j*N_OUTPUT_ROUTES  ].formation.rate;
                          min_rate_index = 0;
                          for (l = 1; l < N_OUTPUT_ROUTES; l++)
                            {
                              if (routes[ j*N_OUTPUT_ROUTES + l ].formation.rate <
                                  min_rate)
                                {
                                  min_rate =
                                   routes[ j*N_OUTPUT_ROUTES + l ].formation.rate;
                                  min_rate_index = (unsigned int) l;
                                }
                            }
                          if (formation_route.rate > min_rate)
                            {
                              routes[ j*N_OUTPUT_ROUTES + min_rate_index ].formation.rate = formation_route.rate;
                              routes[ j*N_OUTPUT_ROUTES + min_rate_index ].formation.reaction_no = formation_route.reaction_no;
                            }
                        }

                      /* If the species is reactant of the reaction
                         then compute the destruction rate. */
                      bool species_in_reactants = false;
                      int r;
                      for ( r = 0; r < MAX_REACTANTS; r++ )
                        {
                          if ( network->reactions[k].reactants[r] == input_params->output.output_species_idx[j])
                            {
                              species_in_reactants = true;
                              break;
                            }
                        }
                      if( species_in_reactants )
                        {
                          r_t destruction_route;
                          double min_rate;
                          unsigned int min_rate_index;

                          if (network->reactions[k].reaction_type == 0)
                            {
                              destruction_route.rate = astrochem_mem.params.reac_rates[k];
                              destruction_route.rate *=
                               NV_Ith_S (astrochem_mem.y, network->reactions[k].reactants[0]);
                            }
                          else if (network->reactions[k].reaction_type == 23)
                            {
                              destruction_route.rate = astrochem_mem.params.reac_rates[k];
                            }
                          else
                            {
                              destruction_route.rate = astrochem_mem.params.reac_rates[k];
                              for ( r = 0; r < MAX_REACTANTS; r++ )
                                {
                                  if (network->reactions[k].reactants[r] != -1)
                                    {
                                      destruction_route.rate *=
                                       NV_Ith_S (astrochem_mem.y, network->reactions[k].reactants[r]);
                                    }
                                }
                            }
                          destruction_route.reaction_no =
                           network->reactions[k].reaction_no;

                          min_rate = routes[ j*N_OUTPUT_ROUTES  ].destruction.rate;
                          min_rate_index = 0;
                          for (l = 1; l < N_OUTPUT_ROUTES; l++)
                            {
                              if (routes[ j*N_OUTPUT_ROUTES + l ].destruction.rate <
                                  min_rate)
                                {
                                  min_rate =
                                   routes[ j*N_OUTPUT_ROUTES + l ].destruction.rate;
                                  min_rate_index = (unsigned int) l;
                                }
                            }
                          if (destruction_route.rate > min_rate)
                            {
                              routes[ j*N_OUTPUT_ROUTES + min_rate_index ].destruction.rate = destruction_route.rate;
                              routes[ j*N_OUTPUT_ROUTES + min_rate_index ].destruction.reaction_no = destruction_route.reaction_no;
                            }
                        }
                    }
                }
#ifdef HAVE_OPENMP
              omp_set_lock(&lock);
#endif
              // Selecting a chunk of the file
              hsize_t     routeStart[4]={  cell_index, i, 0, 0 };
              H5Sselect_hyperslab( routeFileDataspace, H5S_SELECT_SET, routeStart, NULL, routeCount , NULL );

              int spec_idx;
              for( spec_idx = 0; spec_idx < input_params->output.n_output_species; spec_idx++ )
                {
                  // Writing in each route datasets
                  H5Dwrite( routeDatasets[ spec_idx ], routeDatatype, routeMemDataspace, routeFileDataspace, H5P_DEFAULT,
                            routes );
                }

#ifdef HAVE_OPENMP
              omp_unset_lock(&lock);
#endif
            }
        }

    }
  // Cleaning up hdf5
  H5Sclose(memDataspace);
  H5Sclose(fileDataspace);
  if (input_params->output.trace_routes)
    {
      H5Sclose(routeMemDataspace);
      H5Sclose(routeFileDataspace);
    }

  // Free
  free( output_abundances );
  free( routes );
  free_abundances( abundances );
  solver_close( &astrochem_mem );
  return EXIT_SUCCESS;
}
