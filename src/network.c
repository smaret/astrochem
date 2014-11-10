/*
   network.c - Read the chemical network file

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libastrochem.h"
#include "network.h"



int add_species (char *new_species, net_t * network);
int alloc_network (net_t * network, int n_species, int n_reactions);
int realloc_network_species (net_t * network, int n_species);
int get_element_mass( const char* element );

/**
 * @brief Allocate the network structure.
 *
 * Allocate the network structure. We get the number of reactions
 * from the number of lines in the network file. For the number of
 * species, we assume a number equal to the number of reactions
 * divided by 10, and we reallocate the array if needed.
 *
 * @param chem_file chem file to use
 * @param network network to fill
 * @param verbose quiet if 0, verbose if 1
 * @return EXIT_SUCCESS if sucessfull
 */
int
read_network (const char *chem_file, net_t * network, const int verbose)
{
  FILE *f;
  char line[MAX_LINE];

  /*   Find the input file. We first look in the current directory, and
       then in the PKGDATADIR directory. Exit if we can't find it. */

  char chem_file1[MAX_LINE];
  strcpy (chem_file1, chem_file);
  f = fopen (chem_file1, "r");
  if (!f)
    {
      strncpy (chem_file1, PKGDATADIR, sizeof (chem_file1) - 1);
      strncat (chem_file1, "/",
               sizeof (chem_file1) - strlen (chem_file1) - 1);
      strncat (chem_file1, chem_file,
               sizeof (chem_file1) - strlen (chem_file1) - 1);
      f = fopen (chem_file1, "r");
      if (!f)
        {
          fprintf (stderr, "aaastrochem: error: can't find %s.\n", chem_file);
          return EXIT_FAILURE;
        }
    }
  fclose (f);

  // Get the number of reaction and estimates number of species to alloc
  int n_reactions = get_nb_active_line (chem_file1);
  if( n_reactions == -1 )
    {
      fprintf (stdout, "done1.\n");
      return EXIT_FAILURE;
    }
  if (n_reactions == 0)
    {
      fprintf (stderr,
               "astrochem: error: the number of reactions is zero.\n");
      return EXIT_FAILURE;
    }
  int n_alloc_species = n_reactions / 10;
  if (n_alloc_species == 0)
    {
      n_alloc_species = 1;
    }

  //Allocate network dynamic array
  if( alloc_network (network, n_alloc_species, n_reactions) != EXIT_SUCCESS )
    {
      fprintf (stdout, "done22.\n");
      return EXIT_FAILURE;
    }

  network->n_species = 0;
  if (verbose >= 1)
    {
      fprintf (stdout, "Reading reactions network from %s... ", chem_file1);
      fflush (stdout);
    }


  f = fopen (chem_file1, "r");
  /* Loop over the lines, and look for the reactants and products, and
     the parameters of the reactions. */
  int n = 0;
  while (fgets (line, MAX_LINE, f) != NULL)
    {
      if (line[0] == '#')
        continue;               /* Skip comments. */

      if (n >= network->n_reactions)
        {
          fprintf (stderr,
                   "astrochem: error: incorect number of reactions exceed %i,"
                   "file %s may be corrupt.\n", network->n_reactions,
                   chem_file);
          return EXIT_FAILURE;
        }

      char *tmpLine = line;

      // Separating reaction line ( reactants ) -> ( products params )
      char* reaction_arrow = strchr( tmpLine, '>' );
      *(reaction_arrow-2) = '\0';
      char* reactants = tmpLine;
      char* products = reaction_arrow+1;

      // Parsing reactants
      char* specie;
      int nreactants = 0;

      // First non-whitespace char
      specie = strtok( reactants, " " );
      while( specie != NULL )
        {
          // Not a '+', it is a specie
          if( strcmp( specie, "+" ) != 0 )
            {
              // Check number of reactant
              if( nreactants == MAX_REACTANTS )
                {
                  fprintf (stderr,
                           "astrochem: error: number of reactant %i, is greater than %i,"
                           "file %s may be corrupt.\n", nreactants+1, MAX_REACTANTS ,
                           chem_file1);
                  return EXIT_FAILURE;
                }

              // Add specie in network and in reaction
              // Not storing some elements
              if ((strcmp ( specie, "cosmic-ray") != 0) &&
                  (strcmp ( specie, "uv-photon") != 0) &&
                  (strcmp ( specie, "photon") != 0))
                {
                  network->reactions[n].reactants[ nreactants ] =  add_species (specie, network);
                  if( network->reactions[n].reactants[ nreactants ] == -1 )
                    {
                      fprintf (stdout, "done3.\n");
                      return EXIT_FAILURE;
                    }

                  // Increment the number of reactants
                  nreactants++;
                }
            }
          specie = strtok( NULL, " " );
        }

      int nproducts = 0;
      bool specie_ready = true;
      // First non-whitespace char
      specie = strtok( products, " " );
      while( specie != NULL )
        {
          // Found a '+' , be ready for next specie
          if( strcmp( specie, "+" ) == 0 )
            {
              specie_ready = true;
            }
          // Found a specie
          else if( specie_ready )
            {
              specie_ready = false;
              // Check number of products
              if( nproducts == MAX_PRODUCTS )
                {
                  fprintf (stderr,
                           "astrochem: error: number of products %i, is greater than %i,"
                           "file %s may be corrupt.\n", nproducts+1, MAX_PRODUCTS ,
                           chem_file1);
                  return EXIT_FAILURE;
                }

              // Add specie in network and in reaction
              // Not storing some elements
              if ((strcmp ( specie, "cosmic-ray") != 0) &&
                  (strcmp ( specie, "uv-photon") != 0) &&
                  (strcmp ( specie, "photon") != 0))
                {

                  network->reactions[n].products[ nproducts ] =  add_species (specie, network);
                  if( network->reactions[n].products[ nproducts ] == -1 )
                    {
                      fprintf (stdout, "done4.\n");
                      return EXIT_FAILURE;
                    }

                  // Increment the number of reactants
                  nproducts++;
                }
            }
          // Found a char, wich is not a '+', without being ready for a specie : end of products
          else
            {
              break;
            }
          // Point to next non-whitespace char
          specie = strtok( NULL, " " );
        }

      // After the products, there is reaction params, wich need to be stored also

      // Alpha has already been read as the last products
      if (sscanf ( specie, "%lf",
                   &network->reactions[n].alpha) != 1)
        {
          fprintf (stderr, "astrochem: error: incorrect network file in %s line %i.\n",
                   chem_file1, n+1 );
          return EXIT_FAILURE;
        }

      // The remaining params
      char* params =  strtok( NULL, "" );
      if (sscanf ( params, "%lf %lf %d %d",
                   &network->reactions[n].beta,
                   &network->reactions[n].gamma,
                   &network->reactions[n].reaction_type,
                   &network->reactions[n].reaction_no) != 4)
        {
          fprintf (stderr, "astrochem: error: incorrect network file in %s line %i.\n",
                   chem_file1, n+1 );
          return EXIT_FAILURE;
        }
      // Next reaction
      n++;
    }
  // Check number of reactions
  if (n != network->n_reactions)
    {
      fprintf (stderr,
               "astrochem: error: incorect number of reactions %i, different from %i,"
               "file %s may be corrupt.\n", n, network->n_reactions,
               chem_file1);
      return EXIT_FAILURE;
    }
  /* Free non used memory space */
  if( realloc_network_species (network, network->n_species) != EXIT_SUCCESS )
    {
      fprintf (stdout, "done5.\n");
      return EXIT_FAILURE;
    }
  if (verbose >= 1)
    {
      fprintf (stdout, "done.\n");
      fprintf (stdout, "Found %d reactions involving %d species.\n",
               network->n_reactions, network->n_species);
    }
  /* Close the file. */
  fclose (f);
  return EXIT_SUCCESS;
}

/**
 * @brief Add a species in the network
 *
 * Add a species in the network and return it's index
 * If already present, just return the index
 *
 * @param new_species name of new specie
 * @param network network to add specie to
 * @return index of added specie or -1 in case of failure
 */
int
add_species (char *new_species, net_t * network)
{
  // Check specie is not empty
  if (strcmp (new_species, "") == 0)
    {
      fprintf (stderr, "astrochem: error: Adding empty specie\n", new_species );
      return -1;
    }
  //Look for the species, if it exists, return index
  int i;
  for (i = 0; i < network->n_species; i++)
    {
      if (strcmp (network->species[i].name, new_species) == 0)
        return i;
    }
  // Check species array size and reallocate if necessary
  i = network->n_species;
  while (i >= network->n_alloc_species)
    {
      int new_alloc_size = network->n_alloc_species * 2;
      if( realloc_network_species (network, new_alloc_size) != EXIT_SUCCESS )
        {
          return EXIT_FAILURE;
        }
    }
  //Check length of species name
  if (strlen (new_species) < MAX_CHAR_SPECIES - 1)
    {
      // Initialize element
      char element[ MAX_CHAR_ELEMENT ];
      int element_idx = 0;
      char* specie_pt = new_species;
      int specie_mass = 0;
      int specie_charge = 0;
      int mass_multiplier = 1;

      // Parse specie string
      while( *specie_pt != '\0' )
        {
          // A-Z : new element, stack precedent element if any
          if( *specie_pt >= 65 && *specie_pt <= 90 )
            {
              // If there is already a completed element
              if( element_idx > 0 )
                {
                  // Null terminate it
                  element[element_idx] = '\0';
                  // Compute it's mass
                  specie_mass += mass_multiplier * get_element_mass( element );
                  // Reinitialize element
                  element_idx = 0;
                  mass_multiplier = 1;
                }
              // Store the char in element
              element[ element_idx ] = *specie_pt;
              // Next char
              element_idx++;
              specie_pt++;
            }
          // a-z element char part or electron (e)
          else if( *specie_pt >= 97 && *specie_pt <= 122 )
            {
              // Check element is not too big
              if( element_idx >= MAX_CHAR_ELEMENT-1 )
                {
                  fprintf (stderr, "astrochem: error: this specie: %s contain an element longer than MAX_CHAR_ELEMENT\n", new_species );
                  return -1;
                }
              // Store the char in element
              element[ element_idx ] = *specie_pt;
              // Next char
              element_idx++;
              specie_pt++;
            }
          // 1->9 element multiplier
          else if( *specie_pt >= 49 && *specie_pt <= 57 )
            {
              // check there is an element for this multiplier and there is not already a multiplier
              if( element_idx == 0 || mass_multiplier != 1 )
                {
                  fprintf (stderr, "astrochem: error: this specie: %s is not correctly written\n", new_species );
                  return -1;
                }
              // Store multiplier and position to next non numeric char
              mass_multiplier = strtol( specie_pt, &specie_pt, 10 );
            }
          // (*) Charge
          else if( *specie_pt == 40 && *(specie_pt+2) == 41 )
            {
              // jump to charge
              specie_pt++;
              switch( *specie_pt )
                {
                  // '+'
                case(43):
                  specie_charge = 1;
                  break;
                  // '-'
                case(45):
                  specie_charge = -1;
                  break;
                default:
                  fprintf (stderr, "astrochem: error: this specie: %s is not correctly written\n", new_species );
                  return -1;
                  break;
                }
              break;
            }
          else
            {
              fprintf (stderr, "astrochem: error: this specie: %s is not correctly written\n", new_species );
              return -1;
            }

        }
      // Specie string have been parsed, but last element still have to be taken in account
      if( element_idx > 0 )
        {
          // Null terminate element
          element[element_idx] = '\0';
          // Compute it's mass
          specie_mass += mass_multiplier * get_element_mass( element );
        }

      // Check mass of specie is not null
      if( specie_mass == 0 )
        {
          fprintf (stderr, "astrochem: error: Problam parsing this specie: %s \n", new_species );
          return -1;
        }

      // Store specie name charge and mass in network
      strcpy (network->species[i].name, new_species);
      network->species[i].charge = specie_charge;
      network->species[i].mass = specie_mass;
      network->n_species++;
    }
  else
    {
      fprintf (stderr, "astrochem: error: the number of characters of some "
               "species of the chemical network file exceeds %i.\n",
               MAX_CHAR_SPECIES);
      return -1;
    }
  return i;
}

/**
 * @brief Look up the index of a given species in the species array.
 * @param specie specie name to find
 * @param network network to search for specie
 * @return index of specie, -1 if specie name is empty, -2 if specie does not exist.
 */
int
find_species (const species_name_t specie, const net_t * network)
{
  int i;
  /* Return -1 if the specie name is empty. */
  if (strcmp (specie, "") == 0)
    {
      return -1;
    }
  for (i = 0; i < network->n_species; i++)
    {
      if (strncmp (network->species[i].name, specie,
                   sizeof (species_name_t)) == 0)
        {
          return i;
        }
    }

  /* Return -2 if we can not find the specie */
  return -2;
}

/**
 * @brief Alloc the network structure.
 * @param network network to allocate
 * @param n_species number of species
 * @param n_reactions numer of reactions
 * @return EXIT_SUCCESS if sucessfull
 */
int
alloc_network (net_t * network, int n_species, int n_reactions)
{
  network->n_alloc_species = n_species;
  if ((network->species =
       malloc (sizeof (species_t) * n_species)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__,
               "array allocation failed.\n");
      return EXIT_FAILURE;
    }

  network->n_reactions = n_reactions;
  if ((network->reactions = malloc (sizeof (react_t) * n_reactions)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__,
               "array allocation failed.\n");
      return EXIT_FAILURE;
    }
  int i,j;
  for (i = 0; i < n_reactions; i++)
    {
      for( j = 0; j < MAX_REACTANTS; j++ )
        {
          network->reactions[i].reactants[j] = -1;
        }
      for( j = 0; j < MAX_PRODUCTS; j++ )
        {
          network->reactions[i].products[j] = -1;
        }
    }
  return EXIT_SUCCESS;
}

/**
 * @brief Reallocate species array, without overwriting existing species
 * @param network network to reallocate specie in
 * @param n_species new nimber of species
 * @return EXIT_SUCCESS if sucessfull
 */
int
realloc_network_species (net_t * network, int n_species)
{
  if (n_species < network->n_species)
    {
      fprintf (stderr,
               "astrochem: %s:%d: cannot realloc over existing species : "
               "n_species : %i, new_size: %i\n", __FILE__, __LINE__,
               network->n_species, n_species);
      return EXIT_FAILURE;
    }
  network->n_alloc_species = n_species;
  if ((network->species =
       realloc (network->species,
                sizeof (species_t) * n_species)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__,
               "array allocation failed.\n");
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

/**
 * @brief Free the network structure.
 * @param network network to free
 */
void
free_network (net_t * network)
{
  free (network->reactions);
  free (network->species);
}

/**
 * @brief get the mass of an element
 * @param element element to get the mass of
 * @todo TODO actually use a periodic table instead of returning 1!
 * @return mess of elemen
 */
int
get_element_mass( const char* element )
{
  return 1;
}
