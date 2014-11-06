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

void realloc_network_species (net_t * network, int n_species);

/* Allocate the network structure. We get the number of reactions
   from the number of lines in the network file. For the number of
   species, we assume a number equal to the number of reactions
   divided by 10, and we reallocate the array if needed. */
void
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
          exit (1);
        }
    }
  fclose (f);

  // Get the number of reaction and estimates number of species to alloc
  int n_reactions = get_nb_active_line (chem_file1);
  if (n_reactions == 0)
    {
      fprintf (stderr,
               "astrochem: error: the number of reactions is zero.\n");
      exit (1);
    }
  int n_alloc_species = n_reactions / 10;
  if (n_alloc_species == 0)
    {
      n_alloc_species = 1;
    }

  //Allocate network dynamic array
  alloc_network (network, n_alloc_species, n_reactions);
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
          exit (1);
        }

      char *localLine = line;
      char *localLine2 = strchr (line, ' ');
      char str[MAX_CHAR_SPECIES];
      unsigned int mode = 0;
      int ending = 0;

      /* Read the reactants, products, and reaction parameters. one after another */
      while (localLine != NULL)
        {
          // Analyse current char
          if (localLine[0] == '-')
            {
              // We read a "->" time to read products
              localLine = localLine2 + 1;
              localLine2 = strchr (localLine, ' ');
              mode = MAX_REACTANTS;
              ending = 0;
              continue;
            }
          else if (localLine[0] == '+')
            {
              // We read a "+", let's read another species
              localLine = localLine2 + 1;
              localLine2 = strchr (localLine, ' ');
              ending = 0;
              continue;
            }
          else if (localLine[0] == ' ')
            {
              // We read a " ", lets read another char
              localLine = localLine2 + 1;
              localLine2 = strchr (localLine, ' ');
              continue;
            }
          // Last part of the line
          else if (ending == 1)
            {
              //No more species to read, let's read the ending of reactions.
              if (sscanf (localLine, "%lf %lf %lf %d %d",
                          &network->reactions[n].alpha,
                          &network->reactions[n].beta,
                          &network->reactions[n].gamma,
                          &network->reactions[n].reaction_type,
                          &network->reactions[n].reaction_no) != 5)
                {
                  network_file_error (chem_file1, n + 1);
                }
              break;
            }
          else
            {
              // Yet another specie to add
              strncpy (str, localLine, localLine2 - localLine);
              str[localLine2 - localLine] = '\0';
              //Some species are not to read
              if ((strcmp (str, "cosmic-ray") != 0) &&
                  (strcmp (str, "uv-photon") != 0) &&
                  (strcmp (str, "photon") != 0))
                {
                  /* Add the species in list
                     Find the correct place for this species in the reactions */
                  if ( mode < MAX_REACTANTS )
                    {
                      network->reactions[n].reactants[ mode ] = add_species (str, network);
                      mode++;
                    }
                  else if ( mode < MAX_REACTANTS + MAX_PRODUCTS )
                    {
                      network->reactions[n].products[ mode - MAX_REACTANTS ] =  add_species (str, network);
                      mode++;
                    }
                  else
                    {
                      network_file_error (chem_file1, n + 1);
                      break;
                    }
                }
              // Move to next space
              localLine = localLine2 + 1;
              localLine2 = strchr (localLine, ' ');
              ending = 1;
            }
        }
      n++;
    }
  // Check number of reactions
  if (n != network->n_reactions)
    {
      fprintf (stderr,
               "astrochem: error: incorect number of reactions %i, different from %i,"
               "file %s may be corrupt.\n", n, network->n_reactions,
               chem_file1);
      exit (1);
    }
  /* Free non used memory space */
  realloc_network_species (network, network->n_species);
  if (verbose >= 1)
    {
      fprintf (stdout, "done.\n");
      fprintf (stdout, "Found %d reactions involving %d species.\n",
               network->n_reactions, network->n_species);
    }
  /* Close the file. */
  fclose (f);
}

/*
   Add a species in the network and return it's index
   If already present, just return the index
   */
int
add_species (char *new_species, net_t * network)
{
  int i;
  if (strcmp (new_species, "") == 0)
    return -1;
  //Look for the species, if it exists, return index
  for (i = 0; i < network->n_species; i++)
    {
      if (strcmp (network->species_names[i], new_species) == 0)
        return i;
    }
  // Check species array size and reallocate if necessary
  i = network->n_species;
  while (i >= network->n_alloc_species)
    {
      int new_alloc_size = network->n_alloc_species * 2;
      realloc_network_species (network, new_alloc_size);
    }
  //Check length of species name
  if (strlen (new_species) < MAX_CHAR_SPECIES - 1)
    {
      //Add species in network
      strcpy (network->species_names[i], new_species);
      (network->n_species)++;
    }
  else
    {
      fprintf (stderr, "astrochem: error: the number of characters of some "
               "species of the chemical network file exceeds %i.\n",
               MAX_CHAR_SPECIES);
      exit (1);
    }
  return i;
}

/*
   Look up the index of a given species in the species array.
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
      if (strncmp (network->species_names[i], specie,
                   sizeof (species_name_t)) == 0)
        {
          return i;
        }
    }

  /* Return -2 if we can not find the specie */
  return -2;
}

/*
   Alloc the network structure.
 */
void
alloc_network (net_t * network, int n_species, int n_reactions)
{
  network->n_alloc_species = n_species;
  if ((network->species_names =
       malloc (sizeof (species_name_t) * n_species)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__,
               "array allocation failed.\n");
      exit (1);
    }

  network->n_reactions = n_reactions;
  if ((network->reactions = malloc (sizeof (react_t) * n_reactions)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__,
               "array allocation failed.\n");
      exit (1);
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
}

/*
   Reallocate species array, without overwriting existing species
   */
void
realloc_network_species (net_t * network, int n_species)
{
  if (n_species < network->n_species)
    {
      fprintf (stderr,
               "astrochem: %s:%d: cannot realloc over existing species : "
               "n_species : %i, new_size: %i\n", __FILE__, __LINE__,
               network->n_species, n_species);
      exit (1);
    }
  network->n_alloc_species = n_species;
  if ((network->species_names =
       realloc (network->species_names,
                sizeof (species_name_t) * n_species)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__,
               "array allocation failed.\n");
      exit (1);
    }
}

/*
   Free the network structure.
 */
void
free_network (net_t * network)
{
  free (network->reactions);
  free (network->species_names);
}

/**
 * @brief Display an error while reading input
 *
 * Display an error message and exit of an error is encountered while
 * reading the input file.
 *
 * @param input_file file from wich the error is from
 * @param line_number line number where the error occured
 */
void
network_file_error (const char *chem_file, int line_number)
{
  fprintf (stderr, "astrochem: error: incorrect network file in %s line %i.\n",
           chem_file, line_number);
  exit (1);
}
