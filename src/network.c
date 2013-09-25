/* 
   network.c - Read the chemical network file

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

int add_species (char *new_specie, net_t * network);

void realloc_network_species (net_t * network, int n_species);

void
read_network (const char *chem_file, net_t * network, const int verbose)
{
  FILE *f;
  char line[MAX_LINE];
/* Allocate the network structure. We get the number of reactions
   from the number of lines in the network file. For the number of
   species, we assume a number equal to the number of reactions
   divided by 10, and we reallocate the array if needed. 
   
   Fixme: get_nb_active_line() does not look for chem_file in
   PKGDATADIR if it isn't found in the current directory (see
   below). */ 
  // Get size of dynamic arrays from file
  int n_reactions = get_nb_active_line (chem_file);
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
      fprintf (stdout, "Reading reactions network from %s... ", chem_file);
      fflush (stdout);
    }

  /* Open the input file. We first look in the current directory, and 
     then in the PKGDATADIR directory. Exit if we can't find it. */
/* FixMe: this part of the code should be executed *before*
          +     get_nb_active_line() is called (see above). */ 
  f = fopen (chem_file, "r");
  if (!f)
    {
      char chem_file1[MAX_LINE];

      strncpy (chem_file1, PKGDATADIR, sizeof (chem_file1) - 1);
      strncat (chem_file1, "/",
               sizeof (chem_file1) - strlen (chem_file1) - 1);
      strncat (chem_file1, chem_file,
               sizeof (chem_file1) - strlen (chem_file1) - 1);
      f = fopen (chem_file1, "r");
      if (!f)
        {
          fprintf (stderr, "astrochem: error: can't find %s.\n", chem_file);
          exit (1);
        }
    }

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

      /* Read the reactants, products, and reaction parameters. */
      char *localLine = line;
      char *localLine2 = strchr (line, ' ');
      char str[MAX_CHAR_SPECIES];
      int mode = 1;
      int ending = 0;
      while (localLine != NULL)
        {
          // Analyse current char
          if (localLine[0] == '-')
            {
              localLine = localLine2 + 1;
              localLine2 = strchr (localLine, ' ');
              mode = 4;
              ending = 0;
              continue;
            }
          else if (localLine[0] == '+')
            {
              localLine = localLine2 + 1;
              localLine2 = strchr (localLine, ' ');
              ending = 0;
              continue;
            }
          else if (localLine[0] == ' ')
            {
              localLine = localLine2 + 1;
              localLine2 = strchr (localLine, ' ');
              continue;
            }
          // Last part of the line
          else if (ending == 1)
            {
              if (sscanf (localLine, "%lf %lf %lf %d %d",
                          &network->reactions[n].alpha,
                          &network->reactions[n].beta,
                          &network->reactions[n].gamma,
                          &network->reactions[n].reaction_type,
                          &network->reactions[n].reaction_no) != 5)
                {
                  input_error (chem_file, n + 1);
                }
              break;
            }
          // Yet another specie to add
          else
            {
              strncpy (str, localLine, localLine2 - localLine);
              str[localLine2 - localLine] = '\0';
              if ((strcmp (str, "cosmic-ray") != 0) &&
                  (strcmp (str, "uv-photon") != 0) &&
                  (strcmp (str, "photon") != 0))
                {
                  if (mode == 1)
                    {
                      network->reactions[n].reactant1 =
                        add_species (str, network);
                      mode++;
                    }
                  else if (mode == 2)
                    {
                      network->reactions[n].reactant2 =
                        add_species (str, network);
                      mode++;
                    }
                  else if (mode == 3)
                    {
                      network->reactions[n].reactant3 =
                        add_species (str, network);
                      mode++;
                    }
                  else if (mode == 4)
                    {
                      network->reactions[n].product1 =
                        add_species (str, network);
                      mode++;
                    }
                  else if (mode == 5)
                    {
                      network->reactions[n].product2 =
                        add_species (str, network);
                      mode++;
                    }
                  else if (mode == 6)
                    {
                      network->reactions[n].product3 =
                        add_species (str, network);
                      mode++;
                    }
                  else if (mode == 7)
                    {
                      network->reactions[n].product4 =
                        add_species (str, network);
                      mode++;
                    }
                  else
                    {
                      input_error (chem_file, n + 1);
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
  if (n != network->n_reactions)
    {
      fprintf (stderr,
               "astrochem: error: incorect number of reactions %i, different from %i,"
               "file %s may be corrupt.\n", n, network->n_reactions,
               chem_file);
      exit (1);
    }
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

int
add_species (char *new_specie, net_t * network)
{
  int i;
  if (strcmp (new_specie, "") == 0)
    return -1;
  for (i = 0; i < network->n_species; i++)
    {
      if (strcmp (network->species_names[i], new_specie) == 0)
        return i;
    }
  i = network->n_species;
  while (i >= network->n_alloc_species)
    {
      int new_alloc_size = network->n_alloc_species * 2;
      realloc_network_species (network, new_alloc_size);
    }
  if (strlen (new_specie) < MAX_CHAR_SPECIES - 1)
    {
      strcpy (network->species_names[i], new_specie);
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
   Look up the index of a given specie in the species array.
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
  int i;
  for (i = 0; i < n_reactions; i++)
    {
      network->reactions[i].reactant1 = -1;
      network->reactions[i].reactant2 = -1;
      network->reactions[i].reactant3 = -1;
      network->reactions[i].product1 = -1;
      network->reactions[i].product2 = -1;
      network->reactions[i].product3 = -1;
      network->reactions[i].product4 = -1;
    }
}

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
