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

/*void add_specie (char *new_specie, char *species[], 
		 int *n_species);*/
int add_specie (char *new_specie, net_t* network);

void alloc_network ( net_t * network, int n_species, int n_reactions );

void realloc_network_species ( net_t * network, int n_species );
/*
  Add a specie in the species array, if not already present.
*/

void
read_network (const char *chem_file, net_t *network, const int verbose)
{
  FILE *f;
  char line[MAX_LINE];
  int  line_number = 0;
  char reactant1[MAX_CHAR_SPECIES];
  char reactant2[MAX_CHAR_SPECIES];
  char reactant3[MAX_CHAR_SPECIES]; 
  char product1[MAX_CHAR_SPECIES];
  char product2[MAX_CHAR_SPECIES];
  char product3[MAX_CHAR_SPECIES];
  char product4[MAX_CHAR_SPECIES];
  double alpha;
  double beta;
  double gamma;
  int reaction_type;
  int reaction_no;
  int n=0;
  int n_reactions = get_nb_active_line(chem_file);
  if( n_reactions > MAX_REACTIONS )
  {
	  fprintf (stderr, "astrochem: error: the number of reactions exceed %i.\n", 
		   MAX_REACTIONS);
	  exit(1);
  }
  if( n_reactions == 0)
  {
	  fprintf (stderr, "astrochem: error: the number of reactions is zero.\n", 
		   MAX_REACTIONS);
	  exit(1);
  }
  int n_alloc_species = n_reactions/10;
  if(  n_alloc_species == 0 )
  {
      n_alloc_species=1;
  }
  alloc_network(network,n_alloc_species,n_reactions);

  /*network->n_species = 0;
  network->n_reactions = 0;
*/
  if (verbose >= 1)
    {
      fprintf (stdout, "Reading reactions network from %s... ", chem_file);
      fflush (stdout);
    }

  /* Open the input file. We first look in the current directory, and 
     then in the PKGDATADIR directory. Exit if we can't find it. */
  
  f = fopen (chem_file, "r");
  if ( !f )
    {
      char chem_file1[MAX_LINE];
      
      strncpy (chem_file1, PKGDATADIR, sizeof (chem_file1) - 1);
      strncat (chem_file1, "/", sizeof (chem_file1) - strlen (chem_file1) - 1);
      strncat (chem_file1, chem_file, sizeof (chem_file1) - strlen (chem_file1) - 1);
      f = fopen (chem_file1, "r");
      if ( !f )
	{
	  fprintf (stderr, "astrochem: error: can't find %s.\n", chem_file);
	  exit (1);
	}
    }
  
  /* Loop over the lines, and look for the reactants and products, and
     the parameters of the reactions. */
  
  while (fgets (line, MAX_LINE, f) != NULL)
    {
      line_number++;
      if (line[0] == '#') continue; /* Skip comments. */
       
      /* Initialize reactants, products, and reaction parameters. */
    
      strcpy (reactant1, "");
      strcpy (reactant2, "");
      strcpy (reactant3, "");
      strcpy (product1, "");
      strcpy (product2, "");
      strcpy (product3, "");
      strcpy (product4, "");
      alpha = 0;
      beta = 0;
      gamma = 0;
      reaction_type = 0;
      reaction_no = 0;
      
      /* Read the reactants, products, and reaction parameters. */
      
      char * localLine=line;
      char * localLine2=strchr(line,' ');
      char str [MAX_CHAR_SPECIES];
      int mode=1;
      int ending=0;
      while(localLine!=NULL)
      {
        if(localLine[0]=='-')
        {
          localLine=localLine2+1;
          localLine2=strchr(localLine,' ');
          mode=4;
          ending=0;
          continue;
        }
        else if(localLine[0]=='+')
        {
           localLine=localLine2+1;
           localLine2=strchr(localLine,' ');
           ending=0;
           continue;
        }
        else if(localLine[0]==' ')
        {
            localLine=localLine2+1;
            localLine2=strchr(localLine,' ');
            continue;
        }
        else if(ending == 1)
        {
          if( sscanf (localLine, "%lf %lf %lf %d %d",
		   &alpha, &beta, &gamma, &reaction_type, &reaction_no) != 5) 
           {
	            input_error (chem_file, line_number);
           }
          break;
        }
        else
        {
          strncpy(str,localLine,localLine2-localLine);
          str[localLine2-localLine]='\0';
          if ((strcmp (str, "cosmic-ray") != 0) &&
	          (strcmp (str, "uv-photon") != 0) &&
	          (strcmp (str, "photon") != 0))
          {
            if(mode==1)
            {
                strcpy(reactant1,str);
                mode++;
            }
            else if(mode==2)
            {
                strcpy(reactant2,str);
                mode++;
            }
             else if(mode==3)
            {
                strcpy(reactant3,str);
                mode++;
            }
            else if(mode==4)
            {
                strcpy(product1,str);
                mode++;
            }
             else if(mode==5)
            {
                strcpy(product2,str);
                mode++;
            }
            else if(mode==6)
            {
                strcpy(product3,str);
                mode++;
            }
             else if(mode==7)
            {
                strcpy(product4,str);
                mode++;
            }
            else
            {
	            input_error (chem_file, line_number);
                break;
            }
            add_specie (str, network);
          }
          localLine=localLine2+1;
          localLine2=strchr(localLine,' ');
          ending = 1;
        }

      }
      /* Fill the array of reactions. Exit if of the reactant and
	 product is not in the specie array. */
      
      if (n < network->n_reactions) 
	{
	  if (((network->reactions[n].reactant1 = 
		specie_index (reactant1,  (const char * const *) network->species, network->n_species)) == -2) ||
	      ((network->reactions[n].reactant2 = 
		specie_index (reactant2,  (const char *const *) network->species, network->n_species)) == -2) ||
	      ((network->reactions[n].reactant3 = 
		specie_index (reactant3,  (const char *const *) network->species, network->n_species)) == -2) ||
	      ((network->reactions[n].product1 = 
		specie_index (product1,  (const char *const *) network->species, network->n_species)) == -2)  ||
	      ((network->reactions[n].product2 = 
		specie_index (product2,  (const char *const *) network->species, network->n_species)) == -2)  ||
	      ((network->reactions[n].product3 = 
		specie_index (product3,  (const char *const *) network->species, network->n_species)) == -2)  ||
	      ((network->reactions[n].product4 = 
		specie_index (product4,  (const char *const *) network->species, network->n_species)) == -2))
	    {
	      fprintf (stderr, "astrochem: %s:%d: can't find specie index.\n",
		       __FILE__, __LINE__); 
	      exit(1);
	    }
	  network->reactions[n].alpha = alpha;
	  network->reactions[n].beta = beta;
	  network->reactions[n].gamma = gamma;
	  network->reactions[n].reaction_type = reaction_type;
	  network->reactions[n].reaction_no = reaction_no;
	  n++;
	}
      else
	{
	  fprintf (stderr, "astrochem: error: incorect number of reactions exceed %i," 
              "file %s may be corrupt.\n", MAX_REACTIONS,chem_file);
	  exit(1);
	}
    }
    if(n!=network->n_reactions)
    {
      fprintf (stderr, "astrochem: error: incorect number of reactions, different from %i," 
              "file %s may be corrupt.\n", MAX_REACTIONS,chem_file);
	   exit(1);
    }
    realloc_network_species(network,network->n_species);
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
add_specie (char *new_specie, net_t* network)
{
  int i;
  if (strcmp (new_specie, "") == 0)
    return;
  for (i = 0; i < network->n_species; i++)
    {
      if (strcmp (network->species[i], new_specie) == 0)
	return i;
    }
  i = network->n_species;
  while( i >= network->n_alloc_species )
  {
      if(i==MAX_SPECIES)
      {
         fprintf (stderr, "astrochem: error: the number of species in the chemical"
	       "network file exceeds %i.\n", MAX_SPECIES);
        exit (1); 
      }
      int new_alloc_size = network->n_alloc_species*2;
      if(new_alloc_size > MAX_SPECIES )
      {
          new_alloc_size = MAX_SPECIES;
      }
      realloc_network_species(network,new_alloc_size);
  }
  if ((network->species[i] = malloc (sizeof (char) * MAX_CHAR_SPECIES)) == NULL)
	{
	  fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__, 
		   "array allocation failed.\n");
	  exit (1);
	}
      if (strlen (new_specie) < MAX_CHAR_SPECIES - 1)
	{
	  strcpy (network->species[i], new_specie);
	  (network->n_species)++;
	}
      else
	{
	  fprintf (stderr, "astrochem: error: the number of characters of some "
		   "species of the chemical network file exceeds %i.\n",
		   MAX_CHAR_SPECIES);
	  exit (1);
	}
  return;
}

/* 
  Look up the index of a given specie in the species array.
*/

int 
specie_index (const char *specie, const char * const species[], int n_species)
{
  int i;
  
  /* Return -1 if the specie name is empty. */
  if (strcmp (specie, "") == 0)
    {
      return -1;
    }
  for (i = 0; i < n_species; i++)
    {
      if (strncmp (species[i], specie, 
		   sizeof (char) * MAX_CHAR_SPECIES) == 0)
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
alloc_network ( net_t * network, int n_species, int n_reactions )
{
  network->n_alloc_species = n_species;
  network->n_species = 0;
  if ((network->species = malloc (sizeof(char*) * n_species )) == NULL )
  {
      fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__, 
		   "array allocation failed.\n");
	  exit (1);
  }

  network->n_reactions = n_reactions;
  if ((network->reactions = malloc (sizeof( react_t) * n_reactions )) == NULL )
  {
     fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__, 
		   "array allocation failed.\n");
	 exit (1);
  }
  int i;
  for (i=0; i<n_species; i++)
  {
    network->species[i] = NULL;
  }
}

void
realloc_network_species ( net_t * network, int n_species )
{
  if(n_species < network->n_species)
  {
 	  fprintf (stderr, "astrochem: %s:%d: %s cannot realloc over existing species : "
              "n_species : %i, new_size: %i\n", __FILE__, __LINE__,  network->n_species, n_species );
	  exit (1);
  }
  network->n_alloc_species = n_species;
  if(( network->species = realloc (network->species,sizeof(char*) * n_species )) == NULL )
  {
      fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__, 
		   "array allocation failed.\n");
	  exit (1);
  }
  int i;
  for (i=network->n_species; i<n_species; i++)
  {
    network->species[i] = NULL;
  }
}

/*
  Free the network structure.
*/
void
free_network (net_t * network)
{
  int i;
  for (i=0; i<network->n_species; i++)
  {
    if(network->species[i] != NULL )
    {
      free (network->species[i]);
    }
  }
  free(network->reactions);
  free(network->species);
}
