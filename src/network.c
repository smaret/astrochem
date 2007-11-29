/*
  Read the chemical network file.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "astrochem.h"

void add_specie (char *new_specie, char *species[], 
		 int *n_species);

void 
read_network (char *chem_file, struct react reactions[],
	      int *n_reactions, char *species[],
	      int *n_species, int verbose)
{
  FILE *f, *f1;
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
  
  char s1[512], s2[512];

  *n_species = 0;
  *n_reactions = 0;

  if (verbose >= 1)
    {
      fprintf (stdout, "Reading reactions network from %s... ", chem_file);
      fflush (stdout);
    }

  /* Open the input file. We first look in the current directory, and 
     then in the DATADIR directory. Exit if we can't find it. */
  
  f = fopen (chem_file, "r");
  if ( !f )
    {
      char chem_file1[MAX_LINE];
      
      strncpy (chem_file1, DATADIR, sizeof (chem_file1));
      strncat (chem_file1, "/", sizeof (chem_file1));
      strncat (chem_file1, chem_file, sizeof (chem_file1));
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
      
      if ((sscanf (line, "%s -> %s %lf %lf %lf %d %d",
		   reactant1, product1, 
		   &alpha, &beta, &gamma, &reaction_type, &reaction_no) == 7) ||
	  (sscanf (line, "%s + %s -> %s %lf %lf %lf %d %d",
		   reactant1, reactant2, product1, 
		   &alpha, &beta, &gamma, &reaction_type, &reaction_no) == 8) || 
	  (sscanf (line, "%s + %s -> %s + %s %lf %lf %lf %d %d",
		   reactant1, reactant2, product1, product2,
		   &alpha, &beta, &gamma, &reaction_type, &reaction_no) == 9) ||
	  (sscanf (line, "%s + %s -> %s + %s + %s %lf %lf %lf %d %d",
		   reactant1, reactant2, product1, product2, product3,
		   &alpha, &beta, &gamma, &reaction_type, &reaction_no) == 10) ||
	  (sscanf (line, "%s + %s -> %s + %s + %s + %s %lf %lf %lf %d %d",
		   reactant1, reactant2, product1, product2, product3, product4,
		   &alpha, &beta, &gamma, &reaction_type, &reaction_no) == 11) ||
	  (sscanf (line, "%s + %s + %s -> %s + %s %lf %lf %lf %d %d",
		   reactant1, reactant2, reactant3, product1, product2,
		   &alpha, &beta, &gamma, &reaction_type, &reaction_no) == 10))
	;
      else 
	{
	  input_error (chem_file, line_number);
	}

      /* Ignore the following species: cosmic-ray, uv-photon,
	 photon. Replace them by an empty string, and re-sort
	 species. */
      
      if ((strcmp (reactant1, "cosmic-ray") == 0) ||
	  (strcmp (reactant1, "uv-photon") == 0) ||
	  (strcmp (reactant1, "photon") == 0))
	{
	  strcpy (reactant1, reactant2);
	  strcpy (reactant2, reactant3);
	  strcpy (reactant3, "");
	}
      if ((strcmp (reactant2, "cosmic-ray") == 0) ||
	  (strcmp (reactant2, "uv-photon") == 0) ||
	  (strcmp (reactant2, "photon") == 0))
	{
	  strcpy (reactant2, reactant3);
	  strcpy (reactant3, "");
	}
      if ((strcmp (reactant3, "cosmic-ray") == 0) ||
	  (strcmp (reactant3, "uv-photon") == 0) ||
	  (strcmp (reactant3, "photon") == 0))
	{
	  strcpy (reactant3, "");
	}
      if ((strcmp (product1, "cosmic-ray") == 0) ||
	  (strcmp (product1, "uv-photon") == 0) ||
	  (strcmp (product1, "photon") == 0))
	{
	  strcpy (product1, product2);
	  strcpy (product2, product3);
	  strcpy (product3, "");
	}
      if ((strcmp (product2, "cosmic-ray") == 0) ||
	  (strcmp (product2, "uv-photon") == 0) ||
	  (strcmp (product2, "photon") == 0))
	{
	  strcpy (product2, product3);
	  strcpy (product3, "");
	}
      if ((strcmp (product3, "cosmic-ray") == 0) ||
	  (strcmp (product3, "uv-photon") == 0) ||
	  (strcmp (product3, "photon") == 0))
	{
	  strcpy (product3, "");
	}
      if ((strcmp (product4, "cosmic-ray") == 0) ||
	  (strcmp (product4, "uv-photon") == 0) ||
	  (strcmp (product4, "photon") == 0))
	{
	  strcpy (product4, "");
	}

      /* Fill the array of species. */
    
      add_specie (reactant1, species, n_species);
      add_specie (reactant2, species, n_species);
      add_specie (reactant3, species, n_species);
      add_specie (product1, species, n_species);
      add_specie (product2, species, n_species);
      add_specie (product3, species, n_species);
      add_specie (product4, species, n_species);

      /* Fill the array of reactions. Exit if of the reactant and
	 product is not in the specie array. */
      
      if (*n_reactions < MAX_REACTIONS) 
	{
	  if (((reactions[*n_reactions].reactant1 = 
		specie_index (reactant1, species, *n_species)) == -2) ||
	      ((reactions[*n_reactions].reactant2 = 
		specie_index (reactant2, species, *n_species)) == -2) ||
	      ((reactions[*n_reactions].reactant3 = 
		specie_index (reactant3, species, *n_species)) == -2) ||
	      ((reactions[*n_reactions].product1 = 
		specie_index (product1, species, *n_species)) == -2)  ||
	      ((reactions[*n_reactions].product2 = 
		specie_index (product2, species, *n_species)) == -2)  ||
	      ((reactions[*n_reactions].product3 = 
		specie_index (product3, species, *n_species)) == -2)  ||
	      ((reactions[*n_reactions].product4 = 
		specie_index (product4, species, *n_species)) == -2))
	    {
	      fprintf (stderr, "astrochem: %s:%d: can't find specie index.\n",
		       __FILE__, __LINE__); 
	      exit(1);
	    }
	  reactions[*n_reactions].alpha = alpha;
	  reactions[*n_reactions].beta = beta;
	  reactions[*n_reactions].gamma = gamma;
	  reactions[*n_reactions].reaction_type = reaction_type;
	  reactions[*n_reactions].reaction_no = reaction_no;
	  (*n_reactions)++;
	}
      else
	{
	  fprintf (stderr, "astrochem: error: the number of reactions exceed %i.\n", 
		   MAX_REACTIONS);
	  exit(1);
	}
    }
  
  if (verbose >= 1)
    {
      fprintf (stdout, "done.\n");  
      fprintf (stdout, "Found %d reactions involving %d species.\n", 
	       *n_reactions, *n_species);
    }
 
  /* Close the file. */

  fclose (f);
}

/*
  Add a specie in the species array, if not already present.
*/

void  
add_specie (char *new_specie, char *species[],
	    int *n_species)
{
  int i;
     
  if (strcmp (new_specie, "") == 0)
    return;
  for (i = 0; i < *n_species; i++)
    {
      if (strcmp (species[i], new_specie) == 0)
	return;
    }
  i = *n_species;
  if (i < MAX_SPECIES)
    {
      if ((species[i] = malloc (sizeof (char) * MAX_CHAR_SPECIES)) == NULL)
	{
	  fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__, 
		   "array allocation failed.\n");
	  exit (1);
	}
      if (strlen (new_specie) < MAX_CHAR_SPECIES - 1)
	{
	  strcpy (species[i], new_specie);
	  (*n_species)++;
	}
      else
	{
	  fprintf (stderr, "astrochem: error: the number of characters of some "
		   "species of the chemical network file exceeds %i.\n",
		   MAX_CHAR_SPECIES);
	  exit (1);
	}
    } 
  else
    {
      fprintf (stderr, "astrochem: error: the number of species in the chemical"
	       "network file exceeds %i.\n", MAX_SPECIES);
      exit (1); 
    }
  
  return;
}

/* 
  Look up the index of a given specie in the species array.
*/

int 
specie_index (char *specie, char *species[], int n_species)
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
