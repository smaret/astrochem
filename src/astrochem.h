/* 
   astrochem.h - Function prototypes, various constant and data
   structures for Astrochem.

   Copyright (c) 2006-2009 Sebastien Maret

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

#define MAX_LINE 512                 /* Maximum number of characters in each input file
					line */

#define CHI_DEFAULT 1
#define COSMIC_DEFAULT 1.3e-17
#define GRAIN_SIZE_DEFAULT 1e-5     /* Grain radius, in cm */
#define TI_DEFAULT 1e-7
#define TF_DEFAULT 1e7
#define ABS_ERR_DEFAULT 1e-15
#define REL_ERR_DEFAULT 1e-3
#define TIME_STEPS_DEFAULT 32
#define TRACE_ROUTES_DEFAULT 0
#define N_OUTPUT_ROUTES 16

#define MAX_INITIAL_ABUNDANCES 128  /* Maximum initial abundances in the input file */
#define MAX_CHAR_SPECIES 32         /* Maximum number of characters in a specie name */
#define MAX_OUTPUT_ABUNDANCES 32    /* Maximum output abundances in the output file */
#define MAX_TIME_STEPS 128          /* Maximum number of time steps the output file */

#define CONST_MKSA_YEAR 3.1536e7
#define CONST_CGSM_BOLTZMANN (1.3806503e-16)
#define CONST_CGSM_MASS_PROTON (1.67262158e-24)

#define MIN_ABUNDANCE 1e-20        /* Minimum abundance to write in output files */

struct abund {
  char specie[MAX_CHAR_SPECIES];
  double abundance;
};

#define MAX_SHELLS 256

struct spec {
  char specie[MAX_CHAR_SPECIES];
  int index;
};

struct react {
  int reactant1;
  int reactant2;
  int reactant3;
  int product1;
  int product2;
  int product3;
  int product4;
  double alpha;
  double beta;
  double gamma;
  int reaction_type;
  int reaction_no;
};

struct r {
  int reaction_no;
  double rate;
};

struct rout {
  struct r destruction;
  struct r formation;
};

#define MAX_REACTIONS 8192
#define MAX_SPECIES 1024

void read_input (char *input_file, char *chem_file, char *source_file,
		 double *chi, double *cosmic, double *grain_size,
		 double *grain_abundance, double *ti, double *tf, double *abs_err,
		 double *rel_err, struct abund initial_abundances[],
		 int *n_initial_abundances, char *output_species[],
		 int *n_output_abundances, int *time_steps, 
		 int *trace_routes, char *suffix, int verbose); 

void read_source (char *source_f, int shell[], int *n_shells,
		  double av[], double nh[], double tgas[],
		  double tdust[], int verbose);

void input_error (char *input_f, int line_number);

void check_species (struct abund initial_abundances[], int
		    n_initial_abundances, char *output_species[], int
		    n_output_species, char *species[], int n_species);

void read_network (char *chem_file, struct react reactions[],
		   int *n_reactions, char *species[],
		   int *n_species, int verbose);

int specie_index (char specie[], char *species[], int n_species);

double rate(double alpha, double beta, double gamm, int reaction_type,
	    int reaction_no, double nh, double av, double tgas, double tdust,
	    double chi, double cosmic, double grain_size,
	    double grain_abundance);
  
int solve (double chi, double cosmic, double grain_size, double grain_abundance,
	   double abs_err, double rel_err,
	   struct abund initial_abundances[],
	   int n_initial_abundances, char *output_species[],
	   int n_output_species, double av, double nh,
	   double tgas, double tdust,
	   struct react reactions[], int n_reactions, 
	   char *species[], int n_species,
	   int shell_index, double tim[], int time_steps,
	   double abundances[MAX_SHELLS][MAX_TIME_STEPS][MAX_OUTPUT_ABUNDANCES],
	   int trace_routes, 
	   struct rout routes[MAX_SHELLS][MAX_TIME_STEPS][MAX_OUTPUT_ABUNDANCES][N_OUTPUT_ROUTES],
	   int verbose);

void output (int n_shells, double tim[], int time_steps,
	     char *output_species[], int n_output_species,
	     double abundances[MAX_SHELLS][MAX_TIME_STEPS][MAX_OUTPUT_ABUNDANCES],
	     char *species[], int n_species, int trace_routes, 
	     struct rout routes[MAX_SHELLS][MAX_TIME_STEPS][MAX_OUTPUT_ABUNDANCES][N_OUTPUT_ROUTES],
	     char *suffix, int verbose);
