/* 
   astrochem.h - Function prototypes, various constant and data
   structures for Astrochem.

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

/* Various definitions and constants */

#define MAX_LINE 512                 /* Maximum number of characters in each input file
					line */

#define CHI_DEFAULT 1
#define COSMIC_DEFAULT 1.3e-17
#define GRAIN_SIZE_DEFAULT 1e-5     /* Grain radius, in cm */
#define TI_DEFAULT 1e-6
#define TF_DEFAULT 1e7
#define ABS_ERR_DEFAULT 1e-20
#define REL_ERR_DEFAULT 1e-3
#define TIME_STEPS_DEFAULT 32
#define TRACE_ROUTES_DEFAULT 0
#define N_OUTPUT_ROUTES 16

#define MAX_INITIAL_ABUNDANCES 128  /* Maximum initial abundances in the input file */
#define MAX_CHAR_SPECIES 32         /* Maximum number of characters in a specie name */
#define MAX_SHELLS 256              /* Maximum number of shells in the model file */
#define MAX_OUTPUT_ABUNDANCES 32    /* Maximum output abundances in the output file */
#define MAX_TIME_STEPS 128          /* Maximum number of time steps the output file */
#define MAX_REACTIONS 32768         /* Maximum number of reactions in the network file */
#define MAX_SPECIES 4096            /* Maximum number of species in the network file */

#define CONST_MKSA_YEAR 3.1536e7
#define CONST_CGSM_BOLTZMANN (1.3806503e-16)
#define CONST_CGSM_MASS_PROTON (1.67262158e-24)

#define MIN_ABUNDANCE 1e-20        /* Minimum abundance to write in output files */

#define FRACTION_TIME_GRAIN_70K 3.16e-19 
#define GAS_DUST_NUMBER_RATIO 7.57e+11
#define GRAIN_SITES_PER_CM2 3.00e+15 /* cm-2 */
#define AVERAGE_UV_IRSF 1e8          /* photons cm-2 */

/* Data structures */

typedef struct {
  int specie_idx;
  double abundance;
} abund_t;

typedef struct {
  char chem_file[MAX_LINE];
  char source_file[MAX_LINE];
} files_t;

typedef struct {
  double chi;
  double cosmic;
  double grain_size;
  double grain_abundance;
} phys_t;

typedef struct {
  double ti;
  double tf;
  double abs_err;
  double rel_err;
} solver_t;

typedef struct {
  abund_t * initial_abundances;
  int n_initial_abundances;
} abundances_t;

typedef struct {
  int * output_species_idx;
  int n_output_species;
  int time_steps;
  int trace_routes;
  char suffix[MAX_LINE];
} output_t;

typedef struct {
  files_t files;
  phys_t phys;  
  solver_t solver;
  abundances_t abundances;
  output_t output;
} inp_t;

typedef struct {
    double av;
    double nh;
    double tgas;
    double tdust;
  } shell_t;

typedef struct {
  shell_t * shell;
  int n_shells;
} mdl_t;


typedef struct {
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
} react_t;

typedef struct {
  int n_species;
  int n_alloc_species;
  char ** specie_names;
  int n_reactions;
  react_t * reactions;
} net_t;

typedef struct {
  int reaction_no;
  double rate;
} r_t;

typedef struct {
  r_t destruction;
  r_t formation;
} rout_t;

typedef struct {
  double * abundances;
  rout_t * routes;
  double * tim;
  int n_shells;
  int n_time_steps;
  int n_output_abundances;
} res_t;

/* Fonction prototypes

   Fontions names with a _new suffix correspond to the new prototypes
   that uses the new "input_params", "source_mdl", "network" and
   "results" data structures. */

void read_input (const char *input_file, inp_t *input_params, const net_t * network, int verbose);
void read_input_file_names (const char *input_file, files_t *files, int verbose);

void free_input ( inp_t * input_params );
int get_nb_active_line(const char * file);

void read_source (const char *source_file, mdl_t *source_mdl,
		      const int verbose);
void free_mdl( mdl_t * source_mdl );


void input_error (const char *input_f, int line_number);

void check_species ( abund_t initial_abundances[], int
		    n_initial_abundances, char *output_species[], int
		    n_output_species, char *species[], int n_species);

void read_network (const char *chem_file, net_t *network, const int verbose);
void free_network ( net_t *network );

void alloc_results( res_t * results, int n_time_steps, int n_shells, int n_output_abundances);
void free_results ( res_t * results );

int specie_index (const char specie[],const char * const species[], int n_species);

double rate(double alpha, double beta, double gamm, int reaction_type,
	    int reaction_no, double nh, double av, double tgas, double tdust,
	    double chi, double cosmic, double grain_size,
	    double grain_abundance, double ice_abundance);

int solve (int shell_index, const inp_t *input_params, const shell_t *shell,
	   const net_t *network, res_t *results, int verbose);
int get_abundance_idx( const res_t * results,int shell_idx, int ts_idx, int abund_idx);
int get_route_idx( const res_t * results, int shell_idx, int ts_idx, int abund_idx, int route_idx);
void output (int n_shells, const inp_t *input_params, const net_t *network,
	     const res_t *results, int verbose);
