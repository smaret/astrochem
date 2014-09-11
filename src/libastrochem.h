/*
   astrochem.h - Function prototypes, various constant and data
   structures for Astrochem.

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

/* Various definitions and constants */

#ifndef _LIBASTROCHEM_H_
#define _LIBASTROCHEM_H_

#include <nvector/nvector_serial.h>

#define MAX_LINE 512            /* Maximum number of characters in each input file
                                   line */
#define CHI_DEFAULT 1
#define COSMIC_DEFAULT 1.3e-17
#define GRAIN_SIZE_DEFAULT 1e-5 /* Grain radius, in cm */
#define TI_DEFAULT 1e-6
#define TF_DEFAULT 1e7
#define ABS_ERR_DEFAULT 1e-20
#define REL_ERR_DEFAULT 1e-3
#define TIME_STEPS_DEFAULT 32
#define TRACE_ROUTES_DEFAULT 0
#define N_OUTPUT_ROUTES 16

#ifndef M_PI
#define M_PI  3.14159265358979323846264338327950288
#endif

#define MAX_CHAR_SPECIES 32     /* Maximum number of characters in a specie name */

#define CONST_MKSA_YEAR 3.1536e7
#define CONST_CGSM_BOLTZMANN (1.3806503e-16)
#define CONST_CGSM_MASS_PROTON (1.67262158e-24)

#define MIN_ABUNDANCE 1e-20     /* Minimum abundance to write in output files */

#define FRACTION_TIME_GRAIN_70K 3.16e-19
#define GAS_DUST_NUMBER_RATIO 7.57e+11
#define GRAIN_SITES_PER_CM2 3.00e+15    /* cm-2 */
#define DRAINE_STANDARD_ISRF_FUV 1.7e8  /* photons cm-2 */

/* Data structures */

typedef enum
{ STATIC = 0, DYNAMIC = 1 } SOURCE_MODE;

typedef char species_name_t[MAX_CHAR_SPECIES];

typedef struct
{
  int species_idx;
  double abundance;
} abund_t;

typedef struct
{
  char chem_file[MAX_LINE];
  char source_file[MAX_LINE];
} files_t;

/**
 * @brief struct containing physics parameters
 */

typedef struct
{
  double chi;
  double cosmic;
  double grain_size;
  double grain_abundance;
} phys_t;

typedef struct
{
  double ti;
  double tf;
  double abs_err;
  double rel_err;
} solver_t;

typedef struct
{
  abund_t *initial_abundances;
  int n_initial_abundances;
} abundances_t;

typedef struct
{
  int *output_species_idx;
  int n_output_species;
  int time_steps;
  int trace_routes;
  char suffix[MAX_LINE];
} output_t;

typedef struct
{
  files_t files;
  phys_t phys;
  solver_t solver;
  abundances_t abundances;
  output_t output;
} inp_t;

typedef struct
{
  double *av;    /*!< av */
  double *nh;    /*!< density */
  double *tgas;  /*!< gas temperature */
  double *tdust; /*!< dust temperature */
} cell_table_t;

/**
 * @brief struct containing cell parameters
 */
typedef struct
{
  double av;    /*!< av */
  double nh;    /*!< density */
  double tgas;  /*!< gas temperature */
  double tdust; /*!< dust temperature */
} cell_t;


/**
 * @brief struct containing array of time steps
 */
typedef struct
{
  double *time_steps;
  int n_time_steps;
} time_steps_t;

typedef struct
{
  cell_table_t *cell;    /*!< Array of cells */
  time_steps_t ts; /*!< Time steps */
  int n_cells;     /*!< Number of cells */
  SOURCE_MODE mode; /*!< Source mode */
} mdl_t;


typedef struct
{
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

typedef struct
{
  int n_species;
  int n_alloc_species;
  species_name_t *species_names;
  int n_reactions;
  react_t *reactions;
} net_t;

typedef struct
{
  int reaction_no;
  double rate;
} r_t;

typedef struct
{
  r_t destruction;
  r_t formation;
} rout_t;

typedef struct
{
  double *abundances;
  rout_t *routes;
  int n_cells;
  int n_time_steps;
  int n_output_abundances;
} res_t;

typedef enum { false, true } bool;

typedef struct
{
  double *reac_rates;
  const react_t *reactions;
  int n_reactions;
  int n_species;
  double nh;
  double av;
  double tgas;
  double tdust;
  double chi;
  double cosmic;
  double grain_size;
  double grain_abundance;
} params_t;

typedef struct
{
  void* cvode_mem;
  N_Vector y;
  params_t params;
  double density;
} astrochem_mem_t;


/* Fonction prototypes */

int alloc_abundances( const net_t* network, double** abundances );

void free_abundances( double* abundances );

int set_initial_abundances( const char** species, int n_initialized_abundances,
                            const double* initial_abundances, const net_t* network, double* abundances );

int solver_init( const cell_t* cell, const net_t* network, const phys_t* phys,
                 const double* abundances , double density, double abs_err, double rel_err,
                 astrochem_mem_t* astrochem_mem );

int solve( astrochem_mem_t* astrochem_mem, const net_t* network,
           double* abundances, double time , const cell_t* new_cell, int verbose );

void solver_close( astrochem_mem_t* astrochem_mem );


void read_input (const char *input_file, inp_t * input_params,
                 const net_t * network, int verbose);

void read_input_file_names (const char *input_file, files_t * files,
                            int verbose);

void free_input (inp_t * input_params);

void read_source (const char *source_file, mdl_t * source_mdl,
                  const inp_t * input_params, const int verbose);

void free_mdl (mdl_t * source_mdl);

/*void check_species (abund_t initial_abundances[], int
                    n_initial_abundances, char *output_species[], int
                    n_output_species, char *species[], int n_species);*/

void read_network (const char *chem_file, net_t * network, const int verbose);

void free_network (net_t * network);

void alloc_results (res_t * results, int n_time_steps, int n_cells,
                    int n_output_abundances);

void free_results (res_t * results);

int full_solve (int cell_index, const inp_t * input_params, SOURCE_MODE mode,
                const cell_table_t * cell, const net_t * network,
                const time_steps_t * ts, res_t * results, int verbose);


void output (int n_cells, const inp_t * input_params,
             const mdl_t * source_mdl, const net_t * network,
             const res_t * results, int verbose);


#endif // _LIBASTROCHEM_H_
