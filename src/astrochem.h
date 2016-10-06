/*
   astrochem.h - Function prototypes, various constant and data
   structures for Astrochem.

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

#ifndef _ASTROCHEM_H_
#define _ASTROCHEM_H_

#include <hdf5.h>
#include "libastrochem.h"

/* Various definitions and constants */

#define MAX_CHAR_FILENAME 64            /* Maximum number of characters in file names */
#define MAX_LINE 512                    /* Maximum number of characters in each input file line */
#define CHI_DEFAULT 1                   /* Default chi value */
#define COSMIC_DEFAULT 1.3e-17          /* Default cosmic value */
#define GRAIN_SIZE_DEFAULT 1e-5         /* Default grain radius, in cm */
#define GRAIN_GAS_MASS_RATIO_DEFAULT 0  /* Default grain mass ratio */
#define GRAIN_MASS_DENSITY_DEFAULT 3    /* Default grain mass density, Olivine grains, g/cm3  */
#define TI_DEFAULT 1e-6                 /* Default initial time */
#define TF_DEFAULT 1e7                  /* Default final time */
#define ABS_ERR_DEFAULT 1e-20           /* Default absolute error */
#define REL_ERR_DEFAULT 1e-3            /* Default relative error */
#define TIME_STEPS_DEFAULT 32           /* Default number of times steps */
#define TRACE_ROUTES_DEFAULT 0          /* Deactivate route tracing by default */
#define N_OUTPUT_ROUTES 16              /* Defaults number of output routes */
#define CVODE_MXSTEPS   1e6             /* Maximum number of time steps in CVODE */

#ifndef M_PI
#define M_PI  3.14159265358979323846264338327950288
#endif

#define CONST_MKSA_YEAR 3.1536e7                /* Number of seconds in a year */
#define CONST_CGSM_BOLTZMANN (1.3806503e-16)    /* Boltzmann constant */
#define CONST_CGSM_MASS_PROTON (1.67262158e-24) /* Proton Mass */
#define MASS_PROTON       1.672621777e-27       /* Proton Mass */

#define MIN_ABUNDANCE 1e-20   /* Minimum abundance to write in output files */

#define FRACTION_TIME_GRAIN_70K 3.16e-19
#define GAS_DUST_NUMBER_RATIO 7.57e+11
#define GRAIN_SITES_PER_CM2 3.00e+15    /* cm-2 */
#define DRAINE_STANDARD_ISRF_FUV 1.7e8  /* photons cm-2 */

/* Data structures */

typedef enum
{ STATIC = 0, DYNAMIC = 1 } SOURCE_MODE;

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
  int *output_species_idx;
  int n_output_species;
  int time_steps;
  int trace_routes;
  char suffix[MAX_LINE];
} output_t;

typedef struct
{
  char chem_file[MAX_LINE];
  char source_file[MAX_LINE];
} files_t;

typedef struct
{
  int species_idx;
  double abundance;
} abund_t;

typedef struct
{
  abund_t *initial_abundances;
  int n_initial_abundances;
} abundances_t;

typedef struct
{
  double ti;
  double tf;
  double abs_err;
  double rel_err;
} solver_t;

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
  double *av;
  double *nh;
  double *tgas;
  double *tdust;
} cell_table_t;

typedef struct
{
  double *time_steps;
  int n_time_steps;
} time_steps_t;

typedef struct
{
  cell_table_t *cell;
  time_steps_t ts;
  int n_cells;
  SOURCE_MODE mode;
} mdl_t;

/* Fonction prototypes */

int full_solve (hid_t fid, hid_t dataset, hid_t* routeDatasets, hid_t dataspace, hid_t routeDataspace, hid_t datatype, hid_t routeDatatype,
                int cell_index, const inp_t * input_params, SOURCE_MODE mode,
                const cell_table_t * cell, const net_t * network, const time_steps_t * ts, int verbose);

int get_nb_active_line (const char *file);

#endif // _ASTROCHEM_H_

