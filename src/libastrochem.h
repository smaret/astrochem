/*
   libastrochem.h - Function prototypes, various constant and data
   structures for Astrochem library.

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

/* Various definitions and constants */

#ifndef _LIBASTROCHEM_H_
#define _LIBASTROCHEM_H_

#include <nvector/nvector_serial.h>
#include <string.h>

#define MAX_CHAR_SPECIES 32   /* Maximum number of characters in a specie name */
#define MAX_REACTANTS 3
#define MAX_PRODUCTS 4

/* Data structures */

typedef char species_name_t[MAX_CHAR_SPECIES];

typedef struct
{
  double chi;
  double cosmic;
  double grain_size;
  double grain_abundance;
  double grain_gas_mass_ratio;
  double grain_mass_density;
} phys_t;

typedef struct
{
  double av;
  double nh;
  double tgas;
  double tdust;
} cell_t;

typedef struct
{
  int reactants[MAX_REACTANTS];
  int products[MAX_PRODUCTS];
  double alpha;
  double beta;
  double gamma;
  int reaction_type;
  int reaction_no;
} react_t;

typedef struct
{
  species_name_t name;
  double mass;
  int charge;
} species_t;

typedef struct
{
  int n_species;
  int n_alloc_species;
  species_t *species;
  int n_reactions;
  react_t *reactions;
} net_t;

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

int alloc_abundances (const net_t* network, double** abundances);

void free_abundances (double* abundances);

int set_initial_abundances (const char** species, int n_initialized_abundances,
                            const double* initial_abundances, const net_t* network,
			    double* abundances );

int solver_init (const cell_t* cell, const net_t* network, const phys_t* phys,
                 const double* abundances , double density, double abs_err, double rel_err,
                 astrochem_mem_t* astrochem_mem);

int solve (astrochem_mem_t* astrochem_mem, const net_t* network,
           double* abundances, double time , const cell_t* new_cell, int verbose);

void solver_close (astrochem_mem_t* astrochem_mem);

int read_network (const char *chem_file, net_t * network, const int verbose);

void free_network (net_t * network);

#endif // _LIBASTROCHEM_H_
