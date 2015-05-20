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
#ifndef _ASTROCHEM_H_
#define _ASTROCHEM_H_

#include <hdf5.h>
#include "libastrochem.h"

#define MAX_CHAR_FILENAME 64

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
  cell_table_t *cell;
  time_steps_t ts;
  int n_cells;
  SOURCE_MODE mode;
} mdl_t;


int full_solve (hid_t fid, hid_t dataset, hid_t* routeDatasets, hid_t dataspace, hid_t routeDataspace, hid_t datatype, hid_t routeDatatype,
                int cell_index, const inp_t * input_params, SOURCE_MODE mode,
                const cell_table_t * cell, const net_t * network, const time_steps_t * ts, int verbose);

#endif // _ASTROCHEM_H_
