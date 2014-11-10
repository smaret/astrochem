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

/**
 * @brief struct containing a rates of a reaction
 */
typedef struct
{
  int reaction_no; /*!< number of concerned reaction */
  double rate;     /*!< rate value */
} r_t;

/**
 * @brief struct containing a route
 */
typedef struct
{
  r_t destruction; /*!< rate of destruction */
  r_t formation;   /*!< rate of formation */
} rout_t;


/**
 * @brief struct containing output parameters
 */
typedef struct
{
  int *output_species_idx; /*!< array of output species idx */
  int n_output_species;    /*!< number of output species */
  int time_steps;          /*!< time steps used */
  int trace_routes;        /*!< If routes have been traced */
  char suffix[MAX_LINE];   /*!< Suffix to sue in output files */
} output_t;

/**
 * @brief struct containing file names of chem file and source file
 */
typedef struct
{
  char chem_file[MAX_LINE]; /*!< Path to chem file */
  char source_file[MAX_LINE]; /*!< Path to source file */
} files_t;

/**
 * @brief struct containing input parametrs
 */
typedef struct
{
  files_t files;    /*!< Input files */
  phys_t phys;      /*!< Physics parameters */
  solver_t solver;  /*!< Solver parameters */
  abundances_t abundances; /*!< abundances */
  output_t output; /*!< Output parameters */
} inp_t;


/**
 * @brief struct containing cell parameters
 */
typedef struct
{
  double *av;    /*!< av */
  double *nh;    /*!< density */
  double *tgas;  /*!< gas temperature */
  double *tdust; /*!< dust temperature */
} cell_table_t;

/**
 * @brief struct containing a source model
 */
typedef struct
{
  cell_table_t *cell;    /*!< Array of cells */
  time_steps_t ts; /*!< Time steps */
  int n_cells;     /*!< Number of cells */
  SOURCE_MODE mode; /*!< Source mode */
} mdl_t;


int full_solve (hid_t fid, hid_t dataset, hid_t* routeDatasets, hid_t dataspace, hid_t routeDataspace, hid_t datatype, hid_t routeDatatype,
                int cell_index, const inp_t * input_params, SOURCE_MODE mode,
                const cell_table_t * cell, const net_t * network, const time_steps_t * ts, int verbose);

#endif // _ASTROCHEM_H_
