/*
   input.h - private input Function prototypes, various constant and data
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

/* Various definitions and constants */

#ifndef _INPUT_H_
#define _INPUT_H_

#include "astrochem.h"

int read_input (const char *input_file, inp_t * input_params,
                 const net_t * network, int verbose);

int read_input_file_names (const char *input_file, files_t * files,
                            int verbose);

void free_input (inp_t * input_params);


int alloc_input (inp_t * input_params, int n_initial_abundances,
                  int n_output_abundances);

int read_source (const char *source_file, mdl_t * source_mdl,
                  const inp_t * input_params, const int verbose);

void free_mdl (mdl_t * source_mdl);

int get_nb_active_line (const char *file);

#endif // _INPUT_H_
