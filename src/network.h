/*
   network.h - private network Function prototypes, various constant and data
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

#ifndef _NETWORK_H_
#define _NETWORK_H_

#define MAX_CHAR_ELEMENT 8                /* Maximum number of characters in a element name */
#define ELECTRON_MASS    9.10938291e-31   /* Mass of an electron in kg */
#define UMA              1.660538921e-27  /* Atomic mass unity in kg */

typedef char elem_name_t[MAX_CHAR_ELEMENT];

typedef struct
{
  elem_name_t name; /* Name of the element */
  int mass;         /* Mass of the element in UMA */
} elem_t;


int find_species (const species_name_t specie, const net_t * network);
int add_species (char *new_species, net_t * network);
int alloc_network (net_t * network, int n_species, int n_reactions);
int realloc_network_species (net_t * network, int n_species);
bool get_species_mass_and_charge( char* species, double* mass, int* charge );
int get_element_mass( const char* element );

#endif // _NETWORK_H_
