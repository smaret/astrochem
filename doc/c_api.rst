=========================
Astrochem C API reference
=========================

.. default-domain:: c

This document gives a comprehensive list of the :ref:`types <sec-c-types>` and
:ref:`functions <sec-c-functions>` that are defined in Astrochem C
API. These types and functions are defined in the `libastrochem.h`
header file. For an example of how to use the C API, see
:ref:`calling Astrochem from C <sec-calling-astrochem-from-c>`.

.. _sec-c-types:

Types
=====

.. c:type:: astrochem_mem_t

   Astrochem memory structure

.. c:type:: cell_t

   Structure containing the parameters of a gas cell

   .. c:member:: double av

      Visual extinction in the cell, in magnitudes
      
   .. c:member:: double nh
		 
      H nuclei density in the cell, in :math:`\mathrm{cm^{-3}}`
      
   .. c:member:: double tgas
		 
      Gas temperature in the cell, in K
      
   .. c:member:: double tdust

      Dust temperature in the cell, in K

.. c:type:: net_t

   Structure containing a chemical network

   .. c:member:: int n_species
		 
      Number of species in the network
      
   .. c:member:: int n_alloc_species
		 
      Number of allocated species in the network structure
      
   .. c:member:: species_t *species

      Structure containing the species in the network
      
   .. c:member:: int n_reactions
		 
      Number of reactions
      
   .. c:member:: react_t *reactions
		 
      Structure containing the reactions in the network

.. c:type:: phys_t
   
   Structure containing the physical parameters

  .. c:member:: double chi
		
     External UV radiation field, in Draine units
     
  .. c:member:: double cosmic
		
     Cosmic ray ionization rate of molecular hydrogen, in :math:`\mathrm{s^{-1}}`
     
  .. c:member:: double grain_size
		
     The grain radius, in micron
     
  .. c:member:: double grain_abundance
		
     The grain abundance
     
  .. c:member:: double grain_gas_mass_ratio
		
     The grain-to-gas mass ratio
     
  .. c:member:: double grain_mass_density
     
     The grain mass density in :math:`\mathrm{kg \, m^{-3}}`

.. _sec-c-functions:

Functions
=========

.. c:function:: int alloc_abundances (const net_t* network, double** abundances)

   Allocate an array to store the abundances for all species in a
   network

   :param net_t* network: Network structure
   :param double** abundances: Pointer on the abundance array
   :return: EXIT_SUCCESS if the allocation was successful, EXIT_FAILURE otherwise

.. c:function:: void free_abundances (double* abundances)

   Free the array containing the abundances

   :param double** abundances: Pointer on the abundance array
   
.. c:function:: int set_initial_abundances(const char** species, int n_initialized_abundances, const double* initial_abundances, const net_t* network, double* abundances)

   Set the initial abundances

   :param char** species: Array containing the species name
   :param int n_initialized_abundances: Number of initial abundances
   :param double* initial_abundances: Array containing the initial abundances
   :param net_t* network: Network structure
   :param double* abundances: Array containing the abundances of all species
   :return: 0

.. c:function:: int solver_init (const cell_t* cell, const net_t* network, const phys_t* phys, const double* abundances , double density, double abs_err, double rel_err, astrochem_mem_t* astrochem_mem)

   Initialize the solver

   :param cell_t* cell: Cell structure
   :param net_t* network: Network structure
   :param phys_t* phys: Physical parameters structure
   :param double* abundances: Array containing the abundances of all species
   :param double density: Initial density, in :math:`\mathrm{cm^{-3}}`
   :param double abs_err: Solver absolute error on the abundances
   :param double rel_err: Solver relative error on the abundances
   :param astrochem_mem_t* astrochem_mem: Astrochem memory structure
   :return: 0 if the initialization was successful, -1 otherwise.	 

.. c:function:: int solve(astrochem_mem_t* astrochem_mem, const net_t* network, double* abundances, double time , const cell_t* new_cell, int verbose)

   Solve the system of ODE

   This function solve the system of ODE up to a given time, and
   update the abundance array. If the physical parameters in the gas
   cell have changed since the last call, a pointer to cell structure
   must be passed to the function; if not, a null pointer must be
   passed instead.

   :param astrochem_mem_t* astrochem_mem: Astrochem memory structure
   :param net_t* network: Network structure
   :param double* abundances: Array containing the abundances of all species
   :param double time: Time, in seconds
   :param cell_t* new_cell: New cell structure if the physical parameters have changed since the last call
   :param int verbose: Verbosity (0 for quiet, 1 for verbose)
   :return: 0

.. c:function:: void solver_close(astrochem_mem_t* astrochem_mem)

   Close the solver

   :param astrochem_mem_t* astrochem_mem: Astrochem memory structure

.. c:function:: int read_network (const char* chem_file, net_t* network, const int verbose)

   Read a chemical network

   This function reads a chemical network in Astrochem format (.chm)
   and creates a network structure containing all the reactions.

   :param char* chem_file: Network filename
   :param net_t* network: Network structure
   :param int verbose: Verbosity (0 for quiet, 1 for verbose)
   :return: EXIT_SUCCESS after a successful call, EXIT_FAILURE otherwise
		
.. c:function:: void free_network (net_t * network)

   Free a chemical network structure 

   :param net_t* network: Network structure
