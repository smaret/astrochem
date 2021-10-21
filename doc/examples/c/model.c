/*
  Example model for Astrochem, C version.

  Compile with:

    % cc -I${prefix}/include -L${prefix}/lib -lastrochem model.c -o model

  where ${prefix} is the Astrochem installation directory (/usr/local
  by default).

  See the "From C" section in the documentation manual for details.
*/

#include <libastrochem.h>
#include <math.h>

int
main ()
{
  // Select the chemical network
  int verbose = 1;
  char *chem_file = "osu2009.chm";
  net_t network;
  read_network (chem_file, &network, verbose);

  // Set the physical parameters
  phys_t phys;
  phys.cosmic = 1.3e-17;
  phys.chi = 1.0;

  // Allocate an array to store the abundances at the set the initial abundance
  double *abundances;
  const char *species[]
      = { "H2",    "He",    "N",     "O",    "C(+)",  "S(+)", "Si(+)",
          "Fe(+)", "Na(+)", "Mg(+)", "P(+)", "Cl(+)", "F",    "e(-)" };
  const double initial_abundances[]
      = { 1.0,     0.14,    2.14e-5, 1.76e-4,  7.30e-5, 8.00e-8, 8.00e-9,
          3.00e-9, 2.00e-9, 7.00e-9, 2.00e-10, 1.00e-9, 6.68e-9, 7.31012e-5 };
  alloc_abundances (&network, &abundances);
  set_initial_abundances (species, 14, initial_abundances, &network,
                          abundances);

  // Set the initial density, visual extinction and temperature of the gas cell
  double density = 1e4;
  double av = 20;
  double temperature = 10;
  cell_t cell;
  cell.nh = density;
  cell.av = av;
  cell.tgas = temperature;
  cell.tdust = temperature;

  // Initialize the solver
  double abs_err = 1e-15;
  double rel_err = 1e-6;
  astrochem_mem_t astrochem_mem;
  solver_init (&cell, &network, &phys, abundances, density, abs_err, rel_err,
               &astrochem_mem);

  // Set the time
  int time_steps = 32;
  double ti = 1e-6; // years
  double tf = 1e7;  // years
  double dt = exp ((log (tf) - log (ti)) / (time_steps - 1));
  double time_year = ti;

  // Find the indexes of the species for which we want to print the abundance
  int spec_index_co = find_species ("CO", &network);
  int spec_index_h3p = find_species ("H3(+)", &network);
  int spec_index_e = find_species ("e(-)", &network);
  int spec_index_hcop = find_species ("HCO(+)", &network);

  printf ("-------------------------------------------------------------\n");
  printf ("| Time (yr) | X(CO)     | X(H3(+))  | X(e(-))   | X(HCO(+)) |\n");
  printf ("-------------------------------------------------------------\n");

  // Solve the chemistry at each time step
  for (int i = 0; i < time_steps; i++)
    {
      double time_sec = time_year * 365.25 * 24 * 3600;
      solve (&astrochem_mem, &network, abundances, time_sec, NULL, verbose);
      printf ("| %8.2e  | %8.2e  | %8.2e  | %8.2e  | %8.2e  |\n", time_year,
              abundances[spec_index_co], abundances[spec_index_h3p],
              abundances[spec_index_e], abundances[spec_index_hcop]);
      time_year *= dt;
    }

  printf ("-------------------------------------------------------------\n");

  return 0;
}
