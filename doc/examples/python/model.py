"""
Example model for Astrochem, Python version
See the "From Python" section in the documentation manual for details.

"""

import astrochem.wrapper as ac
import numpy as np

# Set the physical parameters
p = ac.Phys()
p.cosmic = 1.3e-17
p.chi = 1.0

# Set the initial abundances
initial_abundances = {
    "H2": 1.0,
    "He": 0.14,
    "N": 2.14e-5,
    "O": 1.76e-4,
    "C(+)": 7.30e-5,
    "S(+)": 8.00e-8,
    "Si(+)": 8.00e-9,
    "Fe(+)": 3.00e-9,
    "Na(+)": 2.00e-9,
    "Mg(+)": 7.00e-9,
    "P(+)": 2.00e-10,
    "Cl(+)": 1.00e-9,
    "F": 6.68e-9,
    "e(-)": 7.31012e-5,
}

# Set the density, visual extinction and temperature in the cell
density = 1e4
av = 20
tgas = 10
c = ac.Cell(av, density, tgas, tgas)

# Initialize the solver
verbose = 0
abs_err = 1e-15
rel_err = 1e-6
s = ac.Solver(
    c, "osu2009.chm", p, abs_err, rel_err, initial_abundances, density, verbose
)

# Set the time
time_steps = 32
ti = 1e-6  # years
tf = 1e7  # years
time_yr = np.logspace(np.log10(ti), np.log10(tf), time_steps)  # years
time_sec = time_yr * 365.25 * 24 * 3600  # seconds

print("-------------------------------------------------------------")
print("| Time (yr) | X(CO)     | X(H3(+))  | X(e(-))   | X(HCO(+)) |")
print("-------------------------------------------------------------")

# Solve the chemistry at each time step
for i in range(time_steps):
    abundances = s.solve(time_sec[i])
    print(
        "| %8.2e  | %8.2e  | %8.2e  | %8.2e  | %8.2e  |"
        % (
            time_yr[i],
            abundances["CO"],
            abundances["H3(+)"],
            abundances["e(-)"],
            abundances["HCO(+)"],
        )
    )

print("-------------------------------------------------------------")
