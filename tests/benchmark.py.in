#!@PYTHON@
#
#  benchmark - Benckmark Astrochem
#
#  Copyright (c) 2006-2014 Sebastien Maret
# 
#  This file is part of Astrochem.
#
#  Astrochem is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  Astrochem is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with Astrochem.  If not, see <http://www.gnu.org/licenses/>.

from subprocess import *

input_ini = """[files]
source = source.mdl
chem = ../networks/osu2009.chm
# Physical paramaters
[phys]
chi = 1.0
cosmic = 1.3e-17
# Solver parameters
[solver]
ti = 1e-6
tf = 1e7
abs_err = 1e-20
rel_err = 1e-6
# Initial abundances
[abundances]
H2      = 0.5
He      = 0.14
N       = 2.14e-5
O       = 1.76e-4
C(+)    = 7.30e-5
S(+)    = 8.00e-8
Si(+)   = 8.00e-8
Fe(+)   = 3.00e-9
Na(+)   = 2.00e-9
Mg(+)   = 7.00e-9
P(+)    = 3.00e-9
Cl(+)   = 4.00e-9
e(-)    = 7.32e-5
# Output
[output]
time_steps = 32
abundances = CO,C(+),C,e(-),OH,H3O(+),H,H2,HCO(+),CO(+),C4H,HCO(+),CH(+),CH
trace_routes = 0
"""

source_mdl = """# Source model file example
# shell number, Av [mag], n(H) [cm^-3], Tgas [K], Tdust [K]
#
0	20.0	1e+04	10.0	10.0
"""

def main():
    """
    Benchmark Astrochem

    """

    # Create the input files

    f = open("input.ini", 'w')
    f.write(input_ini)
    f.close()

    f = open("source.mdl", 'w')
    f.write(source_mdl)
    f.close()

    # Run Astrochem several times

    nrun = 10
    real_avg, user_avg, sys_avg = 0., 0., 0.

    for i in range(nrun):

        # Start time as a subprocess and parse the output
        p = Popen(["time", "-p", "../src/astrochem", "-q", "input.ini"], stderr = PIPE)
        output = p.communicate()[1]
        real, user, sys = map(float, output.split()[1::2])
        print "[run %2i/%2i] real: %6.3fs  user: %6.3fs  sys: %6.3fs" % (i+1, nrun, real, user, sys) 

        # Compute the average of the run execution times
        real_avg += real
        user_avg += user
        sys_avg += sys

    real_avg /= nrun
    user_avg /= nrun
    sys_avg /= nrun

    print "=========================================="
    print "Benchmark results:"
    print "real: %6.3fs  user: %6.3fs  sys: %6.3fs" % (real_avg, user_avg, sys_avg)
    print "=========================================="    

if __name__ == "__main__":
    main()

