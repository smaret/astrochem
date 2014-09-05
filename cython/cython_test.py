import astrochem
import numpy

density = numpy.logspace(4, 6, 128) # typically this would come from an analytical formulae
temperature = numpy.linspace(10, 100, 128)
times = numpy.zeros(128)
times[1:] = numpy.logspace(-6, 8, 127) # first time is zero
av = 0

initial_abundances = {"CO": 1e-4, "HCO(+)": 1e-9, "e(-)": 1e-9}

p = astrochem.Phys()
p.cosmic = 1e-17
p.chi = 0

c = astrochem.Cell( av , density[0] , temperature[0]  )

verbose = 1

s = astrochem.Solver( c,  "network.chm", p , astrochem._ABS_ERR_DEFAULT, astrochem._REL_ERR_DEFAULT, initial_abundances, density[0], verbose )

for i in range(1, len(times)):
    s.density = density[i]
    s.temperature = temperature[i]
    try:
        abundances = s.solve(times[i])
    except ArithmeticError as e:
        raise "something went wrong: %s" % e
    CO_abundance = abundances["CO"] # abundances is a dict of all species
    print CO_abundance
