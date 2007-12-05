# test - test the results of astrochem

import math
import commands
import sys

def read_abund(specie, time):
    """Get the abundance at a given time."""

    # Read the abundance file
    file = open(specie + ".abun")
    file.readline()
    file.readline()
    file.readline()
    file.readline()
    abund = []
    tim = []
    for line in file.readlines():
        t, a = line.split()
        abund.append(float(a))
        tim.append(abs(float(t) - time))
    
    # Return the abundance that correspond to the given time
    index = tim.index(min (tim))
    
    return abund[index]


def main():
    """Compare the abundances computed by Wakelam et al. (2006) in their
    model 1 with those computed by astrochem."""

    # Run astrochem
    
    log = open("test.log", 'w')

    log.write("#-----------------------#\n")
    log.write("#   Running astrochem   #\n")
    log.write("#-----------------------#\n")
    log.write("\n")

    sys.stdout.write("Running astrochem for test model... ")
    sys.stdout.flush()

    exit_status, output = commands.getstatusoutput("../src/astrochem -v input.ini")

    log.write(output)
    sys.stdout.flush()

    if exit_status != 0:
        sys.stdout.write("error\n")
        sys.stderr.write("test: error: astrochem failed to run for test model. "
                         "See 'tests/test.log' for more details \n")
        sys.exit(1)

    sys.stderr.write("done\n")

    # Check results. Got those by eye from their Fig. 6, so it might
    # be not very accurate...

    err_tol = 3.5
    times = [1e3, 1e5, 1e7]
    species = ['CO', 'H2O', 'CS']
    abuns_fiducial = [[-5.1, -3.9, -3.8],
                      [-10.4, -6.6, -6.6],
                      [-7.4, -8.3, -8.8]]

    log.write("\n\n")
    log.write("#-----------------------#\n")
    log.write("#   Checking  results   #\n")
    log.write("#-----------------------#\n")
    log.write("\n")

    sys.stdout.write("Checking results for test model... ")

    for i in range(len(species)):
        for j in range(len(times)):
            abun_astrochem = read_abund(species[i], times[j])
            # Abundance in Wakelam et al. are relative to H2, so we divide
            # them by 2 to get abundances relative to H nuclei
            abun_fiducial = math.pow(10, abuns_fiducial[i][j]) / 2.
            if abun_fiducial > abun_astrochem:
                err = abun_fiducial / abun_astrochem
            else:
                err = abun_astrochem / abun_fiducial
            log.write("specie = %5s   t = %8.2e   astrochem = %8.2e   "
                      "fiducial = %8.2e   rel_err = %6.2f\n" %
                      (species[i], times[j], abun_astrochem, 
                       abun_fiducial, err))
            if err > err_tol:
                log.write("\n")
                log.write("error: computed abundance differ from fiducial " 
                          "by more than a factor %.1i\n" % err_tol)
                sys.stdout.write("failed\n")
                sys.exit(1)
    sys.stdout.write("ok\n")
    sys.exit(0)

if __name__ == "__main__":
    main()


