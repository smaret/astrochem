#!@PYTHON@
#
#  chmconvert - Convert chemical networks into .chm format
#
#  Copyright (c) 2006-2011 Sebastien Maret
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

import sys
import getopt

VERSION = "0.4"

class reaction:
    """A chemical reaction."""
    pass

def usage():
    """Display usage."""

    print """Usage: chmconvert [option] [file]

Common options:
   -h, --help         Display this help
   -V, --version      Display chmconvert version information
   -o, --output       Write the edited network in a file

See chmconvert(1) man page for a complete list of commands and options.
Report bugs to <sebastien.maret@obs.ujf-grenoble.fr>."""

def version():
    """Display version number."""

    print "This is chmconvert, version %s" % VERSION
    print """Copyright (c) 2006-2011 Sebastien Maret

This is free software. You may redistribute copies of it under the terms
of the GNU General Public License. There is NO WARRANTY, to the extent
permitted by law."""

def format_specie(specie):
    """Format a specie name to astrochem format."""

    specie = specie.strip()

    # In astrochem format, the charge of ions is in parenthesis.
    specie = specie.replace("+", "(+)")
    specie = specie.replace("-", "(-)")
    
    # For elements that have two letters in their names, the
    # second one is in lower case.
    specie = specie.replace("HE", "He")
    specie = specie.replace("NA", "Na")
    specie = specie.replace("MG", "Mg")
    specie = specie.replace("CL", "Cl")
    specie = specie.replace("SI", "Si")
    specie = specie.replace("FE", "Fe")
    
    # Grains are in lowercase. Neutral grains are simply noted
    # 'grain'
    specie = specie.replace("GRAIN0", "grain")
    specie = specie.replace("GRAIN", "grain")
    
    # Electrons are noted as e(-)
    if specie == "E": specie = "e(-)"

    return specie

def format_react(react):
    """Format a reaction to astrochem format."""
    
    # In .chm format, we have one, two or three reactants,
    # and one, two, three or four products
    
    f_react= "%-12s" % react.reactant1
    if react.reactant2:
        f_react= f_react+ " + "
    else:
        f_react= f_react+ "   "
    f_react= f_react+ "%-12s" % react.reactant2
    if react.reactant3:
        f_react= f_react+ " + "
    else:
        f_react= f_react+ "   "
    f_react= f_react+ "%-12s" % react.reactant3
    f_react= f_react+ " -> %-12s" % react.product1
    if react.product2:
        f_react= f_react+ " + "
    else:
        f_react= f_react+ "   "
    f_react= f_react+ "%-12s" % react.product2
    if react.product3:
        f_react= f_react+ " + "
    else:
        f_react= f_react+ "   "
    f_react= f_react+ "%-12s" % react.product3
    if react.product4:
        f_react= f_react+ " + "
    else:
        f_react= f_react+ "   "
    f_react= f_react+ "%-12s" % react.product4
    f_react= f_react+ "   %9.2e %9.2e %9.2e %2i %4i\n" % \
        (react.alpha, react.beta, react.gamma, react.type, 
         react.number)

    return f_react

def convert(filein, format, fileout, ignore_unknown = True):
    """Convert a chemical network file to astrochem native format."""

    # Read the input file and convert it
    
    if format == "osu":
	
	# OSU file are fix format text files. The first part of the
	# file lists the species, and the second part the reaction. We
	# disantangle them by checking the line length. Networks older
	# than 2007 have 113 characters for each reaction lines.
	# Newer networks have 119 characters, with an additional
	# column that list the error on the rate.
        
	for line in filein:
	    if len(line) == 113 or len(line) == 119:
                react = reaction()
                try:
                    react.reactant1 = format_specie(line[0:8])
                    react.reactant2 = format_specie(line[8:16])
                    react.reactant3 = format_specie(line[16:24])
                    react.product1 = format_specie(line[24:32])
                    react.product2 = format_specie(line[32:40])
                    react.product3 = format_specie(line[40:48])
                    react.product4 = format_specie(line[48:56])
                    react.alpha = float(line[64:73])
                    react.beta = float(line[73:82])
                    react.gamma = float(line[82:91])
                    react.type = int(line[91:93])
                    react.number = int(line[107:111])
                except:
                    sys.stderr.write("chmconvert: error while reading network file, line %s.\n%"
                                      % line_number)
                    exit(1)

                # In OSU format, cosmic rays and photons are implicit
                # reactants/products for cosmic-ray ionization,
                # photo-ionization, photo-dissociation, radiative
                # association and radiative recombination.
		
		if react.type == 1:
		    # Be sure to not overwrite anything
		    if not (react.reactant2):
			react.reactant2 = "cosmic-ray"
                    else:
                        sys.stderr.write("chmconvert: error while reading network file, line %s.\n%"
                                          % line_number)
                        exit(1)
		if react.type == 13:
		    if not (react.reactant2):
			react.reactant2 = "uv-photon"
                    else:
                        sys.stderr.write("chmconvert: error while reading network file, line %s.\n%"
                                          % line_number)
                        exit(1)
		if react.type == 8 or react.type == 10:
		    if not (react.product2):
			react.product2 = "photon"
		    elif not (react.product3):
			react.product3 = "photon"
                    elif not (react.product4):
                        react.product4 = "photon"
                    else:
                        sys.stderr.write("chmconvert: error while reading network file, line %s.\n%"
                                          % line_number)
                        exit(1)

                # H2 formation, electron attachement and ion
                # recombination on grains have the same type in
                # OSU. However, the rate are not computed in the same
                # fashion. In Astrochem format, only H2 formation has
                # type 0, other reactions have type -1.

                if react.type == 0:
                    if react.reactant1 == "H" and react.reactant2 == "H" and react.product1 == "H2":
                        pass
                    else:
                        react.type = -1

		fileout.write(format_react(react))

    elif format == "udfa":

        sys.stderr.write("chmconvert: conversion of UDFA networks is not yet implemented.\n")
        sys.exit(1)

def main():

    # Parse options and check arguments. Display help message if
    # an unknown option is given.
    try:
	opts, args = getopt.getopt(sys.argv[1:], "hVo:",
				   ["help", "version", "output"])
    except getopt.GetoptError:
	usage()
	sys.exit(1)

    output = None
    for opt, arg in opts:
	if opt in ("-h", "--help"):
	    usage()
	    sys.exit()
	if opt in ("-V", "--version"):
	    version()
	    sys.exit()
	if opt in ("-o", "--output"):
	    output = arg

    # Check that we have a least a filename
    if len(args) == 1:
        network_file = args[0]
        
        # Guess the format from the file extension
        if len(network_file.rsplit('.', 1)) == 2:
            network_file_base, network_file_ext = network_file.rsplit('.', 1)
        else:
            sys.stderr.write("chmconvert: file has no extension.\n")
            sys.exit(1)  
        if network_file_ext in ("osu", "udfa"):
            format = network_file_ext
        else:
            sys.stderr.write("chmconvert: unknown network format \"%s\".\n" 
                              % network_file_ext)
            sys.exit(1)
		
        # Open input and output files. We set fileout to stdout if
        # the --output option was not set.
        try:
            filein = open(network_file)
        except:
            sys.stderr.write("chmconvert: can't open %s.\n" % network_file)
            sys.exit(1)
        if not output:
            fileout = sys.stdout
        else:		
            try:
                fileout = open(output, 'w')
            except:
                sys.stderr.write("chmconvert: can't open %s.\n" % output)
                sys.exit(1)
	    
        # Convert the file
        convert(filein, format, fileout)
    
    else:
        sys.stderr.write("chmconvert: no file to convert.\n")
        sys.stderr.write("Type \'chmconvert --help\' for more information.\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
