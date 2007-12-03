#!@PYTHON@
# chmconvert - tool to convert chemical networks into .chm format

import sys
import getopt

VERSION = "0.1"



def usage():
    """Display usage."""

    print """Usage: chmconvert [option] [file]

Common options:
   -h, --help         Display this help
   -V, --version      Display chemutil version information
   -o, --output       Write the edited network in a file

See chmconvert(1) man page for a complete list of commands and options.
Report bugs to <smaret@umich.edu>."""

def version():
    """Display version number."""

    print "This is chmconvert, version %s" % VERSION
    print """Copyright (c) 2006-2007 Sebastien Maret

This is free software. You may redistribute copies of it under the terms
of the GNU General Public License. There is NO WARRANTY, to the extent
permitted by law."""

def format_specie(specie, format):
    """Format a specie name to astrochem format."""

    specie = specie.strip()

    if format == "osu":

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

    elif format == "udfa":

        specie = specie.replace("+", "(+)")
        specie = specie.replace("-", "(-)")

        # Both cosmic-rays and cosmic-ray induced photons
        # cosmic-ray-photons are noted 'cosmic-ray'
        # FixMe: should we distinguish between the two?
        specie = specie.replace("CRPHOT", "cosmic-ray")
        specie = specie.replace("CRP", "cosmic-ray")

        # Photons are noted 'photon'
        specie = specie.replace("PHOTON", "photon") 

    return specie

def format_react(reactant1, reactant2, reactant3,
                 product1, product2, product3, product4,
                 alpha, beta, gamma, reaction_type,
                 reaction_number):
    """Format a reaction to astrochem format."""
    
    # In .chm format, we have one, two or three reactants,
    # and one, two, three or four products
    
    reaction = "%-12s" % reactant1
    if reactant2:
        reaction = reaction + " + "
    else:
        reaction = reaction + "   "
    reaction = reaction + "%-12s" % reactant2
    if reactant3:
        reaction = reaction + " + "
    else:
        reaction = reaction + "   "
    reaction = reaction + "%-12s" % reactant3
    reaction = reaction + " -> %-12s" % product1
    if product2:
        reaction = reaction + " + "
    else:
        reaction = reaction + "   "
    reaction = reaction + "%-12s" % product2
    if product3:
        reaction = reaction + " + "
    else:
        reaction = reaction + "   "
    reaction = reaction + "%-12s" % product3
    if product4:
        reaction = reaction + " + "
    else:
        reaction = reaction + "   "
    reaction = reaction + "%-12s" % product4
    reaction = reaction + "   %9.2e %9.2e %9.2e %2i %4i\n" % \
        (alpha, beta, gamma, reaction_type, reaction_number)

    return reaction

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
                try:
                    reactant1 = format_specie(line[0:8], "osu")
                    reactant2 = format_specie(line[8:16], "osu")
                    reactant3 = format_specie(line[16:24], "osu")
                    product1 = format_specie(line[24:32], "osu")
                    product2 = format_specie(line[32:40], "osu")
                    product3 = format_specie(line[40:48], "osu")
                    product4 = format_specie(line[48:56], "osu")
                    alpha = float(line[64:73])
                    beta = float(line[73:82])
                    gamma = float(line[82:91])
                    reaction_type = int(line[91:93])
                    reaction_number = int(line[107:111])
                except:
                    sys.stderr.write("chmconvert: error while reading network file, line %s.\n%"
                                      % line_number)
                    exit(1)

                # In OSU format, cosmic rays and photons are implicit
                # reactants/products for cosmic-ray ionization,
                # photo-ionization, photo-dissociation, radiative
                # association and radiative recombination.
		
		if reaction_type == 1:
		    # Be sure to not overwrite anything
		    if not (reactant2):
			reactant2 = "cosmic-ray"
                    else:
                        sys.stderr.write("chmconvert: error while reading network file, line %s.\n%"
                                          % line_number)
                        exit(1)
		if reaction_type == 13:
		    if not (reactant2):
			reactant2 = "uv-photon"
                    else:
                        sys.stderr.write("chmconvert: error while reading network file, line %s.\n%"
                                          % line_number)
                        exit(1)
		if reaction_type == 8 or reaction_type == 10:
		    if not (product2):
			product2 = "photon"
		    elif not (product3):
			product3 = "photon"
                    elif not (product4):
                        product4 = "photon"
                    else:
                        sys.stderr.write("chmconvert: error while reading network file, line %s.\n%"
                                          % line_number)
                        exit(1)

		reaction = format_react(reactant1, reactant2, reactant3,
					product1, product2, product3, product4,
					alpha, beta, gamma, reaction_type,
					reaction_number)
		fileout.write(reaction)

    elif format == "udfa":

        # UDFA files are comma separated values files. Comments lines
        # start with a pound sign.
        
        line_number = 0
        for line in filein:
            line_number = line_number + 1
            if line[0] == '#':
                continue
            line = line.split(",")
            try:
                reaction_number = int(line[0])
                reaction_code = line[1]
                reactant1 = format_specie(line[2], "udfa")
                reactant2 = format_specie(line[3], "udfa")
                reactant3 = format_specie(line[4], "udfa")
                product1 = format_specie(line[5], "udfa")
                product2 = format_specie(line[6], "udfa")
                product3 = format_specie(line[7], "udfa")
                product4 = format_specie(line[8], "udfa")
                alpha = float(line[9])
                beta = float(line[10])
                gamma = float(line[11])
                cleam = line[12]
                tmin = float(line[13])
                tmax = float(line[14])
                accuracy = line[15]
                source = line[16]
            except:
                sys.stderr.write("chmconvert: error while reading network file, line %s.\n"
                                  % line_number)
                print line
                exit(1)
            
            # In the UDFA format, the reaction type is given by a
            # string code (see Woodall et al. 2007 paper for the
            # meaning of these codes). Map these to .chm reaction
            # scheme. CD and CL have no equivalent in OU reaction
            # scheme, so we label these as "other" (14)
    
            format_code = {"NN": 7, "IN": 2, "CE": 2, "II": 11, "DR": 9, 
                           "RR": 10, "AD": 5, "RA": 4, "PH": 13, "CP": 1,
                           "CR": 1, "CD": 14, "CI": 6,  "IM": 11, "CL": 14}
            try:
                reaction_type = format_code [reaction_code]
            except KeyError:
                sys.stderr.write("chmconvert: error: unknown reaction type in network file, line %s.\n"
                                  % line_number)
                exit(1)

            # In astrochem reaction scheme, the rate of type 1
            # reactions is k = \alpha * \chi, where \chi is the H2
            # direct cosmic-ray ionization rate. In the OSU, it is
            # simply k = \alpha. Therefore \alpha should be divided by
            # the H2 cosmic-ray ionization value adopted in the UDFA
            # network, i.e. 1.2e-17 s^-1, for type 1 reactions.

            if reaction_type == 1:
                alpha = alpha / 1.2e-17

            reaction = format_react(reactant1, reactant2, reactant3,
                                    product1, product2, product3, product4,
                                    alpha, beta, gamma, reaction_type,
                                    reaction_number)
            fileout.write(reaction)

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
        if network_file_ext in ("udfa", "osu"):
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
	    
main()			  
	
	
    
