#!@PYTHON@
#
#  chmconvert - Convert chemical networks into .chm format
#
#  Copyright (c) 2006-2009 Sebastien Maret
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

VERSION = "0.1"

class reaction:
    """A chemical reaction."""
    pass

def usage():
    """Display usage."""

    print """Usage: chmconvert [option] [file]

Common options:
   -h, --help         Display this help
   -V, --version      Display chemutil version information
   -o, --output       Write the edited network in a file

See chmconvert(1) man page for a complete list of commands and options.
Report bugs to <sebastien.maret@obs.ujf-grenoble.fr>."""

def version():
    """Display version number."""

    print "This is chmconvert, version %s" % VERSION
    print """Copyright (c) 2006-2009 Sebastien Maret

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
                    react.reactant1 = format_specie(line[0:8], "osu")
                    react.reactant2 = format_specie(line[8:16], "osu")
                    react.reactant3 = format_specie(line[16:24], "osu")
                    react.product1 = format_specie(line[24:32], "osu")
                    react.product2 = format_specie(line[32:40], "osu")
                    react.product3 = format_specie(line[40:48], "osu")
                    react.product4 = format_specie(line[48:56], "osu")
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

		fileout.write(format_react(react))

    elif format == "udfa":

        # UDFA files are comma separated values files. Comments lines
        # start with a pound sign.
        
        line_number = 0
        for line in filein:
            line_number = line_number + 1
            if line[0] == '#':
                continue
            line = line.split(",")
            react = reaction()
            try:
                react.number = int(line[0])
                react.type = line[1]
                react.reactant1 = format_specie(line[2], "udfa")
                react.reactant2 = format_specie(line[3], "udfa")
                react.reactant3 = format_specie(line[4], "udfa")
                react.product1 = format_specie(line[5], "udfa")
                react.product2 = format_specie(line[6], "udfa")
                react.product3 = format_specie(line[7], "udfa")
                react.product4 = format_specie(line[8], "udfa")
                react.alpha = float(line[9])
                react.beta = float(line[10])
                react.gamma = float(line[11])
                react.cleam = line[12]
                react.tmin = float(line[13])
                react.tmax = float(line[14])
                react.accuracy = line[15]
                react.source = line[16]
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
                react.type = format_code [react.type]
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

            if react.type == 1:
                react.alpha = react.alpha / 1.2e-17

            # In the UDFA, some reactions are duplicated when their
            # temperature dependence cannot be approximated by one
            # Arrhenius-type formula. Choose the one that correspond
            # to the lowest temperature range, and remove the others.

            duplicate = [6, 93, 140, 171, 316, 660, 666, 667, 949, 1684, 1731, 
                         2939, 2946, 3533, 3534, 3536, 3537, 3552, 3353, 3555,
                         3556, 3558, 3559, 4007, 4008, 4013, 4014, 4016, 4079,
                         4109, 4112, 4115, 4138, 4558, 4590]

            # Some reaction have negative gamma. Extrapolating these
            # outside their temperature range can give unrealistic
            # high rates. In practice, these are neutral-neutral
            # reactions, so it is safe to neglect them at low
            # temperatures

            negative_gamma = [266, 282, 288, 312, 351, 353, 377, 419, 431, 443,
                              446, 448, 466, 488, 493, 501, 502, 520, 533, 540,
                              3581, 3582, 3583, 3584, 3585, 3586, 4111, 4563, 4571,
                              4573, 4574, 4575, 4579, 4582, 4583]

            # Reactions that Woodall et al. recommand to ignore to low
            # temperature. Some of them are valid only at high
            # temperature, and can't probably interpolated at lower
            # temperatures. For other, it is less clear: for example,
            # reaction 384 is valid between 10 and 300 K, yet it is
            # flagged.

            woodall = [6, 93, 141, 174, 270, 288, 319, 323, 358, 360, 384, 426,
                       438, 450, 453, 455, 495, 500, 508, 509, 527, 540, 547, 660,
                       666, 667, 949, 1684, 1731, 2938, 2945, 3532, 3533, 3535, 3536,
                       3551, 3552, 3554, 3555, 3557, 3558, 3581, 3582, 3583, 3584,
                       3585, 4006, 4007, 4012, 4013, 4015, 4078, 4108, 4110, 4111,
                       4114, 4139, 4557]

            if ((react.number not in duplicate) and
                (react.number not in negative_gamma) and
                (react.number not in woodall)):
                fileout.write(format_react(react))

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

if __name__ == "__main__":
    main()
	
	
    
