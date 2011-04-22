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

VERSION = "0.3"

class reaction:
    """
    Chemical reaction class.
    
    """

    def __init__(self, reactants, products, alpha, beta, gamma, rtype, rnumber):
        """
        Creates a chemical reaction

        Arguments:
        reactants          -- list of reactants
        products           -- list of products
        alpha, beta, gamma -- reaction constants
        rtype              -- reaction type
        rnumber            -- reaction number

        """

        self.reactants = reactants
        self.products = products
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.type = rtype
        self.number = rnumber
        
class network:
    """
    Chemical network class

    """

    def __init__(self, f, format):
        """
        Read a network from a file.

        This function reads a chemistry network from a file and
        creates a network object. Supported format are chm, osu and
        kida.

        Arguments:
        f      -- Network file handle
        format -- Network file format ("chm", "osu" or "kida")

        """

        if format == "chm":
            self.data = self._read_chm(f)
        elif format == "osu":
            self.data = self._read_osu(f)
        elif format == "kida":
            self.data = self._read_kida(f)
        else:
            raise ValueError, "Unknown format"
        
    @staticmethod
    def _read_chm(f):
        """
        Read a network in the chm format.

        Arguments:
        f -- Network file handle

        """

        pass

    @staticmethod
    def _read_osu(f):
        """
        Read a network in the osu format.

        Arguments:
        f -- Network file handle

        """

        def _format_species_osu(species):
            """
            Format a species in osu format to chm format.
            
            """

            # In chm format, electrons are noted e(-) and grains are
            # simply noted 'grain'
            convert = {"E": "e(-)", "GRAIN0": "grain", "GRAIN+": "grain(+)",
                       "GRAIN-": "grain(-)"}
            if convert.has_key(species):
                return convert[species]

            # The charge of ions is in parenthesis.
            species = species.replace("+", "(+)")
            species = species.replace("-", "(-)")
    
            # For elements that have two letters in their names, the
            # second one is in lower case.
            species = species.replace("HE", "He")
            species = species.replace("NA", "Na")
            species = species.replace("MG", "Mg")
            species = species.replace("CL", "Cl")
            species = species.replace("SI", "Si")
            species = species.replace("FE", "Fe")
        
            return species

	# OSU file are fix format text files. The first part of the
	# file lists the species, and the second part the reaction. We
	# disantangle them by checking the line length. Networks older
	# than 2007 have 113 characters for each reaction lines.
	# Newer networks have 119 characters, with an additional
	# column that lists the error on the rate.

        net = []
        linenumber = 0

	for line in f:
            linenumber += 1
 	    if len(line) == 113 or len(line) == 119:
                try:
                    reactant1 = line[0:8].strip()
                    reactant2 = line[8:16].strip()
                    reactant3 = line[16:24].strip()
                    product1 = line[24:32].strip()
                    product2 = line[32:40].strip()
                    product3 = line[40:48].strip()
                    product4 = line[48:56].strip()
                    alpha = float(line[64:73])
                    beta = float(line[73:82])
                    gamma = float(line[82:91])
                    rtype = int(line[91:93])
                    rnumber = int(line[107:111])
                except:
                    raise Exception, "incorrect input on line %i" % linenumber
                    
                reactants = []
                for species in [reactant1, reactant2, reactant3]:
                    if species != "":
                        reactants.append(_format_species_osu(species))

                products = []
                for species in [product1, product2, product3, product4]:
                    if species != "":
                        products.append(_format_species_osu(species))

                react = reaction(reactants, products, alpha, beta, gamma,
                                 rtype, rnumber)
                
                # In OSU format, cosmic rays and photons are implicit
                # reactants/products for cosmic-ray ionization,
                # photo-ionization, photo-dissociation, radiative
                # association and radiative recombination.
		
		if react.type == 1:
                    react.reactants.append("cosmic-ray")
                elif react.type == 8 or react.type == 10:
                    react.products.append("photon")
		elif react.type == 13:
                    react.reactants.append("uv-photon")

                # H2 formation, electron attachement and ion
                # recombination on grains have the same type in
                # OSU. However, the rate are not computed in the same
                # fashion. In Astrochem format, only H2 formation has
                # type 0, other reactions have type -1.

                if react.type == 0:
                    if react.reactants == ["H", "H"] and react.products == ["H2"]:
                        pass
                    else:
                        react.type = -1

                net.append(react)

        return net

    @staticmethod
    def _read_kida(self, filename):
        """
        Read a network in the kida format.

        Arguments:
        filename -- Network file name        

        """

        pass

    def write(self, f):
        """
        Write network in a file.
        
        Arguments:
        f -- Network file handle

        """

        # In order to format the file, we first to find out the
        # maximum number of reactants and products in the network, the
        # lenght of the largest species name, the lenght of the
        # largest reaction types, and the lenght of the largest
        # reaction number.
        
        max_reactants = 0
        max_products = 0
        max_species_lenght = 0
        max_type_lenght = 0
        max_number_lenght = 0

        for react in self.data:
            if len(react.reactants) > max_reactants:
                max_reactants = len(react.reactants)
            if len(react.products) > max_products:
                max_products = len(react.products)
            for species in react.reactants + react.products:
                if len(species) > max_species_lenght:
                    max_species_lenght = len(species)
            if len("%s" % react.type) >  max_type_lenght:
                max_type_lenght = len("%s" % react.type)
            if len("%s" % react.number) >  max_number_lenght:
                max_number_lenght = len("%s" % react.number)

        species_format = "%-" + "%i" % max_species_lenght + "s"
        type_format = "%" + "%i" % max_type_lenght + "i"
        number_format = "%" + "%i" % max_number_lenght + "i"

        # Now that we have defined the format, write each reaction in
        # that format.
        
        for react in self.data:

            for i in range(max_reactants + 1):
                if i == 0:
                    string = species_format % react.reactants[i]
                elif i <= len(react.reactants) - 1:
                    string += " + " + species_format % react.reactants[i]
                else:
                    string += "   " + species_format % ""
                    
            for i in range(max_products + 1):
                if i == 0:
                    string += " -> " + species_format % react.products[i]
                elif i <= len(react.products) - 1:
                    string += " + " + species_format % react.products[i]
                else:
                    string += "   " + species_format % ""
                    
            string += ("   %9.2e %9.2e %9.2e " + type_format + " " + number_format + "\n") \
                % (react.alpha, react.beta, react.gamma, react.type, 
                   react.number)
            f.write(string)

    def check(self):
        """
        Check network.

        """        
        
        pass

def usage():
    """
    Display usage.

    """

    print """Usage: chmconvert [option] [file]

Common options:
   -h, --help         Display this help
   -V, --version      Display chmconvert version information
   -o, --output       Write the edited network in a file

See chmconvert(1) man page for a complete list of commands and options.
Report bugs to <sebastien.maret@obs.ujf-grenoble.fr>."""

def version():
    """
    Display version number.

    """

    print "This is chmconvert, version %s" % VERSION
    print """Copyright (c) 2006-2011 Sebastien Maret

This is free software. You may redistribute copies of it under the terms
of the GNU General Public License. There is NO WARRANTY, to the extent
permitted by law."""

def main():
    """
    Main function for chmconvert

    """
    
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
    else:
        sys.stderr.write("chmconvert: no file to convert.\n")
        sys.stderr.write("Type \'chmconvert --help\' for more information.\n")
        sys.exit(1)
        
    # Guess the format from the file extension
    if len(network_file.rsplit('.', 1)) == 2:
        network_file_base, network_file_ext = network_file.rsplit('.', 1)
    else:
        sys.stderr.write("chmconvert: file has no extension.\n")
        sys.exit(1)  
    if network_file_ext in ("osu"):
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
	    
    # Read the network file and convert it
    try:
        net = network(filein, format = format)
    except Exception, err:
        sys.stderr.write("chmconvert: %s.\n" % err)
        sys.exit(1)        
    filein.close()

    # Write the converted file:
    net.write(fileout)
    fileout.close()
    
if __name__ == "__main__":
    main()
