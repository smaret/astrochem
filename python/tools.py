#  tools.py - Python tools for astrochem
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

"""
Various tools for Astrochem.

"""

import string
import numpy

class reaction:
    """
    Chemical reaction class.

    Attributes
    ----------
    reactants : list of str
        List of reactants.
    products : list of str
        List of products.
    alpha : float 
        Reaction constant.
    beta : float 
        Reaction constant.
    gamma : float 
        Reaction constant.
    rtype : int
        Reaction type.
    rnumber: int
        Reaction number.

    """

    def __init__(self, reactants, products, alpha, beta, gamma, rtype, rnumber):
        """
        Create a chemical reaction instance.

        Parameters
        ----------
        reactants : list of str
            List of reactants.
        products : list of str
            List of products.
        alpha : float 
            Reaction constant.
        beta : float 
            Reaction constant.
        gamma : float 
            Reaction constant.
        rtype : int
            Reaction type.
        rnumber: int
            Reaction number.

        """
        self.reactants = reactants
        self.products = products
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.type = rtype
        self.number = rnumber

    def __repr__(self):
        """
        Returns the string representation of a reaction.

        Returns
        -------
        str
            The string representation of the reaction.

        """

        return  "reaction(" \
            + self.reactants.__repr__() + ", " \
            + self.products.__repr__() + ", " \
            + self.alpha.__repr__() + ", " \
            + self.beta.__repr__() + ", " \
            + self.gamma.__repr__() + ", " \
            + self.type.__repr__() + ", " \
            + self.number.__repr__() + ")"

    def __eq__(self, other):
        """
        Compares the reaction with another.

        This method compares two reactions. The reactions are supposed
        to be equal if both the reactants and products are equal,
        regardless of the reaction rates.

        Parameters
        ----------
        other : `reaction`
            Reaction instance to compare with.

        Returns
        -------
        bool
            True if the reactions are equal, False otherwise.

        """

        if not(isinstance(other, reaction)):
            raise ValueError, "Argument should be a network instance"

        if len(self.reactants) != len(other.reactants) or len(self.products) != len(other.products):
            return False

        for r in self.reactants:
            if not r in other.reactants:
                return False

        for p in self.products:
            if not p in other.products:
                return False

        return True

    def totex(self):
        """
        Returns a reaction in TeX format.

        Returns
        -------        
        str
            TeX formated reaction string.

        """

        tex_r = "$"
        for s in self.reactants:
            tex_r += _totex_species(s) + " + "
        tex_r = tex_r[0:-3] + " \\rightarrow  "
        for s in self.products:
            tex_r += _totex_species(s) + " + "
        tex_r = tex_r[0:-3] + "$"

        return tex_r

class network_reader:
    """Chemical network reader class.

    """

    def __init__(self, reactions ):
        """Creates a network reader object.

        Parameters
        ----------
        reactions : list of `reaction`
            List of reactions

        """

        for r in reactions:
            if isinstance(r, reaction):
                continue
            raise ValueError, "Argument should be a list of network instances"

        self.data = reactions

    def __repr__(self):
        """Returns the string representation of a network.

        Returns
        -------
        str
            The string representation of the network.

        """

        string = "network([ "
        for r in self.data:
            string += "          " + r.__repr__() + ",\n"
        string = string[:-3] + "])"

        return string

    @staticmethod
    def fromfile(f, fileformat):
        """Read a network from a file.

        This function reads a chemistry network from a file and
        creates a network instance. Supported formats are chm, osu and
        kida.

        Parameters
        ----------
        f : file
            Network file
        fileformat : str
            Network file format ("chm", "osu" or "kida")

        Returns
        -------
        `network_reader`
            Network

        """

        def _read_chm(f):
            """Read a network in the chm format.

            Parameters
            ----------
            f : file
                Network file

            Returns
            -------
            `network_reader`
                Network

            """

            net = []
            linenumber = 0
            for line in f:

                linenumber += 1

                if line[0] == "#":
                    continue

                try:
                    react = line.rsplit(None, 5)[0]
                    rconsts = line.rsplit(None, 5)[1:]

                    reactants = map(lambda x: x.strip(), react.split("->")[0].split(" + "))
                    products = map(lambda x: x.strip(), react.split("->")[1].split(" + "))

                    alpha = float(rconsts[0])
                    beta = float(rconsts[1])
                    gamma = float(rconsts[2])
                    rtype = int(rconsts[3])
                    rnumber = int(rconsts[4])
                except:
                    raise Exception, "incorrect input on line %i" % linenumber

                react = reaction(reactants, products, alpha, beta, gamma,
                                 rtype, rnumber)
                net.append(react)

            return net

        def _read_osu(f):
            """Read a network in the osu format.

            Parameters
            ----------
            f : file
                Network file

            Returns
            -------
            `network_reader`
                Network

            """

            def _format_species_osu(species):
                """Format a species in osu format to chm format.

                """

                # In chm format, electrons are noted e(-) and grains are
                # simply noted 'grain'
                convert = {"E": "e(-)", "GRAIN0": "grain", "GRAIN+": "grain(+)",
                           "GRAIN-": "grain(-)"}
                if convert.has_key(species):
                    return convert[species]

                # The charge of ions is in parenthesis.  Be carefull with
                # "-" because some species names may start with "c-",
                # "l-", etc.
                species = species.replace("+", "(+)")
                if species[-1] == "-":
                    species = species[:-1] + "(-)"

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
            # column that gives the error on the rate.

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
                    elif react.type == 4 or react.type == 8 or react.type == 10:
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

        def _read_kida(f):
            """Read a network in the kida format.

            Parameters
            ----------
            f : file
                Network file

            Returns
            -------
            `network_reader`
                Network

            """

            def _format_species_kida(species):
                """Format a species in kida format to chm format.

                """

                # In chm format, electrons are noted e(-), cosmic-rays (or
                # cosmic-ray secondary photons) are noted cosmic-ray and
                # UV photons are noted "uv-photons"
                convert = {"e-": "e(-)", "CR": "cosmic-ray", "CRP": "cosmic-ray",
                           "Photon": "photon"}
                if convert.has_key(species):
                    return convert[species]

                # The charge of ions is in parenthesis.  Be carefull with
                # "-" because some species names may start with "c-",
                # "l-", etc.
                species = species.replace("+", "(+)")
                if species[-1] == "-":
                    species = species[:-1] + "(-)"

                return species

            # KIDA uses a different convention for reaction types than OSU
            # and Astrochem.

            kida2chm_type = {1: 1,    # direct cosmic-ray processes
                             2: 1,    # photo-processes induced by cosmic-rays
                             3: 13,   # photo-reactions
                             4: 2,    # bimolecular (ion-neutral or neutral-neutral)
                             5: 2,    # charge exchange reactions
                             6: 8,    # radiative association
                             7: 5,    # associative detachment
                             8: 9}    # electronic recombination

            # KIDA files are fixed format text files, although it is also
            # possible to download CSV files from the database. Here we
            # assume that the file is in fixed format.

            net = []
            linenumber = 0

            for line in f:
                linenumber += 1
                try:
                    reactant1 = line[0:10].strip()
                    reactant2 = line[11:21].strip()
                    reactant3 = line[22:32].strip()
                    product1 = line[34:44].strip()
                    product2 = line[45:55].strip()
                    product3 = line[56:66].strip()
                    product4 = line[67:77].strip()
                    product5 = line[78:88].strip()
                    alpha = float(line[90:100])
                    beta = float(line[101:111])
                    gamma = float(line[112:122])
                    uncertainty_factor = float(line[123:131])
                    uncertainty_factor_temp_dep = float(line[132:140])
                    uncertainty_type = line[141:145]
                    rtype = int(line[146:148])
                    trange = float(line[149:155]), float(line[156:162])
                    rnumber = int(line[163:168])
                    rate_number = int(line[169:170])
                    recommendation = int(line[171:173])
                except:
                    raise Exception, "incorrect input on line %i" % linenumber

                reactants = []
                for species in [reactant1, reactant2, reactant3]:
                    if species != "":
                        reactants.append(_format_species_kida(species))

                products = []
                for species in [product1, product2, product3, product4]:
                    if species != "":
                        products.append(_format_species_kida(species))

                if rtype in kida2chm_type.keys():
                    rtype = kida2chm_type[rtype]
                else:
                    raise Exception, "unknown reaction type on line %i" % linenumber

                react = reaction(reactants, products, alpha, beta, gamma,
                                 rtype, rnumber)

                net.append(react)

            return net

        if fileformat == "chm":
            l = _read_chm(f)
        elif fileformat == "osu":
            l = _read_osu(f)
        elif fileformat == "kida":
            l = _read_kida(f)
        else:
            raise ValueError, "Unknown format"

        return network_reader(l)

    def tofile(self, f, renumber = False):
        """Write network in a file.

        Parameters
        ----------

        f : file
            Network file handle
        renumber : bool, optional
            Renumber reactions (default False)

        """

        # Renumber the reactions, if requested

        if renumber:
            react_number = 1
            for react in self.data:
                react.number = react_number
                react_number += 1

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

    def duplicate_react_numbers(self):
        """Find reactions with the same reaction number.

        Returns
        -------
        list of str
            List of duplicated reaction_numbers.

        """

        react_numbers = []
        count = []
        dups = []
        for react in self.data:
            if not react.number in react_numbers:
                react_numbers.append(react.number)
                count.append(1)
            else:
                index = react_numbers.index(react.number)
                count[index] += 1
        for i in range(len(react_numbers)):
            if count[i] > 1:
                dups.append(react_numbers[i])

        return dups

    def duplicate_reactions(self):
        """Find duplicate reactions.

        Returns
        -------
        list of int
            List of reaction_numbers of duplicated reactions.

        """

        import copy

        dups = []
        reacts = copy.copy(self.data)
        while reacts:
            react = reacts.pop()
            dup = [react.number]
            while True:
                try:
                    index = reacts.index(react)
                except ValueError:
                    break
                dup.append(reacts.pop(index).number)
            if len(dup) > 1:
                dups.append(dup)

        return dups

    def getreact(self, number):
        """Returns the reaction with a given number.

        Parameters
        ----------
        number : int
            The reaction number.
        
        Returns
        -------
        `reaction`
            The reaction found.

        Raises
        ------
        ValueError
            If no reaction with this number if found.

        """
        for react in self.data:
            if react.number == number:
                return react
        raise ValueError("no reaction with number %i." % number)

def readfilesattrs(filename):
    """Read chem_file and source_file attributes from an hdf5 output file
    
    Parameters
    ----------
    filename : str
        Path to output file.

    Returns
    -------
    chemfile : str
        chem_file attribute
    sourcefile : str
        source_file attribute

    """
    import h5py
    f = h5py.File( filename, "r")
    if "chem_file" in f.attrs:
        chemfile = f.attrs.get("chem_file")
    else:
        raise IOError("%s file doest not contain \"chem_file\" attribute" % filename )
    if "source_file" in f.attrs:
        sourcefile = f.attrs.get("source_file")
    else:
        raise IOError("%s file doest not contain \"source_file\" attribute" % filename )

    return chemfile, sourcefile

def listspecies( filename ):
    """Print available species from an hdf5 output file and return it in a array

    Parameters
    ----------
    filename : str
        Path to output file.

    Returns
    -------
    list if str
        Species list

    """
    import h5py
    f = h5py.File( filename, "r")
    if "Species" in f:
        s_d = f.get("Species")
        s = s_d[...]
    else:
        raise IOError("%s file doest not contain \"Species\" dataset" % filename )

    return s;

def readabun(filename, specie):
    """Read abundances for a specific specie from an hdf5 output file and
    return arrays of time and abundances

    Parameters
    ----------
    filename : str 
        Path to the output file
    specie : str
        Name of specie to read abundance of

    Returns
    -------
    timesteps : list of float
        List of timesteps
    abundance : list of float
        List of abundances

    """
    import h5py
    f = h5py.File( filename, "r")
    if "TimeSteps" in f:
        t_d = f.get("TimeSteps")
        t = t_d[...]
    else:
        raise IOError("%s file doest not contain \"TimeSteps\" dataset" % filename )
    if "Species" in f:
        s_d = f.get("Species")
        s = s_d[...]
    else:
        raise IOError("%s file doest not contain \"Species\" dataset" % filename )
    if specie in s:
        specie_index = s.tolist().index(specie)
    else:
        raise ValueError("%s file doest not contain %s specie" % ( filename, specie ) )
    if "Abundances" in f:
        a_d = f.get("Abundances")
        a = a_d[:,:,specie_index]
    else:
        raise IOError("%s file doest not contain \"Abundances\" dataset" % filename )
    return t,a.transpose()

def readrout(filename, specie):
    """Read a rout from a hdf5 output file and return arrays of time, shell number, formation/destruction rates.

    Parameters
    ----------
    filename : str
        Path to the output file.
    specie : str
        Name of specie to read route of.
    
    Returns
    -------
    timesteps : list of float
        List of timesteps.
    shells : list of float
        List of shell numbers.
    formation_reac : list of `reaction`
        List of formation reactions.
    formation_rate : list of float
        List of formation rates.
    destruction_reac : list of `reaction`
        List of destruction reactions.
    destruction_rate : list of float
        List of destruction rates.

    """
    import h5py
    f = h5py.File( filename, "r")
    if "TimeSteps" in f:
        t_d = f.get("TimeSteps")
        t = t_d[...]
    else:
        raise IOError("%s file doest not contain \"TimeSteps\" dataset" % filename )
    if "Species" in f:
        s_d = f.get("Species")
        s = s_d[...]
    else:
        raise IOError("%s file doest not contain \"Species\" dataset" % filename )
    if specie in s:
        specie_index = s.tolist().index(specie)
    else:
        raise ValueError("%s file doest not contain specie %s" % (filename, specie) )
    if "Routes" in f:
        g = f.get("Routes")
    else:
        raise IOError("%s file doest not contain \"Routes\" group" % filename )
    routeDatasetName = "route_" + specie
    if routeDatasetName in g:
        r_d = g.get( routeDatasetName )
    else:
        raise IOError("%s file doest not contain %s dataset", filename, routeDatasetName )
    formation_reac = r_d[:,:,specie_index,:]['formation_rate']['reaction_number']
    formation_rate = r_d[:,:,specie_index,:]['formation_rate']['reaction_rate']
    destruction_reac = r_d[:,:,specie_index,:]['destruction_rate']['reaction_number']
    destruction_rate = r_d[:,:,specie_index,:]['destruction_rate']['reaction_rate']

    shell = numpy.empty( r_d.len() )
    for i in range(  r_d.len() ):
        shell[i] = i
    return t, shell, formation_reac, formation_rate, destruction_reac, destruction_rate

def converttolegacy( filename, specie ):
    """Convert a hdf5 output file specific species to .abun and .rout legacy format.

    Parameters
    ----------
    filename : str
        Path to the output file
    species : str
        Name of specie to read abundance and route of

    """
    import h5py
    f = h5py.File( filename, "r")
    if "Species" in f:
        s_d = f.get("Species")
        s = s_d[...]
    else:
        raise IOError("%s file doest not contain \"Species\" dataset" % filename )
    if specie == "ALL":
        specie_list = s.tolist()
    else:
        specie_list = [specie]
    for specie_val in specie_list:
        if specie_val in s:
            specie_index = s.tolist().index(specie_val)
        else:
            raise ValueError("%s file doest not contain specie %s" % (filename, specie_val) )

        if "TimeSteps" in f:
            t_d = f.get("TimeSteps")
            t = t_d[...]
        else:
            raise IOError("%s file doest not contain \"TimeSteps\" dataset" % filename )

        if "Abundances" in f:
            a_d = f.get("Abundances")
            a = a_d[:,:,specie_index]
        else:
            raise IOError("%s file doest not contain \"Abundances\" dataset" % filename )

        abunfile = open(specie_val+".abun","w");
        abunfile.write("# "+specie_val+" abundance computed by astrochem\n")
        abunfile.write("# time [yr] / shell number\n")
        abunfile.write("#\n")

        abunfile.write("             ")
        for shell in range( 0, a.shape[0] ):
            abunfile.write("{:>13d} ".format(shell) )
        abunfile.write("\n")

        for t_idx in range( 0, t.shape[0] ):
            abunfile.write(  "{:.2e}  ".format( t[t_idx] ) )
            for a_val in  a[:,t_idx]:
                abunfile.write(  "{:.2e}  ".format( a_val ) )
            abunfile.write("\n")
        abunfile.close()

        if "Routes" in f:
            g = f.get("Routes")
            routeDatasetName = "route_" + specie_val
            if routeDatasetName in g:
                r_d = g.get( routeDatasetName )
            else:
                raise IOError("%s file doest not contain %s dataset", filename, routeDatasetName )
            formation = r_d[:,:,specie_index,:]['formation_rate']
            destruction = r_d[:,:,specie_index,:]['destruction_rate']

            routfile = open(specie_val+".rout","w");
            routfile.write("# Main "+specie_val+" formation/destruction routes computed by astrochem\n")
            routfile.write("# shell number  time [yr]  reaction number 1  reaction rate 1 [cm-3/s]... \n")
            routfile.write("#\n")

            for shell in range( 0, a.shape[0] ):
                for t_idx in range( 0, t.shape[0] ):
                    routfile.write("   {:d} ".format(shell) )
                    routfile.write(  "  {:.2e}".format( t[t_idx] ) )
                    for r_val in  formation[shell,t_idx,:]:
                        routfile.write(  "  {:>4d}".format( r_val[0] ) )
                        routfile.write(  "  {:.2e}".format( r_val[1] ) )
                    routfile.write("\n")
                    routfile.write("   {:d} ".format(shell) )
                    routfile.write(  "  {:.2e}".format( t[t_idx] ) )
                    for r_val in  destruction[shell,t_idx,:]:
                        routfile.write(  "  {:>4d}".format( r_val[0] ) )
                        routfile.write(  " {:.2e}".format( -r_val[1] ) )
                    routfile.write("\n")
            routfile.close()

def readabunlegacy(filename):
    """Read an abund file and return arrays of time and abundances.

    Parameters
    ----------
    filename : str
        Path to the abund file.

    Returns
    -------
    timesteps : list of floats
        List of timesteps
    abundances : list of floats
        List of abundances

    """

    f = open(filename)

    # Skip comments and cell number
    f.readline()
    f.readline()
    f.readline()
    f.readline()

    # Read all lines
    lines = map(string.strip, f.readlines())

    # Construct a list of the elements of the array
    a = []
    for line in lines:
        newline = []
        for elem in string.split(line):
	    elem = string.atof(elem)
            newline.append(elem)
        a.append(newline)

    # Create an array from the list
    a = numpy.array(a)

    time = a[:,0]
    abund = a[:,1:]

    return time, abund

def readroutlegacy(filename):
    """Read a rout file and return arrays of time, cell number,
    formation/destruction rates

    Parameters
    ----------
    filename : str
        Path to the output file.
    
    Returns
    -------
    timesteps : list of float
        List of timesteps.
    shells : list of float
        List of shell numbers.
    formation_reac : list of `reaction`
        List of formation reactions.
    formation_rate : list of float
        List of formation rates.
    destruction_reac : list of `reaction`
        List of destruction reactions.
    destruction_rate : list of float
        List of destruction rates.

    """

    f = open(filename)

    # Skip comments
    f.readline()
    f.readline()
    f.readline()

    # Read all lines
    lines = map(string.strip, f.readlines())

    # Construct a list of the elements of the array
    a = []
    for line in lines:
        newline = []
        for elem in string.split(line):
	    elem = string.atof(elem)
            newline.append(elem)
        a.append(newline)

    # Create an array from the list
    a = numpy.array(a)

    # Extract cell number, time, formation/destruction reactions and
    # rates.
    cell = numpy.array(a[::2,0], dtype = int)
    time = a[::2,1]
    formation_reac = numpy.array(a[0::2,2::2], dtype = int)
    formation_rate = a[0::2,3::2]
    destruction_reac = numpy.array(a[1::2,2::2], dtype = int)
    destruction_rate = a[1::2,3::2]

    cell = numpy.unique(cell)
    ncell = len(numpy.unique(cell))
    time = numpy.unique(time)
    ntime = len(numpy.unique(time))

    # Reshape formation/destruction reactions and rates arrays 
    formation_reac = numpy.reshape(formation_reac, (ncell, ntime, -1))
    formation_rate = numpy.reshape(formation_rate, (ncell, ntime, -1))
    destruction_reac = numpy.reshape(destruction_reac, (ncell, ntime, -1))
    destruction_rate = numpy.reshape(destruction_rate, (ncell, ntime, -1))

    return time, cell, formation_reac, formation_rate, destruction_reac, destruction_rate

def _totex_species(spec):
    """Returns a species in TeX format

    Parameters
    ----------
    species : str
        Species name

    Returns
    -------
    str
        TeX formated species string

    """

    if spec in ["cosmic-ray", "uv-photon"]:
        return spec

    tex_s = ""
    for char in spec:
        if char == '(' or char == ')':
            continue
        elif (char == '+' or char == '-'):
            tex_s = tex_s + "^{" + char + "}"
        elif char.isdigit():
            tex_s = tex_s + "_{" + char + "}"
        else:
            tex_s = tex_s + char

    return tex_s
