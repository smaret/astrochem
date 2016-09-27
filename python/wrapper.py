#  wrapper.py - Python wrapper for libpyastrochem
#
#  Copyright (c) 2006-2016 Sebastien Maret
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
Python wrapper for libpyastrochem.

"""

import libpyastrochem

class cell:
    """Cell class.

    Attributes
    ----------
    av : float
        Visual extinction in magnitudes.
    nh : float
        Hydrogen density in cm^-3.
    tgas : float
        Gas temperature in K.
    tdust : float
        Dust temperature in K.

    """

    def __init__(self, av, nh, tgas, tdust):
        """Create a cell class instance.

        Arguments
        ---------
        av : float
            Visual extinction in magnitudes.
        nh : float
            Hydrogen density in cm^-3.
        tgas : float
            Gas temperature in K.
        tdust : float
            Dust temperature in K.

        """

        self.data = libpyastrochem.Cell(av, nh, tgas, tdust)
        self.av = self.data.av
        self.nh = self.data.nh
        self.tgas = self.data.tgas
        self.tgas = self.data.tdust

class network:
    """Network class.

    Attributes
    ----------
    chem_file : str
        File containing a network to load.
    verbose : int
        Verbose if 1, Quiet if 0.

    """

    def __init__(self, chem_file, verbose):
        """Create a network class instance

        Arguments
        ---------
        chem_file : str
            File containing a network to load.
        verbose : int
            Verbose if 1, Quiet if 0.

        """

        self = libpyastrochem.Network(chem_file, verbose)

class phys:
    """Physical parameters to use in chemical reaction solver.

    Attributes
    ----------
    chi : float
        Chi physical property.
    cosmic : float
        Cosmic physical property.
    grain_abundance : float
        Grain Abundance physical property.
    grain_size : float
        Grain Size physical property.

    """

    def __init__(self):
        """Create a network class instance.

        """

        self.data = libpyastrochem.Phys()
        self.chi = self.data.chi
        self.cosmic = self.data.cosmic
        self.grain_abundance = self.data.grain_abundance
        self.grain_size = self.data.grain_size

class solver:
    """Chemical reaction solver.

    Attributes
    ----------
    cell : `cell`
        Chemical cell class to use in solver.
    chem_file : str
        Chemical network file string to load network from and use in solver.
    phys : `phys`
        Physical properties class to use in solver.
    abs_err : float
        Absolute acceptable error to use in solver.
    rel_err : float
        Relative acceptable error to use in solver.
    initial_abundances : dict
        Initial abundances (format {Species:Value}).
    density : float
        Density to use in solver.
    verbose : int
        verbose if 1, quiet if 0.

    """

    def __init__(self, cell, chem_file, phys, abs_err, rel_err, initial_abundances, density, verbose):
        """
        Create a solver instance.

        Arguments
        ----------
        cell : `cell`
            Chemical cell class to use in solver.
        chem_file : str
            Chemical network file string to load network from and use in solver.
        phys : `phys`
            Physical properties class to use in solver.
        abs_err : float
            Absolute acceptable error to use in solver.
        rel_err : float
            Relative acceptable error to use in solver.
        initial_abundances : dict
            Initial abundances (format {Species:Value}).
        density : float
            Density to use in solver.
        verbose : int
            verbose if 1, quiet if 0.

        """

        # Update the values of phys.data, in case phys values have changed
        phys.data.chi = phys.chi
        phys.data.cosmic = phys.cosmic
        phys.data.grain_abundance = phys.grain_abundance
        phys.data.grain_size = phys.grain_size

        self.data = libpyastrochem.Solver(cell.data, chem_file, phys.data, abs_err,
                                          rel_err, initial_abundances,
                                          density, verbose)

    def solve(self, time, new_cell=None):
        """
        Solve chemical reaction for a certain time

        Arguments
        ---------
        time : float
           Time to solve the system at
        new_cell : `cell`
           Cell class to use in solver, optionnal

        """

        if new_cell:
            return self.data.solve(time, new_cell.data)
        else:
            return self.data.solve(time, 0)
