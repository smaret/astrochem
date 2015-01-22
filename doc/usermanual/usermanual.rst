--------------

height 4pt width

| 
| Institut de Planétologie et d’Astrophysique de Grenoble

This manual documents Astrochem, a code to compute the abundance of
chemical species in the interstellar medium. It corresponds to Astrochem
version 0.7-beta. Please report any errors in this manual to
http://github.com/smaret/astrochem/issues.
Copyright © 2006-2014 Sébastien Maret

What is Astrochem?
==================

Astrochem is a code to compute the abundances of chemical species in the
interstellar medium, as function of time. It is designed to study the
chemistry in a variety of astronomical objects, including diffuse
clouds, dense clouds, photodissociation regions, prestellar cores,
protostars, and protostellar disks [1]_. Astrochem reads a network of
chemical reactions from a text file, builds a system of kinetic rates
equations, and solves it using a state-of-the-art stiff ordinary
differential equations (ODE) solver. The Jacobian matrix of the system
is computed implicitly, so the resolution of the system is extremely
fast: large networks containing several thousands of reactions are
usually solved in a few seconds. In addition, Astrochem may be run in
parallel on multi-cores and or multi-CPU computers. A variety of gas
phase process are considered, as well as simple gas-grain interactions,
such as freeze-out and desorption via several mechanisms (thermal
desorption, cosmic-ray desorption and photo-desorption). The computed
abundances are written in text files, and can be plotted in different
ways with the tools provided with Astrochem. Chemical reactions and
their rates are written in a format which is meant to be easy to read
and to edit. A tool to convert chemical networks from the OSU [2]_ and
KIDA [3]_ databases into this format is provided. Astrochem is written
in C, and its source code is distributed under the terms of the GNU
General Public License (GPL).

This manual documents Astrochem. It is organized as follows. The model
is described in section [sec:model]. The various chemical processes
considered in the code and the underlying assumptions are discussed in
this section. Astrochem usage is described in
section [sec:using-astrochem]. This section is written as a tutorial:
the different steps needed to computed abundances with Astrochem for a
simple case are described. Section [sec:input-source-files] gives a
comprehensive description of Astrochem input and source model files,
while section [sec:chemical-networks] presents the file format used for
chemical networks in Astrochem, as well as the different networks that
are provided with the code. Section [sec:runn-astr-parall] explains how
to run Astrochem in parallel. Finally, Appendix [sec:astrochem-authors]
gives a list of Astrochem authors, and Appendix [sec:contr-astr]
explains how you can contribute to the project.

The model
=========

Abundances in the interstellar medium are usually described by a system
of kinetic equations. For example, let us consider the HCO\ :math:`^{+}`
ion. In the dense interstellar medium, this ion is mainly formed through
the reaction:

.. math::

   \mathrm{H_{3}^{+} + CO \rightarrow HCO^{+} + H_{2}}
     \label{eq:hcop-formation}

with a rate :math:`k_{1}` (in cm\ :math:`^{3}` s:math:`^{-1}` units).
HCO\ :math:`^{+}` it is destroyed principally through the following
reaction:

.. math::

   \mathrm{HCO^{+} + e^{-} \rightarrow H + CO}
     \label{eq:hcop-destruction}

with a rate constant :math:`k_{2}`. If we neglect other formation and
destruction channels, the derivative of the HCO\ :math:`^{+}` density –
i.e. the number of HCO\ :math:`^{+}` per volume unit that we note
:math:`{n(\mathrm{HCO^{+}})}` in the following – writes as:

.. math::

   \frac{d{n(\mathrm{HCO^{+}})}}{dt} = k_{1} \, {n(\mathrm{H_{3}^{+}})} \, {n(\mathrm{CO})}
     - k_{2} \, {n(\mathrm{HCO^{+}})} \, {n(\mathrm{e^{-}})}

Similar kinetic equations can be written for all species; therefore to
compute the abundances as a function of time, one need to solve a system
of :math:`n_\mathrm{s}` equations for a set of initial conditions, with
:math:`n_\mathrm{s}` the number of species in the chemical network.
These are ordinary differential equations (ODE), that may be solved by
different techniques (Euler method, Runge-Kutta method, etc.). One
difficulty is that the ODE system is usually *stiff*, because of the
different timescales considered in the code; for example,
neutral-neutral reactions are typically several orders of magnitude
slower than a ion-neutral reactions. Another characteristic of the
system is that it is *sparse*; many species do not react together,
resulting in a large number of zeros in the Jacobian matrix of the
system [4]_.

Astrochem reads a set of reactions from a text file (see section 
[sec:chemical-networks] for a description a chemical network files),
build-up the system of kinetic equation and solve them using the CVODE
solver from the SUNDIALS library . CVODE uses the Backward
Differentiation formula (BDF method) with a variable step-size, which is
suitable for stiff systems. The resulting system of linear equations is
solved using the Newton iteration. The Jacobian matrix of the system is
computed explicitly to speed-up the computations. In addition, Astrochem
may be run in parallel on multi core/CPU computers (see
§ [sec:runn-astr-parall]). Astrochem includes a number of gas-phase
chemical processes as well as gas-grain interactions that we describe in
the following.

Gas-phase processes
-------------------

Gas-phase reactions
~~~~~~~~~~~~~~~~~~~

Most gas-phase reactions have rates that can be described by an Ahrrenus
law:

.. math::

   k = \alpha  \left( \frac{T}{300} \right)^\beta  \mathrm{exp} \left(
       -\frac{\gamma}{T} \right)
     \label{eq:arrhenius}

where :math:`T` is the gas temperature, and :math:`\alpha`,
:math:`\beta` and :math:`\gamma` are the rate constants. Usually
:math:`\gamma` corresponds to the energy barrier of the reaction,
expressed in Kelvins. It is generally equal to zero for ion-neutral
reactions, and equal or greater than zero for neutral-neutral reactions.
The units of :math:`k` depends on the order of the reaction: for a two
body reaction, which is of the second order, these are
cm \ :math:`^{3}` s:math:`^{-1}`.

Astrochem reads the :math:`\alpha`, :math:`\beta` and :math:`\gamma`
constants from the chemical network file (see
section [sec:network-file-format] for a description of the network file
format). The formation rate of each products of a given reaction is then
computed by multiplying the densities of the reactants by :math:`k`.
Similarly the destruction rate of each reactant is computed by
multiplying the densities of the reactants by :math:`k`. Reactions in
Astrochem may have up to three reactants, and four products.

Cosmic-ray ionization
~~~~~~~~~~~~~~~~~~~~~

Cosmic-ray particles can ionize molecules and atoms. This may happen in
a direct or indirect fashion. In the first case, the molecule (in
general H\ :math:`_{2}`) is ionized by a direct interaction with the
cosmic-ray particle. In the second case, the particle first ionizes
H\ :math:`_{2}`, forming H\ :math:`_{2}^{+}` and an electron. The
electron then recombines with H\ :math:`_{2}^{+}` and emit a UV photon.
This *secondary* UV photon may then ionize other molecules or atoms.
Astrochem assumes that the rate for these (either direct of indirect)
cosmic-ray ionization reactions scale with the H\ :math:`_{2}`
ionization rate :math:`\zeta`, such as:

.. math::

   k = \alpha  \, \zeta
     \label{eq:cosmic-ray-ionization}

The value of :math:`\zeta` is read from the input file (see
section [sec:input-file]). Typical values are comprised between
10\ :math:`^{-17}` and 10\ :math:`^{-15}` s:math:`^{-1}`. Note that in
this case the units of :math:`k` are s\ :math:`^{-1}` because cosmic-ray
ionization reactions are of the first order.

Photo-ionization and photo-dissociation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

UV photons from nearby stars may also dissociate and ionize molecules
and atoms. For sources with a plane-parallel or spherical symmetry, the
ionization or dissociation rate may be written as:

.. math::

   k = \alpha \, \mathrm{exp} \left( -\gamma A_{v} \right) \, \chi
     \label{eq:photo-ionization}

where :math:`A_{v}` is the visual extinction in magnitude, and
:math:`\chi` is the external UV flux in units of the standard Draine
interstellar radiation field .

This formulation implicitly assumes that the external radiation field
has the same spectral shape than the the ISRF. In addition the
*self-shielding* of species that dissociate through a line process is
neglected.

Gas-grain interactions
----------------------

H\ :math:`_{2}` formation on grains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the interstellar medium H\ :math:`_{2}` is mainly formed on dust
grains . The process is complex and involves the absorption of an H atom
in a grain site, the tunneling of the H atom from one site to the other,
and the reaction with another H to form H\ :math:`_{2}`. The energy
released during the reactions causes the evaporation of the
H\ :math:`_{2}` molecule, which returns to the gas phase. Astrochem uses
a simple treatment of this process. We assume that each H atom that
strikes a grain forms H\ :math:`_{2}` with a given efficiency. Under
this assumption, the formation rate of H\ :math:`_{2}` on the grains is
given by:

.. math::

   \frac{\mathrm{d} {n(\mathrm{H_{2}})}}{\mathrm{d} t} = k \, {n(\mathrm{H})}
     \label{eq:h2-formation-rate}

with:

.. math::

   k = \alpha  \left( \frac{T}{300} \right)^\beta
     \label{eq:h2-formation-rate-coeff}

The value of :math:`k` may be estimated by assuming that the efficiency
of the process is close to 1 (i.e. that each atom H that strikes a grain
forms an H\ :math:`_{2}` molecule). The rate coefficient is then simply
:math:`1/2` of the collision rate between H atoms and grains. For
0.1 \ :math:`\mu m` olivine grains and gas-to-dust mass ratio of 100, we
obtain a value of :math:`\sim 10^{-17}` s:math:`^{-1}` at 10 K. This is
close to the value of :math:`5
\times 10^{-17}` s:math:`^{-1}` determined observationally by . However,
because of the numerous uncertainties associated with the formation of
H\ :math:`_{2}`, in Astrochem the rate is not computed in this fashion.
Instead we use the :math:`\alpha` and :math:`\beta` values from the
network file, and compute it with Eq.([eq:h2-formation-rate-coeff]).

It is important to note that although the formation of H\ :math:`_{2}`
is a two body reaction – if we forget about the grain that only works as
a catalyst – this reaction has a first order kinetics: the formation
rate of H\ :math:`_{2}` depends on :math:`{n(\mathrm{H})}` and not on
:math:`{n(\mathrm{H})}^{2}`. Because of this, the reaction has its own
type number, 0 (see Table [tab:react-type-numb]). At present the
formation of H\ :math:`_{2}` on grains is the only grain surface
reaction that is considered in Astrochem.

Electron attachment and ion recombination on grains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Electron may hit grains and charge them. Ions may then recombine on
charged grains. For example, let us consider the following reactions:

.. math:: \mathrm{grain} + e^{-} \rightarrow \mathrm{grain}^{-}

.. math:: \mathrm{C^{+}} + \mathrm{grain}^{-} \rightarrow \mathrm{C} + \mathrm{grain}

The formation rate of charged grains writes as:

.. math:: \frac{d {n(\mathrm{grain^{-}})}}{\mathrm{d}t} = k_{1} \, {n(\mathrm{grain})} \, {n(\mathrm{e^{-}})}

while the recombination rate of :math:`\mathrm{C^{+}}` is:

.. math:: \frac{d {n(\mathrm{C^{+}})}}{\mathrm{d}t} = - k_{2} \, {n(\mathrm{grain^{-}})} \, {n(\mathrm{C^{+}})}

Both :math:`k_{1}` and :math:`k_{2}` are computed from the following
expression:

.. math::

   k = \alpha \left( \frac{T}{300} \right)^\beta \, \frac{n_\mathrm{H}}{n_\mathrm{d}}
     \label{eq:grain-attach-neutralization}

where :math:`n_\mathrm{H}` is the total hydrogen nuclei density [5]_ and
:math:`n_\mathrm{d}` is the total (neutral + charged) grain density. The
:math:`\frac{n_\mathrm{H}}{n_\mathrm{d}}` ratio is assumed to be
:math:`7.57 \times 10^{11}`, a value adequate for olivine grains of
0.1 \ :math:`\mu`\ m and a gas-to-dust mass ratio of 100 [6]_

:math:`k_{1}` may be estimated by assuming that each electron that hits
a grain will attach to it. For 0.1 \ :math:`\mu m` olivine grains and
gas-to-dust mass ratio of 100, we obtain a value of :math:`\sim
10^{-3}` cm:math:`^{-3}` s:math:`^{-1}` at 10 K. This process is
extremely fast because electron have large thermal velocities (thanks to
their small masses). In practice, in simulations free electrons almost
immediately stick on the grains, so that grains become negatively
charged very rapidly.

Depletion
~~~~~~~~~

Molecules may accrete on dust grains and freeze-out (a process often
called depletion). The formation rate of e.g. ices CO on the grains
through this process is given by:

.. math:: \frac{\mathrm{d}{n(\mathrm{CO_{ice}})}}{\mathrm{d}t} = k \, {n(\mathrm{CO})}

with:

.. math::

   k = S \, \pi r_{d}^2 \, v_{th} \, n_\mathrm{d}
     \label{eq:freeze-out}

and:

.. math::

   v_{th} = \left( \frac{8 k_{B} T_{d}}{\pi m} \right)^{1/2}
     \label{eq:thermal-veloc}

Here :math:`S` is a sticking probability (comprised between 0 and 1),
:math:`r_{d}` is the grain radius, :math:`v_{th}` is the thermal
velocity, :math:`n_{d}` is the total grain density (neutral + charged)
and :math:`m` is the mass of the accreting species .

Because no grain destruction or formation mechanisms are considered in
Astrochem, :math:`n_{d}` does not varies with time. It is therefore
computed from the initial abundances of neutral and charged grains given
in the input file (see section [sec:initial-abundances]). The grain size
:math:`r_{d}` is also read from this file. Both :math:`S` and :math:`m`
are read from the network file.

Thermal desorption
~~~~~~~~~~~~~~~~~~

Once frozen on the dust grains, molecules may evaporate through thermal
or non-thermal processes. The formation rate of gaseous CO by CO ices
thermal evaporation is:

.. math:: \frac{\mathrm{d}{n(\mathrm{CO})}}{\mathrm{d}t} = k \, {n(\mathrm{CO_{ice}})}

where :math:`k` is given by the Polanyi-Wigner equation:

.. math::

   k = \nu_{0} \, \mathrm{exp} \left( - \frac{E_{B}}{T_{d}} \right)
     \label{eq:thermal-desorption}

with:

.. math::

   \nu_{0} = \left( \frac{2 N_{S} E_{B}}{\pi^2 m} \right)^{1/2}
     \label{eq:vibration-freq}

Here :math:`\nu_{0}` is the characteristic vibrational frequency of the
desorbing species, :math:`E_{B}` is the binding energy of the desorbing
species on the grain surface expressed in Kelvins, :math:`T_{d}` is the
grain temperature and :math:`N_{S}` is the number of sites per unit
surface assumed to be :math:`3 \times 10^{15}` cm:math:`^{-2}` . The
values of :math:`E_{b}` and :math:`m` are both read from the network
file.

Cosmic-ray desorption
~~~~~~~~~~~~~~~~~~~~~

As mentioned above, ices may also evaporate by non-thermal processes.
For example, cosmic-rays may desorb molecules from grains, either by
creating hot-spots on the grain surfaces, or by heating the whole grains
. Because the energy deposited in a grain varies as :math:`Z^{2}`,
cosmic-ray desorption in mainly caused by heavy cosmic-ray ions, such as
Fe. suggested that desorption by spot-heating dominates over desorption
by whole-grain heating for grains smaller than 2.5 :math:`\mu`\ m.
However, recent molecular dynamics simulations indicate that for 0.1
:math:`\mu`\ m grains the whole grain heating contribution is small .

Because of the uncertainties on this process, two different treatments
are implemented in Astrochem. First, cosmic-ray desorption rates can be
computed following , who assume that desorption occurs mostly through
whole-grain heating; when impacting grains, heavy cosmic-ray ions are
assumed to impulsively heat the grains to a peak temperature of 70 K, at
which most of the desorption occurs. The rate is then similar to that of
thermal desorption:

.. math::

   k = f \, \nu_{0} \, \mathrm{exp} \left( -\frac{E_{B}}{70} \right)
     \label{eq:cosmic-ray-desorption}

where :math:`f` is the fraction of the time spent by a grain in the
vicinity of 70 K between two cosmic-ray heating events, assumed to be
:math:`3.16
\times 10^{-19}` .

Alternatively, the cosmic-ray desorption rate of any specie can be given
explicitly in the network file. This allows for the use of the
cosmic-ray desorption rates that have been computed and/or measured for
some species . The user can specify in the network file which treatment
to use for each species (see Table. [tab:rate-const-meaning]). Note that
no scaling of :math:`k` with the cosmic ray ionization rate is
performed.

Photo-desorption
~~~~~~~~~~~~~~~~

Photo-desorption (i.e. desorption by UV photons) is another non-thermal
desorption process. UV photons can originate in the ISRF, or in the
ionization of H\ :math:`_{2}` by cosmic-rays followed by recombination
(secondary UV photons). At present only photo-desorption from ISRF UV
photons is implemented in Astrochem.

The photo-desorption rate of CO is for example :

.. math::

   \frac{\mathrm{d} {n(\mathrm{CO})}}{\mathrm{d} t} = k
     \label{eq:photo-desorption-rate}

with:

.. math::

   k = \chi \, I_\mathrm{ISRF,FUV} \, \mathrm{exp} \left( -2 A_{v} \right)
     \, \pi r_{d}^{2} \, n_{d} \, Y_\mathrm{PD}
     \label{eq:photo-desorption}

Here :math:`I_\mathrm{ISRF,FUV}` is the standard interstellar radiation
field in the FUV , and :math:`Y_\mathrm{PD}` is the photo-desorption
yield, i.e. the number of molecules ejected per incident photon. The
latter is given by:

.. math::

   Y_\mathrm{PD} = Y_{0} \left[ 1 - \mathrm{exp} \left( -x / l \right) \right]
     \label{eq:photo-desorption-yield}

where :math:`x` is the ice thickness of the considered species expressed
in monolayers (ML), :math:`l` is the diffusion length in ML, and
:math:`Y_{0}` is the photo-desorption yield for thick ices (i.e.
:math:`x \gg l`). Typical values for :math:`Y_{0}` and :math:`l` are
:math:`10^{-3}` molecules photon\ :math:`^{-1}` and 2 ML,
respectively [7]_. The density of e.g. CO ices is given by:

.. math::

   x = \frac{{n(\mathrm{CO_{ice}})}}{N_{s} \, \pi r_{d}^2 \, n_{d}}
     \label{eq:ice-thickness}

It is interesting to note that for thick ices, photo-desorption is
zeroth order process: the desorption rate does not depends on the amount
of e.g. CO ices on the grains. This because UV photons can penetrate
only the first ices monolayers; the bulk of ice is not affected. On the
other hand, for thin ices (i.e. :math:`x \ll l`) the desorption rate
become linearly proportional to the ice thickness, and therefore on the
ice abundance. Consequently for thin ices, photo-desorption is a first
order process.

Astrochem follows the ice thickness of each species as a function of
time. The desorption rate is then computed from the above equations,
using the values of :math:`\chi`, :math:`r_{d}^2` and :math:`n_{d}` from
the input file, the :math:`A_{v}` from the source file, and the
:math:`Y_{0}` and :math:`l` from the network file.

Using Astrochem
===============

In this section we present a simple example of Astrochem usage. We
propose to use Astrochem to study the formation of the HCO\ :math:`^{+}`
ion in a dense interstellar cloud. We suppose that the cloud is isodense
and isothermal, and that it is shielded from the ISRF, so that
photo-processes can be ignored. For a sake of simplicity, we also
neglect the freeze-out of molecules on dust grains. In the following, we
describe the various steps needed to solve this problem.

Describing the problem
----------------------

In order to describe our problem, we first need create an input file
that contains the various parameters the code. This file has several
sections, that set the physical parameters (e.g. the cosmic ionization
rate), the solver parameters (e.g. the initial and final time in the
computation), the initial abundances, and which species we want in
output. Some of these parameters are optional; if they are not specified
in the input file, Astrochem will use a default value that should be
suitable for most problems. Here is a what the input file for our
example problem looks like (for a comprehensive description of the
parameters in input files and their default value, see
section [sec:input-file]):

::

    [files]
    source = source.mdl
    chem = osu2009.chm
    # Physical paramaters
    [phys]
    chi = 1.0
    cosmic = 1.3e-17
    # Solver parameters
    [solver]
    ti = 1e-6
    tf = 1e7
    # Initial abundances
    [abundances]
    H2      = 0.5
    He      = 0.14
    N       = 2.14e-5
    O       = 1.76e-4
    C(+)    = 7.30e-5
    S(+)    = 8.00e-8
    Si(+)   = 8.00e-9
    Fe(+)   = 3.00e-9
    Na(+)   = 2.00e-9
    Mg(+)   = 7.00e-9
    P(+)    = 2.00e-10
    Cl(+)   = 1.00e-9
    F       = 6.68e-9
    e(-)    = 7.31012e-5
    # Output
    [output]
    abundances = H3(+),e(-),CO,HCO(+)
    trace_routes = 1

The various sections of the file are indicated by keywords within
brackets. Lines starting with ``#`` are comments. The first section
(``[files]``) indicates the name of the file describing our source
(``source``), and the chemical network to use (``chem``). The following
section (``[phys]``) sets the physical parameters of the source. Here we
set the UV radiation field in Draine units (``chi``) to 1.0, and the
cosmic ionization rate (``cosmic``) to
:math:`1.3 \times 10^{-17} s^{-1}`. The solver parameters are set in
following section (``[solver]``). ``ti`` and ``tf`` are the initial and
final time for the calculation respectively. Both are expressed in
years. The ``[abundance]`` section sets the initial abundances;
abundances that are not specified are set to zero. The last section
(``[output]``) sets parameters relative to the output of the code.
``abundances`` set the name of the species for which we want to create
output containing the abundances as a function of time. ``ALL`` can be
used instead so all species are available in output. ``trace_route`` is
an optional parameter that allow to trace the various formation and
destruction routes of these species.

In addition to the input file, we need to provide a file describing our
source. The file corresponding to our problem looks like this (for more
information on source model files, see section [sec:source-file]):

::

    # Source model file example
    # cell number, Av [mag], nH [cm^-3], Tgas [K], Tdust [K]
    #
    0   20.0    1e+04   10.0    10.0

As for the input file, lines that starts with a ``#`` are comments. The
file contains one line for each *cell* of our source. In this simple
example, our source is isodense and isothermal, and therefore there is
only one cell in the our source file. A more realistic source with a
temperature and density gradient would be sampled in more cells.

Each line corresponding to a cell has five columns. The first column is
the index of the cell, the second one is the visual extinction
:math:`A_{v}` (in magnitudes), the third one is the number density (in
cm\ :math:`^{-3}`). and the fourth and the fifth are the gas and dust
temperature respectively (in Kelvins). Note that we have adopted a large
visual extinction (20 mag.), that is suitable for dense clouds shielded
from the ISRF.

Running Astrochem
-----------------

Astrochem is run from the command line, and takes the name of the input
file as an argument:

::

    % astrochem input.ini
    Reading input from input.ini.
    Reading source model from source.mdl.
    Reading reactions network from osu2009.chm... done.
    Found 6046 reactions involving 468 species.
    Computing abundances in cell 0...
    Done with cell 0.
    Writing abundances in output files... done.
    Writing formation/destruction in output files... done.
    %

Astrochem produces one hdf5 output file ``astrochem_output.h5``
containing different datasets. A list a output species names, a list of
time steps, a matrix of abundances, and if the ``trace_route`` parameter
is set to 1, a group of matrices containing the formation and
destruction routes.

Plotting abundances
-------------------

Astrochem comes with a program that make plots of the abundances
computed by Astrochem. The program, named Plabun, allows to plot the
abundances of one or several species as of function of time in a given
cell. For example, the following command plots the CO,
H\ :math:`_{3}^{+}`, e\ :math:`^{-}` and HCO\ :math:`^{+}` abundances:

::

    % plabun --xrange=1,1e7 --yrange=1e-12,1e-4 astrochem_output.h5 CO
    H3(+) e(-) HCO

In the example above, we have set the x-axis range from 1 to
10\ :math:`^{7}` yr and y-axis range from :math:`10^{-12}` to
:math:`10^{-4}` with the ``--xrange`` and ``--yrange`` options,
respectively. Note that specie name containing special character are
handled smoothly by the program. The command above produces the plot
shown on Fig. [fig:example-abundances]. Plabun has a number or other
options, including a legacy mode to read old output format ``.abun``;
see ``man plabun`` for a complete list.

.. figure:: fig1.pdf
   :alt: Abundances as a function of time for the example problem.

   Abundances as a function of time for the example problem.
[fig:example-abundances]

Identifying the main formation and destruction channels of a species
--------------------------------------------------------------------

It is often useful to identify the main formation and destruction routes
of a given species, for example to check that the rate of these main
reactions have been determined accurately (e.g. from experiments) or are
rather uncertain. This also allows to understand how the various species
of a chemical network are linked together.

As already mentioned, if one turns on the ``trace_route`` option in the
input file, Astrochem creates some datasets in it’s ouput file
containing the main formation and destruction routes of the species
listed with ``output`` option. Of course these routes may change with
time; therefore Astrochem saves the sixteen most important formation
routes as well as the sixteen most important destruction routes at each
time step and position.

A command allows to plot the main formation and destruction rate of a
given species as a function of time or position [8]_. For example, one
can plot the main formation and destruction routes of HCO\ :math:`^{+}`
with the following command (see ``man plroute`` for a complete list of
commands and options, including a legacy mode allowing to read old
``.rout`` format ):

::

    % plroute astrochem_output.h5 HCO(+)

.. figure:: fig2.pdf
   :alt: Main HCO\ :math:`{+}` formation and destruction routes as
   function of time for the example problem.

   Main HCO\ :math:`{+}` formation and destruction routes as function of
   time for the example problem.
[fig:example-routes]

which produces the plot shown on Fig. [fig:example-routes].

The left panel of the plot shows the formation rate of HCO\ :math:`^{+}`
(in cm\ :math:`^{-3}` s:math:`^{-1}`) through the six most important
formation channels [9]_, together with the total formation rate.
Conversely, the right panel shows the destruction rate of the same
species through the six most important destruction channels as well as
the total destruction rate. On the left panel, we see that for
:math:`t > 10^{5}` yr, the formation of HCO\ :math:`^{+}` is dominated
by the reaction of CO with H\ :math:`_{3}^{+}`. On the other hand, at
any time in the simulation the destruction of HCO\ :math:`^{+}` is
dominated by the dissociative recombination with electrons.

Note that because Astrochem only keeps tracks of the sixteen most
important formation or destruction routes at any time, some gaps may
appear in the plot. This can been seen for the reaction in red on the
left panel of the plot for times between roughly 0.1 and 10 years. In
this time range this reaction is not one the sixteen most important
formation reaction, so Astrochem does not keep track of it in this time
range.

Input and source files
======================

Input file
----------

As mentioned already, the various parameters of Astrochem are read from
an input file. Although this is not mandatory, the file usually has the
``.ini`` file extension. The file has several sections that are
delimited by a keyword within brackets (e.g. ``[files]``). Each section
has a number of parameters that we describe in the following.

Files
~~~~~

This section of the file starts with the ``[files]`` keyword, and
specifies which file Astrochem should use for the source description
(``.mdl`` file) and chemical network (``.chm`` file). The parameters
allowed in this section are:

-  ``source``: The name of the file describing the source.

-  ``network``: The name of the chemical network file. First Astrochem
   searches for this file in the current directory. If not found, it
   then searches for the file in Astrochem data installation directory
   (``/usr/local/share/astrochem`` by default).

Physical parameters
~~~~~~~~~~~~~~~~~~~

This section of the file starts with the ``[phys]`` keyword and
specifies the physical parameters of the problem. These parameters are:

-  ``chi``: The external UV radiation field, expressed in units. This
   corresponds to :math:`\chi` in Eq. ([eq:photo-ionization]) and
   Eq. ([eq:photo-desorption]).

-  ``cosmic``: The cosmic ray ionization rate of molecular hydrogen
   expressed in s\ :math:`^{-1}`. The default value is :math:`1.3 \times
     10^{-17}`. This corresponds to :math:`\zeta` in
   Eq. ([eq:cosmic-ray-ionization]).

-  ``grain_size``: The grain radius in microns. The default value is
   0.1. This corresponds to :math:`r_{d}` in Eq. ([eq:freeze-out]) and
   Eq. ([eq:photo-desorption]).

Solver parameters
~~~~~~~~~~~~~~~~~

This section of the file starts with the ``[solver]`` and specifies the
ODE solver parameters. These parameters are:

-  ``ti``: The initial time for the computation, expressed in years. The
   default value is :math:`10^{-6}`.

-  ``tf``: The final time for the computation, expressed in years. The
   default value is :math:`10^{7}`.

-  ``abs_err``: The solver absolute error (or tolerance) on the computed
   abundances. The default value is :math:`10^{-20}`.

-  ``rel_err``: The solver relative error on the computed abundances.
   The default value is :math:`10^{-6}`.

A note on tolerances: Astrochem adjusts the internal time step so that
the relative error on any abundance is always lower that ``rel_err``,
unless the given abundance is lower that ``abs_err``. Because errors on
the abundances at each time step may add-up, we recommend to chose these
errors quite conservatively. The default values should be suitable for
most problems.

Initial abundances
~~~~~~~~~~~~~~~~~~

This section specifies the initial abundances in the computation. Each
line should contain a specie name followed by a equal sign and the
initial abundance with respect to H nuclei. The initial abundances of
species that are not listed in this section are assumed to be zero. Note
that grains must be written as ``grain``, ``grain(-)`` or ``grain(+)``
so that the total grain density :math:`n_{d}` is computed correctly (see
section [sec:depletion]).

Output
~~~~~~

This section specifies what file Astrochem should create at the end of
the computation. These parameters are:

-  ``output``: The name of the species for which Astrochem creates an
   output file containing the abundance as a function of time and
   position. Species names must be separated by a comma. File names are
   formed with the species name, possibly followed by a suffix (see
   below) and the ``.abun`` file extension.

-  ``suffix``: A suffix to append to the name of the species before the
   file extension (``.abun`` or ``.reac``) in abundance files. This is
   useful when you want to run Astrochem for a number of different input
   files all located in the same directory; this way the results of a
   given simulation will not be overwritten by the results of others. A
   leading underscore is added to this suffix.

-  ``time_steps``: The number of time steps in output the files [10]_.
   The default value is 32. If plots of abundances v.s. time appear
   somewhat boxy, you may increase this number. Note that this parameter
   only affects the number of time steps in the output files. The
   internal time step size is set by the ODE solver in order to reach
   the specified absolute and relative errors on the abundances.

-  ``trace_routes``: This parameter is used to toggle the computation of
   the major formation and destruction routes for all species listed
   with the ``output`` parameter. If ``trace_route`` is set to 1,
   Astrochem will create a file containing the formation/destruction
   rate and reaction number of the 16 most important
   formation/destruction reaction for each specie, as a function of time
   and position (i.e. cell number). As for abundance files, file names
   for formation and destruction routes are formed with the species name
   possibly followed by a suffix (see below) and the ``.rout`` file
   extension.

Source file
-----------

Astrochem reads the physical parameters (density, temperature, visual
extinction) of the astronomical source in a source file. This file can
be of two formats, depending on whether the source physical parameters
vary as a function of time or not.

Time-independent physical parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This type of source file is used to describe astronomical source in
which the physical parameters do not change with time, or change on
timescales longer than the chemical timescales. Although the source may
have in principle any dimension, Astrochem is, for the moment, designed
to study 1D spherical sources in an external radiation field only (such
as a dense cloud or a prestellar core in the ISRF). Future versions will
allow to study 2D sources with axisymmetrical geometries, such as
protoplanetary disks.

In order to construct a time-independent source model file, one needs to
sample the astronomical source in a number of spherical cells (or
shells) at different source radius. What constitutes a good sampling
depends on the source. Often density profiles of astronomical sources
are well described by a power laws, so it is usually a good idea to
sample the source in a number of logarithmically spaced cells. Of
course, the larger number of cells, the longer computational time.
However, Astrochem may be run in parallel on multi-core computers in
order to reduce the computational time (see
section [sec:runn-astr-parall]).

Each line of the file corresponds to a different cell, while each column
corresponds to a different parameter. These parameters are, from the
leftmost to the rightmost columns:

#. The cell index. This is an integer that is used to identify each
   cell. The first index should be 0. Other indexes in the file should
   be in increasing order. All cells should have a different index.

#. The visual extinction in the cell, expressed in magnitudes.

#. The H nuclei density in the cell, expressed in cm\ :math:`^{-3}`.

#. The gas temperature in the cell, expressed in K.

#. The dust temperature in the cell, expressed in K.

#. The radius corresponding to the cell, expressed in astronomical units
   (AU). This optional parameter is used for bookkeeping only; Astrochem
   ignores it.

Columns may be separated by any number of white spaces or tabs. Comments
may written in the source file; comment lines must start with a ``#``
sign.

Time-dependent physical parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This type of source file is used to describe an astronomical source in
which the physical parameters vary as a function of time, for example a
protostar that undergoes gravitational collapse. In this case, the
physical parameters (density, temperature), may come from theoretical
prescriptions or from numerical simulations.

In order to construct this kind of model, one need to sample the source
in a number of different cells that follow the dynamical evolution of
the object as a function of time. These could correspond to the *smooth
particles* from a SPH simulation, or the *buoy particles* from an
adaptive mesh refinement (AMR) (magneto-)hydrodynamical simulation. The
source needs to be properly sampled both spatially and temporally.
Indeed, Astrochem assumes that the physical parameters in each cell
remain constant within each time step [11]_. Therefore this time step
should be sufficiently small so that the physical parameters do not vary
significantly between two time steps.

Time-dependent source models have two sections, which are separated by
keywords in brackets. The first one contains the times, expressed in
years, and starts with the ``[times]`` keyword.

::

    # Source model file example for a time-dependant source structure
    [times]
         0     1.00e-06
         1     1.27e-06
         2     1.60e-06
         3     2.03e-06
         4     2.57e-06
    (...)
       126     7.90e+06
       127     1.00e+07

The number on the first column is the time index, which must start at 0.
The second column contains the time, expressed in years. In this
example, the time on the first line is :math:`10^{-6}` years. This
corresponds to the time at the end of the computation of the first time
step; in other words, the time step ``#0`` is between :math:`t = 0` and
:math:`t = 10^{-6}` years. Likewise, the time step ``#127`` is between
:math:`7.9 \times 10^{6}` and :math:`10^{7}` years, at which time the
computation ends.

The second part of the file contains the physical parameters of each
cell at each time step. This section must start with the ``[cells]``
keyword, followed by the physical properties of the first cell.

::

    # cell number, time index, Av [mag], nH [cm^-3], Tgas [K], Tdust [K]
    [cells]
         0        0    20.00     1.00e+04    10.00    10.00
         0        1    20.00     1.04e+04    10.00    10.00
         0        2    20.00     1.08e+04    10.00    10.00
    (...)
         0      127    20.00     1.00e+06    10.00    10.00
         1        0    20.00     1.00e+05    10.00    10.00
    (...)
         1      127    20.00     1.00e+07    10.00    10.00

The columns in this section are:

#. The cell index, starting at 0.

#. The time index.

#. The visual extinction in the cell, expressed in magnitudes.

#. The hydrogen density, in cm\ :math:`^{-3}`.

#. The gas temperature, in Kelvins.

#. The dust temperature, in Kelvins.

Each line correspond to a different time step. In this example, the
first cell (``#0``) has a density of :math:`10^4` cm:math:`^{-3}` and a
temperature of 10 K during the first time step (``#0``), i.e. between
:math:`t = 0` and :math:`t = 10^{-6}` years. During the last time step
(``#127``), this cell has a density of :math:`10^6` cm:math:`^{-3}` and
a temperature of 10 K.

Any number of cells may be given in the file. The physical parameters of
the second cell (``#1``), directly follows that of the first one. For
example, the second cell (``#1``) has a density of
:math:`10^4` cm:math:`^{-3}` and a temperature of 10 K during the same
time step (``#0``). Note that all cells must have the same time steps.

Chemical networks
=================

Astrochem reads the chemical reactions and their rate coefficient from a
chemical network file. This file should have ``.chm`` extension. Several
networks are distributed with Astrochem, but the user may also write its
own network. In the following, we describe the networks files that are
distributed with Astrochem, as well as the format of these files.
Finally, we explain how networks in other formats can be converted to
Astrochem format.

Networks provided with Astrochem
--------------------------------

The following networks are provided with Astrochem:

-  ``osu2009.chm``: This network file contains the reactions and rates
   from the Ohio State University (OSU) astrochemistry database, that is
   maintained by Eric Herbst. It corresponds to the January 2009 version
   of the database. This network contains 6046 reactions and 468
   species, including anions.

-  ``osu2008.chm``: This network file contains the September 2008
   version of OSU database. It includes 4457 reactions and 452 species
   (no anions).

| These networks can be found in the network directory of the source
distribution. When installing Astrochem, these file are copied in the
data installation directory
| (``/usr/local/share/astrochem`` by default).

Network file format
-------------------

Astrochem network file format is meant to be easily read and edited.
Therefore chemical reactions in this files are written as you would
write them on a piece of paper. Here is an example of a (incomplete)
network file:

::

    # A few reactions extracted from osu2008.chm
    H     + H          -> H2                 4.95e-17  5.00e-01  0.00e+00  0    1
    H2    + cosmic-ray -> H2(+)  + e(-)      9.30e-01  0.00e+00  0.00e+00  1   39
    H3(+) + CO         -> HCO(+) + H2        1.61e-09  0.00e+00  0.00e+00  2 1756
    H3(+) + e(-)       -> H      + H    + H  4.36e-08 -5.20e-01  0.00e+00  9 3746
    CO    + uv-photon  -> C      + O         3.10e-11  0.00e+00  2.54e+00 13 4297
    (...)

As for input and model files, lines that starts with the ``#`` character
are comments. Each line of the file corresponds to a different reaction.
Lines have two parts: the chemical equation, and a list of five numbers
that correspond to the rate constants, reaction type and reaction
number.

The chemical equation is composed of one, two or three reactants, and
one, two, three or four products. Each reactants and products are
separated by a white space, a ``+`` sign, and another white space. To
disentangle the ``+`` sign between reactants from the ones corresponding
to ions, ion charge must be put in parenthesis: for example, the
HCO\ :math:`^{+}` ion must be written as ``HCO(+)`` in the network file.
Reactants and products are separated by a white space, an arrow
(``->``), and another white space.

In general, Astrochem computes the formation rate of each product (or
the destruction rate of each reactant) through a given reaction by
multiplying the reaction rate by the product of the reactants. For
example, for the reaction:

::

    H3(+) + e(-)       -> H      + H    + H  4.36e-08 -5.20e-01  0.00e+00  9 3746

the destruction rate of H\ :math:`_{3}^{+}` and e\ :math:`^{-}` is
computed as:

.. math::

   \frac{d {n(\mathrm{H_{3}^{+}})}}{dt} = \frac{d {n(\mathrm{e^{-}})}}{dt} = - k \,
     {n(\mathrm{H_{3}^{+}})} \, {n(\mathrm{e^{-}})}

while the formation rate of H is computed as:

.. math:: \frac{d {n(\mathrm{H})}}{dt} = k \, {n(\mathrm{H_{3}^{+}})} \, {n(\mathrm{e^{-}})}

In other words, one body reactions (e.g. cosmic-ray ionization, or UV
ionization) are assumed to have a first order kinetics, two body
reactions are assumed to have a second order kinetics, etc. However,
there are two exceptions to this rule. First the formation of
H\ :math:`_{2}` on the grains is assumed to be of the first order (see
section [sec:h\ :sub:`2`-formation-grains]). Second, UV photodesorption
is assumed to have a zeroth order kinetics when the ice thickness is
large enough (see § refsec:photo-desorption).

Several reactants and products are not *bona fide* chemical species, but
are just meant to make the reading of the network easier. These
*pseudo-species* are ``cosmic-ray`` (for cosmic-ray direct or indirect
ionization or desorption reactions), ``uv-photon`` (for
photo-ionization, photo-dissociation and photo-desorption reactions) and
``photon`` (for radiative association reactions). All of these are
ignored by Astrochem.

The five numbers following the products are the three rate constants
(that we note :math:`a`, :math:`b` and :math:`c` respectively), the
reaction type, and the reaction number. The reaction type is a signed
integer that identifies the kind of the reaction (e.g. ion-molecule,
dissociative recombination, etc.). Table [tab:react-typea-numb] lists
the various reaction types together corresponding type numbers [12]_.
The reaction number is an integer that identify each reaction in a
unique fashion: every reaction must have a different number. Reaction
numbers must start at 1, but they are not necessarily contiguous. For
example you may want to identify gas-phase reactions by numbers between
1 and 6046, and gas-grain reactions by numbers starting at 10000.

+---------------+---------------------------------------------------------------+
| Type number   | Reaction type                                                 |
+===============+===============================================================+
| -1            | Electron attachment and ion recombination on grains           |
+---------------+---------------------------------------------------------------+
| 0             | H\ :math:`_{2}` formation on grains                           |
+---------------+---------------------------------------------------------------+
| 1             | Cosmic-ray ionization or cosmic-ray induced photo-reactions   |
+---------------+---------------------------------------------------------------+
| 2             | Ion-molecule reactions, Charge exchange reactions             |
+---------------+---------------------------------------------------------------+
| 3             | Negative ion - neutral species reactions                      |
+---------------+---------------------------------------------------------------+
| 4             | Radiative association                                         |
+---------------+---------------------------------------------------------------+
| 5             | Associative ejection                                          |
+---------------+---------------------------------------------------------------+
| 6             | Neutral + Neutral :math:`\rightarrow` ion + electron          |
+---------------+---------------------------------------------------------------+
| 7             | Neutral-Neutral chemical reactions                            |
+---------------+---------------------------------------------------------------+
| 8             | Neutral-Neutral radiative association                         |
+---------------+---------------------------------------------------------------+
| 9             | Dissociative recombination                                    |
+---------------+---------------------------------------------------------------+
| 10            | Radiative recombination                                       |
+---------------+---------------------------------------------------------------+
| 11            | Positive ion - Negative ion recombination                     |
+---------------+---------------------------------------------------------------+
| 12            | Electron attachment                                           |
+---------------+---------------------------------------------------------------+
| 13            | Photo-ionization, Photo-dissociation                          |
+---------------+---------------------------------------------------------------+
| 20            | Freeze-out on grains                                          |
+---------------+---------------------------------------------------------------+
| 21            | Thermal desorption                                            |
+---------------+---------------------------------------------------------------+
| 22            | Cosmic-ray induced desorption                                 |
+---------------+---------------------------------------------------------------+
| 23            | Photo-desorption                                              |
+---------------+---------------------------------------------------------------+

Table: Reaction types and their corresponding type numbers.

[tab:react-type-numb]

Astrochem computes the rate of each reaction from the :math:`a`,
:math:`b` and :math:`c` rate constants. The physical meaning of these
constants depends on the type of the reaction.
Table [tab:rate-const-meaning] gives physical meaning of these for each
reaction type, as well as the equation that is used to computed the
rate. Because reaction rates generally do not vary with time, Astrochem
computes them only once for each cell. The only exception is the rate of
UV photo-desorption reactions, that depends on the ice thickness;
therefore the rate for these reactions is computed for each solver
internal time step.

| llccc Type number & Equation &
|  & & :math:`a` & :math:`b` & :math:`c`
| -1 & Eq. ([eq:grain-attach-neutralization]) & :math:`\alpha` &
:math:`\beta`
| 0 & Eq. ([eq:h2-formation-rate]) & :math:`\alpha` & :math:`\beta`
| 1 & Eq. ([eq:cosmic-ray-ionization]) & :math:`\alpha` & - & -
| 2-12 & Eq. ([eq:arrhenius]) & :math:`\alpha` & :math:`\beta` &
:math:`\gamma`
| 13 & Eq. ([eq:photo-ionization]) & :math:`\alpha` & - & -
| 20 & Eq. ([eq:freeze-out]) & :math:`S` & :math:`m/m_\mathrm{H}` & -
| 21 & Eq. ([eq:thermal-desorption]) & - & :math:`m/m_\mathrm{H}` &
:math:`E_{b}`
| 22 & Eq. ([eq:cosmic-ray-desorption])\ :math:`^\dagger` & :math:`k` &
:math:`m/m_\mathrm{H}` & :math:`E_{b}`
| 23 & Eq. ([eq:photo-desorption]) & :math:`Y_{0}` & - & :math:`l`

| 
|  [tab:rate-const-meaning]

Convert networks to Astrochem format
------------------------------------

Astrochem comes with a tool called Chmconvert that converts network
files to Astrochem native format (``.chm``). Chmconvert supports the OSU
database and the KIDA formats. The file format is determined from the
file extension, which must be ``.osu`` for Ohio State University
database files, and ``.kida`` for KIDA. Networks can be converted as
follows (see ``man chmconvert`` for more information):

::

    % chmconvert -o osu2009.chm osu2009.osu

The ``-o`` option is used to select the name of the output file. If not
specified, the network is copied to the standard output. Two networks
provided with Astrochem (``osu2008.chm`` and ``osu2009.chm``; see
section [sec:netw-prov-with]) are automatically generated from the OSU
files using this tool.

Running Astrochem in parallel
=============================

Astrochem may be run in parallel for sources containing more than one
cell. In this kind of source, if diffusion and advection are neglected,
each cell is independent from others: abundances in a given cell are not
affected by abundances in other cells. Therefore abundances may be
computed in several cells simultaneously.

Parallelism in Astrochem is implemented using the OpenMP standard [13]_.
Basically, OpenMP is a set of compiler directives that allows a program
to be run in parallel on a shared-memory parallel computer, such as a
multi-core and/or multi-CPU machine. Astrochem forks in several threads
that are run simultaneously on several cores and/or CPUs. Each thread
corresponds to a computation in a different cell. When run in parallel,
Astrochem execution time scales almost linearly with the inverse of the
number of cores or CPUs. For example, for a 64 cells source, Astrochem
should run almost eight times faster on a quad-core dual-CPU (i.e. eight
cores in total) than on single-core single-CPU computer.

The compilation of the parallel version of Astrochem can be turned on by
specifying the ``--enable-openmp`` option to the configure script, e.g.
(see the ``INSTALL`` file in the source distribution for more
information):

::

    ./configure --enable-openmp

The configure script will attempt to detect if your compiler supports
the OpenMP standard [14]_. Then one need to set the ``OMP_NUM_THREADS``
environment variable to the number of threads to be run in parallel. It
is recommended to set this variable the number of cores on the machine
(e.g. 8 for a eight core computer). With the bash shell, this is done as
follows:

::

    export OMP_NUM_THREADS=8

It is also a good idea to sample the source in a number of cell that is
a multiple of the number of threads (e.g. 64 cells for 8 threads) so
that all cores 100% of the time during the computation.

Calling Astrochem from another code
===================================

Astrochem can be called from another code through an application
programming interface (API). This allows to compute the chemical
evolution of a gas cell within another code, typically an hydrodynamic
or magneto-hydrodynamic (MHD) simulation. In practice, Astrochem
functions that are part of the API are included in a dynamic library,
against which other codes may be linked. For the moment, only C is
supported.

A typical call of Astrochem library from a C code can be decomposed in
several steps. First, one must select a chemical network file, which is
done as follows:

::

    #include <libastrochem.h>

    int verbose = 0;
    char *chem_file = "osu2009.chm";
    net_t network;
    read_network (chem_file, &network, verbose);

This creates a network structure of type ``net_t``\  [15]_. In
principle, one should check the return code of the ``read_network ()``
function to make sure that no error occured while reading the network
file. This is not done in this example for a sake of clarity.

The second step is to set the physical parameters of the model. This is
done by defining a structure of type ``phys_t``, and by setting the
different elements of this structure:

::

    phys_t phys;
    phys.cosmic = 1e-17;
    phys.chi = 0;
    phys.grain_abundance = 0;

The parameters are the same than in input file. Please refer to
section [sec:input-file] of this manual for details. As for input files,
parameters that are not set explicitly are set to their default value.

The third step is to allocate a work array that will contain the
abundances of all species at a given time. We also need to set the
initial abundances:

::

    const char* species[]  = {"CO", "HCO(+)", "e(-)"};
    const double initial_abundances[] = {1e-4, 1e-9, 1e-9};
    double *abundances;
    alloc_abundances (&network, &abundances);
    set_initial_abundances (species, 3, initial_abundances, &network,
                            abundances);

Obviously, the species in the ``char*`` array must be present in the
network.

The fourth step is to set the initial density, visual extinction and
temperature of the gas cell:

::

    double density = 1000;
    double av = 20;
    double temperature = 10;

    cell_t cell;
    cell.nh = &density;
    cell.av = &av;
    cell.tgas = &temperature;
    cell.tdust = &temperature;

The fifth step is to initialize the solver, and to set the absolute and
relative error:

::

    double abs_err = ABS_ERR_DEFAULT;
    double rel_err = REL_ERR_DEFAULT;
    astrochem_mem_t astrochem_mem;

    solver_init (&cell, &network, &phys, abundances , density,
                 abs_err, rel_err, &astrochem_mem );

The ``astrochem_mem`` variable is a structure that contains the various
parameters and work arrays of the solver.

The sixth step is to call the solver to *advance time,* i.e. to compute
the abundances up to a given time:

::

    double time = 0;
    time += 1e-6;
    solve (&astrochem_mem, &network, abundances, time, verbose);

The ``solve()`` function updates the abundance array; after a
successfull call, the array contains the abundances of all species at a
given time. The choice of the time step is left to the user. It should
be sufficiently small so that the physical properties of the cell do not
change much, and can approximated as being constant.

This step is repeated an number of times, until the dynamical simulation
finishes. Between two calls, the cell properties needs to be updated
with the new values of the density, temperature, etc. that are computed
in the dynamical code:

::

    time += 1e-6;
    density = 2000;
    temperature = 15;
    solve (&astrochem_mem, &network, abundances, time, verbose);

In this exemple, ``cell.nh``, ``cell.tgas``, ``cell.tdust`` and
``cell.av`` are pointers to ``density``, ``temperature``, etc, so we can
simply update the values of ``density`` and ``temperature``.

The seventh and final step is to free the arrays and structures that are
used by Astrochem:

::

    solver_close (&astrochem_mem);
    free_abundances (abundances);
    free_network (&network);

Calling Astrochem from another code using python API
====================================================

Astrochem can be called from another code through an application
programming interface (API). This allows to compute the chemical
evolution of a gas cell within another code, typically an hydrodynamic
or magneto-hydrodynamic (MHD) simulation. In practice, Astrochem
functions that are part of the API are included in a dynamic library,
against which other codes may be linked. For the moment, only C is
supported.

A typical call of Astrochem library from a python code can be decomposed
in several steps. First, one must import the astrochem python api,
considering it have been installed in a directory contained in
PYTHONPATH

::

    import astrochem

Then, an optionnal step concerning only user who want to compute
abundances with dynamic cell parameters. If so define array containing
your dynamic cell parameter, for example using numpy

::

    import numpy

    density = numpy.logspace(4, 6, 128) # typically this would come from an analytical formulae
    tgas = numpy.linspace(10, 100, 128)
    tdust = numpy.linspace(10, 100, 128)
    times = numpy.zeros(128)
    times[1:] = numpy.logspace(-6, 8, 127) # first time is zero
    av = 0 #av is not dynamic here

Then define all the physical parameters, wich are the same than in input
file. Please refer to section [sec:input-file] of this manual for
details. As for input files, parameters that are not set explicitly are
set to their default value.

::

    p = astrochem.Phys()
    p.cosmic = 1e-17
    p.chi = 0

Then define the initial abundances in a dictionnary. Obviously, the
species in the dictionnary must be present in the network wich will be
used in the solver, see below.

::

    initial_abundances = {"CO": 1e-4, "HCO(+)": 1e-9, "e(-)": 1e-9}

Then set the initial density, visual extinction and temperature of the
gas cell. In the dynamic case:

::

    c = astrochem.Cell( av , density[0] , tgas[0],  tdust[0]  )

Or in static:

::

    density = 1000
    av = 20
    tgas = 10
    tdust = 10
    c = astrochem.Cell( av , density , tgas,  tdust  )

Then initialize the solver, and set the absolute and relative error,
also verbose flag, along with the network file to use and initial
density for the solver:

::

    verbose = 0
    abs_err = astrochem._ABS_ERR_DEFAULT
    rel_err = astrochem._REL_ERR_DEFAULT
    s = astrochem.Solver( c,  "network.chm", p , abs_err, rel_err,
        initial_abundances, density[0], verbose )

Then actual solving can be done by *advancing time,* i.e. we call the
solve function in a loop like this. With dynamic cell parameters, wich
are updated before each call to solve:

::

    for i in range(1, len(times)):
        c = astrochem.Cell( av , density[i] , temperature[i],  temperature[i]  )
        try:
            abundances = s.solve(times[i], c )
        except ArithmeticError as e:
            raise "something went wrong: %s" % e

Or static:

::

    time = 0
    for i in range(1, len(times)):
        try:
            abundances = s.solve(times[i] )
            time+=1e-6
        except ArithmeticError as e:
            raise "something went wrong: %s" % e

The ``solve()`` function updates the internal abundance array then
return it as a dictionnary; after a successfull call, the array contains
the abundances of all species at a given time. The choice of the time
step is left to the user. It should be sufficiently small so that the
physical properties of the cell do not change much, and can approximated
as being constant.

This step is repeated an number of times, until the dynamical simulation
finishes. Between two calls, the cell properties can be updated with the
new values of the density, temperature, etc. as seen in the dynamic
code.

Astrochem authors
=================

The following authors have contributed to Astrochem:

**Sébastien Maret**
    | 
    | wrote Astrochem.

**Ted Bergin**
    | 
    | wrote an original version of this code in Fortran. Although the
    present version of Astrochem has been written from scratch, many
    ideas on its implementation were borrowed from the original Fortran
    version.

Contributing to Astrochem
=========================

Astrochem is hosted by GitHub, a web-based hosting service for software
development that uses the Git revision control system. The Astrochem
repository on GitHub is located at the following address:

http://github.com/smaret/astrochem

From this page, you can download Astrochem source code, either as a tar
file, or using Git. Latest releases source are available from this page
as well. Bugs and issues may also be reported there.

You can help the development of Astrochem in a number of ways: by
testing the package and reporting any bugs that you find; by making
improvements that you develop available to others; by working on the
problems or items listed on GitHub.

.. [1]
   For the moment, only spherical objects with an external radiation
   field are supported.

.. [2]
   http://www.physics.ohio-state.edu/~eric/research.html

.. [3]
   http://kida.obs.u-bordeaux1.fr/

.. [4]
   The Jacobian matrix is a :math:`n_\mathrm{s}^2` elements square
   matrix, whose elements are defined as

   .. math::

      J_{i,j} = \frac{\partial
             \dot{n(\mathrm{x_{i}})}}{\partial {n(\mathrm{x_{i}})}}

   with:

   .. math:: \dot{n(\mathrm{x_{i}})} = \frac{\mathrm{d} {n(\mathrm{x_{i}})}}{\mathrm{d}t}

.. [5]
   In other words:

   .. math::

      n_\mathrm{H} = {n(\mathrm{H})} + 2 \, {n(\mathrm{H_{2}})} + 3 \,
          {n(\mathrm{H_{3}^{+}})} + ...

   Abundances computed by Astrochem are always with respect to
   :math:`n_\mathrm{H}`.

.. [6]
   The formulation for this process, and in particular normalization of
   the rate by the :math:`\frac{n_\mathrm{H}}{n_\mathrm{d}}` may sound a
   bit awkward. In fact, Astrochem treats the electron attachment and
   recombination in this fashion essentially for compatibility with the
   OSU database. This may change in future releases of the code.

.. [7]
   For some species, :math:`Y_{0}` depends weakly on the dust
   temperature . This effect is not implemented in Astrochem.

.. [8]
   Plotting the main formation/destruction routes as a function of
   position is not implemented yet.

.. [9]
   The formation and destruction reactions are ordered by the integral
   of their formation/destruction rate over time. i.e. their total
   contribution to the formation/destruction of the species. The
   ``--percent`` option allows to display the relative contribution,
   expressed in percent, to the formation/destruction of the species.

.. [10]
   Output files contain abundances for ``time_steps`` time values that
   are logarithmically sampled between ``ti`` and ``tf``.

.. [11]
   This dynamical time step is different from the chemical time step.
   The later is adjusted internally by Astrochem to ensure a sufficient
   precision on the computed abundances. The former is provided by the
   user.

.. [12]
   The reaction types adopted in Astrochem follows closely that used in
   the OSU database.

.. [13]
   For more information on OpenMP, see http://openmp.org/

.. [14]
   Most compilers (including GCC, starting from version 4.2) support
   OpenMP.

.. [15]
   For a full description of the datatype and functions defined in the
   API, please see the Doxygen documentation in the ``doc/doxygen``
   directory of Astrochem source code.
