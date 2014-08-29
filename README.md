Astrochem
=========

What is Astrochem?
------------------

Astrochem is a code to compute the abundances of chemical species in
the interstellar medium, as function of time. It is designed to study
the chemistry in a variety of astronomical objects, including diffuse
clouds, dense clouds, photodissociation regions, prestellar cores,
protostars, and protostellar disks. Astrochem reads a network of
chemical reactions from a text file, builds up a system of kinetic
rates equations, and solve it using a state-of-the-art stiff ordinary
differential equation (ODE) solver. The Jacobian matrix of the system
is computed implicitly, so the resolution of the system is extremely
fast: large networks containing several thousands of reactions are
usually solved in a few seconds. A variety of gas phase process are
considered, as well as simple gas-grain interactions, such as the
freeze-out and the desorption via several mechanisms (thermal
desorption, cosmic-ray desorption and photo-desorption). The computed
abundances are written in text files, and can be plotted in different
ways with the tools provided with Astrochem. Chemical reactions and
their rates are written in a format which is meant to be easy to read
and to edit. A tool to convert the chemical networks from the OSU and
KIDA databases into this format is also provided. Astrochem is written
in C, and its source code is distributed under the terms of the [GNU
General Public License (GPL)](COPYING.md).

Installation
------------

Astrochem installation follows a standard procedure and is fairly
straightforward. Instructions can be found in the
[INSTALL.md](./INSTALL.md) file.

Documentation
-------------

Astrochem comes with a comprehensive documentation manual. The LaTeX
source of the manual can be found in the [doc/](./doc/)
directory. Each command also has a man page.

Reporting bugs
--------------

A list of known bugs can be found
[here.](http://github.com/smaret/astrochem/issues?labels=Bug) If you
find a bug which is not listed there, please open a new issue.

How to contribute
-----------------

Contributions to Astrochem are welcome. Please see the
[CONTRIBUTING.md](./CONTRIBUTING.md/) file for instructions.

More information
----------------

For more information on Astrochem and its development, please see the
[project page on GitHub.](http://github.com/smaret/astrochem)

[![Build Status](https://travis-ci.org/smaret/astrochem.svg)]
(https://travis-ci.org/smaret/astrochem)
