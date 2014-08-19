Installation
============

Astrochem is written in ANSI C and should compile on most UNIX
platforms including GNU/Linux and Mac OSX.

Pre-requisites
--------------

To compile Astrochem you will need:

1. An ANSI C-compiler (e.g. gcc, icc, etc.)

2. The [SUNDIALS (SUite of Nonlinear and DIfferential/ALgebraic
   equation Solvers) library](http://computation.llnl.gov/casc/sundials/).
   You will need at least version 2.4.0 of the library installed on
   your computer.

3. Python (version 2.4 or later, but earlier than version 3.0).

4. The [NumPy](http://numpy.scipy.org/) and
   [Biggles](http://biggles.sourceforge.net/) Python modules.  Numpy
   is mandatory, but Biggles is optional. However, the later is
   required to use the plotting tools provided with Astrochem. Both
   modules are present in major Linux distributions (Debian, Ubuntu,
   etc.)  and in Fink or MacPorts on Mac OSX.

5. To compile the documentation, you will also need PDFLaTeX and
   BibTeX.
 
Basic Installation
------------------

Astrochem follows the standard GNU installation procedure. A configure
script attempts to find on your system the C compiler, the libraries
and the programs needed to compile Astrochem and creates the
Makefiles. This script is invoked with:

```
./configure
```

By default, the script will look for the gcc compiler. If you want to
use another compiler, you may do so by setting the CC variable, e.g.:

```
./configure CC=icc
```

You may also set some compiler specific optimization options using the
CFLAGS variable:

```
./configure CC=icc CFLAGS=-O3
```

The script will look for the SUNDIALS library in standard directories
(/usr, /usr/local, etc.). If you have it installed in some other
place, e.g. /sw, so you can set the search path as follows:

```
./configure CPPFLAGS=-I/sw/include LDFLAGS=-L/sw/lib
```

You can then build Astrochem by typing:

```
make
```

A test suite is available.  After compiling the library with `make`,
it can be invoked with `make check` at the top level. If you run the
tests and get some failures, please report
[here](http://github.com/smaret/astrochem/issues?labels=Bug).

The documentation manual can built by typing `make pdf` at the top
level.

Astrochem can be installed using the command:

```
make install
```

The default installation directory prefix is /usr/local.  Installing
in this directory will require root privileges on most systems (use
`su` or `sudo`). The installation directory can be changed with the
`--prefix` option to configure.

Optional features
-----------------

Astrochem can use the LAPACK library together with SUNDIALS, which
results in better performance than when using SUNDIALS alone. The use
of LAPACK in Astrochem can be turned on by specifying the
`--with-lapack` option to `configure`, e.g.

```
./configure --with-lapack
```

In addition, Astrochem can be run in parallel on multi CPU (or core)
computers. For this it uses OpenMP compilation directives. If your
compiler supports OpenMP, compilation of the parallel version of
Astrochem can be turn on by specifying the `--enable-openmp` option to
`configure`, e.g.

```
./configure --enable-openmp
```

The configure script accepts some other options.  Run `configure
--help` for more details.
