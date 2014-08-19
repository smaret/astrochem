What's new in Astrochem ?
=========================

astrochem-0.6 (Not released yet)
--------------------------------

* Numpy module is now mandatory (fixes issue #3)

* Suppress autoconf and automake warnings (fixes issue #5)

* Code parameters have been regrouped into data structures, for a sake
  of clarity. These data structures are dynamically allocated, so
  there is no limit on the numbers of shells or reactions anymore.

astrochem-0.5 (October 15, 2012)
--------------------------------

* Fixed a bug in plroute that occured when some of the formation or
  destruction routes were equal to 0.

* Fixed a segmentation fault that ocurred when the number of species
  in output was greater than the maximum.

* Dropped the `time` and `av` commands from `plroute` and `plabun`; both
  commands now plot abundances and formation/destruction routes as a
  function of time.

* Added some information on the OpenMP and LAPACK support in
  `astrochem --version`.

* Fixed a bug in the detection of the LAPACK library. Astrochem now
  uses LAPACK by default.

* Some of the Python code has been factorized into a Python module
  (libastrochem.py). This module defines a network class that allows
  to work on chemical networks in different formats (.osu, .kida, and
  .chm). The module also allows to read the abundances and the
  formation/destruction routes computed by Astrochem.

* chmconvert can now convert networks from the KInetic Database for
  Astrochemistry (KIDA) into Astrochem format (.chm).

astrochem-0.4 (June 3, 2011)
----------------------------

* The main destruction/formation channels of a species plotted by
  plroute are now ordered by their total contribution to the
  formation/destruction of that species over the time of the
  computation. An aption has been added to print the relative
  contribution of each channel.

* Fixed a bug in the computation of the main destruction routes.

astrochem-0.3 (February 28, 2011)
---------------------------------

* The rate for cosmic-ray desorption is now either 1) computed
  following Hasegawa & Herbst (1993) or 2) given explicitly in the
  network file (e.g. using values from Bringa & Johnson 2004).

* If the LAPACK library is installed, Astrochem can now use the
  CVLapackDense solver from SUNDIALS instead of CVDense. This results
  in a performance increase of 30% or more, depending on the LAPACK
  version.

astrochem-0.2 (January 21, 2011)
--------------------------------

* First public release
