#  configure.ac - Autoconf file for Astrochem
#
#  Copyright (c) 2006-2019 Sebastien Maret
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

AC_PREREQ(2.59)
AC_INIT(astrochem, 0.8, http://github.com/smaret/astrochem/issues)
RELEASE_DATE="November 27, 2019"
AC_SUBST([RELEASE_DATE])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Set package options
AC_ARG_ENABLE([openmp],	AC_HELP_STRING([--enable-openmp], [enable openmp (default is no)]))
AC_ARG_ENABLE([lapack],	AC_HELP_STRING([--enable-lapack], [enable LAPACK (default is no)]))

# Check for programs
AC_PROG_CC
LT_INIT
AM_PATH_PYTHON([3.4],, AC_MSG_ERROR(Could not find Python >=3.4. Please set PYTHON to Python interpreter path.))
AX_PYTHON_DEVEL([])

# Checks for libraries.
AC_CHECK_LIB(m, main)
AC_CHECK_LIB(sundials_cvode, main,, AC_MSG_ERROR(Could not find SUNDIALS/CVODE library. Make sure you set CFLAGS to include the include directory and set LDFLAGS to include the library directory.))
AC_CHECK_LIB(sundials_cvode, CVodeSetLinearSolver,, AC_MSG_ERROR(Could not find the CVodeSetLinearSolver function in SUNDIALS/CVODE  library. Make sure that your SUNDIALS version is >= 4.1.0.))
AC_CHECK_LIB(sundials_nvecserial, main,, AC_MSG_ERROR(Could not find SUNDIALS/NVECSERIAL library. Make sure you set CFLAGS to include the include directory and set LDFLAGS to include the library directory.))
if test x$enable_lapack = "xyes"; then
   AC_CHECK_LIB(sundials_sunlinsollapackdense, main,, AC_MSG_WARN(SUNDIALS was compiled without LAPACK.))
   if test x$ac_cv_lib_sundials_sunlinsollapackdense_main = xyes; then
      use_lapack="yes"
      AC_DEFINE([USE_LAPACK], 1, [Define to use LAPACK.])
   else
      use_lapack="no"
   fi
else
   use_lapack="no"
fi

# Check for Python modules
AC_PYTHON_MODULE(numpy)
if test $HAVE_PYMOD_NUMPY = "no"; then
   AC_MSG_ERROR(Could not find Python module Numpy.)
fi

AC_PYTHON_MODULE(matplotlib, )
if test $HAVE_PYMOD_MATPLOTLIB = "no"; then
   AC_MSG_WARN(Could not find Python module Matplotlib.)
fi
AM_CONDITIONAL([HAVE_PYMOD_MATPLOTLIB], [test $HAVE_PYMOD_MATPLOTLIB = "yes"])

# Checks for typedefs, structures, and compiler characteristics.
if test x$enable_openmp = "xyes"; then
   AX_OPENMP(have_openmp="yes", have_openmp="no")
   AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])
   AC_SUBST([OPENMP_CFLAGS])
else
   have_openmp="no"
fi

#check for cython
if test -z $CYTHON; then
    AC_PATH_PROG( CYTHON, cython )
else
    AC_MSG_CHECKING([Checking Cython path $CYTHON is correct])
    AS_IF($CYTHON -V > /dev/null 2>&1, , CYTHON="")
    if test -z $CYTHON; then
        AC_MSG_RESULT([ no ])
    else
        AC_MSG_RESULT([ yes ])
    fi
fi

HAVE_CYTHON="yes"
if test -z $CYTHON; then
HAVE_CYTHON="no"
AC_MSG_WARN([[Cannot find cython, which is required to build libpyastrochem bindings.  They will not be built.]])
fi
AM_CONDITIONAL([HAVE_CYTHON], [test  "x$CYTHON" != "x" ])

#check for hdf5
AX_LIB_HDF5()
if test "$with_hdf5" = "no"; then
    AC_MSG_ERROR([Unable to find HDF5, cannot compile astrochem.])
fi

#check for h5py
AC_PYTHON_MODULE(h5py, )
if test $HAVE_PYMOD_H5PY = "no"; then
   AC_MSG_WARN(Cannot find python module h5py which is required to use python astrochem tools output file readers. They will not be usable.)
fi

AC_CONFIG_FILES([Makefile
                 doc/Makefile
		 doc/examples/Makefile
		 doc/examples/input.ini
		 doc/examples/dynamic/Makefile
		 doc/examples/dynamic/input.ini
		 doc/man/Makefile
		 doc/man/astrochem.1
		 doc/man/plabun.1
		 doc/man/plroute.1
		 doc/man/chmconvert.1
		 doc/man/converttolegacy.1
         	 src/Makefile
		 src/chmconvert.py
		 src/plabun.py
		 src/plroute.py
		 src/converttolegacy.py
		 networks/Makefile
         	 python/Makefile
		 python/libpyastrochem/Makefile
         	 tests/astrochem_test.py
         	 tests/python_lib_test.py
         	 tests/python_tools_test.py
         	 tests/benchmark.py
		 tests/benchmark_openmp.py
		 tests/Makefile])

AC_OUTPUT

# Print configuration summary

echo "=============================================
Configuration summary:
C compiler:         ${CC} ${CFLAGS}
OpenMP support:     ${have_openmp}"
if test ${have_openmp} = "yes"; then
   echo "OpenMP flag:        ${OPENMP_CFLAGS}"
fi
echo "LAPACK support:     ${use_lapack}"
echo "Python interpreter: ${PYTHON}
Python version:     ${PYTHON_VERSION}
Matplotlib module:  ${HAVE_PYMOD_MATPLOTLIB}
H5py module:        ${HAVE_PYMOD_H5PY}
Python API:         ${HAVE_CYTHON}
============================================="
