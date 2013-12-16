#  configure.ac - Autoconf file for Astrochem
#
#  Copyright (c) 2006-2013 Sebastien Maret
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
AC_INIT(astrochem, 0.6-dev, http://github.com/smaret/astrochem/issues)
RELEASE_DATE="Not released yet"
AC_SUBST([RELEASE_DATE])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADER([config.h])

# Set package options
AC_ARG_ENABLE([openmp],	AC_HELP_STRING([--enable-openmp], [enable openmp (default is no)]))

# Check for programs
AC_PROG_CC
AC_PROG_F77
AC_F77_WRAPPERS
AM_PATH_PYTHON([2.4],, AC_MSG_ERROR(Could not find Python >=2.4. Please set PYTHON to Python interpreter path.))
AC_CHECK_PROGS(PDFLATEX, pdflatex, no)
AC_CHECK_PROGS(BIBTEX, bibtex, no)
AM_CONDITIONAL(HAVE_PDFLATEX, test $PDFLATEX != no)
AM_CONDITIONAL(HAVE_BIBTEX, test $BIBTEX != no)

# Checks for libraries.
AC_CHECK_LIB(m, main)
AX_LAPACK(, AC_MSG_WARN(Could not find LAPACK library. Make sure you set CFLAGS to include the include directory and set LDFLAGS to include the library directory.))
if test x$ax_lapack_ok = xyes; then
   LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
fi
AC_CHECK_LIB(sundials_cvode, main,, AC_MSG_ERROR(Could not find SUNDIALS/CVODE library. Make sure you set CFLAGS to include the include directory and set LDFLAGS to include the library directory.))
AC_CHECK_LIB(sundials_nvecserial, main,, AC_MSG_ERROR(Could not find SUNDIALS/NVECSERIAL library. Make sure you set CFLAGS to include the include directory and set LDFLAGS to include the library directory.))
AC_CHECK_LIB(sundials_cvode, CVLapackDense,, AC_MSG_WARN(SUNDIALS was compiled without LAPACK.))
if test x$ac_cv_lib_sundials_cvode_CVLapackDense = xyes && test x$ax_lapack_ok = xyes; then
   use_lapack="yes"
   AC_DEFINE([USE_LAPACK], 1, [Define to use LAPACK.])
else
   use_lapack="no"
fi  

# Check for Python modules
AC_PYTHON_MODULE(numpy)
if test $HAVE_PYMOD_NUMPY = "no"; then
   AC_MSG_ERROR(Could not find Python module Numpy.)
fi
AC_PYTHON_MODULE(biggles, )
if test $HAVE_PYMOD_BIGGLES = "no"; then
   AC_MSG_WARN(Could not find Python module Biggles.)
fi
AM_CONDITIONAL([HAVE_PYMOD_BIGGLES], [test $HAVE_PYMOD_BIGGLES = "yes"])

# Checks for typedefs, structures, and compiler characteristics.
if test x$enable_openmp = "xyes"; then
   AX_OPENMP(have_openmp="yes", have_openmp="no")
   AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])
   AC_SUBST([OPENMP_CFLAGS])
else
   have_openmp="no"
fi

AC_CONFIG_FILES([Makefile
                 doc/Makefile
		 doc/astrochem.tex
		 doc/astrochem.1
		 doc/plabun.1
		 doc/plroute.1
		 doc/chmconvert.1
		 doc/examples/Makefile
		 doc/examples/dynamic/Makefile
		 lib/Makefile
		 lib/libastrochem.py
                 src/Makefile
		 src/chmconvert.py
		 src/plabun.py
		 src/plroute.py
		 networks/Makefile
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
if test ${use_lapack} = "yes"; then
   echo "LAPACK libraries:  ${LAPACK_LIBS} ${BLAS_LIBS}"
fi
echo "Python interpreter: ${PYTHON}
Python version:     ${PYTHON_VERSION}
Biggles module:     ${HAVE_PYMOD_BIGGLES}
============================================="
