#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([dana],
    [m4_esyscmd_s([git describe --always --abbrev=4 --dirty])],
    [pauvsaravia@unc.edu.ar])

AC_CONFIG_SRCDIR(src/)

AM_INIT_AUTOMAKE([subdir-objects foreign])

# Avoid make to rebuild configure
AM_MAINTAINER_MODE

# Set OpenMP flags
# : ${enable_openmp=no}
AC_OPENMP

# Set AR
AM_PROG_AR

# Set FCFLAGS and FC
: ${FCFLAGS=""} # Override default (i.e. "-g -O2")
AC_PROG_FC

AC_ARG_ENABLE([coarray],
  AS_HELP_STRING([--enable-coarray@<:@=single|lib|compiler@:>@],
                 [enable Fortran coarray runtime support; yes is equivalent to lib]),
  [enable_coarray="${enableval}"],
  [enable_coarray=no])

AC_ARG_WITH([coarray-lib],
  AS_HELP_STRING([--with-coarray-lib=LIBS],
                 [libraries needed by the coarray runtime, e.g. -lcaf_mpi, -lcaf_shmem, or -lcaf_single]),
  [COARRAY_LIBS="${withval}"],
  [COARRAY_LIBS=""])

COARRAY_FCFLAGS=
AS_CASE([$enable_coarray],
  [no], [],
  [yes], [enable_coarray=lib],
  [single], [COARRAY_FCFLAGS="-fcoarray=single"],
  [lib], [COARRAY_FCFLAGS="-fcoarray=lib"],
  [compiler], [],
  [AC_MSG_ERROR([bad value "${enable_coarray}" for --enable-coarray; use no, single, lib, or compiler])])

AS_IF([test "x$enable_coarray" != "xno"], [
  AC_DEFINE([DANA_COARRAY], [1], [Define to enable Fortran coarray support])
])

AC_SUBST([COARRAY_FCFLAGS])
AC_SUBST([COARRAY_LIBS])
 
# External software
# =================

# subtrees
# --------	
AC_CONFIG_SUBDIRS([
    lib/fpt
])
   
# other options
# =============

AM_CONDITIONAL([COMPILER_IS_IFORT], [test x"$FC" = xifort])
AM_CONDITIONAL([COMPILER_IS_GFORTRAN], [test x"$FC" = xgfortran])
    
# Debug
# -----
AC_ARG_ENABLE([debug],
  AS_HELP_STRING([--enable-debug],[Turn on debugging]),
  [case "${enableval}" in
    yes) debug=true ;;
    no)  debug=false ;;
    *) AC_MSG_ERROR([bad value "${enableval}" for --enable-debug]) ;;
  esac],enable_debug=no)
AM_CONDITIONAL([DEBUG], [test "x$enable_debug" = "xyes"])


AC_CONFIG_FILES([
 Makefile
 src/Makefile
])


# A precompiler variable
# AC_DEFINE([MYDEFINE],[SOMEVALUE],[Description: Make #define MYDEFINE SOMEVALUE])

AC_OUTPUT  
