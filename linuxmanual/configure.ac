#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([liblip],[1.2],[gleb@deakin.edu.au esteban@v7w.com])
AM_INIT_AUTOMAKE([liblip],[1.2],[])

AC_CONFIG_SRCDIR(config.h.in)
AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([src/forest.h])

AM_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

AC_LANG_CPLUSPLUS

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [cos])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([src/malloc.h src/memory.h src/stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor memset sqrt])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
