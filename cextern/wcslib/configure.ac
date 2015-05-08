#-----------------------------------------------------------------------------
# Process this file with autoconf-2.53 or later to produce a configure script.
#-----------------------------------------------------------------------------
# Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
# http://www.atnf.csiro.au/people/Mark.Calabretta
# $Id: configure.ac,v 5.5 2015/05/05 13:16:33 mcalabre Exp $
#-----------------------------------------------------------------------------

AC_INIT([WCSLIB], [5.5], [mark@calabretta.id.au], [wcslib-5.5])
AC_PREREQ([2.53])
AC_REVISION([$Revision: 5.5 $])
AC_SUBST([PACKAGE_VERSION])
AC_DEFINE_UNQUOTED([WCSLIB_VERSION], [$PACKAGE_VERSION], [Define wcslib version])

# Library version number, same as package version.
LIBVER="$PACKAGE_VERSION"
AC_SUBST([LIBVER])

AC_CONFIG_SRCDIR([C/wcs.h])
AC_CONFIG_AUX_DIR([config])

# Get the system type.
AC_CANONICAL_BUILD
ARCH="${build_cpu}-$build_os"
AC_SUBST([ARCH])


# Look for Flex.
AC_CHECK_PROG([FLEX], [flex], [flex], [], [], [])
if test "x$FLEX" = xflex ; then
  # Version 2.5.9 or later is required.
  V=`flex --version | awk '{print $2}'`
  W=`echo $V | awk -F. '{if ((($1*100 + $2)*100 + $3) < 20509) print "no"}'`
  if test "x$W" != x ; then
    AC_MSG_NOTICE([Flex version $V is too old, ignored.])
    FLEX=
  else
    AC_MSG_NOTICE([Using Flex version $V.])
  fi
fi

if test "x$FLEX" = x ; then
  AC_MSG_WARN([Flex version 2.5.9 or later does not appear to be
           available, will use pre-generated sources.])
fi


# Look for an ANSI C compiler.
AC_PROG_CPP
AC_PROG_CC
AC_PROG_CC_STDC
AC_C_CONST
AC_TYPE_SIZE_T
if test "x$ac_cv_prog_cc_stdc" = xno -o \
        "x$ac_cv_c_const"      = xno -o \
        "x$ac_cv_type_size_t"  = xno; then
  AC_MSG_ERROR([
    -------------------------------------------------------
    An ANSI standard C library is required to build WCSLIB.

    ERROR: WCSLIB configuration failure.
    -------------------------------------------------------], [1])
fi

# Check for data types (suggested by autoscan - off_t is only required by
# fitshdr).
AC_TYPE_OFF_T
AC_TYPE_INT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T

# Check for ANSI C headers.
AC_HEADER_STDC
AC_CHECK_HEADERS([ctype.h errno.h limits.h locale.h math.h setjmp.h stdarg.h \
                  stdio.h stdlib.h string.h])
if test "x$ac_cv_header_stdc" = xno; then
  AC_MSG_ERROR([
    -------------------------------------------------------------------
    An ANSI standard C library is required to build WCSLIB.  One of the
    ANSI C header files it requires is missing or unusable.

    ERROR: WCSLIB configuration failure.
    -------------------------------------------------------------------], [1])
fi

# Checks for ANSI C library functions (suggested by autoscan - fseeko and
# stat are only required by fitshdr).
AC_CHECK_LIB([m], [floor])
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_STAT
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor memset pow setlocale sqrt strchr strstr])

# System libraries that may be required by WCSLIB itself.
# SunOS, extra maths functions.
AC_CHECK_LIB([sunmath], [cosd], [LIBS="-lsunmath $LIBS"], [], [])

# See if we can find sincos().
AC_CHECK_FUNCS([sincos])

# Check the size and availability of integer data types.
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long int])
AC_CHECK_SIZEOF([long long int])

# 64-bit integer data type; use long long int preferentially since that
# accords with "%lld" formatting used in fitshdr.l, e.g.
#                int   size_t  long int  long long int
#                ---   ------  --------  -------------
#   gcc x86:      32     32       32          64
#   gcc x86_64:   32     64       64          64
if test "x$ac_cv_sizeof_long_long_int" = x8; then
  AC_DEFINE([WCSLIB_INT64], [long long int], [64-bit integer data type.])
elif test "x$ac_cv_sizeof_long_int" = x8; then
  AC_DEFINE([WCSLIB_INT64], [long int], [64-bit integer data type.])
elif test "x$ac_cv_sizeof_int" = x8; then
  AC_DEFINE([WCSLIB_INT64], [int], [64-bit integer data type.])
fi

# Does printf() have the z modifier for size_t type?  Important for 64-bit.
AC_MSG_CHECKING([for printf z format modifier for size_t type])
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([AC_INCLUDES_DEFAULT],
                   [[char buf[64];
                     if (sprintf(buf, "%zu", (size_t)1) != 1)
                       return 1;
                     else if (strcmp(buf, "1"))
                       return 2;]])],
  AC_DEFINE([MODZ], ["z"], [printf format modifier for size_t type.])
    AC_MSG_RESULT(yes),
  AC_DEFINE([MODZ], [""],  [printf format modifier for size_t type.])
    AC_MSG_RESULT(no),
  AC_DEFINE([MODZ], [""],  [printf format modifier for size_t type.])
    AC_MSG_RESULT(assumed not)
)


# Starting values, may be augmented later.
SUBDIRS="C"
TSTDIRS="C"
INSTDIR="C"


# Ways of specifying the Fortran compiler, in order of precedence:
#   configure --enable-fortran=<compiler>
#   F77=<compiler> configure    ...bash
#
# Ways of disabling Fortran:
#   configure --disable-fortran
#   configure --enable-fortran=no
#   F77=no configure            ...bash
AC_ARG_ENABLE([fortran], [AS_HELP_STRING([--enable-fortran=ARG],
            [Fortran compiler to use])], [])
AC_ARG_ENABLE([fortran], [AS_HELP_STRING([--disable-fortran],
            [don't build the Fortran wrappers or PGSBOX])], [])
if test "x$enable_fortran" != x -a "x$enable_fortran" != xyes ; then
  F77="$enable_fortran"
fi

if test "x$F77" = xno ; then
  F77=

  AC_MSG_WARN([Compilation of Fortran wrappers and PGSBOX disabled])

else
  if test "x$F77" = x ; then
    # Look for a Fortran compiler.
    AC_PROG_F77([gfortran g77 f77 ifort xlf frt pgf77 fl32 af77 fort77 f90 \
                 xlf90 pgf90 epcf90 f95 fort xlf95 lf95 g95])
  fi

  if test "x$F77" = x; then
    AC_MSG_WARN([
      ------------------------------------------------------------------
      Fortran compiler not found, will skip Fortran wrappers and PGSBOX.
      ------------------------------------------------------------------])

    # Best guess at Fortran name mangling for use if a compiler does ever
    # become available.
    AC_DEFINE([F77_FUNC(name,NAME)], [name ## _])

  else
    if test "x$ac_cv_f77_compiler_gnu" = xyes ; then
      if test "x$F77" = xg77 -o "x$F77" = xf77 ; then
        # Not recognized by gfortran.
        FFLAGS="$FFLAGS -Wno-globals"
      fi
    fi

    AC_MSG_CHECKING(whether $F77 accepts -I)
    AC_LANG_PUSH(Fortran 77)
    FFLAGS_save=$FFLAGS
    FFLAGS=-I.
    AC_COMPILE_IFELSE(AC_LANG_PROGRAM([], []),
      [FFLAGS="$FFLAGS_save -I."; AC_MSG_RESULT(yes)],
      [FFLAGS="$FFLAGS_save"; AC_MSG_RESULT(no)])
    AC_LANG_POP()

    # Libraries required by the Fortran compiler itself (sets FLIBS).
    # Required by utilities and test programs written in C that link to
    # Fortran object modules such as pgsbox.
    AC_F77_LIBRARY_LDFLAGS

    # F77 name mangling (defines the F77_FUNC preprocessor macro).
    AC_F77_WRAPPERS

    SUBDIRS="C Fortran"
    TSTDIRS="C Fortran"
    INSTDIR="Fortran"
  fi
fi


# System-dependent system libraries (for building the sharable library).
#-----------------------------------------------------------------------
# Darwin (contains stubs for long double).
AC_CHECK_LIB([SystemStubs], [printf\$LDBLStub], [LIBS="$LIBS -lSystemStubs"],
             [], [])


# Library and installation utilities.
#------------------------------------
# Static library generation.
AC_PROG_RANLIB

# Shared library generation.
if test "x$ac_cv_c_compiler_gnu" = xyes ; then
  SHVER=`echo "$LIBVER" | sed -e 's/\..*$//'`

  # Note that -fPIC is on by default for Macs, this just makes it obvious.
  SHRFLAGS="-fPIC"
  SHRLD="\$(CC) \$(SHRFLAGS)"

  case "$build_os" in
  darwin*)
    SHRLIB="libwcs.$LIBVER.dylib"
    SONAME="libwcs.$SHVER.dylib"
    SHRLD="$SHRLD -dynamiclib -single_module"
    SHRLD="$SHRLD -compatibility_version $SHVER -current_version $LIBVER"
    SHRLN=

    case "$build_cpu" in
    powerpc*)
      # Switch off -fPIC (not applicable for PowerPC Macs).
      CFLAGS="$CFLAGS -mdynamic-no-pic"
      ;;
    esac
    ;;
  *)
    # Covers Linux and Solaris at least.
    SHRLIB="libwcs.so.$LIBVER"
    SONAME="libwcs.so.$SHVER"
    SHRLD="$SHRLD -shared -Wl,-h\$(SONAME)"
    SHRLN="libwcs.so"
    ;;
  esac

else
  SHRLIB=
  SONAME=
  SHRFLAGS=
  SHRLD=
  SHRSFX=
  SHRLN=
fi

AC_SUBST([SHRLIB])
AC_SUBST([SONAME])
AC_SUBST([SHRFLAGS])
AC_SUBST([SHRLD])
AC_SUBST([SHRSFX])
AC_SUBST([SHRLN])

# Installation utilities.
AC_PROG_LN_S
AC_PROG_INSTALL

AC_MSG_NOTICE([End of primary configuration.
])


# The following are required to build utilities and test programs.
# ----------------------------------------------------------------
AC_MSG_NOTICE([Looking for libraries etc. for utilities and test suite...])

# Check for other quasi-standard header files.
AC_CHECK_HEADERS([unistd.h])

# Large file support.
AC_FUNC_FSEEKO
AC_SYS_LARGEFILE


# Extra places to look for third-party libraries and header files.
LIBDIRS=

AC_ARG_WITH([cfitsio], [AS_HELP_STRING([--without-cfitsio],
            [eschew CFITSIO])], [])
if test "x$with_cfitsio" = xno ; then
  AC_MSG_WARN([CFITSIO disabled])
else
  AC_ARG_WITH([cfitsiolib], [AS_HELP_STRING([--with-cfitsiolib=DIR],
              [directory containing cfitsio library])], [])
  if test "x$with_cfitsiolib" != x ; then
    LIBDIRS="$LIBDIRS $with_cfitsiolib"
  fi

  AC_ARG_WITH([cfitsioinc], [AS_HELP_STRING([--with-cfitsioinc=DIR],
              [directory containing cfitsio header files])], [])
  if test "x$with_cfitsioinc" != x ; then
    CFITSIO_INCDIRS="$with_cfitsioinc"
  fi

  CFITSIO_INCDIRS="$CFITSIO_INCDIRS   \
           /usr/local/cfitsio/include \
           /local/cfitsio/include"

  LIBDIRS="$LIBDIRS           \
           /usr/local/cfitsio/lib \
           /local/cfitsio/lib"
fi

AC_ARG_WITH([pgplot], [AS_HELP_STRING([--without-pgplot],
            [eschew PGPLOT])], [])
if test "x$with_pgplot" = xno ; then
  AC_MSG_WARN([PGPLOT disabled])
else
  AC_ARG_WITH([pgplotlib], [AS_HELP_STRING([--with-pgplotlib=DIR],
              [directory containing pgplot library])], [])
  if test "x$with_pgplotlib" != x ; then
    LIBDIRS="$LIBDIRS $with_pgplotlib"
  fi

  AC_ARG_WITH([pgplotinc], [AS_HELP_STRING([--with-pgplotinc=DIR],
              [directory containing pgplot header files])], [])
  if test "x$with_pgplotinc" != x ; then
    PGPLOT_INCDIRS="$with_pgplotinc"
  fi

  PGPLOT_INCDIRS="$PGPLOT_INCDIRS    \
           /usr/local/pgplot/include \
           /local/pgplot/include"

  LIBDIRS="$LIBDIRS           \
           /usr/local/pgplot/lib  \
           /local/pgplot/lib"
fi


if test "x$with_cfitsio" != xno -o \
        "x$with_pgplot"  != xno ; then
  LIBDIRS="$LIBDIRS           \
           /usr/local/lib     \
           /local/lib         \
           /opt/local/lib     \
           /opt/SUNWspro/lib  \
           /sw/lib"

  for LIBDIR in $LIBDIRS ; do
    AC_CHECK_FILE([$LIBDIR], [LDFLAGS="$LDFLAGS -L$LIBDIR"], [continue])
  done

  # Generic include directories.
  INCDIRS="/usr/local/include \
           /local/include     \
           /opt/local/include \
           /sw/include        \
           /local             \
           /usr/include"


  # CFITSIO.
  if test "x$with_cfitsio" != xno ; then
    # Search for CFITSIO.
    for INCDIR in $CFITSIO_INCDIRS $INCDIRS ; do
      AC_CHECK_FILE([$INCDIR/cfitsio/fitsio.h],
                    [CFITSIOINC="-I$INCDIR/cfitsio"; break])
      AC_CHECK_FILE([$INCDIR/fitsio.h], [CFITSIOINC="-I$INCDIR"; break])
    done

    AC_CHECK_LIB([socket],  [recv],   [CFITSIOLIB="-lsocket"], [], [$LIBS])
    AC_CHECK_LIB([cfitsio], [ffopen], [CFITSIOLIB="-lcfitsio $CFITSIOLIB"], [],
                 [$CFITSIOLIB $LIBS])

    if test "x$CFITSIOINC" = x -o "x$CFITSIOLIB" = x; then
      AC_MSG_WARN([CFITSIO not found, skipping CFITSIO-dependent tests.])
    else
      AC_MSG_NOTICE([CFITSIO appears to be available.])
      AC_DEFINE([HAVE_CFITSIO], [1], [Define to 1 if CFITSIO is available.])

      # Check for fits_read_wcstab, present in CFITSIO 3.004beta and later.
      AC_CHECK_LIB([cfitsio], [fits_read_wcstab], [GETWCSTAB=],
                   [GETWCSTAB=getwcstab.o], [$CFITSIOLIB $LIBS])
      if test "x$GETWCSTAB" != x ; then
        AC_MSG_WARN([fits_read_wcstab not found in CFITSIO, will use
                        getwcstab.c to compile test programs.])
      fi
    fi
  fi

  # PGPLOT.
  if test "x$F77" != x -a "x$with_pgplot" != xno ; then
    # Search for PGPLOT.
    for INCDIR in $PGPLOT_INCDIRS $INCDIRS ; do
      AC_CHECK_FILE([$INCDIR/pgplot/cpgplot.h],
                    [PGPLOTINC="-I$INCDIR/pgplot"; break])
      AC_CHECK_FILE([$INCDIR/cpgplot.h], [PGPLOTINC="-I$INCDIR"; break])
    done

    # FLIBS (found above via AC_F77_LIBRARY_LDFLAGS) only helps if PGPLOT was
    # built using the same Fortran compiler that we are using here.

    # PGPLOT compiled by the SUN Fortran compiler but linked with something
    # else.
    AC_CHECK_LIB([M77],     [iand_],     [PGPLOTLIB="-lM77 $PGPLOTLIB"],
                 [], [$PGPLOTLIB $LIBS])
    AC_CHECK_LIB([F77],     [f77_init],  [PGPLOTLIB="-lF77 $PGPLOTLIB"],
                 [], [$PGPLOTLIB $LIBS])

    if test "x$F77" != xg77; then
      # For PGPLOT compiled with g77 but linked with something else.
      AC_CHECK_LIB([frtbegin], [main],     [PGPLOTLIB="-lfrtbegin $PGPLOTLIB"],
                   [], [$PGPLOTLIB $LIBS])
      AC_CHECK_LIB([g2c],      [gerror_],  [PGPLOTLIB="-lg2c $PGPLOTLIB"],
                   [], [$PGPLOTLIB $LIBS])
    fi

    if test "x$F77" != xgfortran; then
      # For PGPLOT compiled with gfortran but linked with something else.
      # Note that if gfortran itself is driving the linker it can be harmful
      # to add -lgfortran to the link list without also adding -lgfortranbegin.
      # Doing so stops gfortran from adding -lgfortranbegin which is needed to
      # resolve "main".
      AC_CHECK_LIB([gfortran], [_gfortran_abort],
                   [PGPLOTLIB="-lgfortran $PGPLOTLIB"], [],
                   [$PGPLOTLIB $LIBS])
    fi

    # Search for X11 includes and libraries.
    AC_PATH_X
    if test "x$no_x" = x; then
      if test "x$ac_x_libraries" != x ; then
        # Not needed for systems that keep the X11 libraries in /usr/lib.
        LDFLAGS="$LDFLAGS -L$ac_x_libraries"
      fi
      PGPLOTLIB="-lX11 $PGPLOTLIB"
    fi

    # It is possible that other libraries may be required depending on what
    # graphics drivers were installed with PGPLOT.
    AC_CHECK_LIB([z],       [deflate],   [PGPLOTLIB="-lz $PGPLOTLIB"],
                 [], [$PGPLOTLIB $LIBS])
    AC_CHECK_LIB([png],     [png_error], [PGPLOTLIB="-lpng $PGPLOTLIB"],
                 [], [$PGPLOTLIB $LIBS])
    AC_CHECK_LIB([pgplot],  [pgbeg_],    [PGPLOTLIB="-lpgplot $PGPLOTLIB"],
                 [], [$PGPLOTLIB $FLIBS $LIBS])
    AC_CHECK_LIB([cpgplot], [cpgbeg],    [PGPLOTLIB="-lcpgplot $PGPLOTLIB"],
                 [PGPLOTLIB=], [$PGPLOTLIB $FLIBS $LIBS])

    # Only need the PGPLOT include file to build PGSBOX.
    if test "x$PGPLOTINC" != x; then
      SUBDIRS="$SUBDIRS pgsbox"
      INSTDIR="pgsbox"
    fi

    # Also need the PGPLOT library to build pgtest and cpgtest.
    if test "x$PGPLOTLIB" = x; then
      AC_MSG_WARN([PGPLOT not found, skipping PGPLOT-dependent tests.])
    else
      AC_MSG_NOTICE([PGPLOT appears to be available.])

      TSTDIRS="$TSTDIRS pgsbox"
    fi
  fi
fi


# Utilities are compiled last since they need the libraries.
# Ways of disabling them:
#   configure --disable-utils
#   configure --enable-utils=no
AC_ARG_ENABLE([utils], [AS_HELP_STRING([--disable-utils],
            [don't build the WCS utilities])], [])
if test "x$enable_utils" != xno ; then
  SUBDIRS="$SUBDIRS utils"
  INSTDIR="$INSTDIR utils"
else
  AC_MSG_WARN([Compilation of WCS utilities disabled])
fi


AC_SUBST([CFITSIOINC])
AC_SUBST([CFITSIOLIB])
AC_SUBST([GETWCSTAB])

AC_SUBST([PGPLOTINC])
AC_SUBST([PGPLOTLIB])

AC_SUBST([SUBDIRS])
AC_SUBST([TSTDIRS])
AC_SUBST([INSTDIR])

AC_MSG_NOTICE([End of auxiliary configuration.
])


# Do it.
AC_MSG_NOTICE([Configuring files...])
AC_CONFIG_FILES([makedefs wcslib.pc])
AC_CONFIG_HEADERS([wcsconfig.h wcsconfig_f77.h wcsconfig_tests.h wcsconfig_utils.h])
AC_OUTPUT
