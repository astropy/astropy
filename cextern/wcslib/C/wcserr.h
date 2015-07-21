/*============================================================================

  WCSLIB 5.9 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2015, Mark Calabretta

  This file is part of WCSLIB.

  WCSLIB is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.

  WCSLIB is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with WCSLIB.  If not, see http://www.gnu.org/licenses.

  Direct correspondence concerning WCSLIB to mark@calabretta.id.au

  Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
  Module author: Michael Droettboom
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: wcserr.h,v 5.9 2015/07/21 09:20:01 mcalabre Exp $
*=============================================================================
*
* WCSLIB 5.9 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
* Summary of the wcserr routines
* ------------------------------
* Most of the structs in WCSLIB contain a pointer to a wcserr struct as a
* member.  Functions in WCSLIB that return an error status code can also
* allocate and set a detailed error message in this struct which also
* identifies the function, source file, and line number where the error
* occurred.
*
* For example:
*
=     struct prjprm prj;
=     wcserr_enable(1);
=     if (prjini(&prj)) {
=       // Print the error message to stderr.
=       wcsprintf_set(stderr);
=       wcserr_prt(prj.err, 0x0);
=     }
*
* A number of utility functions used in managing the wcserr struct are for
* internal use only.  They are documented here solely as an aid to
* understanding the code.  They are not intended for external use - the API
* may change without notice!
*
*
* wcserr struct - Error message handling
* --------------------------------------
* The wcserr struct contains the numeric error code, a textual description of
* the error, and information about the function, source file, and line number
* where the error was generated.
*
*   int status
*     Numeric status code associated with the error, the meaning of which
*     depends on the function that generated it.  See the documentation for
*     the particular function.
*
*   int line_no
*     Line number where the error occurred as given by the __LINE__
*     preprocessor macro.
*
*   const char *function
*     Name of the function where the error occurred.
*
*   const char *file
*     Name of the source file where the error occurred as given by the
*     __FILE__ preprocessor macro.
*
*   char msg[WCSERR_MSG_LENGTH]
*     Informative error message.
*
*
* wcserr_enable() - Enable/disable error messaging
* ------------------------------------------------
* wcserr_enable() enables or disables wcserr error messaging.  By default it
* is disabled.
*
* PLEASE NOTE: This function is not thread-safe.
*
* Given:
*   enable    int       If true (non-zero), enable error messaging, else
*                       disable it.
*
* Function return value:
*             int       Status return value:
*                         0: Error messaging is disabled.
*                         1: Error messaging is enabled.
*
*
* wcserr_prt() - Print a wcserr struct
* ------------------------------------
* wcserr_prt() prints the error message (if any) contained in a wcserr struct.
* It uses the wcsprintf() functions.
*
* Given:
*   err       const struct wcserr*
*                       The error object.  If NULL, nothing is printed.
*
*   prefix    const char *
*                       If non-NULL, each output line will be prefixed with
*                       this string.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         2: Error messaging is not enabled.
*
*
* wcserr_clear() - Clear a wcserr struct
* --------------------------------------
* wcserr_clear() clears the error (if any) contained in a wcserr struct.
*
* Given and returned:
*   err       struct wcserr**
*                       The error object.  If NULL, nothing is done.  Set to
*                       NULL on return.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*
*
* wcserr_set() - Fill in the contents of an error object
* ------------------------------------------------------
* INTERNAL USE ONLY.
*
* wcserr_set() fills a wcserr struct with information about an error.
*
* A convenience macro, WCSERR_SET, provides the source file and line number
* information automatically.
*
* Given and returned:
*   err       struct wcserr**
*                       Error object.
*
*                       If err is NULL, returns the status code given without
*                       setting an error message.
*
*                       If *err is NULL, allocates memory for a wcserr struct
*                       (provided that status is non-zero).
*
* Given:
*   status    int       Numeric status code to set.  If 0, then *err will be
*                       deleted and *err will be returned as NULL.
*
*   function  const char *
*                       Name of the function generating the error.  This
*                       must point to a constant string, i.e. in the
*                       initialized read-only data section ("data") of the
*                       executable.
*
*   file      const char *
*                       Name of the source file generating the error.  This
*                       must point to a constant string, i.e. in the
*                       initialized read-only data section ("data") of the
*                       executable such as given by the __FILE__ preprocessor
*                       macro.
*
*   line_no   int       Line number in the source file generating the error
*                       such as given by the __LINE__ preprocessor macro.
*
*   format    const char *
*                       Format string of the error message.  May contain
*                       printf-style %-formatting codes.
*
*   ...       mixed     The remaining variable arguments are applied (like
*                       printf) to the format string to generate the error
*                       message.
*
* Function return value:
*             int       The status return code passed in.
*
*
* wcserr_copy() - Copy an error object
* ------------------------------------
* INTERNAL USE ONLY.
*
* wcserr_copy() copies one error object to another.  Use of this function
* should be avoided in general since the function, source file, and line
* number information copied to the destination may lose its context.
*
* Given:
*   src       const struct wcserr*
*                       Source error object.  If src is NULL, dst is cleared.
*
* Returned:
*   dst       struct wcserr*
*                       Destination error object.  If NULL, no copy is made.
*
* Function return value:
*             int       Numeric status code of the source error object.
*
*
* WCSERR_SET() macro - Fill in the contents of an error object
* ------------------------------------------------------------
* INTERNAL USE ONLY.
*
* WCSERR_SET() is a preprocessor macro that helps to fill in the argument list
* of wcserr_set().  It takes status as an argument of its own and provides the
* name of the source file and the line number at the point where invoked.  It
* assumes that the err and function arguments of wcserr_set() will be provided
* by variables of the same names.
*
*===========================================================================*/

#ifndef WCSLIB_WCSERR
#define WCSLIB_WCSERR

#ifdef __cplusplus
extern "C" {
#endif

#define WCSERR_MSG_LENGTH 160

struct wcserr {
  int  status;			/* Status code for the error.               */
  int  line_no;			/* Line number where the error occurred.    */
  const char *function;		/* Function name.                           */
  const char *file;		/* Source file name.                        */
  char msg[WCSERR_MSG_LENGTH];	/* Informative error message.               */
};

/* Size of the wcserr struct in int units, used by the Fortran wrappers. */
#define ERRLEN (sizeof(struct wcserr)/sizeof(int))

int wcserr_enable(int enable);

int wcserr_prt(const struct wcserr *err, const char *prefix);

int wcserr_clear(struct wcserr **err);


/* INTERNAL USE ONLY -------------------------------------------------------*/

int wcserr_set(struct wcserr **err, int status, const char *function,
  const char *file, int line_no, const char *format, ...);

int wcserr_copy(const struct wcserr *src, struct wcserr *dst);

/* Convenience macro for invoking wcserr_set(). */
#define WCSERR_SET(status) err, status, function, __FILE__, __LINE__

#ifdef __cplusplus
}
#endif

#endif /* WSCLIB_WCSERR */
