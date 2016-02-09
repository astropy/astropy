/*============================================================================

  WCSLIB 5.11 - an implementation of the FITS WCS standard.
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
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: wcsprintf.h,v 5.11 2015/10/18 09:13:05 mcalabre Exp $
*=============================================================================
*
* WCSLIB 5.11 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of the wcsprintf routines
* ---------------------------------
* Routines in this suite allow diagnostic output from celprt(), linprt(),
* prjprt(), spcprt(), tabprt(), wcsprt(), and wcserr_prt() to be redirected to
* a file or captured in a string buffer.  Those routines all use wcsprintf()
* for output.  Likewise wcsfprintf() is used by wcsbth() and wcspih().  Both
* functions may be used by application programmers to have other output go to
* the same place.
*
*
* wcsprintf() - Print function used by WCSLIB diagnostic routines
* ---------------------------------------------------------------
* wcsprintf() is used by celprt(), linprt(), prjprt(), spcprt(), tabprt(),
* wcsprt(), and wcserr_prt() for diagnostic output which by default goes to
* stdout.  However, it may be redirected to a file or string buffer via
* wcsprintf_set().
*
* Given:
*   format    char*     Format string, passed to one of the printf(3) family
*                       of stdio library functions.
*
*   ...       mixed     Argument list matching format, as per printf(3).
*
* Function return value:
*             int       Number of bytes written.
*
*
* wcsfprintf() - Print function used by WCSLIB diagnostic routines
* ----------------------------------------------------------------
* wcsfprintf() is used by wcsbth(), and wcspih() for diagnostic output which
* they send to stderr.  However, it may be redirected to a file or string
* buffer via wcsprintf_set().
*
* Given:
*   stream    FILE*     The output stream if not overridden by a call to
*                       wcsprintf_set().
*
*   format    char*     Format string, passed to one of the printf(3) family
*                       of stdio library functions.
*
*   ...       mixed     Argument list matching format, as per printf(3).
*
* Function return value:
*             int       Number of bytes written.
*
*
* wcsprintf_set() - Set output disposition for wcsprintf() and wcsfprintf()
* -------------------------------------------------------------------------
* wcsprintf_set() sets the output disposition for wcsprintf() which is used by
* the celprt(), linprt(), prjprt(), spcprt(), tabprt(), wcsprt(), and
* wcserr_prt() routines, and for wcsfprintf() which is used by wcsbth() and
* wcspih().
*
* Given:
*   wcsout    FILE*     Pointer to an output stream that has been opened for
*                       writing, e.g. by the fopen() stdio library function,
*                       or one of the predefined stdio output streams - stdout
*                       and stderr.  If zero (NULL), output is written to an
*                       internally-allocated string buffer, the address of
*                       which may be obtained by wcsprintf_buf().
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*
*
* wcsprintf_buf() - Get the address of the internal string buffer
* ---------------------------------------------------------------
* wcsprintf_buf() returns the address of the internal string buffer created
* when wcsprintf_set() is invoked with its FILE* argument set to zero.
*
* Function return value:
*             const char *
*                       Address of the internal string buffer.  The user may
*                       free this buffer by calling wcsprintf_set() with a
*                       valid FILE*, e.g. stdout.  The free() stdlib library
*                       function must NOT be invoked on this const pointer.
*
*
* WCSPRINTF_PTR() macro - Print addresses in a consistent way
* -----------------------------------------------------------
* WCSPRINTF_PTR() is a preprocessor macro used to print addresses in a
* consistent way.
*
* On some systems the "%p" format descriptor renders a NULL pointer as the
* string "0x0".  On others, however, it produces "0" or even "(nil)".  On
* some systems a non-zero address is prefixed with "0x", on others, not.
*
* The WCSPRINTF_PTR() macro ensures that a NULL pointer is always rendered as
* "0x0" and that non-zero addresses are prefixed with "0x" thus providing
* consistency, for example, for comparing the output of test programs.
*
*===========================================================================*/

#ifndef WCSLIB_WCSPRINTF
#define WCSLIB_WCSPRINTF

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define WCSPRINTF_PTR(str1, ptr, str2) \
  if (ptr) { \
    wcsprintf("%s%#lx%s", (str1), (unsigned long)(ptr), (str2)); \
  } else { \
    wcsprintf("%s0x0%s", (str1), (str2)); \
  }

int wcsprintf_set(FILE *wcsout);
int wcsprintf(const char *format, ...);
int wcsfprintf(FILE *stream, const char *format, ...);
const char *wcsprintf_buf(void);

#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_WCSPRINTF */
