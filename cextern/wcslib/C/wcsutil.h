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
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: wcsutil.h,v 5.9 2015/07/21 09:20:01 mcalabre Exp $
*=============================================================================
*
* WCSLIB 5.9 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of the wcsutil routines
* -------------------------------
* Simple utility functions for internal use only by WCSLIB.  They are
* documented here solely as an aid to understanding the code.  They are not
* intended for external use - the API may change without notice!
*
*
* wcsutil_blank_fill() - Fill a character string with blanks
* ----------------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_blank_fill() pads a character string with blanks starting with the
* terminating NULL character.
*
* Used by the Fortran wrapper functions in translating C character strings
* into Fortran CHARACTER variables.
*
* Given:
*   n         int       Length of the character array, c[].
*
* Given and returned:
*   c         char[]    The character string.  It will not be null-terminated
*                       on return.
*
* Function return value:
*             void
*
*
* wcsutil_null_fill() - Fill a character string with NULLs
* --------------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_null_fill() strips off trailing blanks and pads the character array
* holding the string with NULL characters.
*
* Used mainly to make character strings intelligible in the GNU debugger which
* prints the rubbish following the terminating NULL, obscuring the valid part
* of the string.
*
* Given:
*   n         int       Number of characters.
*
* Given and returned:
*   c         char[]    The character string.
*
* Function return value:
*             void
*
*
* wcsutil_allEq() - Test for equality of a particular vector element
* ------------------------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_allEq() tests for equality of a particular element in a set of
* vectors.
*
* Given:
*   nvec      int       The number of vectors.
*
*   nelem     int       The length of each vector.
*
*   first     const double*
*                       Pointer to the first element to test in the array.
*                       The elements tested for equality are
*
=                         *first == *(first + nelem)
=                                == *(first + nelem*2)
=                                           :
=                                == *(first + nelem*(nvec-1));
*
*                       The array might be dimensioned as
*
=                         double v[nvec][nelem];
*
* Function return value:
*             int       Status return value:
*                         0: Not all equal.
*                         1: All equal.
*
*
* wcsutil_Eq() - Test for equality of two double arrays
* -----------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_Eq() tests for equality of two double-precision arrays.
*
* Given:
*   nelem     int       The number of elements in each array.
*
*   tol       double    Tolerance for comparison of the floating-point values.
*                       For example, for tol == 1e-6, all floating-point
*                       values in the arrays must be equal to the first 6
*                       decimal places.  A value of 0 implies exact equality.
*
*   arr1      const double*
*                       The first array.
*
*   arr2      const double*
*                       The second array
*
* Function return value:
*             int       Status return value:
*                         0: Not equal.
*                         1: Equal.
*
*
* wcsutil_intEq() - Test for equality of two int arrays
* -----------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_intEq() tests for equality of two int arrays.
*
* Given:
*   nelem     int       The number of elements in each array.
*
*   arr1      const int*
*                       The first array.
*
*   arr2      const int*
*                       The second array
*
* Function return value:
*             int       Status return value:
*                         0: Not equal.
*                         1: Equal.
*
*
* wcsutil_strEq() - Test for equality of two string arrays
* --------------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_strEq() tests for equality of two string arrays.
*
* Given:
*   nelem     int       The number of elements in each array.
*
*   arr1      const char**
*                       The first array.
*
*   arr2      const char**
*                       The second array
*
* Function return value:
*             int       Status return value:
*                         0: Not equal.
*                         1: Equal.
*
*
* wcsutil_setAll() - Set a particular vector element
* --------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_setAll() sets the value of a particular element in a set of vectors.
*
* Given:
*   nvec      int       The number of vectors.
*
*   nelem     int       The length of each vector.
*
* Given and returned:
*   first     double*   Pointer to the first element in the array, the value
*                       of which is used to set the others
*
=                         *(first + nelem) = *first;
=                         *(first + nelem*2) = *first;
=                                 :
=                         *(first + nelem*(nvec-1)) = *first;
*
*                       The array might be dimensioned as
*
=                         double v[nvec][nelem];
*
* Function return value:
*             void
*
*
* wcsutil_setAli() - Set a particular vector element
* --------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_setAli() sets the value of a particular element in a set of vectors.
*
* Given:
*   nvec      int       The number of vectors.
*
*   nelem     int       The length of each vector.
*
* Given and returned:
*   first     int*      Pointer to the first element in the array, the value
*                       of which is used to set the others
*
=                         *(first + nelem) = *first;
=                         *(first + nelem*2) = *first;
=                                 :
=                         *(first + nelem*(nvec-1)) = *first;
*
*                       The array might be dimensioned as
*
=                         int v[nvec][nelem];
*
* Function return value:
*             void
*
*
* wcsutil_setBit() - Set bits in selected elements of an array
* ------------------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_setBit() sets bits in selected elements of an array.
*
* Given:
*   nelem     int       Number of elements in the array.
*
*   sel       const int*
*                       Address of a selection array of length nelem.  May
*                       be specified as the null pointer in which case all
*                       elements are selected.
*
*   bits      int       Bit mask.
*
* Given and returned:
*   array     int*      Address of the array of length nelem.
*
* Function return value:
*             void
*
*
* wcsutil_fptr2str() - Translate pointer-to-function to string
* ------------------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_fptr2str() translates a pointer-to-function to hexadecimal string
* representation for output.  It is used by the various routines that print
* the contents of WCSLIB structs, noting that it is not strictly legal to
* type-pun a function pointer to void*.  See
* http://stackoverflow.com/questions/2741683/how-to-format-a-function-pointer
*
* Given:
*   fptr      int(*)()  Pointer to function.
*
* Returned:
*   hext      char[19]  Null-terminated string.  Should be at least 19 bytes
*                       in size to accomodate a 64-bit address (16 bytes in
*                       hex), plus the leading "0x" and trailing '\0'.
*
* Function return value:
*             char *    The address of hext.
*
*
* wcsutil_double2str() - Translate double to string ignoring the locale
* ---------------------------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_double2str() converts a double to a string, but unlike sprintf() it
* ignores the locale and always uses a '.' as the decimal separator.  Also,
* unless it includes an exponent, the formatted value will always have a
* fractional part, ".0" being appended if necessary.
*
* Returned:
*   buf       char *    The buffer to write the string into.
*
* Given:
*   format    char *    The formatting directive, such as "%f".  This
*                       may be any of the forms accepted by sprintf(), but
*                       should only include a formatting directive and
*                       nothing else.  For "%g" and "%G" formats, unless it
*                       includes an exponent, the formatted value will always
*                       have a fractional part, ".0" being appended if
*                       necessary.
*
*   value     double    The value to convert to a string.
*
*
* wcsutil_str2double() - Translate string to a double, ignoring the locale
* ------------------------------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_str2double() converts a string to a double, but unlike sscanf() it
* ignores the locale and always expects a '.' as the decimal separator.
*
* Given:
*   buf       char *    The string containing the value
*
*   format    char *    The formatting directive, such as "%lf".  This
*                       may be any of the forms accepted by sscanf(), but
*                       should only include a single formatting directive.
*
* Returned:
*   value     double *  The double value parsed from the string.
*
*
* wcsutil_dpkey_int() - Get the data value in a dpkey struct as int
* -----------------------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_dpkey_int() returns the data value in a dpkey struct as an integer
* value.
*
* Given and returned:
*   dp        const struct dpkey *
*                       Parsed contents of a DPja or DQia keyrecords.
*
* Function return value:
*             int       The record's value as int.
*
*
* wcsutil_dpkey_double() - Get the data value in a dpkey struct as double
* -----------------------------------------------------------------------
* INTERNAL USE ONLY.
*
* wcsutil_dpkey_double() returns the data value in a dpkey struct as a
* floating point value.
*
* Given and returned:
*   dp        const struct dpkey *
*                       Parsed contents of a DPja or DQia keyrecords.
*
* Function return value:
*             double    The record's value as double.
*
*===========================================================================*/

#ifndef WCSLIB_WCSUTIL
#define WCSLIB_WCSUTIL

#include "dis.h"

#ifdef __cplusplus
extern "C" {
#endif

void wcsutil_blank_fill(int n, char c[]);
void wcsutil_null_fill (int n, char c[]);

int  wcsutil_allEq (int nvec, int nelem, const double *first);
int  wcsutil_Eq(int nelem, double tol, const double *arr1,
                const double *arr2);
int  wcsutil_intEq(int nelem, const int *arr1, const int *arr2);
int  wcsutil_strEq(int nelem, char (*arr1)[72], char (*arr2)[72]);
void wcsutil_setAll(int nvec, int nelem, double *first);
void wcsutil_setAli(int nvec, int nelem, int *first);
void wcsutil_setBit(int nelem, const int *sel, int bits, int *array);
char *wcsutil_fptr2str(int (*func)(void), char hext[19]);
int  wcsutil_str2double(const char *buf, const char *format, double *value);
void wcsutil_double2str(char *buf, const char *format, double value);
int    wcsutil_dpkey_int(const struct dpkey *dp);
double wcsutil_dpkey_double(const struct dpkey *dp);

#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_WCSUTIL */
