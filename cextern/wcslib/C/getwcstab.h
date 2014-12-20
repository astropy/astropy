/*============================================================================

  WCSLIB 4.25 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2014, Mark Calabretta

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
  $Id: getwcstab.h,v 4.25 2014/12/14 14:29:36 mcalabre Exp $
*=============================================================================
*
* Summary of the getwcstab routines
* ---------------------------------
* fits_read_wcstab(), an implementation of a FITS table reading routine for
* 'TAB' coordinates, is provided for CFITSIO programmers.  It has been
* incorporated into CFITSIO as of v3.006 with the definitions in this file,
* getwcstab.h, moved into fitsio.h.
*
* fits_read_wcstab() is not included in the WCSLIB object library but the
* source code is presented here as it may be useful for programmers using an
* older version of CFITSIO than 3.006, or as a programming template for
* non-CFITSIO programmers.
*
*
* fits_read_wcstab() - FITS 'TAB' table reading routine
* ----------------------------------------------------
* fits_read_wcstab() extracts arrays from a binary table required in
* constructing 'TAB' coordinates.
*
* Given:
*   fptr      fitsfile *
*                       Pointer to the file handle returned, for example, by
*                       the fits_open_file() routine in CFITSIO.
*
*   nwtb      int       Number of arrays to be read from the binary table(s).
*
* Given and returned:
*   wtb       wtbarr *  Address of the first element of an array of wtbarr
*                       typedefs.  This wtbarr typedef is defined to match the
*                       wtbarr struct defined in WCSLIB.  An array of such
*                       structs returned by the WCSLIB function wcstab() as
*                       discussed in the notes below.
*
* Returned:
*   status    int *     CFITSIO status value.
*
* Function return value:
*             int       CFITSIO status value.
*
* Notes:
*   In order to maintain WCSLIB and CFITSIO as independent libraries it is not
*   permissible for any CFITSIO library code to include WCSLIB header files,
*   or vice versa.  However, the CFITSIO function fits_read_wcstab() accepts
*   an array of wtbarr structs defined in wcs.h within WCSLIB.
*
*   The problem therefore is to define the wtbarr struct within fitsio.h
*   without including wcs.h, especially noting that wcs.h will often (but not
*   always) be included together with fitsio.h in an applications program that
*   uses fits_read_wcstab().
*
*   The solution adopted is for WCSLIB to define "struct wtbarr" while
*   fitsio.h defines "typedef wtbarr" as an untagged struct with identical
*   members.  This allows both wcs.h and fitsio.h to define a wtbarr data type
*   without conflict by virtue of the fact that structure tags and typedef
*   names share different name spaces in C; Appendix A, Sect. A11.1 (p227) of
*   the K&R ANSI edition states that:
*
*     Identifiers fall into several name spaces that do not interfere with one
*     another; the same identifier may be used for different purposes, even in
*     the same scope, if the uses are in different name spaces.  These classes
*     are: objects, functions, typedef names, and enum constants; labels; tags
*     of structures, unions, and enumerations; and members of each structure
*     or union individually.
*
*   Therefore, declarations within WCSLIB look like
*
=     struct wtbarr *w;
*
*   while within CFITSIO they are simply
*
=     wtbarr *w;
*
*   As suggested by the commonality of the names, these are really the same
*   aggregate data type.  However, in passing a (struct wtbarr *) to
*   fits_read_wcstab() a cast to (wtbarr *) is formally required.
*
*   When using WCSLIB and CFITSIO together in C++ the situation is complicated
*   by the fact that typedefs and structs share the same namespace; C++
*   Annotated Reference Manual, Sect. 7.1.3 (p105).  In that case the wtbarr
*   struct in wcs.h is renamed by preprocessor macro substitution to wtbarr_s
*   to distinguish it from the typedef defined in fitsio.h.  However, the
*   scope of this macro substitution is limited to wcs.h itself and CFITSIO
*   programmer code, whether in C++ or C, should always use the wtbarr
*   typedef.
*
*
* wtbarr typedef
* --------------
* The wtbarr typedef is defined as a struct containing the following members:
*
*   int i
*     Image axis number.
*
*   int m
*     Array axis number for index vectors.
*
*   int kind
*     Character identifying the array type:
*       - c: coordinate array,
*       - i: index vector.
*
*   char extnam[72]
*     EXTNAME identifying the binary table extension.
*
*   int extver
*     EXTVER identifying the binary table extension.
*
*   int extlev
*     EXTLEV identifying the binary table extension.
*
*   char ttype[72]
*     TTYPEn identifying the column of the binary table that contains the
*     array.
*
*   long row
*     Table row number.
*
*   int ndim
*     Expected dimensionality of the array.
*
*   int *dimlen
*     Address of the first element of an array of int of length ndim into
*     which the array axis lengths are to be written.
*
*   double **arrayp
*     Pointer to an array of double which is to be allocated by the user
*     and into which the array is to be written.
*
*===========================================================================*/

#ifndef WCSLIB_GETWCSTAB
#define WCSLIB_GETWCSTAB

#ifdef __cplusplus
extern "C" {
#endif

#include <fitsio.h>

typedef struct {
  int  i;			/* Image axis number.                       */
  int  m;			/* Array axis number for index vectors.     */
  int  kind;			/* Array type, 'c' (coord) or 'i' (index).  */
  char extnam[72];		/* EXTNAME of binary table extension.       */
  int  extver;			/* EXTVER  of binary table extension.       */
  int  extlev;			/* EXTLEV  of binary table extension.       */
  char ttype[72];		/* TTYPEn of column containing the array.   */
  long row;			/* Table row number.                        */
  int  ndim;			/* Expected array dimensionality.           */
  int  *dimlen;			/* Where to write the array axis lengths.   */
  double **arrayp;		/* Where to write the address of the array  */
				/* allocated to store the array.            */
} wtbarr;


int fits_read_wcstab(fitsfile *fptr, int nwtb, wtbarr *wtb, int *status);


#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_GETWCSTAB */
