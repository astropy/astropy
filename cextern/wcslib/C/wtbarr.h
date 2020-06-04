/*============================================================================

  WCSLIB 7.3 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2020, Mark Calabretta

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
  $Id: wtbarr.h,v 7.3 2020/06/03 03:37:02 mcalabre Exp $
*=============================================================================
*
* WCSLIB 7.3 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of the wtbarr struct
* ----------------------------
* The wtbarr struct is used by wcstab() in extracting coordinate lookup tables
* from a binary table extension (BINTABLE) and copying them into the tabprm
* structs stored in wcsprm.
*
*
* wtbarr struct - Extraction of coordinate lookup tables from BINTABLE
* --------------------------------------------------------------------
* Function wcstab(), which is invoked automatically by wcspih(), sets up an
* array of wtbarr structs to assist in extracting coordinate lookup tables
* from a binary table extension (BINTABLE) and copying them into the tabprm
* structs stored in wcsprm.  Refer to the usage notes for wcspih() and
* wcstab() in wcshdr.h, and also the prologue to tab.h.
*
* For C++ usage, because of a name space conflict with the wtbarr typedef
* defined in CFITSIO header fitsio.h, the wtbarr struct is renamed to wtbarr_s
* by preprocessor macro substitution with scope limited to wtbarr.h itself,
* and similarly in wcs.h.
*
*   int i
*     (Given) Image axis number.
*
*   int m
*     (Given) wcstab array axis number for index vectors.
*
*   int kind
*     (Given) Character identifying the wcstab array type:
*       - c: coordinate array,
*       - i: index vector.
*
*   char extnam[72]
*     (Given) EXTNAME identifying the binary table extension.
*
*   int extver
*     (Given) EXTVER identifying the binary table extension.
*
*   int extlev
*     (Given) EXTLEV identifying the binary table extension.
*
*   char ttype[72]
*     (Given) TTYPEn identifying the column of the binary table that contains
*     the wcstab array.
*
*   long row
*     (Given) Table row number.
*
*   int ndim
*     (Given) Expected dimensionality of the wcstab array.
*
*   int *dimlen
*     (Given) Address of the first element of an array of int of length ndim
*     into which the wcstab array axis lengths are to be written.
*
*   double **arrayp
*     (Given) Pointer to an array of double which is to be allocated by the
*     user and into which the wcstab array is to be written.
*
*===========================================================================*/

#ifndef WCSLIB_WTBARR
#define WCSLIB_WTBARR

#ifdef __cplusplus
extern "C" {
#define wtbarr wtbarr_s		/* See prologue above.                      */
#endif
				/* For extracting wcstab arrays.  Matches   */
				/* the wtbarr typedef defined in CFITSIO    */
				/* header fitsio.h.                         */
struct wtbarr {
  int  i;			/* Image axis number.                       */
  int  m;			/* Array axis number for index vectors.     */
  int  kind;			/* wcstab array type.                       */
  char extnam[72];		/* EXTNAME of binary table extension.       */
  int  extver;			/* EXTVER  of binary table extension.       */
  int  extlev;			/* EXTLEV  of binary table extension.       */
  char ttype[72];		/* TTYPEn of column containing the array.   */
  long row;			/* Table row number.                        */
  int  ndim;			/* Expected wcstab array dimensionality.    */
  int  *dimlen;			/* Where to write the array axis lengths.   */
  double **arrayp;		/* Where to write the address of the array  */
				/* allocated to store the wcstab array.     */
};

#ifdef __cplusplus
#undef wtbarr
}
#endif

#endif /* WCSLIB_WTBARR */
