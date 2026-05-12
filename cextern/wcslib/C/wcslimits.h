/*============================================================================
  WCSLIB 8.7 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2026, Mark Calabretta

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

  Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
  http://www.atnf.csiro.au/computing/software/wcs
  $Id: wcslimits.h,v 8.7 2026/05/11 12:01:10 mcalabre Exp $
*=============================================================================
*
* WCSLIB 8.7 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of wcslimits.h
* ----------------------
* Declaration of global variables that set limits for arrays used by WCSLIB.
* The external variables are defined in wcs.c.
*
*===========================================================================*/

#ifndef WCSLIB_WCSLIMITS
#define WCSLIB_WCSLIMITS

// Maximum number of image axes, NAXIS.  Limited to 31 by cylfix().
#define NAXMAX 31

// Maximum number of PVi_ma and PSi_ma keywords.  May be changed by wcsnpv().
extern int NPVMAX;
extern int NPSMAX;

// Maximum number of DPja or DQia keywords.  May be changed by disndp().
extern int NDPMAX;

#endif // WCSLIB_WCSLIMITS
