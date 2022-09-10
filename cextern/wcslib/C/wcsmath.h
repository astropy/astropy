/*============================================================================
  WCSLIB 7.12 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2022, Mark Calabretta

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
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: wcsmath.h,v 7.12 2022/09/09 04:57:58 mcalabre Exp $
*=============================================================================
*
* WCSLIB 7.12 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of wcsmath.h
* --------------------
* Definition of mathematical constants used by WCSLIB.
*
*===========================================================================*/

#ifndef WCSLIB_WCSMATH
#define WCSLIB_WCSMATH

#ifdef PI
#undef PI
#endif

#ifdef D2R
#undef D2R
#endif

#ifdef R2D
#undef R2D
#endif

#ifdef SQRT2
#undef SQRT2
#endif

#ifdef SQRT2INV
#undef SQRT2INV
#endif

#define PI 3.141592653589793238462643
#define D2R PI/180.0
#define R2D 180.0/PI
#define SQRT2 1.4142135623730950488
#define SQRT2INV 1.0/SQRT2

#ifdef UNDEFINED
#undef UNDEFINED
#endif

#define UNDEFINED 987654321.0e99
#define undefined(value) (value == UNDEFINED)

#endif // WCSLIB_WCSMATH
