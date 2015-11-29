/*============================================================================

  WCSLIB 5.10 - an implementation of the FITS WCS standard.
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
  $Id: wcstrig.h,v 5.10 2015/10/09 08:19:15 mcalabre Exp $
*=============================================================================
*
* WCSLIB 5.10 - C routines that implement the FITS World Coordinate System
* (WCS) standard.  Refer to the README file provided with WCSLIB for an
* overview of the library.
*
*
* Summary of the wcstrig routines
* -------------------------------
* When dealing with celestial coordinate systems and spherical projections
* (some moreso than others) it is often desirable to use an angular measure
* that provides an exact representation of the latitude of the north or south
* pole.  The WCSLIB routines use the following trigonometric functions that
* take or return angles in degrees:
*
*   - cosd()
*   - sind()
*   - tand()
*   - acosd()
*   - asind()
*   - atand()
*   - atan2d()
*   - sincosd()
*
* These "trigd" routines are expected to handle angles that are a multiple of
* 90 degrees returning an exact result.  Some C implementations provide these
* as part of a system library and in such cases it may (or may not!) be
* preferable to use them.  WCSLIB provides wrappers on the standard trig
* functions based on radian measure, adding tests for multiples of 90 degrees.
*
* However, wcstrig.h also provides the choice of using preprocessor macro
* implementations of the trigd functions that don't test for multiples of
* 90 degrees (compile with -DWCSTRIG_MACRO).  These are typically 20% faster
* but may lead to problems near the poles.
*
*
* cosd() - Cosine of an angle in degrees
* --------------------------------------
* cosd() returns the cosine of an angle given in degrees.
*
* Given:
*   angle     double    [deg].
*
* Function return value:
*             double    Cosine of the angle.
*
*
* sind() - Sine of an angle in degrees
* ------------------------------------
* sind() returns the sine of an angle given in degrees.
*
* Given:
*   angle     double    [deg].
*
* Function return value:
*             double    Sine of the angle.
*
*
* sincosd() - Sine and cosine of an angle in degrees
* --------------------------------------------------
* sincosd() returns the sine and cosine of an angle given in degrees.
*
* Given:
*   angle     double    [deg].
*
* Returned:
*   sin       *double   Sine of the angle.
*
*   cos       *double   Cosine of the angle.
*
* Function return value:
*             void
*
*
* tand() - Tangent of an angle in degrees
* ---------------------------------------
* tand() returns the tangent of an angle given in degrees.
*
* Given:
*   angle     double    [deg].
*
* Function return value:
*             double    Tangent of the angle.
*
*
* acosd() - Inverse cosine, returning angle in degrees
* ----------------------------------------------------
* acosd() returns the inverse cosine in degrees.
*
* Given:
*   x         double    in the range [-1,1].
*
* Function return value:
*             double    Inverse cosine of x [deg].
*
*
* asind() - Inverse sine, returning angle in degrees
* --------------------------------------------------
* asind() returns the inverse sine in degrees.
*
* Given:
*   y         double    in the range [-1,1].
*
* Function return value:
*             double    Inverse sine of y [deg].
*
*
* atand() - Inverse tangent, returning angle in degrees
* -----------------------------------------------------
* atand() returns the inverse tangent in degrees.
*
* Given:
*   s         double
*
* Function return value:
*             double    Inverse tangent of s [deg].
*
*
* atan2d() - Polar angle of (x,y), in degrees
* -------------------------------------------
* atan2d() returns the polar angle, beta, in degrees, of polar coordinates
* (rho,beta) corresponding to Cartesian coordinates (x,y).  It is equivalent
* to the arg(x,y) function of WCS Paper II, though with transposed arguments.
*
* Given:
*   y         double    Cartesian y-coordinate.
*
*   x         double    Cartesian x-coordinate.
*
* Function return value:
*             double    Polar angle of (x,y) [deg].
*
*===========================================================================*/

#ifndef WCSLIB_WCSTRIG
#define WCSLIB_WCSTRIG

#include <math.h>

#include "wcsconfig.h"

#ifdef HAVE_SINCOS
  void sincos(double angle, double *sin, double *cos);
#endif

#ifdef __cplusplus
extern "C" {
#endif


#ifdef WCSTRIG_MACRO

/* Macro implementation of the trigd functions. */
#include "wcsmath.h"

#define cosd(X) cos((X)*D2R)
#define sind(X) sin((X)*D2R)
#define tand(X) tan((X)*D2R)
#define acosd(X) acos(X)*R2D
#define asind(X) asin(X)*R2D
#define atand(X) atan(X)*R2D
#define atan2d(Y,X) atan2(Y,X)*R2D
#ifdef HAVE_SINCOS
  #define sincosd(X,S,C) sincos((X)*D2R,(S),(C))
#else
  #define sincosd(X,S,C) *(S) = sin((X)*D2R); *(C) = cos((X)*D2R);
#endif

#else

/* Use WCSLIB wrappers or native trigd functions. */

double cosd(double angle);
double sind(double angle);
void sincosd(double angle, double *sin, double *cos);
double tand(double angle);
double acosd(double x);
double asind(double y);
double atand(double s);
double atan2d(double y, double x);

/* Domain tolerance for asin() and acos() functions. */
#define WCSTRIG_TOL 1e-10

#endif /* WCSTRIG_MACRO */


#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_WCSTRIG */
