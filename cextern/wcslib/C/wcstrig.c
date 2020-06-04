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
  $Id: wcstrig.c,v 7.3 2020/06/03 03:37:02 mcalabre Exp $
*===========================================================================*/

#include <math.h>
#include <stdlib.h>
#include "wcsmath.h"
#include "wcstrig.h"

double cosd(double angle)

{
  int i;

  if (fmod(angle,90.0) == 0.0) {
    i = abs((int)floor(angle/90.0 + 0.5))%4;
    switch (i) {
    case 0:
      return 1.0;
    case 1:
      return 0.0;
    case 2:
      return -1.0;
    case 3:
      return 0.0;
    }
  }

  return cos(angle*D2R);
}

/*--------------------------------------------------------------------------*/

double sind(double angle)

{
  int i;

  if (fmod(angle,90.0) == 0.0) {
    i = abs((int)floor(angle/90.0 - 0.5))%4;
    switch (i) {
    case 0:
      return 1.0;
    case 1:
      return 0.0;
    case 2:
      return -1.0;
    case 3:
      return 0.0;
    }
  }

  return sin(angle*D2R);
}

/*--------------------------------------------------------------------------*/

void sincosd(double angle, double *s, double *c)

{
  int i;

  if (fmod(angle,90.0) == 0.0) {
    i = abs((int)floor(angle/90.0 + 0.5))%4;
    switch (i) {
    case 0:
      *s = 0.0;
      *c = 1.0;
      return;
    case 1:
      *s = (angle > 0.0) ? 1.0 : -1.0;
      *c = 0.0;
      return;
    case 2:
      *s =  0.0;
      *c = -1.0;
      return;
    case 3:
      *s = (angle > 0.0) ? -1.0 : 1.0;
      *c = 0.0;
      return;
    }
  }

#ifdef HAVE_SINCOS
  sincos(angle*D2R, s, c);
#else
  *s = sin(angle*D2R);
  *c = cos(angle*D2R);
#endif

  return;
}

/*--------------------------------------------------------------------------*/

double tand(double angle)

{
  double resid;

  resid = fmod(angle,360.0);
  if (resid == 0.0 || fabs(resid) == 180.0) {
    return 0.0;
  } else if (resid == 45.0 || resid == 225.0) {
    return 1.0;
  } else if (resid == -135.0 || resid == -315.0) {
    return -1.0;
  }

  return tan(angle*D2R);
}

/*--------------------------------------------------------------------------*/

double acosd(double v)

{
  if (v >= 1.0) {
    if (v-1.0 <  WCSTRIG_TOL) return 0.0;
  } else if (v == 0.0) {
    return 90.0;
  } else if (v <= -1.0) {
    if (v+1.0 > -WCSTRIG_TOL) return 180.0;
  }

  return acos(v)*R2D;
}

/*--------------------------------------------------------------------------*/

double asind(double v)

{
  if (v <= -1.0) {
    if (v+1.0 > -WCSTRIG_TOL) return -90.0;
  } else if (v == 0.0) {
    return 0.0;
  } else if (v >= 1.0) {
    if (v-1.0 <  WCSTRIG_TOL) return 90.0;
  }

  return asin(v)*R2D;
}

/*--------------------------------------------------------------------------*/

double atand(double v)

{
  if (v == -1.0) {
    return -45.0;
  } else if (v == 0.0) {
    return 0.0;
  } else if (v == 1.0) {
    return 45.0;
  }

  return atan(v)*R2D;
}

/*--------------------------------------------------------------------------*/

double atan2d(double y, double x)

{
  if (y == 0.0) {
    if (x >= 0.0) {
      return 0.0;
    } else if (x < 0.0) {
      return 180.0;
    }
  } else if (x == 0.0) {
    if (y > 0.0) {
      return 90.0;
    } else if (y < 0.0) {
      return -90.0;
    }
   }

   return atan2(y,x)*R2D;
}
