/*============================================================================

  WCSLIB 4.10 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2012, Mark Calabretta

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
  along with WCSLIB.  If not, see <http://www.gnu.org/licenses/>.

  Correspondence concerning WCSLIB may be directed to:
    Internet email: mcalabre@atnf.csiro.au
    Postal address: Dr. Mark Calabretta
                    Australia Telescope National Facility, CSIRO
                    PO Box 76
                    Epping NSW 1710
                    AUSTRALIA

  Author: Mark Calabretta, Australia Telescope National Facility
  http://www.atnf.csiro.au/~mcalabre/index.html
  $Id: log.c,v 4.10 2012/02/05 23:41:44 cal103 Exp $
*===========================================================================*/

#include <math.h>

#include "log.h"

/* Map status return value to message. */
const char *log_errmsg[] = {
  "Success",
  "",
  "Invalid log-coordinate reference value",
  "One or more of the x coordinates were invalid",
  "One or more of the world coordinates were invalid"};


/*--------------------------------------------------------------------------*/

int logx2s(
  double crval,
  int nx,
  int sx,
  int slogc,
  const double x[],
  double logc[],
  int stat[])

{
  register int ix;
  register int *statp;
  register const double *xp;
  register double *logcp;


  if (crval <= 0.0) {
    return LOGERR_BAD_LOG_REF_VAL;
  }

  xp = x;
  logcp = logc;
  statp = stat;
  for (ix = 0; ix < nx; ix++, xp += sx, logcp += slogc) {
    *logcp = crval * exp((*xp) / crval);
    *(statp++) = 0;
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int logs2x(
  double crval,
  int nlogc,
  int slogc,
  int sx,
  const double logc[],
  double x[],
  int stat[])

{
  int status;
  register int ilogc;
  register int *statp;
  register const double *logcp;
  register double *xp;


  if (crval <= 0.0) {
    return LOGERR_BAD_LOG_REF_VAL;
  }

  xp = x;
  logcp = logc;
  statp = stat;
  status = 0;
  for (ilogc = 0; ilogc < nlogc; ilogc++, logcp += slogc, xp += sx) {
    if (*logcp > 0.0) {
      *xp = crval * log(*logcp / crval);
      *(statp++) = 0;
    } else {
      *(statp++) = 1;
      status = LOGERR_BAD_WORLD;
    }
  }

  return status;
}
