/*============================================================================
  WCSLIB 8.4 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2024, Mark Calabretta

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
  $Id: cel.c,v 8.4 2024/10/28 13:56:16 mcalabre Exp $
*===========================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "wcserr.h"
#include "wcsmath.h"
#include "wcsprintf.h"
#include "wcstrig.h"
#include "sph.h"
#include "cel.h"

// Map status return value to message.
const char *cel_errmsg[] = {
  "Success",
  "Null celprm pointer passed",
  "Invalid projection parameters",
  "Invalid coordinate transformation parameters",
  "Ill-conditioned coordinate transformation parameters",
  "One or more of the (x,y) coordinates were invalid",
  "One or more of the (lng,lat) coordinates were invalid"};

// Map error returns for lower-level routines.
const int cel_prjerr[] = {
  CELERR_SUCCESS,		//  0: PRJERR_SUCCESS
  CELERR_NULL_POINTER,		//  1: PRJERR_NULL_POINTER
  CELERR_BAD_PARAM,		//  2: PRJERR_BAD_PARAM
  CELERR_BAD_PIX,		//  3: PRJERR_BAD_PIX
  CELERR_BAD_WORLD		//  4: PRJERR_BAD_WORLD
};

static const int CELSET = 137;

// Convenience macro for invoking wcserr_set().
#define CEL_ERRMSG(status) WCSERR_SET(status), cel_errmsg[status]

//----------------------------------------------------------------------------

int celini(struct celprm *cel)

{
  if (cel == 0x0) return CELERR_NULL_POINTER;

  cel->offset = 0;
  cel->phi0   = UNDEFINED;
  cel->theta0 = UNDEFINED;
  cel->ref[0] =   0.0;
  cel->ref[1] =   0.0;
  cel->ref[2] = UNDEFINED;
  cel->ref[3] = +90.0;

  for (int k = 0; k < 5; cel->euler[k++] = 0.0);
  cel->latpreq = -1;
  cel->isolat  =  0;

  cel->err = 0x0;

  cel->flag = 0;

  return cel_prjerr[prjini(&(cel->prj))];
}

//----------------------------------------------------------------------------

int celfree(struct celprm *cel)

{
  if (cel == 0x0) return CELERR_NULL_POINTER;

  wcserr_clear(&(cel->err));

  return cel_prjerr[prjfree(&(cel->prj))];
}

//----------------------------------------------------------------------------

int celsize(const struct celprm *cel, int sizes[2])

{
  if (cel == 0x0) {
    sizes[0] = sizes[1] = 0;
    return 0;
  }

  // Base size, in bytes.
  sizes[0] = sizeof(struct celprm);

  // Total size of allocated memory, in bytes.
  sizes[1] = 0;

  // celprm::prj.
  int exsizes[2];
  prjsize(&(cel->prj), exsizes);
  sizes[1] += exsizes[1];

  // celprm::err.
  wcserr_size(cel->err, exsizes);
  sizes[1] += exsizes[0] + exsizes[1];

  return 0;
}

//----------------------------------------------------------------------------

int celenq(const struct celprm *cel, int enquiry)

{
  // Initialize.
  if (cel == 0x0) return CELERR_NULL_POINTER;

  int answer = 0;

  if (enquiry & CELENQ_SET) {
    if (abs(cel->flag) != CELSET) return 0;
    answer = 1;
  }

  if (enquiry & CELENQ_BYP) {
    if (cel->flag != 1 && cel->flag != -CELSET) return 0;
    answer = 1;
  }

  return answer;
}

//----------------------------------------------------------------------------

int celprt(const struct celprm *cel)

{
  if (cel == 0x0) return CELERR_NULL_POINTER;

  // Parameters supplied.
  wcsprintf("       flag: %d\n", cel->flag);
  wcsprintf("     offset: %d\n", cel->offset);
  if (undefined(cel->phi0)) {
    wcsprintf("       phi0: UNDEFINED\n");
  } else {
    wcsprintf("       phi0: %9f\n", cel->phi0);
  }
  if (undefined(cel->theta0)) {
    wcsprintf("     theta0: UNDEFINED\n");
  } else {
    wcsprintf("     theta0: %9f\n", cel->theta0);
  }
  wcsprintf("        ref:");
  for (int i = 0; i < 4; i++) {
    wcsprintf("  %#- 11.5g", cel->ref[i]);
  }
  wcsprintf("\n");
  wcsprintf("        prj: (see below)\n");

  // Derived values.
  wcsprintf("      euler:");
  for (int i = 0; i < 5; i++) {
    wcsprintf("  %#- 11.5g", cel->euler[i]);
  }
  wcsprintf("\n");
  wcsprintf("    latpreq: %d", cel->latpreq);
  if (cel->latpreq == 0) {
    wcsprintf(" (not required)\n");
  } else if (cel->latpreq == 1) {
    wcsprintf(" (disambiguation)\n");
  } else if (cel->latpreq == 2) {
    wcsprintf(" (specification)\n");
  } else {
    wcsprintf(" (UNDEFINED)\n");
  }
  wcsprintf("     isolat: %d\n", cel->isolat);

  // Error handling.
  WCSPRINTF_PTR("        err: ", cel->err, "\n");
  if (cel->err) {
    wcserr_prt(cel->err, "             ");
  }

  // Projection parameters (from above).
  wcsprintf("\n");
  wcsprintf("   prj.*\n");
  prjprt(&(cel->prj));

  return 0;
}

//----------------------------------------------------------------------------

int celperr(const struct celprm *cel, const char *prefix)

{
  if (cel == 0x0) return CELERR_NULL_POINTER;

  if (cel->err && wcserr_prt(cel->err, prefix) == 0) {
    wcserr_prt(cel->prj.err, prefix);
  }

  return 0;
}

//----------------------------------------------------------------------------

int celset(struct celprm *cel)

{
  static const char *function = "celset";

  const double tol = 1.0e-10;

  if (cel == 0x0) return CELERR_NULL_POINTER;
  if (cel->flag == -CELSET) return 0;
  struct wcserr **err = &(cel->err);

  // Initialize the projection driver routines.
  struct prjprm *celprj = &(cel->prj);
  if (cel->offset) {
    celprj->phi0   = cel->phi0;
    celprj->theta0 = cel->theta0;
  } else {
    // Ensure that these are undefined - no fiducial offset.
    celprj->phi0   = UNDEFINED;
    celprj->theta0 = UNDEFINED;
  }

  celprj->flag = 0;
  int status;
  if ((status = prjset(celprj))) {
    return wcserr_set(CEL_ERRMSG(cel_prjerr[status]));
  }

  // Defaults set by the projection routines.
  if (undefined(cel->phi0)) {
    cel->phi0 = celprj->phi0;
  }

  if (undefined(cel->theta0)) {
    cel->theta0 = celprj->theta0;

  } else if (fabs(cel->theta0) > 90.0) {
    if (fabs(cel->theta0) > 90.0 + tol) {
      return wcserr_set(WCSERR_SET(CELERR_BAD_COORD_TRANS),
        "Invalid coordinate transformation parameters: theta0 > 90");
    }

    if (cel->theta0 > 90.0) {
      cel->theta0 =  90.0;
    } else {
      cel->theta0 = -90.0;
    }
  }


  double lng0 = cel->ref[0];
  double lat0 = cel->ref[1];
  double phip = cel->ref[2];
  double lngp;
  double latp = cel->ref[3];

  // Set default for native longitude of the celestial pole?
  if (undefined(phip) || phip == 999.0) {
    phip = (lat0 < cel->theta0) ? 180.0 : 0.0;
    phip += cel->phi0;

    if (phip < -180.0) {
      phip += 360.0;
    } else if (phip > 180.0) {
      phip -= 360.0;
    }

    cel->ref[2] = phip;
  }


  // Compute celestial coordinates of the native pole.
  cel->latpreq = 0;
  if (cel->theta0 == 90.0) {
    // Fiducial point at the native pole.
    lngp = lng0;
    latp = lat0;

  } else {
    // Fiducial point away from the native pole.
    double slat0, clat0, sthe0, cthe0;
    sincosd(lat0, &slat0, &clat0);
    sincosd(cel->theta0, &sthe0, &cthe0);

    double cphip, sphip, u, v;
    if (phip == cel->phi0) {
      sphip = 0.0;
      cphip = 1.0;

      u = cel->theta0;
      v = 90.0 - lat0;

    } else {
      sincosd(phip - cel->phi0, &sphip, &cphip);

      double x = cthe0*cphip;
      double y = sthe0;
      double z = sqrt(x*x + y*y);
      if (z == 0.0) {
        if (slat0 != 0.0) {
          return wcserr_set(WCSERR_SET(CELERR_BAD_COORD_TRANS),
            "Invalid coordinate description:\n"
            "lat0 == 0 is required for |phip - phi0| = 90 and theta0 == 0");
        }

        // latp determined solely by LATPOLEa in this case.
        cel->latpreq = 2;
        if (latp > 90.0) {
          latp = 90.0;
        } else if (latp < -90.0) {
          latp = -90.0;
        }

        // Avert a spurious compiler warning.
        u = v = 0.0;

      } else {
        double slz = slat0/z;
        if (fabs(slz) > 1.0) {
          if ((fabs(slz) - 1.0) < tol) {
            if (slz > 0.0) {
              slz = 1.0;
            } else {
              slz = -1.0;
            }
          } else {
            return wcserr_set(WCSERR_SET(CELERR_BAD_COORD_TRANS),
              "Invalid coordinate description:\n|lat0| <= %.3f is required "
              "for these values of phip, phi0, and theta0", asind(z));
          }
        }

        u = atan2d(y,x);
        v = acosd(slz);
      }
    }

    if (cel->latpreq == 0) {
      double latp1 = u + v;
      if (latp1 > 180.0) {
        latp1 -= 360.0;
      } else if (latp1 < -180.0) {
        latp1 += 360.0;
      }

      double latp2 = u - v;
      if (latp2 > 180.0) {
        latp2 -= 360.0;
      } else if (latp2 < -180.0) {
        latp2 += 360.0;
      }

      if (fabs(latp1) < 90.0+tol &&
          fabs(latp2) < 90.0+tol) {
        // There are two valid solutions for latp.
        cel->latpreq = 1;
      }

      if (fabs(latp-latp1) < fabs(latp-latp2)) {
        if (fabs(latp1) < 90.0+tol) {
          latp = latp1;
        } else {
          latp = latp2;
        }
      } else {
        if (fabs(latp2) < 90.0+tol) {
          latp = latp2;
        } else {
          latp = latp1;
        }
      }

      // Account for rounding error.
      if (fabs(latp) < 90.0+tol) {
        if (latp > 90.0) {
          latp =  90.0;
        } else if (latp < -90.0) {
          latp = -90.0;
        }
      }
    }

    double z = cosd(latp)*clat0;
    if (fabs(z) < tol) {
      if (fabs(clat0) < tol) {
        // Celestial pole at the fiducial point.
        lngp = lng0;

      } else if (latp > 0.0) {
        // Celestial north pole at the native pole.
        lngp = lng0 + phip - cel->phi0 - 180.0;

      } else {
        // Celestial south pole at the native pole.
        lngp = lng0 - phip + cel->phi0;
      }

    } else {
      double x = (sthe0 - sind(latp)*slat0)/z;
      double y =  sphip*cthe0/clat0;
      if (x == 0.0 && y == 0.0) {
        // Sanity check (shouldn't be possible).
        return wcserr_set(WCSERR_SET(CELERR_BAD_COORD_TRANS),
          "Invalid coordinate transformation parameters, internal error");
      }
      lngp = lng0 - atan2d(y,x);
    }

    // Make celestial longitude of the native pole the same sign as at the
    // fiducial point.
    if (lng0 >= 0.0) {
      if (lngp < 0.0) {
        lngp += 360.0;
      } else if (lngp > 360.0) {
        lngp -= 360.0;
      }
    } else {
      if (lngp > 0.0) {
        lngp -= 360.0;
      } else if (lngp < -360.0) {
        lngp += 360.0;
      }
    }
  }

  // Reset LATPOLEa.
  cel->ref[3] = latp;

  // Set the Euler angles.
  cel->euler[0] = lngp;
  cel->euler[1] = 90.0 - latp;
  cel->euler[2] = phip;
  sincosd(cel->euler[1], &cel->euler[4], &cel->euler[3]);
  cel->isolat = (cel->euler[4] == 0.0);

  // Check for ill-conditioned parameters.
  if (fabs(latp) > 90.0+tol) {
    return wcserr_set(WCSERR_SET(CELERR_ILL_COORD_TRANS),
      "Ill-conditioned coordinate transformation parameters\nNo valid "
      "solution for latp for these values of phip, phi0, and theta0");
  }

  cel->flag = (cel->flag == 1) ? -CELSET : CELSET;

  return 0;
}

//----------------------------------------------------------------------------

int celx2s(
  struct celprm *cel,
  int nx,
  int ny,
  int sxy,
  int sll,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  double lng[],
  double lat[],
  int    stat[])

{
  static const char *function = "celx2s";

  // Initialize.
  if (cel == 0x0) return CELERR_NULL_POINTER;
  struct wcserr **err = &(cel->err);

  int status = 0;
  if (abs(cel->flag) != CELSET) {
    if ((status = celset(cel))) return status;
  }

  // Apply spherical deprojection.
  int istat;
  struct prjprm *celprj = &(cel->prj);
  if ((istat = celprj->prjx2s(celprj, nx, ny, sxy, 1, x, y, phi, theta,
                               stat))) {
    if (istat) {
      status = wcserr_set(CEL_ERRMSG(cel_prjerr[istat]));
      if (status != CELERR_BAD_PIX) {
        return status;
      }
    }
  }

  int nphi = (ny > 0) ? (nx*ny) : nx;

  // Compute celestial coordinates.
  sphx2s(cel->euler, nphi, 0, 1, sll, phi, theta, lng, lat);

  return status;
}

//----------------------------------------------------------------------------

int cels2x(
  struct celprm *cel,
  int nlng,
  int nlat,
  int sll,
  int sxy,
  const double lng[],
  const double lat[],
  double phi[],
  double theta[],
  double x[],
  double y[],
  int    stat[])

{
  static const char *function = "cels2x";

  // Initialize.
  if (cel == 0x0) return CELERR_NULL_POINTER;
  struct wcserr **err = &(cel->err);

  int status = 0;
  if (abs(cel->flag) != CELSET) {
    if ((status = celset(cel))) return status;
  }

  // Compute native coordinates.
  sphs2x(cel->euler, nlng, nlat, sll, 1, lng, lat, phi, theta);

  int nphi, ntheta;
  if (cel->isolat) {
    // Constant celestial latitude -> constant native latitude.
    nphi   = nlng;
    ntheta = nlat;
  } else {
    nphi   = (nlat > 0) ? (nlng*nlat) : nlng;
    ntheta = 0;
  }

  // Apply the spherical projection.
  int istat;
  struct prjprm *celprj = &(cel->prj);
  if ((istat = celprj->prjs2x(celprj, nphi, ntheta, 1, sxy, phi, theta, x, y,
                               stat))) {
    if (istat) {
      status = wcserr_set(CEL_ERRMSG(cel_prjerr[istat]));
      if (status != CELERR_BAD_WORLD) {
        return status;
      }
    }
  }

  return status;
}
