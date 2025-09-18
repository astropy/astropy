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
  $Id: prj.c,v 8.4 2024/10/28 13:56:16 mcalabre Exp $
*===========================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcserr.h"
#include "wcsmath.h"
#include "wcsprintf.h"
#include "wcstrig.h"
#include "wcsutil.h"
#include "prj.h"


// Projection categories.
const int ZENITHAL          = 1;
const int CYLINDRICAL       = 2;
const int PSEUDOCYLINDRICAL = 3;
const int CONVENTIONAL      = 4;
const int CONIC             = 5;
const int POLYCONIC         = 6;
const int QUADCUBE          = 7;
const int HEALPIX           = 8;

const char prj_categories[9][32] =
  {"undefined", "zenithal", "cylindrical", "pseudocylindrical",
  "conventional", "conic", "polyconic", "quadcube", "HEALPix"};


// Projection codes.
const int  prj_ncode = 28;
const char prj_codes[28][4] =
  {"AZP", "SZP", "TAN", "STG", "SIN", "ARC", "ZPN", "ZEA", "AIR", "CYP",
   "CEA", "CAR", "MER", "COP", "COE", "COD", "COO", "SFL", "PAR", "MOL",
   "AIT", "BON", "PCO", "TSC", "CSC", "QSC", "HPX", "XPH"};

// Map status return value to message.
const char *prj_errmsg[] = {
  "Success",
  "Null prjprm pointer passed",
  "Invalid projection parameters",
  "One or more of the (x,y) coordinates were invalid",
  "One or more of the (phi,theta) coordinates were invalid"};

static const int AZP = 101;
static const int SZP = 102;
static const int TAN = 103;
static const int STG = 104;
static const int SIN = 105;
static const int ARC = 106;
static const int ZPN = 107;
static const int ZEA = 108;
static const int AIR = 109;
static const int CYP = 201;
static const int CEA = 202;
static const int CAR = 203;
static const int MER = 204;
static const int SFL = 301;
static const int PAR = 302;
static const int MOL = 303;
static const int AIT = 401;
static const int COP = 501;
static const int COE = 502;
static const int COD = 503;
static const int COO = 504;
static const int BON = 601;
static const int PCO = 602;
static const int TSC = 701;
static const int CSC = 702;
static const int QSC = 703;
static const int HPX = 801;
static const int XPH = 802;

// Convenience macros for generating common error messages.
#define PRJERR_BAD_PARAM_SET(function) \
  wcserr_set(&(prj->err), PRJERR_BAD_PARAM, function, __FILE__, __LINE__, \
    "Invalid parameters for %s projection", prj->name);

#define PRJERR_BAD_PIX_SET(function) \
  wcserr_set(&(prj->err), PRJERR_BAD_PIX, function, __FILE__, __LINE__, \
    "One or more of the (x, y) coordinates were invalid for %s projection", \
    prj->name);

#define PRJERR_BAD_WORLD_SET(function) \
  wcserr_set(&(prj->err), PRJERR_BAD_WORLD, function, __FILE__, __LINE__, \
    "One or more of the (lat, lng) coordinates were invalid for " \
    "%s projection", prj->name);

#define copysign(X, Y) ((Y) < 0.0 ? -fabs(X) : fabs(X))


/*============================================================================
* Generic routines:
*
* prjini initializes a prjprm struct to default values.
*
* prjfree frees any memory that may have been allocated to store an error
*        message in the prjprm struct.
*
* prjsize computes the size of a prjprm struct.
*
* prjprt prints the contents of a prjprm struct.
*
* prjbchk performs bounds checking on the native coordinates returned by the
*        *x2s() routines.
*
* prjset invokes the specific initialization routine based on the projection
*        code in the prjprm struct.
*
* prjx2s invokes the specific deprojection routine based on the pointer-to-
*        function stored in the prjprm struct.
*
* prjs2x invokes the specific projection routine based on the pointer-to-
*        function stored in the prjprm struct.
*
*---------------------------------------------------------------------------*/

int prjini(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  strcpy(prj->code, "   ");
  prj->pv[0]  = 0.0;
  prj->pv[1]  = UNDEFINED;
  prj->pv[2]  = UNDEFINED;
  prj->pv[3]  = UNDEFINED;
  for (int k = 4; k < PVN; prj->pv[k++] = 0.0);
  prj->r0     = 0.0;
  prj->phi0   = UNDEFINED;
  prj->theta0 = UNDEFINED;
  prj->bounds = 7;

  strcpy(prj->name, "undefined");
  for (int k = 9; k < 40; prj->name[k++] = '\0');
  prj->category  = 0;
  prj->pvrange   = 0;
  prj->simplezen = 0;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 0;
  prj->divergent = 0;
  prj->x0 = 0.0;
  prj->y0 = 0.0;

  prj->err = 0x0;

  prj->padding = 0x0;
  for (int k = 0; k < 10; prj->w[k++] = 0.0);
  prj->m = 0;
  prj->n = 0;
  prj->prjx2s = 0x0;
  prj->prjs2x = 0x0;

  prj->flag = 0;

  return 0;
}

//----------------------------------------------------------------------------

int prjfree(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  wcserr_clear(&(prj->err));

  return 0;
}

//----------------------------------------------------------------------------

int prjsize(const struct prjprm *prj, int sizes[2])

{
  if (prj == 0x0) {
    sizes[0] = sizes[1] = 0;
    return 0;
  }

  // Base size, in bytes.
  sizes[0] = sizeof(struct prjprm);

  // Total size of allocated memory, in bytes.
  sizes[1] = 0;

  // prjprm::err.
  int exsizes[2];
  wcserr_size(prj->err, exsizes);
  sizes[1] += exsizes[0] + exsizes[1];

  return 0;
}

//----------------------------------------------------------------------------

int prjenq(const struct prjprm *prj, int enquiry)

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int answer = 0;
  int absflag = abs(prj->flag);

  if (enquiry & PRJENQ_SET) {
    if (absflag < 100 || 1000 < absflag) return 0;
    answer = 1;
  }

  if (enquiry & PRJENQ_BYP) {
    if (prj->flag != 1 && !(-1000 < prj->flag && prj->flag < -100)) return 0;
    answer = 1;
  }

  return answer;
}

//----------------------------------------------------------------------------

int prjprt(const struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  // Parameters supplied.
  wcsprintf("       flag: %d\n",  prj->flag);
  wcsprintf("       code: \"%s\"\n", prj->code);
  wcsprintf("         r0: %9f\n", prj->r0);
  wcsprintf("         pv:");
  if (prj->pvrange) {
    int n = (prj->pvrange)%100;

    if (prj->pvrange/100) {
      wcsprintf(" (0)");
    } else {
      wcsprintf(" %#- 11.5g", prj->pv[0]);
      n--;
    }

    for (int i = 1; i <= n; i++) {
      if (i%5 == 1) {
        wcsprintf("\n           ");
      }

      if (undefined(prj->pv[i])) {
        wcsprintf("  UNDEFINED   ");
      } else {
        wcsprintf("  %#- 11.5g", prj->pv[i]);
      }
    }
    wcsprintf("\n");
  } else {
    wcsprintf(" (not used)\n");
  }
  if (undefined(prj->phi0)) {
    wcsprintf("       phi0: UNDEFINED\n");
  } else {
    wcsprintf("       phi0: %9f\n", prj->phi0);
  }
  if (undefined(prj->theta0)) {
    wcsprintf("     theta0: UNDEFINED\n");
  } else {
    wcsprintf("     theta0: %9f\n", prj->theta0);
  }
  wcsprintf("     bounds: %d\n",  prj->bounds);

  // Derived values.
  wcsprintf("\n");
  wcsprintf("       name: \"%s\"\n", prj->name);
  wcsprintf("   category: %d (%s)\n", prj->category,
                                      prj_categories[prj->category]);
  wcsprintf("    pvrange: %d\n", prj->pvrange);
  wcsprintf("  simplezen: %d\n", prj->simplezen);
  wcsprintf("  equiareal: %d\n", prj->equiareal);
  wcsprintf("  conformal: %d\n", prj->conformal);
  wcsprintf("     global: %d\n", prj->global);
  wcsprintf("  divergent: %d\n", prj->divergent);
  wcsprintf("         x0: %f\n", prj->x0);
  wcsprintf("         y0: %f\n", prj->y0);

  // Error handling.
  WCSPRINTF_PTR("        err: ", prj->err, "\n");
  if (prj->err) {
    wcserr_prt(prj->err, "             ");
  }

  // Intermediate values set by prjset().
  wcsprintf("        w[]:");
  for (int i = 0; i < 5; i++) {
    wcsprintf("  %#- 11.5g", prj->w[i]);
  }
  wcsprintf("\n            ");
  for (int i = 5; i < 10; i++) {
    wcsprintf("  %#- 11.5g", prj->w[i]);
  }
  wcsprintf("\n");
  wcsprintf("          m: %d\n", prj->m);
  wcsprintf("          n: %d\n", prj->n);

  // Pointers to projection and deprojection functions.
  char hext[32];
  wcsprintf("     prjx2s: %s\n",
    wcsutil_fptr2str((void (*)(void))prj->prjx2s, hext));
  wcsprintf("     prjs2x: %s\n",
    wcsutil_fptr2str((void (*)(void))prj->prjs2x, hext));

  return 0;
}

//----------------------------------------------------------------------------

int prjperr(const struct prjprm *prj, const char *prefix)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  if (prj->err) {
    wcserr_prt(prj->err, prefix);
  }

  return 0;
}

//----------------------------------------------------------------------------

int prjbchk(
  double tol,
  int nphi,
  int ntheta,
  int spt,
  double phi[],
  double theta[],
  int stat[])

{
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  int    status  = 0;
  for (int itheta = 0; itheta < ntheta; itheta++) {
    for (int iphi = 0; iphi < nphi; iphi++, phip += spt, thetap += spt,
         statp++) {
      // Skip values already marked as illegal.
      if (*statp == 0) {
        if (*phip < -180.0) {
          if (*phip < -180.0-tol) {
            *statp = 1;
            status = 1;
          } else {
            *phip = -180.0;
          }
        } else if (180.0 < *phip) {
          if (180.0+tol < *phip) {
            *statp = 1;
            status = 1;
          } else {
            *phip = 180.0;
          }
        }

        if (*thetap < -90.0) {
          if (*thetap < -90.0-tol) {
            *statp = 1;
            status = 1;
          } else {
            *thetap = -90.0;
          }
        } else if (90.0 < *thetap) {
          if (90.0+tol < *thetap) {
            *statp = 1;
            status = 1;
          } else {
            *thetap = 90.0;
          }
        }
      }
    }
  }

  return status;
}

//----------------------------------------------------------------------------

int prjset(struct prjprm *prj)

{
  static const char *function = "prjset";

  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag < 0) return 0;
  struct wcserr **err = &(prj->err);

  // Invoke the relevant initialization routine.
  prj->code[3] = '\0';

  int status;
  if (strcmp(prj->code, "AZP") == 0) {
    status = azpset(prj);
  } else if (strcmp(prj->code, "SZP") == 0) {
    status = szpset(prj);
  } else if (strcmp(prj->code, "TAN") == 0) {
    status = tanset(prj);
  } else if (strcmp(prj->code, "STG") == 0) {
    status = stgset(prj);
  } else if (strcmp(prj->code, "SIN") == 0) {
    status = sinset(prj);
  } else if (strcmp(prj->code, "ARC") == 0) {
    status = arcset(prj);
  } else if (strcmp(prj->code, "ZPN") == 0) {
    status = zpnset(prj);
  } else if (strcmp(prj->code, "ZEA") == 0) {
    status = zeaset(prj);
  } else if (strcmp(prj->code, "AIR") == 0) {
    status = airset(prj);
  } else if (strcmp(prj->code, "CYP") == 0) {
    status = cypset(prj);
  } else if (strcmp(prj->code, "CEA") == 0) {
    status = ceaset(prj);
  } else if (strcmp(prj->code, "CAR") == 0) {
    status = carset(prj);
  } else if (strcmp(prj->code, "MER") == 0) {
    status = merset(prj);
  } else if (strcmp(prj->code, "SFL") == 0) {
    status = sflset(prj);
  } else if (strcmp(prj->code, "PAR") == 0) {
    status = parset(prj);
  } else if (strcmp(prj->code, "MOL") == 0) {
    status = molset(prj);
  } else if (strcmp(prj->code, "AIT") == 0) {
    status = aitset(prj);
  } else if (strcmp(prj->code, "COP") == 0) {
    status = copset(prj);
  } else if (strcmp(prj->code, "COE") == 0) {
    status = coeset(prj);
  } else if (strcmp(prj->code, "COD") == 0) {
    status = codset(prj);
  } else if (strcmp(prj->code, "COO") == 0) {
    status = cooset(prj);
  } else if (strcmp(prj->code, "BON") == 0) {
    status = bonset(prj);
  } else if (strcmp(prj->code, "PCO") == 0) {
    status = pcoset(prj);
  } else if (strcmp(prj->code, "TSC") == 0) {
    status = tscset(prj);
  } else if (strcmp(prj->code, "CSC") == 0) {
    status = cscset(prj);
  } else if (strcmp(prj->code, "QSC") == 0) {
    status = qscset(prj);
  } else if (strcmp(prj->code, "HPX") == 0) {
    status = hpxset(prj);
  } else if (strcmp(prj->code, "XPH") == 0) {
    status = xphset(prj);
  } else {
    // Unrecognized projection code.
    status = wcserr_set(WCSERR_SET(PRJERR_BAD_PARAM),
               "Unrecognized projection code '%s'", prj->code);
  }

  return status;
}

//----------------------------------------------------------------------------

int prjx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) < 100) {
    if ((status = prjset(prj))) return status;
  }

  return prj->prjx2s(prj, nx, ny, sxy, spt, x, y, phi, theta, stat);
}

//----------------------------------------------------------------------------

int prjs2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) < 100) {
    if ((status = prjset(prj))) return status;
  }

  return prj->prjs2x(prj, nphi, ntheta, spt, sxy, phi, theta, x, y, stat);
}

/*============================================================================
* Internal helper routine used by the *set() routines - not intended for
* outside use.  It forces (x,y) = (0,0) at (phi0,theta0).
*---------------------------------------------------------------------------*/

static int prjoff(
  struct prjprm *prj,
  const double phi0,
  const double theta0)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  prj->x0 = 0.0;
  prj->y0 = 0.0;

  if (undefined(prj->phi0) || undefined(prj->theta0)) {
    // Set both to the projection-specific default if either undefined.
    prj->phi0   = phi0;
    prj->theta0 = theta0;

  } else {
    double x0, y0;
    int    stat;
    if (prj->prjs2x(prj, 1, 1, 1, 1, &(prj->phi0), &(prj->theta0), &x0, &y0,
                    &stat)) {
      return PRJERR_BAD_PARAM_SET("prjoff");
    }

    prj->x0 = x0;
    prj->y0 = y0;
  }

  return 0;
}

/*============================================================================
*   AZP: zenithal/azimuthal perspective projection.
*
*   Given:
*      prj->pv[1]   Distance parameter, mu in units of r0.
*      prj->pv[2]   Tilt angle, gamma in degrees.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to  0.0 if undefined.
*      prj->theta0  Reset to 90.0 if undefined.
*
*   Returned:
*      prj->flag     AZP
*      prj->code    "AZP"
*      prj->x0      Offset in x.
*      prj->y0      Offset in y.
*      prj->w[0]    r0*(mu+1)
*      prj->w[1]    tan(gamma)
*      prj->w[2]    sec(gamma)
*      prj->w[3]    cos(gamma)
*      prj->w[4]    sin(gamma)
*      prj->w[5]    asin(-1/mu) for |mu| >= 1, -90 otherwise
*      prj->w[6]    mu*cos(gamma)
*      prj->w[7]    1 if |mu*cos(gamma)| < 1, 0 otherwise
*      prj->prjx2s  Pointer to azpx2s().
*      prj->prjs2x  Pointer to azps2x().
*===========================================================================*/

int azpset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -AZP) return 0;

  strcpy(prj->code, "AZP");

  if (undefined(prj->pv[1])) prj->pv[1] = 0.0;
  if (undefined(prj->pv[2])) prj->pv[2] = 0.0;
  if (prj->r0 == 0.0) prj->r0 = R2D;

  strcpy(prj->name, "zenithal/azimuthal perspective");
  prj->category  = ZENITHAL;
  prj->pvrange   = 102;
  prj->simplezen = prj->pv[2] == 0.0;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 0;
  prj->divergent = prj->pv[1] <= 1.0;

  prj->w[0] = prj->r0*(prj->pv[1] + 1.0);
  if (prj->w[0] == 0.0) {
    return PRJERR_BAD_PARAM_SET("azpset");
  }

  prj->w[3] = cosd(prj->pv[2]);
  if (prj->w[3] == 0.0) {
    return PRJERR_BAD_PARAM_SET("azpset");
  }

  prj->w[2] = 1.0/prj->w[3];
  prj->w[4] = sind(prj->pv[2]);
  prj->w[1] = prj->w[4] / prj->w[3];

  if (fabs(prj->pv[1]) > 1.0) {
    prj->w[5] = asind(-1.0/prj->pv[1]);
  } else {
    prj->w[5] = -90.0;
  }

  prj->w[6] = prj->pv[1] * prj->w[3];
  prj->w[7] = (fabs(prj->w[6]) < 1.0) ? 1.0 : 0.0;

  prj->prjx2s = azpx2s;
  prj->prjs2x = azps2x;

  prj->flag = (prj->flag == 1) ? -AZP : AZP;

  return prjoff(prj, 0.0, 90.0);
}

//----------------------------------------------------------------------------

int azpx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-13;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != AZP) {
    if ((status = azpset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  double *phip, *thetap;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;

    phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xj;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  int *statp  = stat;

  phip   = phi;
  thetap = theta;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yj = *yp + prj->y0;

    double yc  = yj*prj->w[3];
    double yc2 = yc*yc;

    double q = prj->w[0] + yj*prj->w[4];

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xj = *phip;

      double r = sqrt(xj*xj + yc2);
      if (r == 0.0) {
        *phip = 0.0;
        *thetap = 90.0;
        *statp  = 0;

      } else {
        *phip = atan2d(xj, -yc);

        double s = r / q;
        double t = s*prj->pv[1]/sqrt(s*s + 1.0);

        s = atan2d(1.0, s);

        if (fabs(t) > 1.0) {
          if (fabs(t) > 1.0+tol) {
            *thetap = 0.0;
            *statp  = 1;
            if (!status) status = PRJERR_BAD_PIX_SET("azpx2s");
            continue;
          }
          t = copysign(90.0, t);
        } else {
          t = asind(t);
        }

        double a = s - t;
        double b = s + t + 180.0;

        if (a > 90.0) a -= 360.0;
        if (b > 90.0) b -= 360.0;

        *thetap = (a > b) ? a : b;
        *statp  = 0;
      }
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("azpx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int azps2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != AZP) {
    if ((status = azpset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double sinphi, cosphi;
    sincosd(*phip, &sinphi, &cosphi);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinphi;
      *yp = cosphi;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double sinthe, costhe;
    sincosd(*thetap, &sinthe, &costhe);

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      double s = prj->w[1]*(*yp);
      double t = (prj->pv[1] + sinthe) + costhe*s;

      if (t == 0.0) {
        *xp = 0.0;
        *yp = 0.0;
        *statp = 1;
        if (!status) status = PRJERR_BAD_WORLD_SET("azps2x");

      } else {
        double r = prj->w[0]*costhe/t;

        // Bounds checking.
        int istat = 0;
        if (prj->bounds&1) {
          if (*thetap < prj->w[5]) {
            // Overlap.
            istat = 1;
            if (!status) status = PRJERR_BAD_WORLD_SET("azps2x");

          } else if (prj->w[7] > 0.0) {
            // Divergence.
            t = prj->pv[1] / sqrt(1.0 + s*s);

            if (fabs(t) <= 1.0) {
              s = atand(-s);
              t = asind(t);

              double a = s - t;
              double b = s + t + 180.0;

              if (a > 90.0) a -= 360.0;
              if (b > 90.0) b -= 360.0;

              if (*thetap < ((a > b) ? a : b)) {
                istat = 1;
                if (!status) status = PRJERR_BAD_WORLD_SET("azps2x");
              }
            }
          }
        }

        *xp =  r*(*xp) - prj->x0;
        *yp = -r*(*yp)*prj->w[2] - prj->y0;
        *statp = istat;
      }
    }
  }

  return status;
}

/*============================================================================
*   SZP: slant zenithal perspective projection.
*
*   Given:
*      prj->pv[1]   Distance of the point of projection from the centre of the
*                   generating sphere, mu in units of r0.
*      prj->pv[2]   Native longitude, phi_c, and ...
*      prj->pv[3]   Native latitude, theta_c, on the planewards side of the
*                   intersection of the line through the point of projection
*                   and the centre of the generating sphere, phi_c in degrees.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to  0.0 if undefined.
*      prj->theta0  Reset to 90.0 if undefined.
*
*   Returned:
*      prj->flag     SZP
*      prj->code    "SZP"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    1/r0
*      prj->w[1]    xp = -mu*cos(theta_c)*sin(phi_c)
*      prj->w[2]    yp =  mu*cos(theta_c)*cos(phi_c)
*      prj->w[3]    zp =  mu*sin(theta_c) + 1
*      prj->w[4]    r0*xp
*      prj->w[5]    r0*yp
*      prj->w[6]    r0*zp
*      prj->w[7]    (zp - 1)^2
*      prj->w[8]    asin(1-zp) if |1 - zp| < 1, -90 otherwise
*      prj->prjx2s  Pointer to szpx2s().
*      prj->prjs2x  Pointer to szps2x().
*===========================================================================*/

int szpset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -SZP) return 0;

  strcpy(prj->code, "SZP");

  if (undefined(prj->pv[1])) prj->pv[1] =  0.0;
  if (undefined(prj->pv[2])) prj->pv[2] =  0.0;
  if (undefined(prj->pv[3])) prj->pv[3] = 90.0;
  if (prj->r0 == 0.0) prj->r0 = R2D;

  strcpy(prj->name, "slant zenithal perspective");
  prj->category  = ZENITHAL;
  prj->pvrange   = 103;
  prj->simplezen = prj->pv[3] == 90.0;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 0;
  prj->divergent = prj->pv[1] <= 1.0;

  prj->w[0] = 1.0/prj->r0;

  prj->w[3] = prj->pv[1] * sind(prj->pv[3]) + 1.0;
  if (prj->w[3] == 0.0) {
    return PRJERR_BAD_PARAM_SET("szpset");
  }

  prj->w[1] = -prj->pv[1] * cosd(prj->pv[3]) * sind(prj->pv[2]);
  prj->w[2] =  prj->pv[1] * cosd(prj->pv[3]) * cosd(prj->pv[2]);
  prj->w[4] =  prj->r0 * prj->w[1];
  prj->w[5] =  prj->r0 * prj->w[2];
  prj->w[6] =  prj->r0 * prj->w[3];
  prj->w[7] =  (prj->w[3] - 1.0) * prj->w[3] - 1.0;

  if (fabs(prj->w[3] - 1.0) < 1.0) {
    prj->w[8] = asind(1.0 - prj->w[3]);
  } else {
    prj->w[8] = -90.0;
  }

  prj->prjx2s = szpx2s;
  prj->prjs2x = szps2x;

  prj->flag = (prj->flag == 1) ? -SZP : SZP;

  return prjoff(prj, 0.0, 90.0);
}

//----------------------------------------------------------------------------

int szpx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-13;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != SZP) {
    if ((status = szpset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xr = (*xp + prj->x0)*prj->w[0];

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xr;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int     *statp = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yr = (*yp + prj->y0)*prj->w[0];

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xr = *phip;
      double r2 = xr*xr + yr*yr;

      double x1 = (xr - prj->w[1])/prj->w[3];
      double y1 = (yr - prj->w[2])/prj->w[3];
      double xy = xr*x1 + yr*y1;

      double z;
      if (r2 < 1.0e-10) {
        // Use small angle formula.
        z = r2/2.0;
        *thetap = 90.0 - R2D*sqrt(r2/(1.0 + xy));

      } else {
        double t = x1*x1 + y1*y1;
        double a = t + 1.0;
        double b = xy - t;
        double c = r2 - xy - xy + t - 1.0;
        double d = b*b - a*c;

        // Check for a solution.
        if (d < 0.0) {
          *phip = 0.0;
          *thetap = 0.0;
          *statp  = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("szpx2s");
          continue;
        }
        d = sqrt(d);

        // Choose solution closest to pole.
        double sinth1 = (-b + d)/a;
        double sinth2 = (-b - d)/a;
        double sinthe = (sinth1 > sinth2) ? sinth1 : sinth2;
        if (sinthe > 1.0) {
          if (sinthe-1.0 < tol) {
            sinthe = 1.0;
          } else {
            sinthe = (sinth1 < sinth2) ? sinth1 : sinth2;
          }
        }

        if (sinthe < -1.0) {
          if (sinthe+1.0 > -tol) {
            sinthe = -1.0;
          }
        }

        if (sinthe > 1.0 || sinthe < -1.0) {
          *phip   = 0.0;
          *thetap = 0.0;
          *statp  = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("szpx2s");
          continue;
        }

        *thetap = asind(sinthe);

        z = 1.0 - sinthe;
      }

      *phip = atan2d(xr - x1*z, -(yr - y1*z));
      *statp = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("szpx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int szps2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != SZP) {
    if ((status = szpset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double sinphi, cosphi;
    sincosd(*phip, &sinphi, &cosphi);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinphi;
      *yp = cosphi;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double s = 1.0 - sind(*thetap);
    double t = prj->w[3] - s;

    if (t == 0.0) {
      for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
        *xp = 0.0;
        *yp = 0.0;
        *statp = 1;
      }

      if (!status) status = PRJERR_BAD_WORLD_SET("szps2x");

    } else {
      double r = prj->w[6]*cosd(*thetap)/t;
      double u = prj->w[4]*s/t + prj->x0;
      double v = prj->w[5]*s/t + prj->y0;

      for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
        // Bounds checking.
        int istat = 0;
        if (prj->bounds&1) {
          if (*thetap < prj->w[8]) {
            // Divergence.
            istat = 1;
            if (!status) status = PRJERR_BAD_WORLD_SET("szps2x");

          } else if (fabs(prj->pv[1]) > 1.0) {
            // Overlap.
            s = prj->w[1]*(*xp) - prj->w[2]*(*yp);
            t = 1.0/sqrt(prj->w[7] + s*s);

            if (fabs(t) <= 1.0) {
              s = atan2d(s, prj->w[3] - 1.0);
              t = asind(t);

              double a = s - t;
              double b = s + t + 180.0;

              if (a > 90.0) a -= 360.0;
              if (b > 90.0) b -= 360.0;

              if (*thetap < ((a > b) ? a : b)) {
                istat = 1;
                if (!status) status = PRJERR_BAD_WORLD_SET("szps2x");
              }
            }
          }
        }

        *xp =  r*(*xp) - u;
        *yp = -r*(*yp) - v;
        *statp = istat;
      }
    }
  }

  return status;
}


/*============================================================================
*   TAN: gnomonic projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to  0.0 if undefined.
*      prj->theta0  Reset to 90.0 if undefined.
*
*   Returned:
*      prj->flag     TAN
*      prj->code    "TAN"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->prjx2s  Pointer to tanx2s().
*      prj->prjs2x  Pointer to tans2x().
*===========================================================================*/

int tanset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -TAN) return 0;

  strcpy(prj->code, "TAN");

  if (prj->r0 == 0.0) prj->r0 = R2D;

  strcpy(prj->name, "gnomonic");
  prj->category  = ZENITHAL;
  prj->pvrange   = 0;
  prj->simplezen = 1;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 0;
  prj->divergent = 1;

  prj->prjx2s = tanx2s;
  prj->prjs2x = tans2x;

  prj->flag = (prj->flag == 1) ? -TAN : TAN;

  return prjoff(prj, 0.0, 90.0);
}

//----------------------------------------------------------------------------

int tanx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != TAN) {
    if ((status = tanset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xj;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yj  = *yp + prj->y0;
    double yj2 = yj*yj;

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xj = *phip;

      double r = sqrt(xj*xj + yj2);
      if (r == 0.0) {
        *phip = 0.0;
      } else {
        *phip = atan2d(xj, -yj);
      }

      *thetap = atan2d(prj->r0, r);
      *statp  = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("tanx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int tans2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != TAN) {
    if ((status = tanset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double sinphi, cosphi;
    sincosd(*phip, &sinphi, &cosphi);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinphi;
      *yp = cosphi;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double s = sind(*thetap);
    if (s == 0.0) {
      for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
        *xp = 0.0;
        *yp = 0.0;
        *statp = 1;
      }
      if (!status) status = PRJERR_BAD_WORLD_SET("tans2x");

    } else {
      double r =  prj->r0*cosd(*thetap)/s;

      // Bounds checking.
      int istat = 0;
      if (prj->bounds&1) {
        if (s < 0.0) {
          istat = 1;
          if (!status) status = PRJERR_BAD_WORLD_SET("tans2x");
        }
      }

      for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
        *xp =  r*(*xp) - prj->x0;
        *yp = -r*(*yp) - prj->y0;
        *statp = istat;
      }
    }
  }

  return status;
}

/*============================================================================
*   STG: stereographic projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to  0.0 if undefined.
*      prj->theta0  Reset to 90.0 if undefined.
*
*   Returned:
*      prj->flag     STG
*      prj->code    "STG"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    2*r0
*      prj->w[1]    1/(2*r0)
*      prj->prjx2s  Pointer to stgx2s().
*      prj->prjs2x  Pointer to stgs2x().
*===========================================================================*/

int stgset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -STG) return 0;

  strcpy(prj->code, "STG");

  strcpy(prj->name, "stereographic");
  prj->category  = ZENITHAL;
  prj->pvrange   = 0;
  prj->simplezen = 1;
  prj->equiareal = 0;
  prj->conformal = 1;
  prj->global    = 0;
  prj->divergent = 1;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 360.0/PI;
    prj->w[1] = PI/360.0;
  } else {
    prj->w[0] = 2.0*prj->r0;
    prj->w[1] = 1.0/prj->w[0];
  }

  prj->prjx2s = stgx2s;
  prj->prjs2x = stgs2x;

  prj->flag = (prj->flag == 1) ? -STG : STG;

  return prjoff(prj, 0.0, 90.0);
}

//----------------------------------------------------------------------------

int stgx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != STG) {
    if ((status = stgset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xj;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yj  = *yp + prj->y0;
    double yj2 = yj*yj;

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xj  = *phip;

      double r = sqrt(xj*xj + yj2);
      if (r == 0.0) {
        *phip = 0.0;
      } else {
        *phip = atan2d(xj, -yj);
      }

      *thetap = 90.0 - 2.0*atand(r*prj->w[1]);
      *statp  = 0;
    }
  }

  return 0;
}

//----------------------------------------------------------------------------

int stgs2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != STG) {
    if ((status = stgset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double sinphi, cosphi;
    sincosd(*phip, &sinphi, &cosphi);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinphi;
      *yp = cosphi;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double s = 1.0 + sind(*thetap);
    if (s == 0.0) {
      for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
        *xp = 0.0;
        *yp = 0.0;
        *statp = 1;
      }
      if (!status) status = PRJERR_BAD_WORLD_SET("stgs2x");

    } else {
      double r = prj->w[0]*cosd(*thetap)/s;

      for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
        *xp =  r*(*xp) - prj->x0;
        *yp = -r*(*yp) - prj->y0;
        *statp = 0;
      }
    }
  }

  return status;
}

/*============================================================================
*   SIN: orthographic/synthesis projection.
*
*   Given:
*      prj->pv[1:2] Obliqueness parameters, xi and eta.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to  0.0 if undefined.
*      prj->theta0  Reset to 90.0 if undefined.
*
*   Returned:
*      prj->flag     SIN
*      prj->code    "SIN"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    1/r0
*      prj->w[1]    xi**2 + eta**2
*      prj->w[2]    xi**2 + eta**2 + 1
*      prj->w[3]    xi**2 + eta**2 - 1
*      prj->prjx2s  Pointer to sinx2s().
*      prj->prjs2x  Pointer to sins2x().
*===========================================================================*/

int sinset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -SIN) return 0;

  strcpy(prj->code, "SIN");

  if (undefined(prj->pv[1])) prj->pv[1] = 0.0;
  if (undefined(prj->pv[2])) prj->pv[2] = 0.0;
  if (prj->r0 == 0.0) prj->r0 = R2D;

  strcpy(prj->name, "orthographic/synthesis");
  prj->category  = ZENITHAL;
  prj->pvrange   = 102;
  prj->simplezen = (prj->pv[1] == 0.0 && prj->pv[2] == 0.0);
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 0;
  prj->divergent = 0;

  prj->w[0] = 1.0/prj->r0;
  prj->w[1] = prj->pv[1]*prj->pv[1] + prj->pv[2]*prj->pv[2];
  prj->w[2] = prj->w[1] + 1.0;
  prj->w[3] = prj->w[1] - 1.0;

  prj->prjx2s = sinx2s;
  prj->prjs2x = sins2x;

  prj->flag = (prj->flag == 1) ? -SIN : SIN;

  return prjoff(prj, 0.0, 90.0);
}

//----------------------------------------------------------------------------

int sinx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-13;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != SIN) {
    if ((status = sinset(prj))) return status;
  }

  double xi  = prj->pv[1];
  double eta = prj->pv[2];

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double x0 = (*xp + prj->x0)*prj->w[0];

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = x0;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double y0 = (*yp + prj->y0)*prj->w[0];
    double y02 = y0*y0;

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      // Compute intermediaries.
      double x0 = *phip;
      double r2 = x0*x0 + y02;

      if (prj->w[1] == 0.0) {
        // Orthographic projection.
        if (r2 != 0.0) {
          *phip = atan2d(x0, -y0);
        } else {
          *phip = 0.0;
        }

        if (r2 < 0.5) {
          *thetap = acosd(sqrt(r2));
        } else if (r2 <= 1.0) {
          *thetap = asind(sqrt(1.0 - r2));
        } else {
          *statp = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("sinx2s")
          continue;
        }

      } else {
        // "Synthesis" projection.
        double xy = x0*xi + y0*eta;

        double z;
        if (r2 < 1.0e-10) {
          // Use small angle formula.
          z = r2/2.0;
          *thetap = 90.0 - R2D*sqrt(r2/(1.0 + xy));

        } else {
          double a = prj->w[2];
          double b = xy - prj->w[1];
          double c = r2 - xy - xy + prj->w[3];
          double d = b*b - a*c;

          // Check for a solution.
          if (d < 0.0) {
            *phip = 0.0;
            *thetap = 0.0;
            *statp  = 1;
            if (!status) status = PRJERR_BAD_PIX_SET("sinx2s")
            continue;
          }
          d = sqrt(d);

          // Choose solution closest to pole.
          double sinth1 = (-b + d)/a;
          double sinth2 = (-b - d)/a;
          double sinthe = (sinth1 > sinth2) ? sinth1 : sinth2;
          if (sinthe > 1.0) {
            if (sinthe-1.0 < tol) {
              sinthe = 1.0;
            } else {
              sinthe = (sinth1 < sinth2) ? sinth1 : sinth2;
            }
          }

          if (sinthe < -1.0) {
            if (sinthe+1.0 > -tol) {
              sinthe = -1.0;
            }
          }

          if (sinthe > 1.0 || sinthe < -1.0) {
            *phip = 0.0;
            *thetap = 0.0;
            *statp  = 1;
            if (!status) status = PRJERR_BAD_PIX_SET("sinx2s")
            continue;
          }

          *thetap = asind(sinthe);
          z = 1.0 - sinthe;
        }

        double x1 = -y0 + eta*z;
        double y1 =  x0 -  xi*z;
        if (x1 == 0.0 && y1 == 0.0) {
          *phip = 0.0;
        } else {
          *phip = atan2d(y1,x1);
        }
      }

      *statp = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("sinx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int sins2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != SIN) {
    if ((status = sinset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double sinphi, cosphi;
    sincosd(*phip, &sinphi, &cosphi);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinphi;
      *yp = cosphi;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double costhe, z;

    double t = (90.0 - fabs(*thetap))*D2R;
    if (t < 1.0e-5) {
      if (*thetap > 0.0) {
         z = t*t/2.0;
      } else {
         z = 2.0 - t*t/2.0;
      }
      costhe = t;
    } else {
      z = 1.0 - sind(*thetap);
      costhe = cosd(*thetap);
    }
    double r = prj->r0*costhe;

    if (prj->w[1] == 0.0) {
      // Orthographic projection.
      int istat = 0;
      if (prj->bounds&1) {
        if (*thetap < 0.0) {
          istat = 1;
          if (!status) status = PRJERR_BAD_WORLD_SET("sins2x");
        }
      }

      for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
        *xp =  r*(*xp) - prj->x0;
        *yp = -r*(*yp) - prj->y0;
        *statp = istat;
      }

    } else {
      // "Synthesis" projection.
      z *= prj->r0;
      double z1 = prj->pv[1]*z - prj->x0;
      double z2 = prj->pv[2]*z - prj->y0;

      for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
        int istat = 0;
        if (prj->bounds&1) {
          t = -atand(prj->pv[1]*(*xp) - prj->pv[2]*(*yp));
          if (*thetap < t) {
            istat = 1;
            if (!status) status = PRJERR_BAD_WORLD_SET("sins2x");
          }
        }

        *xp =  r*(*xp) + z1;
        *yp = -r*(*yp) + z2;
        *statp = istat;
      }
    }
  }

  return status;
}

/*============================================================================
*   ARC: zenithal/azimuthal equidistant projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to  0.0 if undefined.
*      prj->theta0  Reset to 90.0 if undefined.
*
*   Returned:
*      prj->flag     ARC
*      prj->code    "ARC"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    r0*(pi/180)
*      prj->w[1]    (180/pi)/r0
*      prj->prjx2s  Pointer to arcx2s().
*      prj->prjs2x  Pointer to arcs2x().
*===========================================================================*/

int arcset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -ARC) return 0;

  strcpy(prj->code, "ARC");

  strcpy(prj->name, "zenithal/azimuthal equidistant");
  prj->category  = ZENITHAL;
  prj->pvrange   = 0;
  prj->simplezen = 1;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 1.0;
    prj->w[1] = 1.0;
  } else {
    prj->w[0] = prj->r0*D2R;
    prj->w[1] = 1.0/prj->w[0];
  }

  prj->prjx2s = arcx2s;
  prj->prjs2x = arcs2x;

  prj->flag = (prj->flag == 1) ? -ARC : ARC;

  return prjoff(prj, 0.0, 90.0);
}

//----------------------------------------------------------------------------

int arcx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != ARC) {
    if ((status = arcset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xj;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yj  = *yp + prj->y0;
    double yj2 = yj*yj;

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xj = *phip;

      double r = sqrt(xj*xj + yj2);
      if (r == 0.0) {
        *phip = 0.0;
        *thetap = 90.0;
      } else {
        *phip = atan2d(xj, -yj);
        *thetap = 90.0 - r*prj->w[1];
      }

      *statp = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("arcx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int arcs2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != ARC) {
    if ((status = arcset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double sinphi, cosphi;
    sincosd(*phip, &sinphi, &cosphi);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinphi;
      *yp = cosphi;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double r =  prj->w[0]*(90.0 - *thetap);

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      *xp =  r*(*xp) - prj->x0;
      *yp = -r*(*yp) - prj->y0;
      *statp = 0;
    }
  }

  return 0;
}

/*============================================================================
*   ZPN: zenithal/azimuthal polynomial projection.
*
*   Given:
*      prj->pv[]    Polynomial coefficients.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to  0.0 if undefined.
*      prj->theta0  Reset to 90.0 if undefined.
*
*   Returned:
*      prj->flag     ZPN
*      prj->code    "ZPN"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->n       Degree of the polynomial, N.
*      prj->w[0]    Co-latitude of the first point of inflection, radian.
*      prj->w[1]    Radius of the first point of inflection (N > 1), radian.
*      prj->prjx2s  Pointer to zpnx2s().
*      prj->prjs2x  Pointer to zpns2x().
*===========================================================================*/

int zpnset(struct prjprm *prj)

{
  const double tol = 1.0e-13;

  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -ZPN) return 0;

  strcpy(prj->code, "ZPN");

  if (undefined(prj->pv[1])) prj->pv[1] = 0.0;
  if (undefined(prj->pv[2])) prj->pv[2] = 0.0;
  if (undefined(prj->pv[3])) prj->pv[3] = 0.0;
  if (prj->r0 == 0.0) prj->r0 = R2D;

  strcpy(prj->name, "zenithal/azimuthal polynomial");
  prj->category  = ZENITHAL;
  prj->pvrange   = 30;
  prj->simplezen = 1;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 0;
  prj->divergent = 0;

  // Find the highest non-zero coefficient.
  int k;
  for (k = PVN-1; k >= 0 && prj->pv[k] == 0.0; k--);
  if (k < 0) {
    return PRJERR_BAD_PARAM_SET("zpnset");
  }

  prj->n = k;

  if (k < 2) {
    // No point of inflection.
    prj->w[0] = PI;

  } else {
    // Find the point of inflection closest to the pole.
    double d, d1, d2;

    d1 = prj->pv[1];
    if (d1 <= 0.0) {
      return PRJERR_BAD_PARAM_SET("zpnset");
    }

    // Find the point where the derivative first goes negative.
    int j;
    double zd, zd1, zd2;

    zd1 = 0.0;
    for (j = 0; j < 180; j++) {
      zd2 = j*D2R;
      d2  = 0.0;
      for (int m = k; m > 0; m--) {
        d2 = d2*zd2 + m*prj->pv[m];
      }

      if (d2 <= 0.0) break;
      zd1 = zd2;
      d1  = d2;
    }

    if (j == 180) {
      // No negative derivative -> no point of inflection.
      zd = PI;
      prj->global = 1;
    } else {
      // Find where the derivative is zero.
      for (j = 1; j <= 10; j++) {
        zd = zd1 - d1*(zd2-zd1)/(d2-d1);

        d = 0.0;
        for (int m = k; m > 0; m--) {
          d = d*zd + m*prj->pv[m];
        }

        if (fabs(d) < tol) break;

        if (d < 0.0) {
          zd2 = zd;
          d2  = d;
        } else {
          zd1 = zd;
          d1  = d;
        }
      }
    }

    double r = 0.0;
    for (int m = k; m >= 0; m--) {
      r = r*zd + prj->pv[m];
    }
    prj->w[0] = zd;
    prj->w[1] = r;
  }

  prj->prjx2s = zpnx2s;
  prj->prjs2x = zpns2x;

  prj->flag = (prj->flag == 1) ? -ZPN : ZPN;

  return prjoff(prj, 0.0, 90.0);
}

//----------------------------------------------------------------------------

int zpnx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-13;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != ZPN) {
    if ((status = zpnset(prj))) return status;
  }

  int k = prj->n;

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xj;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yj  = *yp + prj->y0;
    double yj2 = yj*yj;

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xj = *phip;

      double r = sqrt(xj*xj + yj2)/prj->r0;
      if (r == 0.0) {
        *phip = 0.0;
      } else {
        *phip = atan2d(xj, -yj);
      }

      double zd;
      if (k < 1) {
        // Constant - no solution.
        return PRJERR_BAD_PARAM_SET("zpnx2s");

      } else if (k == 1) {
        // Linear.
        zd = (r - prj->pv[0])/prj->pv[1];

      } else if (k == 2) {
        // Quadratic.
        double a = prj->pv[2];
        double b = prj->pv[1];
        double c = prj->pv[0] - r;

        double d = b*b - 4.0*a*c;
        if (d < 0.0) {
          *thetap = 0.0;
          *statp  = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("zpnx2s");
          continue;
        }
        d = sqrt(d);

        // Choose solution closest to pole.
        double zd1 = (-b + d)/(2.0*a);
        double zd2 = (-b - d)/(2.0*a);

        zd  = (zd1<zd2) ? zd1 : zd2;
        if (zd < -tol) zd = (zd1>zd2) ? zd1 : zd2;
        if (zd < 0.0) {
          if (zd < -tol) {
            *thetap = 0.0;
            *statp  = 1;
            if (!status) status = PRJERR_BAD_PIX_SET("zpnx2s");
            continue;
          }
          zd = 0.0;
        } else if (zd > PI) {
          if (zd > PI+tol) {
            *thetap = 0.0;
            *statp  = 1;
            if (!status) status = PRJERR_BAD_PIX_SET("zpnx2s");
            continue;
          }
          zd = PI;
        }
      } else {
        // Higher order - solve iteratively.
        double zd1 = 0.0;
        double r1  = prj->pv[0];
        double zd2 = prj->w[0];
        double r2  = prj->w[1];

        if (r < r1) {
          if (r < r1-tol) {
            *thetap = 0.0;
            *statp  = 1;
            if (!status) status = PRJERR_BAD_PIX_SET("zpnx2s");
            continue;
          }
          zd = zd1;
        } else if (r > r2) {
          if (r > r2+tol) {
            *thetap = 0.0;
            *statp  = 1;
            if (!status) status = PRJERR_BAD_PIX_SET("zpnx2s");
            continue;
          }
          zd = zd2;
        } else {
          // Dissect the interval.
          for (int j = 0; j < 100; j++) {
            double lambda = (r2 - r)/(r2 - r1);
            if (lambda < 0.1) {
              lambda = 0.1;
            } else if (lambda > 0.9) {
              lambda = 0.9;
            }

            zd = zd2 - lambda*(zd2 - zd1);

            double rt = 0.0;
            for (int m = k; m >= 0; m--) {
              rt = (rt * zd) + prj->pv[m];
            }

            if (rt < r) {
              if (r-rt < tol) break;
              r1 = rt;
              zd1 = zd;
            } else {
              if (rt-r < tol) break;
              r2 = rt;
              zd2 = zd;
            }

            if (fabs(zd2-zd1) < tol) break;
          }
        }
      }

      *thetap = 90.0 - zd*R2D;
      *statp  = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("zpnx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int zpns2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != ZPN) {
    if ((status = zpnset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double sinphi, cosphi;
    sincosd(*phip, &sinphi, &cosphi);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinphi;
      *yp = cosphi;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double s = (90.0 - *thetap)*D2R;

    double r = 0.0;
    for (int m = prj->n; m >= 0; m--) {
      r = r*s + prj->pv[m];
    }
    r *= prj->r0;

    // Bounds checking.
    int istat = 0;
    if (prj->bounds&1) {
      if (s > prj->w[0]) {
        istat = 1;
        if (!status) status = PRJERR_BAD_WORLD_SET("zpns2x");
      }
    }

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      *xp =  r*(*xp) - prj->x0;
      *yp = -r*(*yp) - prj->y0;
      *statp = istat;
    }
  }

  return status;
}

/*============================================================================
*   ZEA: zenithal/azimuthal equal area projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to  0.0 if undefined.
*      prj->theta0  Reset to 90.0 if undefined.
*
*   Returned:
*      prj->flag     ZEA
*      prj->code    "ZEA"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    2*r0
*      prj->w[1]    1/(2*r0)
*      prj->prjx2s  Pointer to zeax2s().
*      prj->prjs2x  Pointer to zeas2x().
*===========================================================================*/

int zeaset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -ZEA) return 0;

  strcpy(prj->code, "ZEA");

  strcpy(prj->name, "zenithal/azimuthal equal area");
  prj->category  = ZENITHAL;
  prj->pvrange   = 0;
  prj->simplezen = 1;
  prj->equiareal = 1;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 360.0/PI;
    prj->w[1] = PI/360.0;
  } else {
    prj->w[0] = 2.0*prj->r0;
    prj->w[1] = 1.0/prj->w[0];
  }

  prj->prjx2s = zeax2s;
  prj->prjs2x = zeas2x;

  prj->flag = (prj->flag == 1) ? -ZEA : ZEA;

  return prjoff(prj, 0.0, 90.0);
}

//----------------------------------------------------------------------------

int zeax2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-12;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != ZEA) {
    if ((status = zeaset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xj;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yj  = *yp + prj->y0;
    double yj2 = yj*yj;

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xj  = *phip;

      double r = sqrt(xj*xj + yj2);
      if (r == 0.0) {
        *phip = 0.0;
      } else {
        *phip = atan2d(xj, -yj);
      }

      double s = r*prj->w[1];
      if (fabs(s) > 1.0) {
        if (fabs(r - prj->w[0]) < tol) {
          *thetap = -90.0;
        } else {
          *thetap = 0.0;
          *statp  = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("zeax2s");
          continue;
        }
      } else {
        *thetap = 90.0 - 2.0*asind(s);
      }

      *statp = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("zeax2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int zeas2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != ZEA) {
    if ((status = zeaset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double sinphi, cosphi;
    sincosd(*phip, &sinphi, &cosphi);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinphi;
      *yp = cosphi;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double r =  prj->w[0]*sind((90.0 - *thetap)/2.0);

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      *xp =  r*(*xp) - prj->x0;
      *yp = -r*(*yp) - prj->y0;
      *statp = 0;
    }
  }

  return 0;
}

/*============================================================================
*   AIR: Airy's projection.
*
*   Given:
*      prj->pv[1]   Latitude theta_b within which the error is minimized, in
*                   degrees.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to  0.0 if undefined.
*      prj->theta0  Reset to 90.0 if undefined.
*
*   Returned:
*      prj->flag     AIR
*      prj->code    "AIR"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    2*r0
*      prj->w[1]    ln(cos(xi_b))/tan(xi_b)**2, where xi_b = (90-theta_b)/2
*      prj->w[2]    1/2 - prj->w[1]
*      prj->w[3]    2*r0*prj->w[2]
*      prj->w[4]    tol, cutoff for using small angle approximation, in
*                   radians.
*      prj->w[5]    prj->w[2]*tol
*      prj->w[6]    (180/pi)/prj->w[2]
*      prj->prjx2s  Pointer to airx2s().
*      prj->prjs2x  Pointer to airs2x().
*===========================================================================*/

int airset(struct prjprm *prj)

{
  const double tol = 1.0e-4;

  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -AIR) return 0;

  strcpy(prj->code, "AIR");

  if (undefined(prj->pv[1])) prj->pv[1] = 90.0;
  if (prj->r0 == 0.0) prj->r0 = R2D;

  strcpy(prj->name, "Airy's zenithal");
  prj->category  = ZENITHAL;
  prj->pvrange   = 101;
  prj->simplezen = 1;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 0;
  prj->divergent = 1;

  prj->w[0] = 2.0*prj->r0;
  if (prj->pv[1] == 90.0) {
    prj->w[1] = -0.5;
    prj->w[2] =  1.0;
  } else if (prj->pv[1] > -90.0) {
    double cosxi = cosd((90.0 - prj->pv[1])/2.0);
    prj->w[1] = log(cosxi)*(cosxi*cosxi)/(1.0-cosxi*cosxi);
    prj->w[2] = 0.5 - prj->w[1];
  } else {
    return PRJERR_BAD_PARAM_SET("airset");
  }

  prj->w[3] = prj->w[0] * prj->w[2];
  prj->w[4] = tol;
  prj->w[5] = prj->w[2]*tol;
  prj->w[6] = R2D/prj->w[2];

  prj->prjx2s = airx2s;
  prj->prjs2x = airs2x;

  prj->flag = (prj->flag == 1) ? -AIR : AIR;

  return prjoff(prj, 0.0, 90.0);
}

//----------------------------------------------------------------------------

int airx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-12;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != AIR) {
    if ((status = airset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xj;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yj  = *yp + prj->y0;
    double yj2 = yj*yj;

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xj = *phip;

      double r = sqrt(xj*xj + yj2)/prj->w[0];
      if (r == 0.0) {
        *phip = 0.0;
      } else {
        *phip = atan2d(xj, -yj);
      }

      double xi;
      if (r == 0.0) {
        xi = 0.0;
      } else if (r < prj->w[5]) {
        xi = r*prj->w[6];
      } else {
        // Find a solution interval.
        double x1 = 1.0;
        double x2 = 1.0;
        double r1 = 0.0;
        double r2 = 0.0;

        int k;
        for (k = 0; k < 30; k++) {
          x2 = x1/2.0;
          double tanxi = sqrt(1.0-x2*x2)/x2;
          r2 = -(log(x2)/tanxi + prj->w[1]*tanxi);

          if (r2 >= r) break;
          x1 = x2;
          r1 = r2;
        }
        if (k == 30) {
          *thetap = 0.0;
          *statp  = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("airx2s");
          continue;
        }

        double cosxi;
        for (k = 0; k < 100; k++) {
          // Weighted division of the interval.
          double lambda = (r2-r)/(r2-r1);
          if (lambda < 0.1) {
            lambda = 0.1;
          } else if (lambda > 0.9) {
            lambda = 0.9;
          }
          cosxi = x2 - lambda*(x2-x1);

          double tanxi = sqrt(1.0-cosxi*cosxi)/cosxi;
          double rt = -(log(cosxi)/tanxi + prj->w[1]*tanxi);

          if (rt < r) {
            if (r-rt < tol) break;
            r1 = rt;
            x1 = cosxi;
          } else {
            if (rt-r < tol) break;
            r2 = rt;
            x2 = cosxi;
          }
        }
        if (k == 100) {
          *thetap = 0.0;
          *statp  = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("airx2s");
          continue;
        }

        xi = acosd(cosxi);
      }

      *thetap = 90.0 - 2.0*xi;
      *statp  = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("airx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int airs2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != AIR) {
    if ((status = airset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double sinphi, cosphi;
    sincosd(*phip, &sinphi, &cosphi);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinphi;
      *yp = cosphi;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    int istat = 0;

    double r;
    if (*thetap == 90.0) {
      r = 0.0;
    } else if (*thetap > -90.0) {
      double xi = D2R*(90.0 - *thetap)/2.0;
      if (xi < prj->w[4]) {
        r = xi*prj->w[3];
      } else {
        double cosxi = cosd((90.0 - *thetap)/2.0);
        double tanxi = sqrt(1.0 - cosxi*cosxi)/cosxi;
        r = -prj->w[0]*(log(cosxi)/tanxi + prj->w[1]*tanxi);
      }
    } else {
      r = 0.0;
      istat = 1;
      if (!status) status = PRJERR_BAD_WORLD_SET("airs2x");
    }

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      *xp =  r*(*xp) - prj->x0;
      *yp = -r*(*yp) - prj->y0;
      *statp = istat;
    }
  }

  return status;
}

/*============================================================================
*   CYP: cylindrical perspective projection.
*
*   Given:
*      prj->pv[1]   Distance of point of projection from the centre of the
*                   generating sphere, mu, in units of r0.
*      prj->pv[2]   Radius of the cylinder of projection, lambda, in units of
*                   r0.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     CYP
*      prj->code    "CYP"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    r0*lambda*(pi/180)
*      prj->w[1]    (180/pi)/(r0*lambda)
*      prj->w[2]    r0*(mu + lambda)
*      prj->w[3]    1/(r0*(mu + lambda))
*      prj->prjx2s  Pointer to cypx2s().
*      prj->prjs2x  Pointer to cyps2x().
*===========================================================================*/

int cypset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -CYP) return 0;

  strcpy(prj->code, "CYP");

  if (undefined(prj->pv[1])) prj->pv[1] = 1.0;
  if (undefined(prj->pv[2])) prj->pv[2] = 1.0;

  strcpy(prj->name, "cylindrical perspective");
  prj->category  = CYLINDRICAL;
  prj->pvrange   = 102;
  prj->simplezen = 0;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = prj->pv[1] < -1.0 || 0.0 < prj->pv[1];
  prj->divergent = !prj->global;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;

    prj->w[0] = prj->pv[2];
    if (prj->w[0] == 0.0) {
      return PRJERR_BAD_PARAM_SET("cypset");
    }

    prj->w[1] = 1.0/prj->w[0];

    prj->w[2] = R2D*(prj->pv[1] + prj->pv[2]);
    if (prj->w[2] == 0.0) {
      return PRJERR_BAD_PARAM_SET("cypset");
    }

    prj->w[3] = 1.0/prj->w[2];
  } else {
    prj->w[0] = prj->r0*prj->pv[2]*D2R;
    if (prj->w[0] == 0.0) {
      return PRJERR_BAD_PARAM_SET("cypset");
    }

    prj->w[1] = 1.0/prj->w[0];

    prj->w[2] = prj->r0*(prj->pv[1] + prj->pv[2]);
    if (prj->w[2] == 0.0) {
      return PRJERR_BAD_PARAM_SET("cypset");
    }

    prj->w[3] = 1.0/prj->w[2];
  }

  prj->prjx2s = cypx2s;
  prj->prjs2x = cyps2x;

  prj->flag = (prj->flag == 1) ? -CYP : CYP;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int cypx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != CYP) {
    if ((status = cypset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double s = prj->w[1]*(*xp + prj->x0);

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = s;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double eta = prj->w[3]*(*yp + prj->y0);
    double t = atan2d(eta,1.0) + asind(eta*prj->pv[1]/sqrt(eta*eta+1.0));

    for (int ix = 0; ix < mx; ix++, thetap += spt, statp++) {
      *thetap = t;
      *statp  = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("cypx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int cyps2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != CYP) {
    if ((status = cypset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double xi = prj->w[0]*(*phip) - prj->x0;

    double *xp = x + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = xi;
      xp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double eta = prj->pv[1] + cosd(*thetap);

    int istat = 0;
    if (eta == 0.0) {
      istat = 1;
      if (!status) status = PRJERR_BAD_WORLD_SET("cyps2x");

    } else {
      eta = prj->w[2]*sind(*thetap)/eta;
    }

    eta -= prj->y0;
    for (int iphi = 0; iphi < mphi; iphi++, yp += sxy, statp++) {
      *yp = eta;
      *statp = istat;
    }
  }

  return status;
}

/*============================================================================
*   CEA: cylindrical equal area projection.
*
*   Given:
*      prj->pv[1]   Square of the cosine of the latitude at which the
*                   projection is conformal, lambda.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     CEA
*      prj->code    "CEA"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    r0*(pi/180)
*      prj->w[1]    (180/pi)/r0
*      prj->w[2]    r0/lambda
*      prj->w[3]    lambda/r0
*      prj->prjx2s  Pointer to ceax2s().
*      prj->prjs2x  Pointer to ceas2x().
*===========================================================================*/

int ceaset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -CEA) return 0;

  strcpy(prj->code, "CEA");

  if (undefined(prj->pv[1])) prj->pv[1] = 1.0;

  strcpy(prj->name, "cylindrical equal area");
  prj->category  = CYLINDRICAL;
  prj->pvrange   = 101;
  prj->simplezen = 0;
  prj->equiareal = 1;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 1.0;
    prj->w[1] = 1.0;
    if (prj->pv[1] <= 0.0 || prj->pv[1] > 1.0) {
      return PRJERR_BAD_PARAM_SET("ceaset");
    }
    prj->w[2] = prj->r0/prj->pv[1];
    prj->w[3] = prj->pv[1]/prj->r0;
  } else {
    prj->w[0] = prj->r0*D2R;
    prj->w[1] = R2D/prj->r0;
    if (prj->pv[1] <= 0.0 || prj->pv[1] > 1.0) {
      return PRJERR_BAD_PARAM_SET("ceaset");
    }
    prj->w[2] = prj->r0/prj->pv[1];
    prj->w[3] = prj->pv[1]/prj->r0;
  }

  prj->prjx2s = ceax2s;
  prj->prjs2x = ceas2x;

  prj->flag = (prj->flag == 1) ? -CEA : CEA;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int ceax2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-13;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != CEA) {
    if ((status = ceaset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double s = prj->w[1]*(*xp + prj->x0);

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = s;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double s = prj->w[3]*(*yp + prj->y0);

    int istat = 0;
    if (fabs(s) > 1.0) {
      if (fabs(s) > 1.0+tol) {
        s = 0.0;
        istat = 1;
        if (!status) status = PRJERR_BAD_PIX_SET("ceax2s");
      } else {
        s = copysign(90.0, s);
      }
    } else {
      s = asind(s);
    }

    for (int ix = 0; ix < mx; ix++, thetap += spt, statp++) {
      *thetap = s;
      *statp  = istat;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("ceax2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int ceas2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != CEA) {
    if ((status = ceaset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double xi = prj->w[0]*(*phip) - prj->x0;

    double *xp = x + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = xi;
      xp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double eta = prj->w[2]*sind(*thetap) - prj->y0;

    for (int iphi = 0; iphi < mphi; iphi++, yp += sxy, statp++) {
      *yp = eta;
      *statp = 0;
    }
  }

  return 0;
}

/*============================================================================
*   CAR: Plate carree projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     CAR
*      prj->code    "CAR"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    r0*(pi/180)
*      prj->w[1]    (180/pi)/r0
*      prj->prjx2s  Pointer to carx2s().
*      prj->prjs2x  Pointer to cars2x().
*===========================================================================*/

int carset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -CAR) return 0;

  strcpy(prj->code, "CAR");

  strcpy(prj->name, "plate caree");
  prj->category  = CYLINDRICAL;
  prj->pvrange   = 0;
  prj->simplezen = 0;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 1.0;
    prj->w[1] = 1.0;
  } else {
    prj->w[0] = prj->r0*D2R;
    prj->w[1] = 1.0/prj->w[0];
  }

  prj->prjx2s = carx2s;
  prj->prjs2x = cars2x;

  prj->flag = (prj->flag == 1) ? -CAR : CAR;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int carx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != CAR) {
    if ((status = carset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double s = prj->w[1]*(*xp + prj->x0);

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = s;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double t = prj->w[1]*(*yp + prj->y0);

    for (int ix = 0; ix < mx; ix++, thetap += spt, statp++) {
      *thetap = t;
      *statp  = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("carx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int cars2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != CAR) {
    if ((status = carset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double xi = prj->w[0]*(*phip) - prj->x0;

    double *xp = x + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = xi;
      xp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double eta = prj->w[0]*(*thetap) - prj->y0;

    for (int iphi = 0; iphi < mphi; iphi++, yp += sxy, statp++) {
      *yp = eta;
      *statp = 0;
    }
  }

  return 0;
}

/*============================================================================
*   MER: Mercator's projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     MER
*      prj->code    "MER"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    r0*(pi/180)
*      prj->w[1]    (180/pi)/r0
*      prj->prjx2s  Pointer to merx2s().
*      prj->prjs2x  Pointer to mers2x().
*===========================================================================*/

int merset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -MER) return 0;

  strcpy(prj->code, "MER");

  strcpy(prj->name, "Mercator's");
  prj->category  = CYLINDRICAL;
  prj->pvrange   = 0;
  prj->simplezen = 0;
  prj->equiareal = 0;
  prj->conformal = 1;
  prj->global    = 0;
  prj->divergent = 1;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 1.0;
    prj->w[1] = 1.0;
  } else {
    prj->w[0] = prj->r0*D2R;
    prj->w[1] = 1.0/prj->w[0];
  }

  prj->prjx2s = merx2s;
  prj->prjs2x = mers2x;

  prj->flag = (prj->flag == 1) ? -MER : MER;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int merx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != MER) {
    if ((status = merset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double s = prj->w[1]*(*xp + prj->x0);

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = s;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double t = 2.0*atand(exp((*yp + prj->y0)/prj->r0)) - 90.0;

    for (int ix = 0; ix < mx; ix++, thetap += spt, statp++) {
      *thetap = t;
      *statp  = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("merx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int mers2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != MER) {
    if ((status = merset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double xi  = prj->w[0]*(*phip) - prj->x0;

    double *xp = x + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = xi;
      xp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    int istat = 0;

    double eta;
    if (*thetap <= -90.0 || *thetap >= 90.0) {
      eta = 0.0;
      istat = 1;
      if (!status) status = PRJERR_BAD_WORLD_SET("mers2x");
    } else {
      eta = prj->r0*log(tand((*thetap+90.0)/2.0)) - prj->y0;
    }

    for (int iphi = 0; iphi < mphi; iphi++, yp += sxy, statp++) {
      *yp = eta;
      *statp = istat;
    }
  }

  return status;
}

/*============================================================================
*   SFL: Sanson-Flamsteed ("global sinusoid") projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     SFL
*      prj->code    "SFL"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    r0*(pi/180)
*      prj->w[1]    (180/pi)/r0
*      prj->prjx2s  Pointer to sflx2s().
*      prj->prjs2x  Pointer to sfls2x().
*===========================================================================*/

int sflset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -SFL) return 0;

  strcpy(prj->code, "SFL");

  strcpy(prj->name, "Sanson-Flamsteed");
  prj->category  = PSEUDOCYLINDRICAL;
  prj->pvrange   = 0;
  prj->simplezen = 0;
  prj->equiareal = 1;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 1.0;
    prj->w[1] = 1.0;
  } else {
    prj->w[0] = prj->r0*D2R;
    prj->w[1] = 1.0/prj->w[0];
  }

  prj->prjx2s = sflx2s;
  prj->prjs2x = sfls2x;

  prj->flag = (prj->flag == 1) ? -SFL : SFL;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int sflx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != SFL) {
    if ((status = sflset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double s = prj->w[1]*(*xp + prj->x0);

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = s;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yj = *yp + prj->y0;
    double s = cos(yj/prj->r0);

    int istat = 0;
    if (s == 0.0) {
      istat = 1;
      if (!status) status = PRJERR_BAD_PIX_SET("sflx2s");
    } else {
      s = 1.0/s;
    }

    double t = prj->w[1]*yj;

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      *phip  *= s;
      *thetap = t;
      *statp  = istat;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-12, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("sflx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int sfls2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != SFL) {
    if ((status = sflset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double xi = prj->w[0]*(*phip);

    double *xp = x + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = xi;
      xp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double xi  = cosd(*thetap);
    double eta = prj->w[0]*(*thetap) - prj->y0;

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      *xp = xi*(*xp) - prj->x0;
      *yp = eta;
      *statp = 0;
    }
  }

  return 0;
}

/*============================================================================
*   PAR: parabolic projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     PAR
*      prj->code    "PAR"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    r0*(pi/180)
*      prj->w[1]    (180/pi)/r0
*      prj->w[2]    pi*r0
*      prj->w[3]    1/(pi*r0)
*      prj->prjx2s  Pointer to parx2s().
*      prj->prjs2x  Pointer to pars2x().
*===========================================================================*/

int parset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -PAR) return 0;

  strcpy(prj->code, "PAR");

  strcpy(prj->name, "parabolic");
  prj->category  = PSEUDOCYLINDRICAL;
  prj->pvrange   = 0;
  prj->simplezen = 0;
  prj->equiareal = 1;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 1.0;
    prj->w[1] = 1.0;
    prj->w[2] = 180.0;
    prj->w[3] = 1.0/prj->w[2];
  } else {
    prj->w[0] = prj->r0*D2R;
    prj->w[1] = 1.0/prj->w[0];
    prj->w[2] = PI*prj->r0;
    prj->w[3] = 1.0/prj->w[2];
  }

  prj->prjx2s = parx2s;
  prj->prjs2x = pars2x;

  prj->flag = (prj->flag == 1) ? -PAR : PAR;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int parx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-13;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != PAR) {
    if ((status = parset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;
    double s = prj->w[1]*xj;
    double t = fabs(xj) - tol;

    double *phip   = phi   + rowoff;
    double *thetap = theta + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip   = s;
      *thetap = t;
      phip   += rowlen;
      thetap += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double r = prj->w[3]*(*yp + prj->y0);

    int istat = 0;
    double s, t;
    if (r > 1.0 || r < -1.0) {
      s = 0.0;
      t = 0.0;
      istat = 1;
      if (!status) status = PRJERR_BAD_PIX_SET("parx2s");

    } else {
      s = 1.0 - 4.0*r*r;
      if (s == 0.0) {
        // Deferred test.
        istat = -1;
      } else {
        s = 1.0/s;
      }

      t = 3.0*asind(r);
    }

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      if (istat < 0) {
        if (*thetap < 0.0) {
          *statp = 0;
        } else {
          *statp = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("parx2s");
        }
      } else {
        *statp = istat;
      }

      *phip  *= s;
      *thetap = t;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-12, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("parx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int pars2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != PAR) {
    if ((status = parset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double xi = prj->w[0]*(*phip);

    double *xp = x + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = xi;
      xp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double s = sind((*thetap)/3.0);
    double xi = (1.0 - 4.0*s*s);
    double eta = prj->w[2]*s - prj->y0;

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      *xp = xi*(*xp) - prj->x0;
      *yp = eta;
      *statp = 0;
    }
  }

  return 0;
}

/*============================================================================
*   MOL: Mollweide's projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     MOL
*      prj->code    "MOL"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    sqrt(2)*r0
*      prj->w[1]    sqrt(2)*r0/90
*      prj->w[2]    1/(sqrt(2)*r0)
*      prj->w[3]    90/r0
*      prj->prjx2s  Pointer to molx2s().
*      prj->prjs2x  Pointer to mols2x().
*===========================================================================*/

int molset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -MOL) return 0;

  strcpy(prj->code, "MOL");

  if (prj->r0 == 0.0) prj->r0 = R2D;

  strcpy(prj->name, "Mollweide's");
  prj->category  = PSEUDOCYLINDRICAL;
  prj->pvrange   = 0;
  prj->simplezen = 0;
  prj->equiareal = 1;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  prj->w[0] = SQRT2*prj->r0;
  prj->w[1] = prj->w[0]/90.0;
  prj->w[2] = 1.0/prj->w[0];
  prj->w[3] = 90.0/prj->r0;
  prj->w[4] = 2.0/PI;

  prj->prjx2s = molx2s;
  prj->prjs2x = mols2x;

  prj->flag = (prj->flag == 1) ? -MOL : MOL;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int molx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-12;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != MOL) {
    if ((status = molset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;
    double s = prj->w[3]*xj;
    double t = fabs(xj) - tol;

    double *phip   = phi   + rowoff;
    double *thetap = theta + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip   = s;
      *thetap = t;
      phip   += rowlen;
      thetap += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yj = *yp + prj->y0;
    double y0 = yj/prj->r0;
    double r  = 2.0 - y0*y0;

    int istat = 0;
    double s;
    if (r <= tol) {
      if (r < -tol) {
        istat = 1;
        if (!status) status = PRJERR_BAD_PIX_SET("molx2s");
      } else {
        // OK if fabs(x) < tol whence phi = 0.0.
        istat = -1;
      }

      r = 0.0;
      s = 0.0;

    } else {
      r = sqrt(r);
      s = 1.0/r;
    }

    double z = yj*prj->w[2];
    if (fabs(z) > 1.0) {
      if (fabs(z) > 1.0+tol) {
        z = 0.0;
        istat = 1;
        if (!status) status = PRJERR_BAD_PIX_SET("molx2s");
      } else {
        z = copysign(1.0, z) + y0*r/PI;
      }
    } else {
      z = asin(z)*prj->w[4] + y0*r/PI;
    }

    if (fabs(z) > 1.0) {
      if (fabs(z) > 1.0+tol) {
        z = 0.0;
        istat = 1;
        if (!status) status = PRJERR_BAD_PIX_SET("molx2s");
      } else {
        z = copysign(1.0, z);
      }
    }

    double t = asind(z);

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      if (istat < 0) {
        if (*thetap < 0.0) {
          *statp = 0;
        } else {
          *statp = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("molx2s");
        }
      } else {
        *statp = istat;
      }

      *phip  *= s;
      *thetap = t;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-11, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("molx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int mols2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  const double tol = 1.0e-13;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != MOL) {
    if ((status = molset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double xi = prj->w[1]*(*phip);

    double *xp = x + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = xi;
      xp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double xi, eta;

    if (fabs(*thetap) == 90.0) {
      xi  = 0.0;
      eta = copysign(prj->w[0], *thetap);

    } else if (*thetap == 0.0) {
      xi  = 1.0;
      eta = 0.0;

    } else {
      double u  = PI*sind(*thetap);
      double v0 = -PI;
      double v1 =  PI;
      double v  = u;
      for (int k = 0; k < 100; k++) {
        double resid = (v - u) + sin(v);
        if (resid < 0.0) {
          if (resid > -tol) break;
          v0 = v;
        } else {
          if (resid < tol) break;
          v1 = v;
        }
        v = (v0 + v1)/2.0;
      }

      double gamma = v/2.0;
      xi  = cos(gamma);
      eta = prj->w[0]*sin(gamma);
    }

    eta -= prj->y0;
    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      *xp = xi*(*xp) - prj->x0;
      *yp = eta;
      *statp = 0;
    }
  }

  return 0;
}

/*============================================================================
*   AIT: Hammer-Aitoff projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     AIT
*      prj->code    "AIT"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    2*r0**2
*      prj->w[1]    1/(2*r0)**2
*      prj->w[2]    1/(4*r0)**2
*      prj->w[3]    1/(2*r0)
*      prj->prjx2s  Pointer to aitx2s().
*      prj->prjs2x  Pointer to aits2x().
*===========================================================================*/

int aitset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -AIT) return 0;

  strcpy(prj->code, "AIT");

  if (prj->r0 == 0.0) prj->r0 = R2D;

  strcpy(prj->name, "Hammer-Aitoff");
  prj->category  = CONVENTIONAL;
  prj->pvrange   = 0;
  prj->simplezen = 0;
  prj->equiareal = 1;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  prj->w[0] = 2.0*prj->r0*prj->r0;
  prj->w[1] = 1.0/(2.0*prj->w[0]);
  prj->w[2] = prj->w[1]/4.0;
  prj->w[3] = 1.0/(2.0*prj->r0);

  prj->prjx2s = aitx2s;
  prj->prjs2x = aits2x;

  prj->flag = (prj->flag == 1) ? -AIT : AIT;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int aitx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-13;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != AIT) {
    if ((status = aitset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;
    double s  = 1.0 - xj*xj*prj->w[2];
    double t  = xj*prj->w[3];

    double *phip   = phi   + rowoff;
    double *thetap = theta + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip   = s;
      *thetap = t;
      phip   += rowlen;
      thetap += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yj  = *yp + prj->y0;
    double yj2 = yj*yj*prj->w[1];

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double s = *phip - yj2;

      int istat = 0;
      if (s < 0.5) {
        if (s < 0.5-tol) {
          istat = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("aitx2s");
        }

        s = 0.5;
      }

      double z = sqrt(s);
      double x0 = 2.0*z*z - 1.0;
      double y0 = z*(*thetap);
      if (x0 == 0.0 && y0 == 0.0) {
        *phip = 0.0;
      } else {
        *phip = 2.0*atan2d(y0, x0);
      }

      double t = z*yj/prj->r0;
      if (fabs(t) > 1.0) {
        if (fabs(t) > 1.0+tol) {
          istat = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("aitx2s");
        }
        t = copysign(90.0, t);

      } else {
        t = asind(t);
      }

      *thetap = t;
      *statp  = istat;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("aitx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int aits2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != AIT) {
    if ((status = aitset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double w = (*phip)/2.0;
    double sinphi, cosphi;
    sincosd(w, &sinphi, &cosphi);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinphi;
      *yp = cosphi;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double sinthe, costhe;
    sincosd(*thetap, &sinthe, &costhe);

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      double w = sqrt(prj->w[0]/(1.0 + costhe*(*yp)));
      *xp = 2.0*w*costhe*(*xp) - prj->x0;
      *yp = w*sinthe - prj->y0;
      *statp = 0;
    }
  }

  return 0;
}

/*============================================================================
*   COP: conic perspective projection.
*
*   Given:
*      prj->pv[1]   sigma = (theta2+theta1)/2
*      prj->pv[2]   delta = (theta2-theta1)/2, where theta1 and theta2 are the
*                   latitudes of the standard parallels, in degrees.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to sigma if undefined.
*      prj->theta0  Reset to sigma if undefined.
*
*   Returned:
*      prj->flag     COP
*      prj->code    "COP"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    C  = sin(sigma)
*      prj->w[1]    1/C
*      prj->w[2]    Y0 = r0*cos(delta)*cot(sigma)
*      prj->w[3]    r0*cos(delta)
*      prj->w[4]    1/(r0*cos(delta)
*      prj->w[5]    cot(sigma)
*      prj->prjx2s  Pointer to copx2s().
*      prj->prjs2x  Pointer to cops2x().
*===========================================================================*/

int copset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -COP) return 0;

  strcpy(prj->code, "COP");

  if (undefined(prj->pv[1])) {
    return PRJERR_BAD_PARAM_SET("copset");
  }
  if (undefined(prj->pv[2])) prj->pv[2] = 0.0;
  if (prj->r0 == 0.0) prj->r0 = R2D;

  strcpy(prj->name, "conic perspective");
  prj->category  = CONIC;
  prj->pvrange   = 102;
  prj->simplezen = 0;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 0;
  prj->divergent = 1;

  prj->w[0] = sind(prj->pv[1]);
  if (prj->w[0] == 0.0) {
    return PRJERR_BAD_PARAM_SET("copset");
  }

  prj->w[1] = 1.0/prj->w[0];

  prj->w[3] = prj->r0*cosd(prj->pv[2]);
  if (prj->w[3] == 0.0) {
    return PRJERR_BAD_PARAM_SET("copset");
  }

  prj->w[4] = 1.0/prj->w[3];
  prj->w[5] = 1.0/tand(prj->pv[1]);

  prj->w[2] = prj->w[3]*prj->w[5];

  prj->prjx2s = copx2s;
  prj->prjs2x = cops2x;

  prj->flag = (prj->flag == 1) ? -COP : COP;

  return prjoff(prj, 0.0, prj->pv[1]);
}

//----------------------------------------------------------------------------

int copx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != COP) {
    if ((status = copset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xj;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double dy  = prj->w[2] - (*yp + prj->y0);
    double dy2 = dy*dy;

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xj = *phip;

      double r = sqrt(xj*xj + dy2);
      if (prj->pv[1] < 0.0) r = -r;

      double alpha;
      if (r == 0.0) {
        alpha = 0.0;
      } else {
        alpha = atan2d(xj/r, dy/r);
      }

      *phip = alpha*prj->w[1];
      *thetap = prj->pv[1] + atand(prj->w[5] - r*prj->w[4]);
      *statp  = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("copx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int cops2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != COP) {
    if ((status = copset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double alpha = prj->w[0]*(*phip);
    double sinalpha, cosalpha;
    sincosd(alpha, &sinalpha, &cosalpha);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinalpha;
      *yp = cosalpha;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  double y0 = prj->y0 - prj->w[2];
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double t = *thetap - prj->pv[1];
    double s = cosd(t);

    int istat = 0;
    double r;
    if (s == 0.0) {
      // Latitude of divergence.
      r = 0.0;
      istat = 1;
      if (!status) status = PRJERR_BAD_WORLD_SET("cops2x");

    } else if (fabs(*thetap) == 90.0) {
      // Return an exact value at the poles.
      r = 0.0;

      // Bounds checking.
      if (prj->bounds&1) {
        if ((*thetap < 0.0) != (prj->pv[1] < 0.0)) {
          istat = 1;
          if (!status) status = PRJERR_BAD_WORLD_SET("cops2x");
        }
      }

    } else {
      r = prj->w[2] - prj->w[3]*sind(t)/s;

      // Bounds checking.
      if (prj->bounds&1) {
        if (r*prj->w[0] < 0.0) {
          istat = 1;
          if (!status) status = PRJERR_BAD_WORLD_SET("cops2x");
        }
      }
    }

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      *xp =  r*(*xp) - prj->x0;
      *yp = -r*(*yp) - y0;
      *statp = istat;
    }
  }

  return status;
}

/*============================================================================
*   COE: conic equal area projection.
*
*   Given:
*      prj->pv[1]   sigma = (theta2+theta1)/2
*      prj->pv[2]   delta = (theta2-theta1)/2, where theta1 and theta2 are the
*                   latitudes of the standard parallels, in degrees.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to sigma if undefined.
*      prj->theta0  Reset to sigma if undefined.
*
*   Returned:
*      prj->flag     COE
*      prj->code    "COE"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    C = (sin(theta1) + sin(theta2))/2
*      prj->w[1]    1/C
*      prj->w[2]    Y0 = chi*sqrt(psi - 2C*sind(sigma))
*      prj->w[3]    chi = r0/C
*      prj->w[4]    psi = 1 + sin(theta1)*sin(theta2)
*      prj->w[5]    2C
*      prj->w[6]    (1 + sin(theta1)*sin(theta2))*(r0/C)**2
*      prj->w[7]    C/(2*r0**2)
*      prj->w[8]    chi*sqrt(psi + 2C)
*      prj->prjx2s  Pointer to coex2s().
*      prj->prjs2x  Pointer to coes2x().
*===========================================================================*/

int coeset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -COE) return 0;

  strcpy(prj->code, "COE");

  if (undefined(prj->pv[1])) {
    return PRJERR_BAD_PARAM_SET("coeset");
  }
  if (undefined(prj->pv[2])) prj->pv[2] = 0.0;
  if (prj->r0 == 0.0) prj->r0 = R2D;

  strcpy(prj->name, "conic equal area");
  prj->category  = CONIC;
  prj->pvrange   = 102;
  prj->simplezen = 0;
  prj->equiareal = 1;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  double theta1 = prj->pv[1] - prj->pv[2];
  double theta2 = prj->pv[1] + prj->pv[2];

  prj->w[0] = (sind(theta1) + sind(theta2))/2.0;
  if (prj->w[0] == 0.0) {
    return PRJERR_BAD_PARAM_SET("coeset");
  }

  prj->w[1] = 1.0/prj->w[0];

  prj->w[3] = prj->r0/prj->w[0];
  prj->w[4] = 1.0 + sind(theta1)*sind(theta2);
  prj->w[5] = 2.0*prj->w[0];
  prj->w[6] = prj->w[3]*prj->w[3]*prj->w[4];
  prj->w[7] = 1.0/(2.0*prj->r0*prj->w[3]);
  prj->w[8] = prj->w[3]*sqrt(prj->w[4] + prj->w[5]);

  prj->w[2] = prj->w[3]*sqrt(prj->w[4] - prj->w[5]*sind(prj->pv[1]));

  prj->prjx2s = coex2s;
  prj->prjs2x = coes2x;

  prj->flag = (prj->flag == 1) ? -COE : COE;

  return prjoff(prj, 0.0, prj->pv[1]);
}

//----------------------------------------------------------------------------

int coex2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-12;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != COE) {
    if ((status = coeset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xj;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double dy  = prj->w[2] - (*yp + prj->y0);
    double dy2 = dy*dy;

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xj = *phip;

      double r = sqrt(xj*xj + dy2);
      if (prj->pv[1] < 0.0) r = -r;

      double alpha;
      if (r == 0.0) {
        alpha = 0.0;
      } else {
        alpha = atan2d(xj/r, dy/r);
      }

      int istat = 0;
      double t;
      if (fabs(r - prj->w[8]) < tol) {
        t = -90.0;
      } else {
        double w = (prj->w[6] - r*r)*prj->w[7];
        if (fabs(w) > 1.0) {
          if (fabs(w-1.0) < tol) {
            t = 90.0;
          } else if (fabs(w+1.0) < tol) {
            t = -90.0;
          } else {
            t = 0.0;
            istat = 1;
            if (!status) status = PRJERR_BAD_PIX_SET("coex2s");
          }
        } else {
          t = asind(w);
        }
      }

      *phip = alpha*prj->w[1];
      *thetap = t;
      *statp  = istat;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("coex2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int coes2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != COE) {
    if ((status = coeset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double alpha = prj->w[0]*(*phip);
    double sinalpha, cosalpha;
    sincosd(alpha, &sinalpha, &cosalpha);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinalpha;
      *yp = cosalpha;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  double y0 = prj->y0 - prj->w[2];
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double r;
    if (*thetap == -90.0) {
      r = prj->w[8];
    } else {
      r = prj->w[3]*sqrt(prj->w[4] - prj->w[5]*sind(*thetap));
    }

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      *xp =  r*(*xp) - prj->x0;
      *yp = -r*(*yp) - y0;
      *statp = 0;
    }
  }

  return 0;
}

/*============================================================================
*   COD: conic equidistant projection.
*
*   Given:
*      prj->pv[1]   sigma = (theta2+theta1)/2
*      prj->pv[2]   delta = (theta2-theta1)/2, where theta1 and theta2 are the
*                   latitudes of the standard parallels, in degrees.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to sigma if undefined.
*      prj->theta0  Reset to sigma if undefined.
*
*   Returned:
*      prj->flag     COD
*      prj->code    "COD"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    C = r0*sin(sigma)*sin(delta)/delta
*      prj->w[1]    1/C
*      prj->w[2]    Y0 = delta*cot(delta)*cot(sigma)
*      prj->w[3]    Y0 + sigma
*      prj->prjx2s  Pointer to codx2s().
*      prj->prjs2x  Pointer to cods2x().
*===========================================================================*/

int codset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -COD) return 0;

  strcpy(prj->code, "COD");

  if (undefined(prj->pv[1])) {
    return PRJERR_BAD_PARAM_SET("codset");
  }
  if (undefined(prj->pv[2])) prj->pv[2] = 0.0;
  if (prj->r0 == 0.0) prj->r0 = R2D;

  strcpy(prj->name, "conic equidistant");
  prj->category  = CONIC;
  prj->pvrange   = 102;
  prj->simplezen = 0;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->pv[2] == 0.0) {
    prj->w[0] = prj->r0*sind(prj->pv[1])*D2R;
  } else {
    prj->w[0] = prj->r0*sind(prj->pv[1])*sind(prj->pv[2])/prj->pv[2];
  }

  if (prj->w[0] == 0.0) {
    return PRJERR_BAD_PARAM_SET("codset");
  }

  prj->w[1] = 1.0/prj->w[0];
  prj->w[2] = prj->r0*cosd(prj->pv[2])*cosd(prj->pv[1])/prj->w[0];
  prj->w[3] = prj->w[2] + prj->pv[1];

  prj->prjx2s = codx2s;
  prj->prjs2x = cods2x;

  prj->flag = (prj->flag == 1) ? -COD : COD;

  return prjoff(prj, 0.0, prj->pv[1]);
}

//----------------------------------------------------------------------------

int codx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != COD) {
    if ((status = codset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xj;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double dy  = prj->w[2] - (*yp + prj->y0);
    double dy2 = dy*dy;

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xj = *phip;

      double r = sqrt(xj*xj + dy2);
      if (prj->pv[1] < 0.0) r = -r;

      double alpha;
      if (r == 0.0) {
        alpha = 0.0;
      } else {
        alpha = atan2d(xj/r, dy/r);
      }

      *phip = alpha*prj->w[1];
      *thetap = prj->w[3] - r;
      *statp  = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("codx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int cods2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != COD) {
    if ((status = codset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double alpha = prj->w[0]*(*phip);
    double sinalpha, cosalpha;
    sincosd(alpha, &sinalpha, &cosalpha);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinalpha;
      *yp = cosalpha;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  double y0 = prj->y0 - prj->w[2];
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double r = prj->w[3] - *thetap;

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      *xp =  r*(*xp) - prj->x0;
      *yp = -r*(*yp) - y0;
      *statp = 0;
    }
  }

  return 0;
}

/*============================================================================
*   COO: conic orthomorphic projection.
*
*   Given:
*      prj->pv[1]   sigma = (theta2+theta1)/2
*      prj->pv[2]   delta = (theta2-theta1)/2, where theta1 and theta2 are the
*                   latitudes of the standard parallels, in degrees.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to sigma if undefined.
*      prj->theta0  Reset to sigma if undefined.
*
*   Returned:
*      prj->flag     COO
*      prj->code    "COO"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    C = ln(cos(theta2)/cos(theta1))/ln(tan(tau2)/tan(tau1))
*                       where tau1 = (90 - theta1)/2
*                             tau2 = (90 - theta2)/2
*      prj->w[1]    1/C
*      prj->w[2]    Y0 = psi*tan((90-sigma)/2)**C
*      prj->w[3]    psi = (r0*cos(theta1)/C)/tan(tau1)**C
*      prj->w[4]    1/psi
*      prj->prjx2s  Pointer to coox2s().
*      prj->prjs2x  Pointer to coos2x().
*===========================================================================*/

int cooset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -COO) return 0;

  strcpy(prj->code, "COO");

  if (undefined(prj->pv[1])) {
    return PRJERR_BAD_PARAM_SET("cooset");
  }
  if (undefined(prj->pv[2])) prj->pv[2] = 0.0;
  if (prj->r0 == 0.0) prj->r0 = R2D;

  strcpy(prj->name, "conic orthomorphic");
  prj->category  = CONIC;
  prj->pvrange   = 102;
  prj->simplezen = 0;
  prj->equiareal = 0;
  prj->conformal = 1;
  prj->global    = 0;
  prj->divergent = 1;

  double theta1 = prj->pv[1] - prj->pv[2];
  double theta2 = prj->pv[1] + prj->pv[2];

  double tan1 = tand((90.0 - theta1)/2.0);
  double cos1 = cosd(theta1);

  if (theta1 == theta2) {
    prj->w[0] = sind(theta1);
  } else {
    double tan2 = tand((90.0 - theta2)/2.0);
    double cos2 = cosd(theta2);
    prj->w[0] = log(cos2/cos1)/log(tan2/tan1);
  }
  if (prj->w[0] == 0.0) {
    return PRJERR_BAD_PARAM_SET("cooset");
  }

  prj->w[1] = 1.0/prj->w[0];

  prj->w[3] = prj->r0*(cos1/prj->w[0])/pow(tan1,prj->w[0]);
  if (prj->w[3] == 0.0) {
    return PRJERR_BAD_PARAM_SET("cooset");
  }
  prj->w[2] = prj->w[3]*pow(tand((90.0 - prj->pv[1])/2.0),prj->w[0]);
  prj->w[4] = 1.0/prj->w[3];

  prj->prjx2s = coox2s;
  prj->prjs2x = coos2x;

  prj->flag = (prj->flag == 1) ? -COO : COO;

  return prjoff(prj, 0.0, prj->pv[1]);
}

//----------------------------------------------------------------------------

int coox2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != COO) {
    if ((status = cooset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xj;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double dy  = prj->w[2] - (*yp + prj->y0);
    double dy2 = dy*dy;

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xj = *phip;

      double r = sqrt(xj*xj + dy2);
      if (prj->pv[1] < 0.0) r = -r;

      double alpha;
      if (r == 0.0) {
        alpha = 0.0;
      } else {
        alpha = atan2d(xj/r, dy/r);
      }

      int istat = 0;
      double t;
      if (r == 0.0) {
        if (prj->w[0] < 0.0) {
          t = -90.0;
        } else {
          t = 0.0;
          istat = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("coox2s");
        }
      } else {
        t = 90.0 - 2.0*atand(pow(r*prj->w[4],prj->w[1]));
      }

      *phip = alpha*prj->w[1];
      *thetap = t;
      *statp  = istat;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("coox2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int coos2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != COO) {
    if ((status = cooset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double alpha = prj->w[0]*(*phip);
    double sinalpha, cosalpha;
    sincosd(alpha, &sinalpha, &cosalpha);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = sinalpha;
      *yp = cosalpha;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  double y0 = prj->y0 - prj->w[2];
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    int istat = 0;

    double r;
    if (*thetap == -90.0) {
      r = 0.0;
      if (prj->w[0] >= 0.0) {
        istat = 1;
        if (!status) status = PRJERR_BAD_WORLD_SET("coos2x");
      }
    } else {
      r = prj->w[3]*pow(tand((90.0 - *thetap)/2.0),prj->w[0]);
    }

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      *xp =  r*(*xp) - prj->x0;
      *yp = -r*(*yp) - y0;
      *statp = istat;
    }
  }

  return status;
}

/*============================================================================
*   BON: Bonne's projection.
*
*   Given:
*      prj->pv[1]   Bonne conformal latitude, theta1, in degrees.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     BON
*      prj->code    "BON"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[1]    r0*pi/180
*      prj->w[2]    Y0 = r0*(cot(theta1) + theta1*pi/180)
*      prj->prjx2s  Pointer to bonx2s().
*      prj->prjs2x  Pointer to bons2x().
*===========================================================================*/

int bonset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -BON) return 0;

  strcpy(prj->code, "BON");

  if (undefined(prj->pv[1])) {
    return PRJERR_BAD_PARAM_SET("bonset");
  }

  if (prj->pv[1] == 0.0) {
    // Sanson-Flamsteed.
    return sflset(prj);
  }

  strcpy(prj->name, "Bonne's");
  prj->category  = POLYCONIC;
  prj->pvrange   = 101;
  prj->simplezen = 0;
  prj->equiareal = 1;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[1] = 1.0;
    prj->w[2] = prj->r0*cosd(prj->pv[1])/sind(prj->pv[1]) + prj->pv[1];
  } else {
    prj->w[1] = prj->r0*D2R;
    prj->w[2] = prj->r0*(cosd(prj->pv[1])/sind(prj->pv[1]) + prj->pv[1]*D2R);
  }

  prj->prjx2s = bonx2s;
  prj->prjs2x = bons2x;

  prj->flag = (prj->flag == 1) ? -BON : BON;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int bonx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->pv[1] == 0.0) {
    // Sanson-Flamsteed.
    return sflx2s(prj, nx, ny, sxy, spt, x, y, phi, theta, stat);
  }

  int status;
  if (abs(prj->flag) != BON) {
    if ((status = bonset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xj;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double dy  = prj->w[2] - (*yp + prj->y0);
    double dy2 = dy*dy;

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xj = *phip;

      double r = sqrt(xj*xj + dy2);
      if (prj->pv[1] < 0.0) r = -r;

      double alpha;
      if (r == 0.0) {
        alpha = 0.0;
      } else {
        alpha = atan2d(xj/r, dy/r);
      }

      double s;
      double t = (prj->w[2] - r)/prj->w[1];
      double costhe = cosd(t);
      if (costhe == 0.0) {
        s = 0.0;
      } else {
        s = alpha*(r/prj->r0)/costhe;
      }

      *phip = s;
      *thetap = t;
      *statp  = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-11, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("bonx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int bons2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->pv[1] == 0.0) {
    // Sanson-Flamsteed.
    return sfls2x(prj, nphi, ntheta, spt, sxy, phi, theta, x, y, stat);
  }

  int status;
  if (abs(prj->flag) != BON) {
    if ((status = bonset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double s = prj->r0*(*phip);

    double *xp = x + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = s;
      xp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  double y0 = prj->y0 - prj->w[2];
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double r = prj->w[2] - prj->w[1]*(*thetap);
    double s = cosd(*thetap)/r;

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      double alpha = s*(*xp);
      double sinalpha, cosalpha;
      sincosd(alpha, &sinalpha, &cosalpha);
      *xp =  r*sinalpha - prj->x0;
      *yp = -r*cosalpha - y0;
      *statp = 0;
    }
  }

  return 0;
}

/*============================================================================
*   PCO: polyconic projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     PCO
*      prj->code    "PCO"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    r0*(pi/180)
*      prj->w[1]    (180/pi)/r0
*      prj->w[2]    2*r0
*      prj->w[3]    (pi/180)/(2*r0)
*      prj->prjx2s  Pointer to pcox2s().
*      prj->prjs2x  Pointer to pcos2x().
*===========================================================================*/

int pcoset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -PCO) return 0;

  strcpy(prj->code, "PCO");

  strcpy(prj->name, "polyconic");
  prj->category  = POLYCONIC;
  prj->pvrange   = 0;
  prj->simplezen = 0;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 1.0;
    prj->w[1] = 1.0;
    prj->w[2] = 360.0/PI;
  } else {
    prj->w[0] = prj->r0*D2R;
    prj->w[1] = 1.0/prj->w[0];
    prj->w[2] = 2.0*prj->r0;
  }
  prj->w[3] = D2R/prj->w[2];

  prj->prjx2s = pcox2s;
  prj->prjs2x = pcos2x;

  prj->flag = (prj->flag == 1) ? -PCO : PCO;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int pcox2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-12;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != PCO) {
    if ((status = pcoset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xj = *xp + prj->x0;

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xj;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yj = *yp + prj->y0;
    double w  = fabs(yj*prj->w[1]);

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xj = *phip;

      if (w < tol) {
        *phip = xj*prj->w[1];
        *thetap = 0.0;

      } else if (fabs(w-90.0) < tol) {
        *phip = 0.0;
        *thetap = copysign(90.0, yj);

      } else {
        double the, ymthe, tanthe;
        if (w < 1.0e-4) {
          // To avoid cot(theta) blowing up near theta == 0.
          the    = yj / (prj->w[0] + prj->w[3]*xj*xj);
          ymthe  = yj - prj->w[0]*the;
          tanthe = tand(the);

        } else {
          // Iterative solution using weighted division of the interval.
          double thepos = yj / prj->w[0];
          double theneg = 0.0;

          // Setting fneg = -fpos halves the interval in the first iter.
          double xx = xj*xj;
          double fpos =  xx;
          double fneg = -xx;

          for (int k = 0; k < 64; k++) {
            // Weighted division of the interval.
            double lambda = fpos/(fpos-fneg);
            if (lambda < 0.1) {
              lambda = 0.1;
            } else if (lambda > 0.9) {
              lambda = 0.9;
            }
            the = thepos - lambda*(thepos-theneg);

            // Compute the residue.
            ymthe  = yj - prj->w[0]*the;
            tanthe = tand(the);
            double f = xx + ymthe*(ymthe - prj->w[2]/tanthe);

            // Check for convergence.
            if (fabs(f) < tol) break;
            if (fabs(thepos-theneg) < tol) break;

            // Redefine the interval.
            if (f > 0.0) {
              thepos = the;
              fpos = f;
            } else {
              theneg = the;
              fneg = f;
            }
          }
        }

        double x1 = prj->r0 - ymthe*tanthe;
        double y1 = xj*tanthe;
        if (x1 == 0.0 && y1 == 0.0) {
          *phip = 0.0;
        } else {
          *phip = atan2d(y1, x1)/sind(the);
        }

        *thetap = the;
      }

      *statp = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-12, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("pcox2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int pcos2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != PCO) {
    if ((status = pcoset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double *xp = x + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = *phip;
      xp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    if (*thetap == 0.0) {
      for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
        *xp =  prj->w[0]*(*xp) - prj->x0;
        *yp = -prj->y0;
        *statp = 0;
      }

    } else if (fabs(*thetap) < 1.0e-4) {
      // To avoid cot(theta) blowing up near theta == 0.
      for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
        *xp = prj->w[0]*(*xp)*cosd(*thetap) - prj->x0;
        *yp = (prj->w[0] + prj->w[3]*(*xp)*(*xp))*(*thetap) - prj->y0;
        *statp = 0;
      }

    } else {
      double therad = (*thetap)*D2R;
      double sinthe, costhe;
      sincosd(*thetap, &sinthe, &costhe);

      for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
        double sinpsi, cospsi;
        sincosd((*xp)*sinthe, &sinpsi, &cospsi);
        double cotthe = costhe/sinthe;
        *xp = prj->r0*cotthe*sinpsi - prj->x0;
        *yp = prj->r0*(cotthe*(1.0 - cospsi) + therad) - prj->y0;
        *statp = 0;
      }
    }
  }

  return 0;
}

/*============================================================================
*   TSC: tangential spherical cube projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     TSC
*      prj->code    "TSC"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    r0*(pi/4)
*      prj->w[1]    (4/pi)/r0
*      prj->prjx2s  Pointer to tscx2s().
*      prj->prjs2x  Pointer to tscs2x().
*===========================================================================*/

int tscset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -TSC) return 0;

  strcpy(prj->code, "TSC");

  strcpy(prj->name, "tangential spherical cube");
  prj->category  = QUADCUBE;
  prj->pvrange   = 0;
  prj->simplezen = 0;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 45.0;
    prj->w[1] = 1.0/45.0;
  } else {
    prj->w[0] = prj->r0*PI/4.0;
    prj->w[1] = 1.0/prj->w[0];
  }

  prj->prjx2s = tscx2s;
  prj->prjs2x = tscs2x;

  prj->flag = (prj->flag == 1) ? -TSC : TSC;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int tscx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != TSC) {
    if ((status = tscset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xf = (*xp + prj->x0)*prj->w[1];

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xf;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yf = (*yp + prj->y0)*prj->w[1];

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xf = *phip;

      // Bounds checking.
      if (fabs(xf) <= 1.0) {
        if (fabs(yf) > 3.0) {
          *phip = 0.0;
          *thetap = 0.0;
          *statp  = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("tscx2s");
          continue;
        }
      } else {
        if (fabs(xf) > 7.0 || fabs(yf) > 1.0) {
          *phip = 0.0;
          *thetap = 0.0;
          *statp  = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("tscx2s");
          continue;
        }
      }

      // Map negative faces to the other side.
      if (xf < -1.0) xf += 8.0;

      // Determine the face.
      double l, m, n;
      if (xf > 5.0) {
        // face = 4
        xf = xf - 6.0;
        m  = -1.0/sqrt(1.0 + xf*xf + yf*yf);
        l  = -m*xf;
        n  = -m*yf;
      } else if (xf > 3.0) {
        // face = 3
        xf = xf - 4.0;
        l  = -1.0/sqrt(1.0 + xf*xf + yf*yf);
        m  =  l*xf;
        n  = -l*yf;
      } else if (xf > 1.0) {
        // face = 2
        xf = xf - 2.0;
        m  =  1.0/sqrt(1.0 + xf*xf + yf*yf);
        l  = -m*xf;
        n  =  m*yf;
      } else if (yf > 1.0) {
        // face = 0
        yf = yf - 2.0;
        n  = 1.0/sqrt(1.0 + xf*xf + yf*yf);
        l  = -n*yf;
        m  =  n*xf;
      } else if (yf < -1.0) {
        // face = 5
        yf = yf + 2.0;
        n  = -1.0/sqrt(1.0 + xf*xf + yf*yf);
        l  = -n*yf;
        m  = -n*xf;
      } else {
        // face = 1
        l  =  1.0/sqrt(1.0 + xf*xf + yf*yf);
        m  =  l*xf;
        n  =  l*yf;
      }

      if (l == 0.0 && m == 0.0) {
        *phip = 0.0;
      } else {
        *phip = atan2d(m, l);
      }

      *thetap = asind(n);
      *statp  = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("tscx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int tscs2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  const double tol = 1.0e-12;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != TSC) {
    if ((status = tscset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double sinphi, cosphi;
    sincosd(*phip, &sinphi, &cosphi);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = cosphi;
      *yp = sinphi;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double sinthe, costhe;
    sincosd(*thetap, &sinthe, &costhe);

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      double l = costhe*(*xp);
      double m = costhe*(*yp);
      double n = sinthe;

      int face = 0;
      double zeta = n;
      if (l > zeta) {
        face = 1;
        zeta = l;
      }
      if (m > zeta) {
        face = 2;
        zeta = m;
      }
      if (-l > zeta) {
        face = 3;
        zeta = -l;
      }
      if (-m > zeta) {
        face = 4;
        zeta = -m;
      }
      if (-n > zeta) {
        face = 5;
        zeta = -n;
      }

      double xf, yf, x0, y0;
      switch (face) {
      case 1:
        xf =  m/zeta;
        yf =  n/zeta;
        x0 =  0.0;
        y0 =  0.0;
        break;
      case 2:
        xf = -l/zeta;
        yf =  n/zeta;
        x0 =  2.0;
        y0 =  0.0;
        break;
      case 3:
        xf = -m/zeta;
        yf =  n/zeta;
        x0 =  4.0;
        y0 =  0.0;
        break;
      case 4:
        xf =  l/zeta;
        yf =  n/zeta;
        x0 =  6.0;
        y0 =  0.0;
        break;
      case 5:
        xf =  m/zeta;
        yf =  l/zeta;
        x0 =  0.0;
        y0 = -2.0;
        break;
      default:
        // face == 0
        xf =  m/zeta;
        yf = -l/zeta;
        x0 =  0.0;
        y0 =  2.0;
        break;
      }

      int istat = 0;
      if (fabs(xf) > 1.0) {
        if (fabs(xf) > 1.0+tol) {
          istat = 1;
          if (!status) status = PRJERR_BAD_WORLD_SET("tscs2x");
        }
        xf = copysign(1.0, xf);
      }
      if (fabs(yf) > 1.0) {
        if (fabs(yf) > 1.0+tol) {
          istat = 1;
          if (!status) status = PRJERR_BAD_WORLD_SET("tscs2x");
        }
        yf = copysign(1.0, yf);
      }

      *xp = prj->w[0]*(xf + x0) - prj->x0;
      *yp = prj->w[0]*(yf + y0) - prj->y0;
      *statp = istat;
    }
  }

  return status;
}

/*============================================================================
*   CSC: COBE quadrilateralized spherical cube projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     CSC
*      prj->code    "CSC"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    r0*(pi/4)
*      prj->w[1]    (4/pi)/r0
*      prj->prjx2s  Pointer to cscx2s().
*      prj->prjs2x  Pointer to cscs2x().
*===========================================================================*/

int cscset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -CSC) return 0;

  strcpy(prj->code, "CSC");

  strcpy(prj->name, "COBE quadrilateralized spherical cube");
  prj->category  = QUADCUBE;
  prj->pvrange   = 0;
  prj->simplezen = 0;
  prj->equiareal = 0;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 45.0;
    prj->w[1] = 1.0/45.0;
  } else {
    prj->w[0] = prj->r0*PI/4.0;
    prj->w[1] = 1.0/prj->w[0];
  }

  prj->prjx2s = cscx2s;
  prj->prjs2x = cscs2x;

  prj->flag = (prj->flag == 1) ? -CSC : CSC;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int cscx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const float p00 = -0.27292696f;
  const float p10 = -0.07629969f;
  const float p20 = -0.22797056f;
  const float p30 =  0.54852384f;
  const float p40 = -0.62930065f;
  const float p50 =  0.25795794f;
  const float p60 =  0.02584375f;
  const float p01 = -0.02819452f;
  const float p11 = -0.01471565f;
  const float p21 =  0.48051509f;
  const float p31 = -1.74114454f;
  const float p41 =  1.71547508f;
  const float p51 = -0.53022337f;
  const float p02 =  0.27058160f;
  const float p12 = -0.56800938f;
  const float p22 =  0.30803317f;
  const float p32 =  0.98938102f;
  const float p42 = -0.83180469f;
  const float p03 = -0.60441560f;
  const float p13 =  1.50880086f;
  const float p23 = -0.93678576f;
  const float p33 =  0.08693841f;
  const float p04 =  0.93412077f;
  const float p14 = -1.41601920f;
  const float p24 =  0.33887446f;
  const float p05 = -0.63915306f;
  const float p15 =  0.52032238f;
  const float p06 =  0.14381585f;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != CSC) {
    if ((status = cscset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    float xf = (float)((*xp + prj->x0)*prj->w[1]);

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xf;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    float yf = (float)((*yp + prj->y0)*prj->w[1]);

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      float xf = (float)(*phip);

      // Bounds checking.
      if (fabs((double)xf) <= 1.0) {
        if (fabs((double)yf) > 3.0) {
          *phip = 0.0;
          *thetap = 0.0;
          *statp  = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("cscx2s");
          continue;
        }
      } else {
        if (fabs((double)xf) > 7.0 || fabs((double)yf) > 1.0) {
          *phip = 0.0;
          *thetap = 0.0;
          *statp  = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("cscx2s");
          continue;
        }
      }

      // Map negative faces to the other side.
      if (xf < -1.0f) xf += 8.0f;

      // Determine the face.
      int face = 0;
      if (xf > 5.0f) {
        face = 4;
        xf = xf - 6.0f;
      } else if (xf > 3.0f) {
        face = 3;
        xf = xf - 4.0f;
      } else if (xf > 1.0f) {
        face = 2;
        xf = xf - 2.0f;
      } else if (yf > 1.0f) {
        face = 0;
        yf = yf - 2.0f;
      } else if (yf < -1.0f) {
        face = 5;
        yf = yf + 2.0f;
      } else {
        face = 1;
      }

      float xx  =  xf*xf;
      float yy  =  yf*yf;

      float z0, z1, z2, z3, z4, z5, z6;
      z0 = p00 + xx*(p10 + xx*(p20 + xx*(p30 + xx*(p40 + xx*(p50 +
                 xx*(p60))))));
      z1 = p01 + xx*(p11 + xx*(p21 + xx*(p31 + xx*(p41 + xx*(p51)))));
      z2 = p02 + xx*(p12 + xx*(p22 + xx*(p32 + xx*(p42))));
      z3 = p03 + xx*(p13 + xx*(p23 + xx*(p33)));
      z4 = p04 + xx*(p14 + xx*(p24));
      z5 = p05 + xx*(p15);
      z6 = p06;

      float chi;
      chi = z0 + yy*(z1 + yy*(z2 + yy*(z3 + yy*(z4 + yy*(z5 + yy*z6)))));
      chi = xf + xf*(1.0f - xx)*chi;

      z0 = p00 + yy*(p10 + yy*(p20 + yy*(p30 + yy*(p40 + yy*(p50 +
                 yy*(p60))))));
      z1 = p01 + yy*(p11 + yy*(p21 + yy*(p31 + yy*(p41 + yy*(p51)))));
      z2 = p02 + yy*(p12 + yy*(p22 + yy*(p32 + yy*(p42))));
      z3 = p03 + yy*(p13 + yy*(p23 + yy*(p33)));
      z4 = p04 + yy*(p14 + yy*(p24));
      z5 = p05 + yy*(p15);
      z6 = p06;

      float psi;
      psi = z0 + xx*(z1 + xx*(z2 + xx*(z3 + xx*(z4 + xx*(z5 + xx*z6)))));
      psi = yf + yf*(1.0f - yy)*psi;

      double l, m, n;
      double t = 1.0/sqrt((double)(chi*chi + psi*psi) + 1.0);
      switch (face) {
      case 1:
        l =  t;
        m =  chi*l;
        n =  psi*l;
        break;
      case 2:
        m =  t;
        l = -chi*m;
        n =  psi*m;
        break;
      case 3:
        l = -t;
        m =  chi*l;
        n = -psi*l;
        break;
      case 4:
        m = -t;
        l = -chi*m;
        n = -psi*m;
        break;
      case 5:
        n = -t;
        l = -psi*n;
        m = -chi*n;
        break;
      default:
        // face == 0
        n =  t;
        l = -psi*n;
        m =  chi*n;
        break;
      }

      if (l == 0.0 && m == 0.0) {
        *phip = 0.0;
      } else {
        *phip = atan2d(m, l);
      }

      *thetap = asind(n);
      *statp  = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("cscx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int cscs2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  const double tol = 1.0e-7;

  const float gstar  =  1.37484847732f;
  const float mm     =  0.004869491981f;
  const float gamma  = -0.13161671474f;
  const float omega1 = -0.159596235474f;
  const float d0  =  0.0759196200467f;
  const float d1  = -0.0217762490699f;
  const float c00 =  0.141189631152f;
  const float c10 =  0.0809701286525f;
  const float c01 = -0.281528535557f;
  const float c11 =  0.15384112876f;
  const float c20 = -0.178251207466f;
  const float c02 =  0.106959469314f;


  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != CSC) {
    if ((status = cscset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double sinphi, cosphi;
    sincosd(*phip, &sinphi, &cosphi);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = cosphi;
      *yp = sinphi;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double sinthe, costhe;
    sincosd(*thetap, &sinthe, &costhe);

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      double l = costhe*(*xp);
      double m = costhe*(*yp);
      double n = sinthe;

      int face = 0;
      double zeta = n;
      if (l > zeta) {
        face = 1;
        zeta = l;
      }
      if (m > zeta) {
        face = 2;
        zeta = m;
      }
      if (-l > zeta) {
        face = 3;
        zeta = -l;
      }
      if (-m > zeta) {
        face = 4;
        zeta = -m;
      }
      if (-n > zeta) {
        face = 5;
        zeta = -n;
      }

      double eta, xi;
      float  x0, y0;
      switch (face) {
      case 1:
        xi  =  m;
        eta =  n;
        x0  =  0.0;
        y0  =  0.0;
        break;
      case 2:
        xi  = -l;
        eta =  n;
        x0  =  2.0;
        y0  =  0.0;
        break;
      case 3:
        xi  = -m;
        eta =  n;
        x0  =  4.0;
        y0  =  0.0;
        break;
      case 4:
        xi  =  l;
        eta =  n;
        x0  =  6.0;
        y0  =  0.0;
        break;
      case 5:
        xi  =  m;
        eta =  l;
        x0  =  0.0;
        y0  = -2.0;
        break;
      default:
        // face == 0
        xi  =  m;
        eta = -l;
        x0  =  0.0;
        y0  =  2.0;
        break;
      }

      float chi = (float)( xi/zeta);
      float psi = (float)(eta/zeta);

      float chi2 = chi*chi;
      float psi2 = psi*psi;
      float chi2co = 1.0f - chi2;
      float psi2co = 1.0f - psi2;

      // Avoid floating underflows.
      float chipsi = (float)fabs((double)(chi*psi));
      float chi4   = (chi2 > 1.0e-16f) ? chi2*chi2 : 0.0f;
      float psi4   = (psi2 > 1.0e-16f) ? psi2*psi2 : 0.0f;
      float chi2psi2 = (chipsi > 1.0e-16f) ? chi2*psi2 : 0.0f;

      float xf, yf;
      xf = chi*(chi2 + chi2co*(gstar + psi2*(gamma*chi2co + mm*chi2 +
                psi2co*(c00 + c10*chi2 + c01*psi2 + c11*chi2psi2 + c20*chi4 +
                c02*psi4)) + chi2*(omega1 - chi2co*(d0 + d1*chi2))));
      yf = psi*(psi2 + psi2co*(gstar + chi2*(gamma*psi2co + mm*psi2 +
                chi2co*(c00 + c10*psi2 + c01*chi2 + c11*chi2psi2 + c20*psi4 +
                c02*chi4)) + psi2*(omega1 - psi2co*(d0 + d1*psi2))));

      int istat = 0;
      if (fabs((double)xf) > 1.0) {
        if (fabs((double)xf) > 1.0+tol) {
          istat = 1;
          if (!status) status = PRJERR_BAD_WORLD_SET("cscs2x");
        }
        xf = (float)copysign(1.0, (double)xf);
      }
      if (fabs((double)yf) > 1.0) {
        if (fabs((double)yf) > 1.0+tol) {
          istat = 1;
          if (!status) status = PRJERR_BAD_WORLD_SET("cscs2x");
        }
        yf = (float)copysign(1.0, (double)yf);
      }

      *xp = prj->w[0]*(xf + x0) - prj->x0;
      *yp = prj->w[0]*(yf + y0) - prj->y0;
      *statp = istat;
    }
  }

  return status;
}

/*============================================================================
*   QSC: quadrilaterilized spherical cube projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     QSC
*      prj->code    "QSC"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    r0*(pi/4)
*      prj->w[1]    (4/pi)/r0
*      prj->prjx2s  Pointer to qscx2s().
*      prj->prjs2x  Pointer to qscs2x().
*===========================================================================*/

int qscset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -QSC) return 0;

  strcpy(prj->code, "QSC");

  strcpy(prj->name, "quadrilateralized spherical cube");
  prj->category  = QUADCUBE;
  prj->pvrange   = 0;
  prj->simplezen = 0;
  prj->equiareal = 1;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 45.0;
    prj->w[1] = 1.0/45.0;
  } else {
    prj->w[0] = prj->r0*PI/4.0;
    prj->w[1] = 1.0/prj->w[0];
  }

  prj->prjx2s = qscx2s;
  prj->prjs2x = qscs2x;

  prj->flag = (prj->flag == 1) ? -QSC : QSC;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int qscx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-12;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != QSC) {
    if ((status = qscset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xf = (*xp + prj->x0)*prj->w[1];

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xf;
      phip += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yf = (*yp + prj->y0)*prj->w[1];

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xf = *phip;

      // Bounds checking.
      if (fabs(xf) <= 1.0) {
        if (fabs(yf) > 3.0) {
          *phip = 0.0;
          *thetap = 0.0;
          *statp  = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("qscx2s");
          continue;
        }
      } else {
        if (fabs(xf) > 7.0 || fabs(yf) > 1.0) {
          *phip = 0.0;
          *thetap = 0.0;
          *statp  = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("qscx2s");
          continue;
        }
      }

      // Map negative faces to the other side.
      if (xf < -1.0) xf += 8.0;

      int face = 0;
      // Determine the face.
      if (xf > 5.0) {
        face = 4;
        xf -= 6.0;
      } else if (xf > 3.0) {
        face = 3;
        xf -= 4.0;
      } else if (xf > 1.0) {
        face = 2;
        xf -= 2.0;
      } else if (yf > 1.0) {
        face = 0;
        yf -= 2.0;
      } else if (yf < -1.0) {
        face = 5;
        yf += 2.0;
      } else {
        face = 1;
      }

      double omega, tau, w, zeta, zeco;
      int direct = (fabs(xf) > fabs(yf));
      if (direct) {
        if (xf == 0.0) {
          omega = 0.0;
          tau  = 1.0;
          zeta = 1.0;
          zeco = 0.0;
        } else {
          w = 15.0*yf/xf;
          omega = sind(w)/(cosd(w) - SQRT2INV);
          tau  = 1.0 + omega*omega;
          zeco = xf*xf*(1.0 - 1.0/sqrt(1.0 + tau));
          zeta = 1.0 - zeco;
        }
      } else {
        if (yf == 0.0) {
          omega = 0.0;
          tau  = 1.0;
          zeta = 1.0;
          zeco = 0.0;
        } else {
          w = 15.0*xf/yf;
          double sinw, cosw;
          sincosd(w, &sinw, &cosw);
          omega = sinw/(cosw - SQRT2INV);
          tau  = 1.0 + omega*omega;
          zeco = yf*yf*(1.0 - 1.0/sqrt(1.0 + tau));
          zeta = 1.0 - zeco;
        }
      }

      if (zeta < -1.0) {
        if (zeta < -1.0-tol) {
          *phip = 0.0;
          *thetap = 0.0;
          *statp  = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("qscx2s");
          continue;
        }

        zeta = -1.0;
        zeco =  2.0;
        w    =  0.0;
      } else {
        w = sqrt(zeco*(2.0-zeco)/tau);
      }

      double l, m, n;
      switch (face) {
      case 1:
        l = zeta;
        if (direct) {
          m = w;
          if (xf < 0.0) m = -m;
          n = m*omega;
        } else {
          n = w;
          if (yf < 0.0) n = -n;
          m = n*omega;
        }
        break;
      case 2:
        m = zeta;
        if (direct) {
          l = w;
          if (xf > 0.0) l = -l;
          n = -l*omega;
        } else {
          n = w;
          if (yf < 0.0) n = -n;
          l = -n*omega;
        }
        break;
      case 3:
        l = -zeta;
        if (direct) {
          m = w;
          if (xf > 0.0) m = -m;
          n = -m*omega;
        } else {
          n = w;
          if (yf < 0.0) n = -n;
          m = -n*omega;
        }
        break;
      case 4:
        m = -zeta;
        if (direct) {
          l = w;
          if (xf < 0.0) l = -l;
          n = l*omega;
        } else {
          n = w;
          if (yf < 0.0) n = -n;
          l = n*omega;
        }
        break;
      case 5:
        n = -zeta;
        if (direct) {
          m = w;
          if (xf < 0.0) m = -m;
          l = m*omega;
        } else {
          l = w;
          if (yf < 0.0) l = -l;
          m = l*omega;
        }
        break;
      default:
        // face == 0
        n = zeta;
        if (direct) {
          m = w;
          if (xf < 0.0) m = -m;
          l = -m*omega;
        } else {
          l = w;
          if (yf > 0.0) l = -l;
          m = -l*omega;
        }
        break;
      }

      if (l == 0.0 && m == 0.0) {
        *phip = 0.0;
      } else {
        *phip = atan2d(m, l);
      }

      *thetap = asind(n);
      *statp  = 0;
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-13, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("qscx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int qscs2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  const double tol = 1.0e-12;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != QSC) {
    if ((status = qscset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }

  status = 0;


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double sinphi, cosphi;
    sincosd(*phip, &sinphi, &cosphi);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      *xp = cosphi;
      *yp = sinphi;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double sinthe, costhe;
    sincosd(*thetap, &sinthe, &costhe);

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      if (fabs(*thetap) == 90.0) {
        *xp = -prj->x0;
        *yp = copysign(2.0*prj->w[0], *thetap) - prj->y0;
        *statp = 0;
        continue;
      }

      double l = costhe*(*xp);
      double m = costhe*(*yp);
      double n = sinthe;

      int face = 0;
      double zeta = n;
      if (l > zeta) {
        face = 1;
        zeta = l;
      }
      if (m > zeta) {
        face = 2;
        zeta = m;
      }
      if (-l > zeta) {
        face = 3;
        zeta = -l;
      }
      if (-m > zeta) {
        face = 4;
        zeta = -m;
      }
      if (-n > zeta) {
        face = 5;
        zeta = -n;
      }

      double zeco = 1.0 - zeta;

      double xi, eta, x0, y0;
      switch (face) {
      case 1:
        xi  = m;
        eta = n;
        if (zeco < 1.0e-8) {
          // Small angle formula.
          double t = (*thetap)*D2R;
          double p = atan2(*yp, *xp);
          zeco = (p*p + t*t)/2.0;
        }
        x0 = 0.0;
        y0 = 0.0;
        break;
      case 2:
        xi  = -l;
        eta =  n;
        if (zeco < 1.0e-8) {
          // Small angle formula.
          double t = (*thetap)*D2R;
          double p = atan2(*yp, *xp) - PI/2.0;
          zeco = (p*p + t*t)/2.0;
        }
        x0 = 2.0;
        y0 = 0.0;
        break;
      case 3:
        xi  = -m;
        eta =  n;
        if (zeco < 1.0e-8) {
          // Small angle formula.
          double t = (*thetap)*D2R;
          double p = atan2(*yp, *xp);
          p -= copysign(PI, p);
          zeco = (p*p + t*t)/2.0;
        }
        x0 = 4.0;
        y0 = 0.0;
        break;
      case 4:
        xi  = l;
        eta = n;
        if (zeco < 1.0e-8) {
          // Small angle formula.
          double t = (*thetap)*D2R;
          double p = atan2(*yp, *xp) + PI/2.0;
          zeco = (p*p + t*t)/2.0;
        }
        x0 = 6;
        y0 = 0.0;
        break;
      case 5:
        xi  =  m;
        eta =  l;
        if (zeco < 1.0e-8) {
          // Small angle formula.
          double t = (*thetap + 90.0)*D2R;
          zeco = t*t/2.0;
        }
        x0 =  0.0;
        y0 = -2;
         break;
      default:
        // face == 0
        xi  =  m;
        eta = -l;
        if (zeco < 1.0e-8) {
          // Small angle formula.
          double t = (90.0 - *thetap)*D2R;
          zeco = t*t/2.0;
        }
        x0 = 0.0;
        y0 = 2.0;
        break;
      }

      double omega, tau;
      double xf = 0.0;
      double yf = 0.0;
      if (xi != 0.0 || eta != 0.0) {
        if (-xi > fabs(eta)) {
          omega = eta/xi;
          tau = 1.0 + omega*omega;
          xf  = -sqrt(zeco/(1.0 - 1.0/sqrt(1.0+tau)));
          yf  = (xf/15.0)*(atand(omega) - asind(omega/sqrt(tau+tau)));
        } else if (xi > fabs(eta)) {
          omega = eta/xi;
          tau = 1.0 + omega*omega;
          xf  =  sqrt(zeco/(1.0 - 1.0/sqrt(1.0+tau)));
          yf  = (xf/15.0)*(atand(omega) - asind(omega/sqrt(tau+tau)));
        } else if (-eta >= fabs(xi)) {
          omega = xi/eta;
          tau = 1.0 + omega*omega;
          yf  = -sqrt(zeco/(1.0 - 1.0/sqrt(1.0+tau)));
          xf  = (yf/15.0)*(atand(omega) - asind(omega/sqrt(tau+tau)));
        } else if (eta >= fabs(xi)) {
          omega = xi/eta;
          tau = 1.0 + omega*omega;
          yf  =  sqrt(zeco/(1.0 - 1.0/sqrt(1.0+tau)));
          xf  = (yf/15.0)*(atand(omega) - asind(omega/sqrt(tau+tau)));
        }
      }

      int istat = 0;
      if (fabs(xf) > 1.0) {
        if (fabs(xf) > 1.0+tol) {
          istat = 1;
          if (!status) status = PRJERR_BAD_WORLD_SET("qscs2x");
        }
        xf = copysign(1.0, xf);
      }
      if (fabs(yf) > 1.0) {
        if (fabs(yf) > 1.0+tol) {
          istat = 1;
          if (!status) status = PRJERR_BAD_WORLD_SET("qscs2x");
        }
        yf = copysign(1.0, yf);
      }

      *xp = prj->w[0]*(xf + x0) - prj->x0;
      *yp = prj->w[0]*(yf + y0) - prj->y0;
      *statp = istat;
    }
  }

  return status;
}

/*============================================================================
*   HPX: HEALPix projection.
*
*   Given:
*      prj->pv[1]   H - the number of facets in longitude.
*      prj->pv[2]   K - the number of facets in latitude
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     HPX
*      prj->code    "HPX"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->m       True if H is odd.
*      prj->n       True if K is odd.
*      prj->w[0]    r0*(pi/180)
*      prj->w[1]    (180/pi)/r0
*      prj->w[2]    (K-1)/K
*      prj->w[3]    90*K/H
*      prj->w[4]    (K+1)/2
*      prj->w[5]    90*(K-1)/H
*      prj->w[6]    180/H
*      prj->w[7]    H/360
*      prj->w[8]    r0*(pi/180)*(90*K/H)
*      prj->w[9]    r0*(pi/180)*(180/H)
*      prj->prjx2s  Pointer to hpxx2s().
*      prj->prjs2x  Pointer to hpxs2x().
*===========================================================================*/

int hpxset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -HPX) return 0;

  strcpy(prj->code, "HPX");

  if (undefined(prj->pv[1])) prj->pv[1] = 4.0;
  if (undefined(prj->pv[2])) prj->pv[2] = 3.0;

  strcpy(prj->name, "HEALPix");
  prj->category  = HEALPIX;
  prj->pvrange   = 102;
  prj->simplezen = 0;
  prj->equiareal = 1;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->pv[1] <= 0.0 || prj->pv[2] <= 0.0) {
    return PRJERR_BAD_PARAM_SET("hpxset");
  }

  prj->m = ((int)(prj->pv[1]+0.5))%2;
  prj->n = ((int)(prj->pv[2]+0.5))%2;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 1.0;
    prj->w[1] = 1.0;
  } else {
    prj->w[0] = prj->r0*D2R;
    prj->w[1] = R2D/prj->r0;
  }

  prj->w[2] = (prj->pv[2] - 1.0) / prj->pv[2];
  prj->w[3] = 90.0 * prj->pv[2] / prj->pv[1];
  prj->w[4] = (prj->pv[2] + 1.0) / 2.0;
  prj->w[5] = 90.0 * (prj->pv[2] - 1.0) / prj->pv[1];
  prj->w[6] = 180.0 / prj->pv[1];
  prj->w[7] = prj->pv[1] / 360.0;
  prj->w[8] = prj->w[3] * prj->w[0];
  prj->w[9] = prj->w[6] * prj->w[0];

  prj->prjx2s = hpxx2s;
  prj->prjs2x = hpxs2x;

  prj->flag = (prj->flag == 1) ? -HPX : HPX;

  return prjoff(prj, 0.0, 0.0);
}

//----------------------------------------------------------------------------

int hpxx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != HPX) {
    if ((status = hpxset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double s = prj->w[1] * (*xp + prj->x0);
    // x_c for K odd or theta > 0.
    double t;
    t = -180.0 + (2.0 * floor((*xp + 180.0) * prj->w[7]) + 1.0) * prj->w[6];
    t = prj->w[1] * (*xp - t);

    double *phip   = phi + rowoff;
    double *thetap = theta + rowoff;
    for (int iy = 0; iy < my; iy++) {
      // theta[] is used to hold (x - x_c).
      *phip   = s;
      *thetap = t;
      phip   += rowlen;
      thetap += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;

  double slim = prj->w[6] + 1e-12;
  double ylim = prj->w[9] * prj->w[4];

  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yr = prj->w[1]*(*yp + prj->y0);
    double absy = fabs(yr);

    int istat = 0;
    if (absy <= prj->w[5]) {
      // Equatorial regime.
      double t = asind(yr/prj->w[3]);
      for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
        *thetap = t;
        *statp  = 0;
      }

    } else if (absy <= ylim) {
      // Polar regime.
      int offset = (prj->n || *yp > 0.0) ? 0 : 1;

      double sigma = prj->w[4] - absy / prj->w[6];

      double s, t;
      if (sigma == 0.0) {
        s = 1e9;
        t = 90.0;

      } else {
        t = 1.0 - sigma*sigma/prj->pv[2];
        if (t < -1.0) {
          s = 0.0;
          t = 0.0;
          istat = 1;
          if (!status) status = PRJERR_BAD_PIX_SET("hpxx2s");
        } else {
          s = 1.0/sigma;
          t = asind(t);
        }
      }
      if (*yp < 0.0) t = -t;

      for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
        if (offset) {
          // Offset the southern polar half-facets for even K.
          int h = (int)floor(*phip / prj->w[6]) + prj->m;
          if (h%2) {
            *thetap -= prj->w[6];
          } else {
            *thetap += prj->w[6];
          }
        }

        // Recall that theta[] holds (x - x_c).
        double r = s * *thetap;

        // Bounds checking.
        if (prj->bounds&2) {
          if (slim <= fabs(r)) {
            istat = 1;
            if (!status) status = PRJERR_BAD_PIX_SET("hpxx2s");
          }
        }

        if (r != 0.0) r -= *thetap;
        *phip  += r;
        *thetap = t;
        *statp  = istat;
      }

    } else {
      // Beyond latitude range.
      for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
        *phip   = 0.0;
        *thetap = 0.0;
        *statp  = 1;
      }
      if (!status) status = PRJERR_BAD_PIX_SET("hpxx2s");
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-12, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("hpxx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int hpxs2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != HPX) {
    if ((status = hpxset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    double xi = prj->w[0] * (*phip) - prj->x0;

    // phi_c for K odd or theta > 0.
    double t;
    t = -180.0 + (2.0*floor((*phip+180.0) * prj->w[7]) + 1.0) * prj->w[6];
    t = prj->w[0] * (*phip - t);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      // y[] is used to hold (phi - phi_c).
      *xp = xi;
      *yp = t;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double sinthe = sind(*thetap);
    double abssin = fabs(sinthe);

    double eta;
    if (abssin <= prj->w[2]) {
      // Equatorial regime.
      eta = prj->w[8] * sinthe - prj->y0;
      for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
        *yp = eta;
        *statp = 0;
      }

    } else {
      // Polar regime.
      int offset = (prj->n || *thetap > 0.0) ? 0 : 1;

      double sigma = sqrt(prj->pv[2]*(1.0 - abssin));
      double xi = sigma - 1.0;

      eta = prj->w[9] * (prj->w[4] - sigma);
      if (*thetap < 0) eta = -eta;
      eta -= prj->y0;

      for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
        if (offset) {
          // Offset the southern polar half-facets for even K.
          int h = (int)floor((*xp + prj->x0) / prj->w[9]) + prj->m;
          if (h%2) {
            *yp -= prj->w[9];
          } else {
            *yp += prj->w[9];
          }
        }

        // Recall that y[] holds (phi - phi_c).
        *xp += *yp * xi;
        *yp = eta;
        *statp = 0;

        // Put the phi = 180 meridian in the expected place.
        if (180.0 < *xp) *xp = 360.0 - *xp;
      }
    }
  }

  return 0;
}

/*============================================================================
*   XPH: HEALPix polar, aka "butterfly" projection.
*
*   Given and/or returned:
*      prj->r0      Reset to 180/pi if 0.
*      prj->phi0    Reset to 0.0 if undefined.
*      prj->theta0  Reset to 0.0 if undefined.
*
*   Returned:
*      prj->flag     XPH
*      prj->code    "XPH"
*      prj->x0      Fiducial offset in x.
*      prj->y0      Fiducial offset in y.
*      prj->w[0]    r0*(pi/180)/sqrt(2)
*      prj->w[1]    (180/pi)/r0/sqrt(2)
*      prj->w[2]    2/3
*      prj->w[3]    tol (= 1e-4)
*      prj->w[4]    sqrt(2/3)*(180/pi)
*      prj->w[5]    90 - tol*sqrt(2/3)*(180/pi)
*      prj->w[6]    sqrt(3/2)*(pi/180)
*      prj->prjx2s  Pointer to xphx2s().
*      prj->prjs2x  Pointer to xphs2x().
*===========================================================================*/

int xphset(struct prjprm *prj)

{
  if (prj == 0x0) return PRJERR_NULL_POINTER;
  if (prj->flag == -XPH) return 0;

  strcpy(prj->code, "XPH");

  strcpy(prj->name, "butterfly");
  prj->category  = HEALPIX;
  prj->pvrange   = 0;
  prj->simplezen = 0;
  prj->equiareal = 1;
  prj->conformal = 0;
  prj->global    = 1;
  prj->divergent = 0;

  if (prj->r0 == 0.0) {
    prj->r0 = R2D;
    prj->w[0] = 1.0;
    prj->w[1] = 1.0;
  } else {
    prj->w[0] = prj->r0*D2R;
    prj->w[1] = R2D/prj->r0;
  }

  prj->w[0] /= sqrt(2.0);
  prj->w[1] /= sqrt(2.0);
  prj->w[2]  = 2.0/3.0;
  prj->w[3]  = 1e-4;
  prj->w[4]  = sqrt(prj->w[2])*R2D;
  prj->w[5]  = 90.0 - prj->w[3]*prj->w[4];
  prj->w[6]  = sqrt(1.5)*D2R;

  prj->prjx2s = xphx2s;
  prj->prjs2x = xphs2x;

  prj->flag = (prj->flag == 1) ? -XPH : XPH;

  return prjoff(prj, 0.0, 90.0);
}

//----------------------------------------------------------------------------

int xphx2s(
  struct prjprm *prj,
  int nx,
  int ny,
  int sxy,
  int spt,
  const double x[],
  const double y[],
  double phi[],
  double theta[],
  int stat[])

{
  const double tol = 1.0e-12;

  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != XPH) {
    if ((status = xphset(prj))) return status;
  }

  int mx, my;
  if (ny > 0) {
    mx = nx;
    my = ny;
  } else {
    mx = 1;
    my = 1;
    ny = nx;
  }

  status = 0;


  // Do x dependence.
  const double *xp = x;
  int rowoff = 0;
  int rowlen = nx*spt;
  for (int ix = 0; ix < nx; ix++, rowoff += spt, xp += sxy) {
    double xr = (*xp + prj->x0)*prj->w[1];

    double *phip = phi + rowoff;
    for (int iy = 0; iy < my; iy++) {
      *phip = xr;
      phip  += rowlen;
    }
  }


  // Do y dependence.
  const double *yp = y;
  double *phip   = phi;
  double *thetap = theta;
  int    *statp  = stat;
  for (int iy = 0; iy < ny; iy++, yp += sxy) {
    double yr = (*yp + prj->y0)*prj->w[1];

    for (int ix = 0; ix < mx; ix++, phip += spt, thetap += spt, statp++) {
      double xr = *phip;

      double xi1, eta1;
      if (xr <= 0.0 && 0.0 < yr) {
        xi1  = -xr - yr;
        eta1 =  xr - yr;
        *phip = -180.0;
      } else if (xr < 0.0 && yr <= 0.0) {
        xi1  =  xr - yr;
        eta1 =  xr + yr;
        *phip = -90.0;
      } else if (0.0 <= xr && yr < 0.0) {
        xi1  =  xr + yr;
        eta1 = -xr + yr;
        *phip = 0.0;
      } else {
        xi1  = -xr + yr;
        eta1 = -xr - yr;
        *phip = 90.0;
      }

      double xi  = xi1  + 45.0;
      double eta = eta1 + 90.0;
      double abseta = fabs(eta);

      if (abseta <= 90.0) {
        if (abseta <= 45.0) {
          // Equatorial regime.
          *phip  += xi;
          *thetap = asind(eta/67.5);
          int istat = 0;

          // Bounds checking.
          if (prj->bounds&2) {
            if (45.0+tol < fabs(xi1)) {
              istat = 1;
              if (!status) status = PRJERR_BAD_PIX_SET("xphx2s");
            }
          }

          *statp = istat;

        } else {
          // Polar regime.
          double sigma = (90.0 - abseta) / 45.0;

          // Ensure an exact result for points on the boundary.
          if (xr == 0.0) {
            if (yr <= 0.0) {
              *phip = 0.0;
            } else {
              *phip = 180.0;
            }
          } else if (yr == 0.0) {
            if (xr < 0.0) {
              *phip = -90.0;
            } else {
              *phip =  90.0;
            }
          } else {
            *phip += 45.0 + xi1/sigma;
          }

          if (sigma < prj->w[3]) {
            *thetap = 90.0 - sigma*prj->w[4];
          } else {
            *thetap = asind(1.0 - sigma*sigma/3.0);
          }
          if (eta < 0.0) *thetap = -(*thetap);

          // Bounds checking.
          int istat = 0;
          if (prj->bounds&2) {
            if (eta < -45.0 && eta+90.0+tol < fabs(xi1)) {
              istat = 1;
              if (!status) status = PRJERR_BAD_PIX_SET("xphx2s");
            }
          }

          *statp = istat;
        }

      } else {
        // Beyond latitude range.
        *phip   = 0.0;
        *thetap = 0.0;
        *statp  = 1;
        if (!status) status = PRJERR_BAD_PIX_SET("xphx2s");
      }
    }
  }


  // Do bounds checking on the native coordinates.
  if (prj->bounds&4 && prjbchk(1.0e-12, nx, my, spt, phi, theta, stat)) {
    if (!status) status = PRJERR_BAD_PIX_SET("xphx2s");
  }

  return status;
}

//----------------------------------------------------------------------------

int xphs2x(
  struct prjprm *prj,
  int nphi,
  int ntheta,
  int spt,
  int sxy,
  const double phi[],
  const double theta[],
  double x[],
  double y[],
  int stat[])

{
  // Initialize.
  if (prj == 0x0) return PRJERR_NULL_POINTER;

  int status;
  if (abs(prj->flag) != XPH) {
    if ((status = xphset(prj))) return status;
  }

  int mphi, mtheta;
  if (ntheta > 0) {
    mphi   = nphi;
    mtheta = ntheta;
  } else {
    mphi   = 1;
    mtheta = 1;
    ntheta = nphi;
  }


  // Do phi dependence.
  const double *phip = phi;
  int rowoff = 0;
  int rowlen = nphi*sxy;

  double chi, psi;
  for (int iphi = 0; iphi < nphi; iphi++, rowoff += sxy, phip += spt) {
    chi = *phip;
    if (180.0 <= fabs(chi)) {
      chi = fmod(chi, 360.0);
      if (chi < -180.0) {
        chi += 360.0;
      } else if (180.0 <= chi) {
        chi -= 360.0;
      }
    }

    // phi is also recomputed from chi to avoid rounding problems.
    chi += 180.0;
    psi = fmod(chi, 90.0);

    double *xp = x + rowoff;
    double *yp = y + rowoff;
    for (int itheta = 0; itheta < mtheta; itheta++) {
      // y[] is used to hold phi (rounded).
      *xp = psi;
      *yp = chi - 180.0;
      xp += rowlen;
      yp += rowlen;
    }
  }


  // Do theta dependence.
  const double *thetap = theta;
  double *xp = x;
  double *yp = y;
  int *statp = stat;
  for (int itheta = 0; itheta < ntheta; itheta++, thetap += spt) {
    double sinthe = sind(*thetap);
    double abssin = fabs(sinthe);

    for (int iphi = 0; iphi < mphi; iphi++, xp += sxy, yp += sxy, statp++) {
      double xi, eta;
      if (abssin <= prj->w[2]) {
        // Equatorial regime.
        xi  = *xp;
        eta = 67.5 * sinthe;

      } else {
        // Polar regime.
        double sigma;
        if (*thetap < prj->w[5]) {
          sigma = sqrt(3.0*(1.0 - abssin));
        } else {
          sigma = (90.0 - *thetap)*prj->w[6];
        }

        xi  = 45.0 + (*xp - 45.0)*sigma;
        eta = 45.0 * (2.0 - sigma);
        if (*thetap < 0.0) eta = -eta;
      }

      xi  -= 45.0;
      eta -= 90.0;

      // Recall that y[] holds phi.
      if (*yp < -90.0) {
        *xp = prj->w[0]*(-xi + eta) - prj->x0;
        *yp = prj->w[0]*(-xi - eta) - prj->y0;

      } else if (*yp <  0.0) {
        *xp = prj->w[0]*(+xi + eta) - prj->x0;
        *yp = prj->w[0]*(-xi + eta) - prj->y0;

      } else if (*yp < 90.0) {
        *xp = prj->w[0]*( xi - eta) - prj->x0;
        *yp = prj->w[0]*( xi + eta) - prj->y0;

      } else {
        *xp = prj->w[0]*(-xi - eta) - prj->x0;
        *yp = prj->w[0]*( xi - eta) - prj->y0;
      }

      *statp = 0;
    }
  }

  return 0;
}
