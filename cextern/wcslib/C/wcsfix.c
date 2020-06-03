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
  $Id: wcsfix.c,v 7.3 2020/06/03 03:37:02 mcalabre Exp $
*===========================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcserr.h"
#include "wcsmath.h"
#include "wcstrig.h"
#include "wcsutil.h"
#include "lin.h"
#include "sph.h"
#include "wcs.h"
#include "wcsunits.h"
#include "wcsfix.h"

extern const int WCSSET;

/* Maximum number of coordinate axes that can be handled. */
#define NMAX 16

/* Map status return value to message. */
const char *wcsfix_errmsg[] = {
  "Success",
  "Null wcsprm pointer passed",
  "Memory allocation failed",
  "Linear transformation matrix is singular",
  "Inconsistent or unrecognized coordinate axis types",
  "Invalid parameter value",
  "Invalid coordinate transformation parameters",
  "Ill-conditioned coordinate transformation parameters",
  "All of the corner pixel coordinates are invalid",
  "Could not determine reference pixel coordinate",
  "Could not determine reference pixel value"};

/* Map error returns for lower-level routines. */
const int fix_linerr[] = {
  FIXERR_SUCCESS,		/*  0: LINERR_SUCCESS         */
  FIXERR_NULL_POINTER,		/*  1: LINERR_NULL_POINTER    */
  FIXERR_MEMORY,		/*  2: LINERR_MEMORY          */
  FIXERR_SINGULAR_MTX,		/*  3: LINERR_SINGULAR_MTX    */
  FIXERR_BAD_PARAM,		/*  4: LINERR_DISTORT_INIT    */
  FIXERR_NO_REF_PIX_COORD,	/*  5: LINERR_DISTORT         */
  FIXERR_NO_REF_PIX_VAL		/*  6: LINERR_DEDISTORT       */
};

const int fix_wcserr[] = {
  FIXERR_SUCCESS,		/*  0: WCSERR_SUCCESS         */
  FIXERR_NULL_POINTER,		/*  1: WCSERR_NULL_POINTER    */
  FIXERR_MEMORY,		/*  2: WCSERR_MEMORY          */
  FIXERR_SINGULAR_MTX,		/*  3: WCSERR_SINGULAR_MTX    */
  FIXERR_BAD_CTYPE,		/*  4: WCSERR_BAD_CTYPE       */
  FIXERR_BAD_PARAM,		/*  5: WCSERR_BAD_PARAM       */
  FIXERR_BAD_COORD_TRANS,	/*  6: WCSERR_BAD_COORD_TRANS */
  FIXERR_ILL_COORD_TRANS,	/*  7: WCSERR_ILL_COORD_TRANS */
  FIXERR_BAD_CORNER_PIX,	/*  8: WCSERR_BAD_PIX         */
  FIXERR_NO_REF_PIX_VAL,	/*  9: WCSERR_BAD_WORLD       */
  FIXERR_NO_REF_PIX_VAL 	/* 10: WCSERR_BAD_WORLD_COORD */
				/*     ...others not used     */
};

/* Convenience macro for invoking wcserr_set(). */
#define WCSFIX_ERRMSG(status) WCSERR_SET(status), wcsfix_errmsg[status]

/*--------------------------------------------------------------------------*/

int wcsfix(int ctrl, const int naxis[], struct wcsprm *wcs, int stat[])

{
  int status = 0;

  if ((stat[CDFIX] = cdfix(wcs)) > 0) {
    status = 1;
  }

  if ((stat[DATFIX] = datfix(wcs)) > 0) {
    status = 1;
  }

  if ((stat[OBSFIX] = obsfix(0, wcs)) > 0) {
    status = 1;
  }

  if ((stat[UNITFIX] = unitfix(ctrl, wcs)) > 0) {
    status = 1;
  }

  if ((stat[SPCFIX] = spcfix(wcs)) > 0) {
    status = 1;
  }

  if ((stat[CELFIX] = celfix(wcs)) > 0) {
    status = 1;
  }

  if ((stat[CYLFIX] = cylfix(naxis, wcs)) > 0) {
    status = 1;
  }

  return status;
}

/*--------------------------------------------------------------------------*/

int wcsfixi(
  int ctrl,
  const int naxis[],
  struct wcsprm *wcs,
  int stat[],
  struct wcserr info[])

{
  int ifix, status = 0;
  struct wcserr err;

  /* Handling the status values returned from the sub-fixers is trickier than
  it might seem, especially considering that wcs->err may contain an error
  status on input which should be preserved if no translation errors occur.
  The simplest way seems to be to save a copy of wcs->err and clear it before
  each sub-fixer.  The last real error to occur, excluding informative
  messages, is the one returned.

  To get informative messages from spcfix() it must precede celfix() and
  cylfix().  The latter call wcsset() which also translates AIPS-convention
  spectral axes. */
  wcserr_copy(wcs->err, &err);

  for (ifix = CDFIX; ifix < NWCSFIX; ifix++) {
    /* Clear (delete) wcs->err. */
    wcserr_clear(&(wcs->err));

    switch (ifix) {
    case CDFIX:
      stat[ifix] = cdfix(wcs);
      break;
    case DATFIX:
      stat[ifix] = datfix(wcs);
      break;
    case OBSFIX:
      stat[ifix] = obsfix(0, wcs);
      break;
    case UNITFIX:
      stat[ifix] = unitfix(ctrl, wcs);
      break;
    case SPCFIX:
      stat[ifix] = spcfix(wcs);
      break;
    case CELFIX:
      stat[ifix] = celfix(wcs);
      break;
    case CYLFIX:
      stat[ifix] = cylfix(naxis, wcs);
      break;
    default:
      continue;
    }

    if (stat[ifix] == FIXERR_NO_CHANGE) {
      /* No change => no message. */
      wcserr_copy(0x0, info+ifix);

    } else if (stat[ifix] == 0) {
      /* Successful translation, but there may be an informative message. */
      if (wcs->err && wcs->err->status < 0) {
        wcserr_copy(wcs->err, info+ifix);
      } else {
        wcserr_copy(0x0, info+ifix);
      }

    } else {
      /* An informative message or error message. */
      wcserr_copy(wcs->err, info+ifix);

      if ((status = (stat[ifix] > 0))) {
        /* It was an error, replace the previous one. */
        wcserr_copy(wcs->err, &err);
      }
    }
  }

  /* Restore the last error to occur. */
  if (err.status) {
    wcserr_copy(&err, wcs->err);
  } else {
    wcserr_clear(&(wcs->err));
  }

  return status;
}

/*--------------------------------------------------------------------------*/

int cdfix(struct wcsprm *wcs)

{
  int  i, k, naxis, status = FIXERR_NO_CHANGE;
  double *cd;

  if (wcs == 0x0) return FIXERR_NULL_POINTER;

  if ((wcs->altlin & 1) || !(wcs->altlin & 2)) {
    /* Either we have PCi_ja or there are no CDi_ja. */
    return FIXERR_NO_CHANGE;
  }

  naxis = wcs->naxis;
  status = FIXERR_NO_CHANGE;
  for (i = 0; i < naxis; i++) {
    /* Row of zeros? */
    cd = wcs->cd + i * naxis;
    for (k = 0; k < naxis; k++, cd++) {
      if (*cd != 0.0) goto next;
    }

    /* Column of zeros? */
    cd = wcs->cd + i;
    for (k = 0; k < naxis; k++, cd += naxis) {
      if (*cd != 0.0) goto next;
    }

    cd = wcs->cd + i * (naxis + 1);
    *cd = 1.0;
    status = 0;

next: ;
  }

  return status;
}

/*--------------------------------------------------------------------------*/

static int parse_date(const char *buf, int *hour, int *minute, double *sec)

{
  char ctmp[72];

  if (sscanf(buf, "%2d:%2d:%s", hour, minute, ctmp) < 3 ||
      wcsutil_str2double(ctmp, sec)) {
    return 1;
  }

  return 0;
}


static void write_date(char *buf, int hour, int minute, double sec)

{
  char ctmp[32];

  wcsutil_double2str(ctmp, "%04.1f", sec);
  sprintf(buf, "T%.2d:%.2d:%s", hour, minute, ctmp);
}


static char *newline(char **cp)

{
  size_t k;

  if ((k = strlen(*cp))) {
    *cp += k;
    strcat(*cp, ".\n");
    *cp += 2;
  }

  return *cp;
}


int datfix(struct wcsprm *wcs)

{
  static const char *function = "datfix";

  /* MJD of J2000.0 and B1900.0. */
  const double mjd2000 = 51544.5;
  const double mjd1900 = 15019.81352;

  /* Days per Julian year and per tropical year. */
  const double djy = 365.25;
  const double dty = 365.242198781;

  const char *dateid;
  char *cp, *date, infomsg[512], orig_date[72];
  int  day, dd, hour = 0, i, jd, minute = 0, month, msec, n4, status, year;
  double bepoch, jepoch, mjd[2], mjdsum, mjdtmp, sec = 0.0, t, *wcsmjd;
  struct wcserr **err;

  if (wcs == 0x0) return FIXERR_NULL_POINTER;
  err = &(wcs->err);

  cp = infomsg;
  *cp = '\0';
  status = FIXERR_NO_CHANGE;

  for (i = 0; i < 5; i++) {
    /* MJDREF is split into integer and fractional parts, wheres MJDOBS and
       the rest are a single value. */
    if (i == 0) {
      /* Note, DATEREF and MJDREF, not DATE-REF and MJD-REF (sigh). */
      dateid = "REF";
      date   = wcs->dateref;
      wcsmjd = wcs->mjdref;
    } else if (i == 1) {
      dateid = "-OBS";
      date   = wcs->dateobs;
      wcsmjd = &(wcs->mjdobs);
    } else if (i == 2) {
      dateid = "-BEG";
      date   = wcs->datebeg;
      wcsmjd = &(wcs->mjdbeg);
    } else if (i == 3) {
      dateid = "-AVG";
      date   = wcs->dateavg;
      wcsmjd = &(wcs->mjdavg);
    } else if (i == 4) {
      dateid = "-END";
      date   = wcs->dateend;
      wcsmjd = &(wcs->mjdend);
    }

    strncpy(orig_date, date, 72);

    if (date[0] == '\0') {
      /* Fill in DATE from MJD if possible. */

      if (i == 1 && undefined(*wcsmjd)) {
        /* See if we have jepoch or bepoch. */
        if (!undefined(wcs->jepoch)) {
          *wcsmjd = mjd2000 + (wcs->jepoch - 2000.0)*djy;
          sprintf(newline(&cp), "Set MJD-OBS to %.6f from JEPOCH", *wcsmjd);

        } else if (!undefined(wcs->bepoch)) {
          *wcsmjd = mjd1900 + (wcs->bepoch - 1900.0)*dty;
          sprintf(newline(&cp), "Set MJD-OBS to %.6f from BEPOCH", *wcsmjd);
        }
      }

      if (undefined(*wcsmjd)) {
        /* No date information was provided. */

      } else {
        /* Calendar date from MJD, with allowance for MJD < 0. */
        if (i == 0) {
          /* MJDREF is already split into integer and fractional parts. */
          mjd[0] = wcsmjd[0];
          mjd[1] = wcsmjd[1];
          if (1.0 < mjd[1]) {
            /* Ensure the fractional part lies between 0 and +1. */
            t = floor(mjd[1]);
            mjd[0] += t;
            mjd[1] -= t;
          }
        } else {
          /* Split it into integer and fractional parts. */
          mjd[0] = floor(*wcsmjd);
          mjd[1] = *wcsmjd - mjd[0];
        }

        jd = 2400001 + (int)mjd[0];

        n4 =  4*(jd + ((2*((4*jd - 17918)/146097)*3)/4 + 1)/2 - 37);
        dd = 10*(((n4-237)%1461)/4) + 5;

        year  = n4/1461 - 4712;
        month = (2 + dd/306)%12 + 1;
        day   = (dd%306)/10 + 1;
        sprintf(date, "%.4d-%.2d-%.2d", year, month, day);

        /* Write time part only if non-zero. */
        if (0.0 < (t = mjd[1])) {
          t *= 24.0;
          hour = (int)t;
          t = 60.0 * (t - hour);
          minute = (int)t;
          sec    = 60.0 * (t - minute);

          /* Round to 1ms. */
          dd = 60000*(60*hour + minute) + (int)(1000*(sec+0.0005));
          hour = dd / 3600000;
          dd -= 3600000 * hour;
          minute = dd / 60000;
          msec = dd - 60000 * minute;
          sprintf(date+10, "T%.2d:%.2d:%.2d", hour, minute, msec/1000);

          /* Write fractions of a second only if non-zero. */
          if (msec%1000) {
            sprintf(date+19, ".%.3d", msec%1000);
          }
        }
      }

    } else {
      if (strlen(date) < 8) {
        /* Can't be a valid date. */
        status = FIXERR_BAD_PARAM;
        sprintf(newline(&cp), "Invalid DATE%s format '%s' is too short",
          dateid, date);
        continue;
      }

      /* Identify the date format. */
      if (date[4] == '-' && date[7] == '-') {
        /* Standard year-2000 form: CCYY-MM-DD[Thh:mm:ss[.sss...]] */
        if (sscanf(date, "%4d-%2d-%2d", &year, &month, &day) < 3) {
          status = FIXERR_BAD_PARAM;
          sprintf(newline(&cp), "Invalid DATE%s format '%s'", dateid, date);
          continue;
        }

        if (date[10] == 'T') {
          if (parse_date(date+11, &hour, &minute, &sec)) {
            status = FIXERR_BAD_PARAM;
            sprintf(newline(&cp), "Invalid time in DATE%s '%s'", dateid,
              date+11);
            continue;
          }
        } else if (date[10] == ' ') {
          hour = 0;
          minute = 0;
          sec = 0.0;
          if (parse_date(date+11, &hour, &minute, &sec)) {
            write_date(date+10, hour, minute, sec);
          } else {
            date[10] = 'T';
          }
        }

      } else if (date[4] == '/' && date[7] == '/') {
        /* Also allow CCYY/MM/DD[Thh:mm:ss[.sss...]] */
        if (sscanf(date, "%4d/%2d/%2d", &year, &month, &day) < 3) {
          status = FIXERR_BAD_PARAM;
          sprintf(newline(&cp), "Invalid DATE%s format '%s'", dateid, date);
          continue;
        }

        if (date[10] == 'T') {
          if (parse_date(date+11, &hour, &minute, &sec)) {
            status = FIXERR_BAD_PARAM;
            sprintf(newline(&cp), "Invalid time in DATE%s '%s'", dateid,
              date+11);
            continue;
          }
        } else if (date[10] == ' ') {
          hour = 0;
          minute = 0;
          sec = 0.0;
          if (parse_date(date+11, &hour, &minute, &sec)) {
            write_date(date+10, hour, minute, sec);
          } else {
            date[10] = 'T';
          }
        }

        /* Looks ok, fix it up. */
        date[4]  = '-';
        date[7]  = '-';

      } else {
        if (i == 1 && date[2] == '/' && date[5] == '/') {
          /* Old format DATE-OBS date: DD/MM/YY, also allowing DD/MM/CCYY. */
          if (sscanf(date, "%2d/%2d/%4d", &day, &month, &year) < 3) {
            status = FIXERR_BAD_PARAM;
            sprintf(newline(&cp), "Invalid DATE%s format '%s'", dateid,
              date);
            continue;
          }

        } else if (i == 1 && date[2] == '-' && date[5] == '-') {
          /* Also recognize DD-MM-YY and DD-MM-CCYY */
          if (sscanf(date, "%2d-%2d-%4d", &day, &month, &year) < 3) {
            status = FIXERR_BAD_PARAM;
            sprintf(newline(&cp), "Invalid DATE%s format '%s'", dateid,
              date);
            continue;
          }

        } else {
          /* Not a valid date format. */
          status = FIXERR_BAD_PARAM;
          sprintf(newline(&cp), "Invalid DATE%s format '%s'", dateid, date);
          continue;
        }

        if (year < 100) year += 1900;

        /* Doesn't have a time. */
        sprintf(date, "%.4d-%.2d-%.2d", year, month, day);
      }

      /* Compute MJD. */
      mjd[0] = (double)((1461*(year - (12-month)/10 + 4712))/4
               + (306*((month+9)%12) + 5)/10
               - (3*((year - (12-month)/10 + 4900)/100))/4
               + day - 2399904);
      mjd[1] = (hour + (minute + sec/60.0)/60.0)/24.0;
      mjdsum = mjd[0] + mjd[1];

      if (undefined(*wcsmjd)) {
        if (i == 0) {
          wcsmjd[0] = mjd[0];
          wcsmjd[1] = mjd[1];
        } else {
          *wcsmjd = mjdsum;
        }
        sprintf(newline(&cp), "Set MJD%s to %.6f from DATE%s", dateid,
          mjdsum, dateid);

      } else {
        /* Check for consistency. */
        if (i == 0) {
          mjdtmp = wcsmjd[0] + wcsmjd[1];
        } else {
          mjdtmp = *wcsmjd;
        }

        if (0.001 < fabs(mjdsum - mjdtmp)) {
          status = FIXERR_BAD_PARAM;
          sprintf(newline(&cp),
            "Invalid parameter values: MJD%s and DATE%s are inconsistent",
            dateid, dateid);
        }
      }

      if (i == 1) {
        if (!undefined(wcs->jepoch)) {
          /* Check consistency of JEPOCH. */
          jepoch = 2000.0 + (*wcsmjd - mjd2000) / djy;

          if (0.000002 < fabs(jepoch - wcs->jepoch)) {
            /* Informational only, no error. */
            sprintf(newline(&cp), "JEPOCH is inconsistent with DATE-OBS");
          }
        }

        if (!undefined(wcs->bepoch)) {
          /* Check consistency of BEPOCH. */
          bepoch = 1900.0 + (*wcsmjd - mjd1900) / dty;

          if (0.000002 < fabs(bepoch - wcs->bepoch)) {
            /* Informational only, no error. */
            sprintf(newline(&cp), "BEPOCH is inconsistent with DATE-OBS");
          }
        }
      }
    }

    if (strncmp(orig_date, date, 72)) {
      if (orig_date[0] == '\0') {
        sprintf(newline(&cp), "Set DATE%s to '%s' from MJD%s", dateid, date,
          dateid);
      } else {
        sprintf(newline(&cp), "Changed DATE%s from '%s' to '%s'", dateid,
          orig_date, date);
      }

      if (status == FIXERR_NO_CHANGE) status = 0;
    }
  }

  if (*infomsg) {
    wcserr_set(WCSERR_SET(FIXERR_DATE_FIX), infomsg);
  }

  return status;
}

/*--------------------------------------------------------------------------*/

int obsfix(int ctrl, struct wcsprm *wcs)

{
  static const char *function = "obsfix";

  /* IAU(1976) ellipsoid (as prescribed by WCS Paper VII). */
  const double a = 6378140.0, f = 1.0 / 298.2577;
  const double e2 = (2.0 - f)*f;

  char   *cp, infomsg[256];
  int    havelbh = 7, havexyz = 7, i, status;
  size_t k;
  double coslat, coslng, d, hgt, lat, lng, n, r2, rho, sinlat, sinlng, x, y,
         z, zeta;
  struct wcserr **err;

  if (wcs == 0x0) return FIXERR_NULL_POINTER;
  err = &(wcs->err);

  /* Set masks for checking partially-defined coordinate triplets. */
  havexyz -= 1*undefined(wcs->obsgeo[0]);
  havexyz -= 2*undefined(wcs->obsgeo[1]);
  havexyz -= 4*undefined(wcs->obsgeo[2]);
  havelbh -= 1*undefined(wcs->obsgeo[3]);
  havelbh -= 2*undefined(wcs->obsgeo[4]);
  havelbh -= 4*undefined(wcs->obsgeo[5]);

  if (ctrl == 2) {
    /* Make no changes. */
    if (0 < havexyz && havexyz < 7) {
      return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
        "Partially undefined Cartesian coordinate triplet");
    }

    if (0 < havelbh && havelbh < 7) {
      return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
        "Partially undefined Geodetic coordinate triplet");
    }

    if (havexyz == 0 || havelbh == 0) {
      return FIXERR_NO_CHANGE;
    }
  }

  if (havexyz == 0 && havelbh == 0) {
    return FIXERR_NO_CHANGE;
  }


  infomsg[0] = '\0';
  status = FIXERR_NO_CHANGE;

  if (havelbh == 7) {
    /* Compute (x,y,z) from (lng,lat,hgt). */
    sincosd(wcs->obsgeo[3], &sinlng, &coslng);
    sincosd(wcs->obsgeo[4], &sinlat, &coslat);
    n = a / sqrt(1.0 - e2*sinlat*sinlat);
    rho = n + wcs->obsgeo[5];

    x = rho*coslng*coslat;
    y = rho*sinlng*coslat;
    z = (rho - n*e2)*sinlat;

    if (havexyz < 7) {
      /* One or more of the Cartesian elements was undefined. */
      status = 0;
      cp = infomsg;

      if (ctrl == 1 || !(havexyz & 1)) {
        wcs->obsgeo[0] = x;
        sprintf(cp, "%s OBSGEO-X to %12.3f from OBSGEO-[LBH]",
          (havexyz & 1) ? "Reset" : "Set", x);
      }

      if (ctrl == 1 || !(havexyz & 2)) {
        wcs->obsgeo[1] = y;

        if ((k = strlen(cp))) {
          strcat(cp+k, ".\n");
          cp += k + 2;
        }

        sprintf(cp, "%s OBSGEO-Y to %12.3f from OBSGEO-[LBH]",
          (havexyz & 2) ? "Reset" : "Set", y);
      }

      if (ctrl == 1 || !(havexyz & 4)) {
        wcs->obsgeo[2] = z;
        if ((k = strlen(cp))) {
          strcat(cp+k, ".\n");
          cp += k + 2;
        }

        sprintf(cp, "%s OBSGEO-Z to %12.3f from OBSGEO-[LBH]",
          (havexyz & 4) ? "Reset" : "Set", z);
      }

      wcserr_set(WCSERR_SET(FIXERR_OBSGEO_FIX), infomsg);

      if (havexyz == 0) {
        /* Skip the consistency check. */
        return status;
      }
    }

  } else if (havexyz == 7) {
    /* Compute (lng,lat,hgt) from (x,y,z). */
    x = wcs->obsgeo[0];
    y = wcs->obsgeo[1];
    z = wcs->obsgeo[2];
    r2 = x*x + y*y;

    /* Iterate over the value of zeta. */
    zeta = z;
    for (i = 0; i < 4; i++) {
      rho = sqrt(r2 + zeta*zeta);
      sinlat = zeta / rho;
      n = a / sqrt(1.0 - e2*sinlat*sinlat);

      zeta = z / (1.0 - n*e2/rho);
    }

    lng = atan2d(y, x);
    lat = asind(sinlat);
    hgt = rho - n;

    if (havelbh < 7) {
      /* One or more of the Geodetic elements was undefined. */
      status = 0;
      cp = infomsg;

      if (ctrl == 1 || !(havelbh & 1)) {
        wcs->obsgeo[3] = lng;
        sprintf(cp, "%s OBSGEO-L to %12.6f from OBSGEO-[XYZ]",
          (havelbh & 1) ? "Reset" : "Set", lng);
      }

      if (ctrl == 1 || !(havelbh & 2)) {
        wcs->obsgeo[4] = lat;
        if ((k = strlen(cp))) {
          strcat(cp+k, ".\n");
          cp += k + 2;
        }

        sprintf(cp, "%s OBSGEO-B to %12.6f from OBSGEO-[XYZ]",
          (havelbh & 2) ? "Reset" : "Set", lat);
      }

      if (ctrl == 1 || !(havelbh & 4)) {
        wcs->obsgeo[5] = hgt;
        if ((k = strlen(cp))) {
          strcat(cp+k, ".\n");
          cp += k + 2;
        }

        sprintf(cp, "%s OBSGEO-H to %12.3f from OBSGEO-[XYZ]",
          (havelbh & 4) ? "Reset" : "Set", hgt);
      }

      wcserr_set(WCSERR_SET(FIXERR_OBSGEO_FIX), infomsg);

      if (havelbh == 0) {
        /* Skip the consistency check. */
        return status;
      }
    }

    /* Compute (x,y,z) from (lng,lat,hgt) for consistency checking. */
    sincosd(wcs->obsgeo[3], &sinlng, &coslng);
    sincosd(wcs->obsgeo[4], &sinlat, &coslat);
    n = a / sqrt(1.0 - e2*sinlat*sinlat);
    rho = n + wcs->obsgeo[5];

    x = rho*coslng*coslat;
    y = rho*sinlng*coslat;
    z = (rho - n*e2)*sinlat;

  } else {
    return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
      "Observatory coordinates incomplete");
  }


  /* Check consistency. */
  r2 = 0.0;
  d = wcs->obsgeo[0] - x;
  r2 += d*d;
  d = wcs->obsgeo[1] - y;
  r2 += d*d;
  d = wcs->obsgeo[2] - z;
  r2 += d*d;

  if (1.0 < r2) {
    d = sqrt(r2);
    return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
      "Observatory coordinates inconsistent by %.1f metres", d);
  }

  return status;
}


/*--------------------------------------------------------------------------*/

int unitfix(int ctrl, struct wcsprm *wcs)

{
  const char *function = "unitfix";

  char   orig_unit[72], msg[512], msgtmp[192];
  int    i, result, status = FIXERR_NO_CHANGE;
  size_t msglen;
  struct wcserr **err;

  if (wcs == 0x0) return FIXERR_NULL_POINTER;
  err = &(wcs->err);

  strncpy(msg, "Changed units:", 512);

  for (i = 0; i < wcs->naxis; i++) {
    strncpy(orig_unit, wcs->cunit[i], 71);
    result = wcsutrne(ctrl, wcs->cunit[i], &(wcs->err));
    if (result == 0 || result == 12) {
      msglen = strlen(msg);
      if (msglen < 511) {
        wcsutil_null_fill(72, orig_unit);
        sprintf(msgtmp, "\n  '%s' -> '%s',", orig_unit, wcs->cunit[i]);
        strncpy(msg+msglen, msgtmp, 511-msglen);
        status = FIXERR_UNITS_ALIAS;
      }
    }
  }

  if (status == FIXERR_UNITS_ALIAS) {
    /* Chop off the trailing ", ". */
    msglen = strlen(msg) - 2;
    msg[msglen] = '\0';
    wcserr_set(WCSERR_SET(FIXERR_UNITS_ALIAS), msg);

    status = 0;
  }

  return status;
}

/*--------------------------------------------------------------------------*/

int spcfix(struct wcsprm *wcs)

{
  static const char *function = "spcfix";

  char ctype[9], specsys[9];
  int  i, status;
  struct wcserr **err;

  if (wcs == 0x0) return FIXERR_NULL_POINTER;
  err = &(wcs->err);

  for (i = 0; i < wcs->naxis; i++) {
    /* Translate an AIPS-convention spectral type if present. */
    status = spcaips(wcs->ctype[i], wcs->velref, ctype, specsys);
    if (status == 0) {
      /* An AIPS type was found but it may match what we already have. */
      status = FIXERR_NO_CHANGE;

      /* Was specsys translated? */
      if (wcs->specsys[0] == '\0' && *specsys) {
        strncpy(wcs->specsys, specsys, 9);
        wcserr_set(WCSERR_SET(FIXERR_SPC_UPDATE),
          "Changed SPECSYS to '%s'", specsys);
        status = 0;
      }

      /* Was ctype translated?  Have to null-fill for comparing them. */
      wcsutil_null_fill(9, wcs->ctype[i]);
      if (strncmp(wcs->ctype[i], ctype, 9)) {
        /* ctype was translated... */
        if (status == 0) {
          /* ...and specsys was also. */
          wcserr_set(WCSERR_SET(FIXERR_SPC_UPDATE),
            "Changed CTYPE%d from '%s' to '%s', and SPECSYS to '%s' "
            "(VELREF=%d)", i+1, wcs->ctype[i], ctype, wcs->specsys,
            wcs->velref);
        } else {
          wcserr_set(WCSERR_SET(FIXERR_SPC_UPDATE),
            "Changed CTYPE%d from '%s' to '%s' (VELREF=%d)", i+1,
            wcs->ctype[i], ctype, wcs->velref);
          status = 0;
        }

        strncpy(wcs->ctype[i], ctype, 9);
      }

      /* Tidy up. */
      if (status == 0) {
        wcsutil_null_fill(72, wcs->ctype[i]);
        wcsutil_null_fill(72, wcs->specsys);
      }

      /* No need to check for others, wcsset() will fail if so. */
      return status;

    } else if (status == SPCERR_BAD_SPEC_PARAMS) {
      /* An AIPS spectral type was found but with invalid velref. */
      return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
        "Invalid parameter value: velref = %d", wcs->velref);
    }
  }

  return FIXERR_NO_CHANGE;
}

/*--------------------------------------------------------------------------*/

int celfix(struct wcsprm *wcs)

{
  static const char *function = "celfix";

  int k, status;
  struct celprm *wcscel = &(wcs->cel);
  struct prjprm *wcsprj = &(wcscel->prj);
  struct wcserr **err;

  if (wcs == 0x0) return FIXERR_NULL_POINTER;
  err = &(wcs->err);

  /* Initialize if required. */
  if (wcs->flag != WCSSET) {
    if ((status = wcsset(wcs))) return fix_wcserr[status];
  }

  /* Was an NCP or GLS projection code translated? */
  if (wcs->lat >= 0) {
    /* Check ctype. */
    if (strcmp(wcs->ctype[wcs->lat]+5, "NCP") == 0) {
      strcpy(wcs->ctype[wcs->lng]+5, "SIN");
      strcpy(wcs->ctype[wcs->lat]+5, "SIN");

      if (wcs->npvmax < wcs->npv + 2) {
        /* Allocate space for two more PVi_ma keyvalues. */
        if (wcs->m_flag == WCSSET && wcs->pv == wcs->m_pv) {
          if (!(wcs->pv = calloc(wcs->npv+2, sizeof(struct pvcard)))) {
            wcs->pv = wcs->m_pv;
            return wcserr_set(WCSFIX_ERRMSG(FIXERR_MEMORY));
          }

          wcs->npvmax = wcs->npv + 2;
          wcs->m_flag = WCSSET;

          for (k = 0; k < wcs->npv; k++) {
            wcs->pv[k] = wcs->m_pv[k];
          }

          if (wcs->m_pv) free(wcs->m_pv);
          wcs->m_pv = wcs->pv;

        } else {
          return wcserr_set(WCSFIX_ERRMSG(FIXERR_MEMORY));
        }
      }

      wcs->pv[wcs->npv].i = wcs->lat + 1;
      wcs->pv[wcs->npv].m = 1;
      wcs->pv[wcs->npv].value = wcsprj->pv[1];
      (wcs->npv)++;

      wcs->pv[wcs->npv].i = wcs->lat + 1;
      wcs->pv[wcs->npv].m = 2;
      wcs->pv[wcs->npv].value = wcsprj->pv[2];
      (wcs->npv)++;

      return 0;

    } else if (strcmp(wcs->ctype[wcs->lat]+5, "GLS") == 0) {
      strcpy(wcs->ctype[wcs->lng]+5, "SFL");
      strcpy(wcs->ctype[wcs->lat]+5, "SFL");

      if (wcs->crval[wcs->lng] != 0.0 || wcs->crval[wcs->lat] != 0.0) {
        /* In the AIPS convention, setting the reference longitude and
         * latitude for GLS does not create an oblique graticule.  A non-zero
         * reference longitude introduces an offset in longitude in the normal
         * way, whereas a non-zero reference latitude simply translates the
         * reference point (i.e. the map as a whole) to that latitude.  This
         * might be effected by adjusting CRPIXja but that is complicated by
         * the linear transformation and instead is accomplished here by
         * setting theta_0. */
        if (wcs->npvmax < wcs->npv + 3) {
          /* Allocate space for three more PVi_ma keyvalues. */
          if (wcs->m_flag == WCSSET && wcs->pv == wcs->m_pv) {
            if (!(wcs->pv = calloc(wcs->npv+3, sizeof(struct pvcard)))) {
              wcs->pv = wcs->m_pv;
              return wcserr_set(WCSFIX_ERRMSG(FIXERR_MEMORY));
            }

            wcs->npvmax = wcs->npv + 3;
            wcs->m_flag = WCSSET;

            for (k = 0; k < wcs->npv; k++) {
              wcs->pv[k] = wcs->m_pv[k];
            }

            if (wcs->m_pv) free(wcs->m_pv);
            wcs->m_pv = wcs->pv;

          } else {
            return wcserr_set(WCSFIX_ERRMSG(FIXERR_MEMORY));
          }
        }

        wcs->pv[wcs->npv].i = wcs->lng + 1;
        wcs->pv[wcs->npv].m = 0;
        wcs->pv[wcs->npv].value = 1.0;
        (wcs->npv)++;

        /* Note that the reference longitude is still zero. */
        wcs->pv[wcs->npv].i = wcs->lng + 1;
        wcs->pv[wcs->npv].m = 1;
        wcs->pv[wcs->npv].value = 0.0;
        (wcs->npv)++;

        wcs->pv[wcs->npv].i = wcs->lng + 1;
        wcs->pv[wcs->npv].m = 2;
        wcs->pv[wcs->npv].value = wcs->crval[wcs->lat];
        (wcs->npv)++;
      }

      return 0;
    }
  }

  return FIXERR_NO_CHANGE;
}

/*--------------------------------------------------------------------------*/

int cylfix(const int naxis[], struct wcsprm *wcs)

{
  static const char *function = "cylfix";

  unsigned short icnr, indx[NMAX], ncnr;
  int    j, k, stat[4], status;
  double img[4][NMAX], lat, lng, phi[4], phi0, phimax, phimin, pix[4][NMAX],
         *pixj, theta[4], theta0, world[4][NMAX], x, y;
  struct wcserr **err;

  if (naxis == 0x0) return FIXERR_NO_CHANGE;
  if (wcs == 0x0) return FIXERR_NULL_POINTER;
  err = &(wcs->err);

  /* Initialize if required. */
  if (wcs->flag != WCSSET) {
    if ((status = wcsset(wcs))) return fix_wcserr[status];
  }

  /* Check that we have a cylindrical projection. */
  if (wcs->cel.prj.category != CYLINDRICAL) return FIXERR_NO_CHANGE;
  if (wcs->naxis < 2) return FIXERR_NO_CHANGE;


  /* Compute the native longitude in each corner of the image. */
  ncnr = 1 << wcs->naxis;

  for (k = 0; k < NMAX; k++) {
    indx[k] = 1 << k;
  }

  phimin =  1.0e99;
  phimax = -1.0e99;
  for (icnr = 0; icnr < ncnr;) {
    /* Do four corners at a time. */
    for (j = 0; j < 4; j++, icnr++) {
      pixj = pix[j];

      for (k = 0; k < wcs->naxis; k++) {
        if (icnr & indx[k]) {
          *(pixj++) = naxis[k] + 0.5;
        } else {
          *(pixj++) = 0.5;
        }
      }
    }

    if (!(status = wcsp2s(wcs, 4, NMAX, pix[0], img[0], phi, theta, world[0],
                          stat))) {
      for (j = 0; j < 4; j++) {
        if (phi[j] < phimin) phimin = phi[j];
        if (phi[j] > phimax) phimax = phi[j];
      }
    }
  }

  if (phimin > phimax) return fix_wcserr[status];

  /* Any changes needed? */
  if (phimin >= -180.0 && phimax <= 180.0) return FIXERR_NO_CHANGE;


  /* Compute the new reference pixel coordinates. */
  phi0 = (phimin + phimax) / 2.0;
  theta0 = 0.0;

  if ((status = prjs2x(&(wcs->cel.prj), 1, 1, 1, 1, &phi0, &theta0, &x, &y,
                       stat))) {
    if (status == PRJERR_BAD_PARAM) {
      status = FIXERR_BAD_PARAM;
    } else {
      status = FIXERR_NO_REF_PIX_COORD;
    }
    return wcserr_set(WCSFIX_ERRMSG(status));
  }

  for (k = 0; k < wcs->naxis; k++) {
    img[0][k] = 0.0;
  }
  img[0][wcs->lng] = x;
  img[0][wcs->lat] = y;

  if ((status = linx2p(&(wcs->lin), 1, 0, img[0], pix[0]))) {
    return wcserr_set(WCSFIX_ERRMSG(fix_linerr[status]));
  }


  /* Compute celestial coordinates at the new reference pixel. */
  if ((status = wcsp2s(wcs, 1, 0, pix[0], img[0], phi, theta, world[0],
                       stat))) {
    return fix_wcserr[status];
  }

  /* Compute native coordinates of the celestial pole. */
  lng =  0.0;
  lat = 90.0;
  (void)sphs2x(wcs->cel.euler, 1, 1, 1, 1, &lng, &lat, phi, theta);

  wcs->crpix[wcs->lng] = pix[0][wcs->lng];
  wcs->crpix[wcs->lat] = pix[0][wcs->lat];
  wcs->crval[wcs->lng] = world[0][wcs->lng];
  wcs->crval[wcs->lat] = world[0][wcs->lat];
  wcs->lonpole = phi[0] - phi0;

  return wcsset(wcs);
}
