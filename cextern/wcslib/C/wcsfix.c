/*============================================================================

  WCSLIB 5.11 - an implementation of the FITS WCS standard.
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
  $Id: wcsfix.c,v 5.11 2015/10/18 09:13:05 mcalabre Exp $
*===========================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcserr.h"
#include "wcsmath.h"
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

int wcsfixi(int ctrl, const int naxis[], struct wcsprm *wcs, int stat[],
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
      wcsutil_str2double(ctmp, "%lf", sec)) {
    return 1;
  }

  return 0;
}

static void write_date(char *buf, int hour, int minute, double sec)

{
  char ctmp[72];

  wcsutil_double2str(ctmp, "%04.1f", sec);
  sprintf(buf, "T%.2d:%.2d:%s", hour, minute, ctmp);
}

int datfix(struct wcsprm *wcs)

{
  static const char *function = "datfix";

  char orig_dateobs[72];
  char *dateobs;
  int  day, dd, hour = 0, jd, minute = 0, month, msec, n4, year;
  double mjdobs, sec = 0.0, t;
  struct wcserr **err;

  if (wcs == 0x0) return FIXERR_NULL_POINTER;
  err = &(wcs->err);

  dateobs = wcs->dateobs;
  strncpy(orig_dateobs, dateobs, 72);
  if (dateobs[0] == '\0') {
    if (undefined(wcs->mjdobs)) {
     /* No date information was provided. */
      return FIXERR_NO_CHANGE;

    } else {
      /* Calendar date from MJD. */
      jd = 2400001 + (int)wcs->mjdobs;

      n4 =  4*(jd + ((2*((4*jd - 17918)/146097)*3)/4 + 1)/2 - 37);
      dd = 10*(((n4-237)%1461)/4) + 5;

      year  = n4/1461 - 4712;
      month = (2 + dd/306)%12 + 1;
      day   = (dd%306)/10 + 1;
      sprintf(dateobs, "%.4d-%.2d-%.2d", year, month, day);

      /* Write time part only if non-zero. */
      if ((t = wcs->mjdobs - (int)wcs->mjdobs) > 0.0) {
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
        sprintf(dateobs+10, "T%.2d:%.2d:%.2d", hour, minute, msec/1000);

        /* Write fractions of a second only if non-zero. */
        if (msec%1000) {
          sprintf(dateobs+19, ".%.3d", msec%1000);
        }
      }
    }

  } else {
    if (strlen(dateobs) < 8) {
      /* Can't be a valid date. */
      return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
        "Invalid parameter value: date string too short '%s'", dateobs);
    }

    /* Identify the date format. */
    if (dateobs[4] == '-' && dateobs[7] == '-') {
      /* Standard year-2000 form: CCYY-MM-DD[Thh:mm:ss[.sss...]] */
      if (sscanf(dateobs, "%4d-%2d-%2d", &year, &month, &day) < 3) {
        return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
          "Invalid parameter value: invalid date '%s'", dateobs);
      }

      if (dateobs[10] == 'T') {
        if (parse_date(dateobs+11, &hour, &minute, &sec)) {
          return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
            "Invalid parameter value: invalid time '%s'", dateobs+11);
        }
      } else if (dateobs[10] == ' ') {
        hour = 0;
        minute = 0;
        sec = 0.0;
        if (parse_date(dateobs+11, &hour, &minute, &sec)) {
          write_date(dateobs+10, hour, minute, sec);
        } else {
          dateobs[10] = 'T';
        }
      }

    } else if (dateobs[4] == '/' && dateobs[7] == '/') {
      /* Also allow CCYY/MM/DD[Thh:mm:ss[.sss...]] */
      if (sscanf(dateobs, "%4d/%2d/%2d", &year, &month, &day) < 3) {
        return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
          "Invalid parameter value: invalid date '%s'", dateobs);
      }

      if (dateobs[10] == 'T') {
        if (parse_date(dateobs+11, &hour, &minute, &sec)) {
          return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
            "Invalid parameter value: invalid time '%s'", dateobs+11);
        }
      } else if (dateobs[10] == ' ') {
        hour = 0;
        minute = 0;
        sec = 0.0;
        if (parse_date(dateobs+11, &hour, &minute, &sec)) {
          write_date(dateobs+10, hour, minute, sec);
        } else {
          dateobs[10] = 'T';
        }
      }

      /* Looks ok, fix it up. */
      dateobs[4]  = '-';
      dateobs[7]  = '-';

    } else {
      if (dateobs[2] == '/' && dateobs[5] == '/') {
        /* Old format date: DD/MM/YY, also allowing DD/MM/CCYY. */
        if (sscanf(dateobs, "%2d/%2d/%4d", &day, &month, &year) < 3) {
          return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
            "Invalid parameter value: invalid date '%s'", dateobs);
        }

      } else if (dateobs[2] == '-' && dateobs[5] == '-') {
        /* Also recognize DD-MM-YY and DD-MM-CCYY */
        if (sscanf(dateobs, "%2d-%2d-%4d", &day, &month, &year) < 3) {
          return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
            "Invalid parameter value: invalid date '%s'", dateobs);
        }

      } else {
        /* Not a valid date format. */
        return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
          "Invalid parameter value: invalid date '%s'", dateobs);
      }

      if (year < 100) year += 1900;

      /* Doesn't have a time. */
      sprintf(dateobs, "%.4d-%.2d-%.2d", year, month, day);
    }

    /* Compute MJD. */
    mjdobs = (double)((1461*(year - (12-month)/10 + 4712))/4
             + (306*((month+9)%12) + 5)/10
             - (3*((year - (12-month)/10 + 4900)/100))/4
             + day - 2399904)
             + (hour + (minute + sec/60.0)/60.0)/24.0;

    if (undefined(wcs->mjdobs)) {
      wcs->mjdobs = mjdobs;
    } else {
      /* Check for consistency. */
      if (fabs(mjdobs - wcs->mjdobs) > 0.5) {
        return wcserr_set(WCSERR_SET(FIXERR_BAD_PARAM),
          "Invalid parameter value: inconsistent date '%s'", dateobs);
      }
    }
  }

  if (strncmp(orig_dateobs, dateobs, 72)) {
    wcserr_set(WCSERR_SET(FIXERR_DATE_FIX),
      "Changed '%s' to '%s'", orig_dateobs, dateobs);

    return 0;
  }

  return FIXERR_NO_CHANGE;
}

/*--------------------------------------------------------------------------*/

int unitfix(int ctrl, struct wcsprm *wcs)

{
  int  i, k, result, status = FIXERR_NO_CHANGE;
  char orig_unit[80], msg[WCSERR_MSG_LENGTH];
  const char *function = "unitfix";
  struct wcserr **err;

  if (wcs == 0x0) return FIXERR_NULL_POINTER;
  err = &(wcs->err);

  strcpy(msg, "Changed units: ");
  for (i = 0; i < wcs->naxis; i++) {
    strncpy(orig_unit, wcs->cunit[i], 80);
    result = wcsutrne(ctrl, wcs->cunit[i], &(wcs->err));
    if (result == 0 || result == 12) {
      k = strlen(msg);
      sprintf(msg+k, "'%s' -> '%s', ", orig_unit, wcs->cunit[i]);
      status = FIXERR_UNITS_ALIAS;
    }
  }

  if (status == FIXERR_UNITS_ALIAS) {
    k = strlen(msg) - 2;
    msg[k] = '\0';
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
            "Changed CTYPE%d from '%s' to '%s', and SPECSYS to '%s'",
            i+1, wcs->ctype[i], ctype, wcs->specsys);
        } else {
          wcserr_set(WCSERR_SET(FIXERR_SPC_UPDATE),
            "Changed CTYPE%d from '%s' to '%s'", i+1, wcs->ctype[i], ctype);
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
