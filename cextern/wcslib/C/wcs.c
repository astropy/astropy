/*============================================================================

  WCSLIB 5.5 - an implementation of the FITS WCS standard.
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
  $Id: wcs.c,v 5.5 2015/05/05 13:16:31 mcalabre Exp $
*===========================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcserr.h"
#include "wcsmath.h"
#include "wcsprintf.h"
#include "wcstrig.h"
#include "wcsunits.h"
#include "wcsutil.h"
#include "lin.h"
#include "dis.h"
#include "log.h"
#include "spc.h"
#include "prj.h"
#include "sph.h"
#include "cel.h"
#include "tab.h"
#include "wcs.h"

const int WCSSET = 137;

/* Maximum number of PVi_ma and PSi_ma keywords. */
int NPVMAX = 64;
int NPSMAX =  8;

/* Map status return value to message. */
const char *wcs_errmsg[] = {
  "Success",
  "Null wcsprm pointer passed",
  "Memory allocation failed",
  "Linear transformation matrix is singular",
  "Inconsistent or unrecognized coordinate axis type",
  "Invalid parameter value",
  "Unrecognized coordinate transformation parameter",
  "Ill-conditioned coordinate transformation parameter",
  "One or more of the pixel coordinates were invalid",
  "One or more of the world coordinates were invalid",
  "Invalid world coordinate",
  "No solution found in the specified interval",
  "Invalid subimage specification",
  "Non-separable subimage coordinate system"};

/* Map error returns for lower-level routines. */
const int wcs_linerr[] = {
  WCSERR_SUCCESS,		/*  0: LINERR_SUCCESS         */
  WCSERR_NULL_POINTER,		/*  1: LINERR_NULL_POINTER    */
  WCSERR_MEMORY,		/*  2: LINERR_MEMORY          */
  WCSERR_SINGULAR_MTX,		/*  3: LINERR_SINGULAR_MTX    */
  WCSERR_BAD_PARAM,		/*  4: LINERR_DISTORT_INIT    */
  WCSERR_BAD_PIX,		/*  5: LINERR_DISTORT         */
  WCSERR_BAD_WORLD		/*  6: LINERR_DEDISTORT       */
};

const int wcs_logerr[] = {
  WCSERR_SUCCESS,		/*  0: LOGERR_SUCCESS         */
  WCSERR_NULL_POINTER,		/*  1: LOGERR_NULL_POINTER    */
  WCSERR_BAD_PARAM,		/*  2: LOGERR_BAD_LOG_REF_VAL */
  WCSERR_BAD_PIX,		/*  3: LOGERR_BAD_X           */
  WCSERR_BAD_WORLD		/*  4: LOGERR_BAD_WORLD       */
};

const int wcs_spcerr[] = {
				/* -1: SPCERR_NO_CHANGE       */
  WCSERR_SUCCESS,		/*  0: SPCERR_SUCCESS         */
  WCSERR_NULL_POINTER,		/*  1: SPCERR_NULL_POINTER    */
  WCSERR_BAD_PARAM,		/*  2: SPCERR_BAD_SPEC_PARAMS */
  WCSERR_BAD_PIX,		/*  3: SPCERR_BAD_X           */
  WCSERR_BAD_WORLD		/*  4: SPCERR_BAD_SPEC        */
};

const int wcs_celerr[] = {
  WCSERR_SUCCESS,		/*  0: CELERR_SUCCESS         */
  WCSERR_NULL_POINTER,		/*  1: CELERR_NULL_POINTER    */
  WCSERR_BAD_PARAM,		/*  2: CELERR_BAD_PARAM       */
  WCSERR_BAD_COORD_TRANS,	/*  3: CELERR_BAD_COORD_TRANS */
  WCSERR_ILL_COORD_TRANS,	/*  4: CELERR_ILL_COORD_TRANS */
  WCSERR_BAD_PIX,		/*  5: CELERR_BAD_PIX         */
  WCSERR_BAD_WORLD		/*  6: CELERR_BAD_WORLD       */
};

const int wcs_taberr[] = {
  WCSERR_SUCCESS,		/*  0: TABERR_SUCCESS         */
  WCSERR_NULL_POINTER,		/*  1: TABERR_NULL_POINTER    */
  WCSERR_MEMORY,		/*  2: TABERR_MEMORY          */
  WCSERR_BAD_PARAM,		/*  3: TABERR_BAD_PARAMS      */
  WCSERR_BAD_PIX,		/*  4: TABERR_BAD_X           */
  WCSERR_BAD_WORLD		/*  5: TABERR_BAD_WORLD       */
};

/* Convenience macro for invoking wcserr_set(). */
#define WCS_ERRMSG(status) WCSERR_SET(status), wcs_errmsg[status]

#ifndef signbit
#define signbit(X) ((X) < 0.0 ? 1 : 0)
#endif

/* Internal helper functions, not for general use. */
static int wcs_types(struct wcsprm *);
static int wcs_units(struct wcsprm *);

/*--------------------------------------------------------------------------*/

int wcsnpv(int npvmax) { if (npvmax >= 0) NPVMAX = npvmax; return NPVMAX; }
int wcsnps(int npsmax) { if (npsmax >= 0) NPSMAX = npsmax; return NPSMAX; }

/*--------------------------------------------------------------------------*/

int wcsini(int alloc, int naxis, struct wcsprm *wcs)

{
  static const char *function = "wcsini";

  int i, j, k, status;
  double *cd;
  struct wcserr **err;

  if (wcs == 0x0) return WCSERR_NULL_POINTER;

  /* Initialize error message handling. */
  err = &(wcs->err);
  if (wcs->flag != -1) {
    if (wcs->err) free(wcs->err);
    if (wcs->lin.err) free(wcs->lin.err);
    if (wcs->cel.err) free(wcs->cel.err);
    if (wcs->spc.err) free(wcs->spc.err);
  }
  wcs->err = 0x0;
  wcs->lin.err = 0x0;
  wcs->cel.err = 0x0;
  wcs->spc.err = 0x0;


  /* Initialize pointers. */
  if (wcs->flag == -1 || wcs->m_flag != WCSSET) {
    if (wcs->flag == -1) {
      wcs->tab   = 0x0;
      wcs->types = 0x0;
      wcs->lin.flag = -1;
    }

    /* Initialize memory management. */
    wcs->m_flag  = 0;
    wcs->m_naxis = 0;
    wcs->m_crpix = 0x0;
    wcs->m_pc    = 0x0;
    wcs->m_cdelt = 0x0;
    wcs->m_crval = 0x0;
    wcs->m_cunit = 0x0;
    wcs->m_ctype = 0x0;
    wcs->m_pv    = 0x0;
    wcs->m_ps    = 0x0;
    wcs->m_cd    = 0x0;
    wcs->m_crota = 0x0;
    wcs->m_colax = 0x0;
    wcs->m_cname = 0x0;
    wcs->m_crder = 0x0;
    wcs->m_csyer = 0x0;
    wcs->m_tab   = 0x0;
    wcs->m_wtb   = 0x0;
  }

  if (naxis < 0) {
    return wcserr_set(WCSERR_SET(WCSERR_MEMORY),
      "naxis must not be negative (got %d)", naxis);
  }


  /* Allocate memory for arrays if required. */
  if (alloc ||
     wcs->crpix == 0x0 ||
     wcs->pc    == 0x0 ||
     wcs->cdelt == 0x0 ||
     wcs->crval == 0x0 ||
     wcs->cunit == 0x0 ||
     wcs->ctype == 0x0 ||
     (NPVMAX && wcs->pv == 0x0) ||
     (NPSMAX && wcs->ps == 0x0) ||
     wcs->cd    == 0x0 ||
     wcs->crota == 0x0 ||
     wcs->colax == 0x0 ||
     wcs->cname == 0x0 ||
     wcs->crder == 0x0 ||
     wcs->csyer == 0x0) {

    /* Was sufficient allocated previously? */
    if (wcs->m_flag == WCSSET &&
       (wcs->m_naxis < naxis  ||
        wcs->npvmax  < NPVMAX ||
        wcs->npsmax  < NPSMAX)) {
      /* No, free it. */
      wcsfree(wcs);
    }

    if (alloc || wcs->crpix == 0x0) {
      if (wcs->m_crpix) {
        /* In case the caller fiddled with it. */
        wcs->crpix = wcs->m_crpix;

      } else {
        if ((wcs->crpix = calloc(naxis, sizeof(double))) == 0x0) {
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_crpix = wcs->crpix;
      }
    }

    if (alloc || wcs->pc == 0x0) {
      if (wcs->m_pc) {
        /* In case the caller fiddled with it. */
        wcs->pc = wcs->m_pc;

      } else {
        if ((wcs->pc = calloc(naxis*naxis, sizeof(double))) == 0x0) {
          wcsfree(wcs);
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_pc    = wcs->pc;
      }
    }

    if (alloc || wcs->cdelt == 0x0) {
      if (wcs->m_cdelt) {
        /* In case the caller fiddled with it. */
        wcs->cdelt = wcs->m_cdelt;

      } else {
        if ((wcs->cdelt = calloc(naxis, sizeof(double))) == 0x0) {
          wcsfree(wcs);
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_cdelt = wcs->cdelt;
      }
    }

    if (alloc || wcs->crval == 0x0) {
      if (wcs->m_crval) {
        /* In case the caller fiddled with it. */
        wcs->crval = wcs->m_crval;

      } else {
        if ((wcs->crval = calloc(naxis, sizeof(double))) == 0x0) {
          wcsfree(wcs);
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_crval = wcs->crval;
      }
    }

    if (alloc || wcs->cunit == 0x0) {
      if (wcs->m_cunit) {
        /* In case the caller fiddled with it. */
        wcs->cunit = wcs->m_cunit;

      } else {
        if ((wcs->cunit = calloc(naxis, sizeof(char [72]))) == 0x0) {
          wcsfree(wcs);
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_cunit = wcs->cunit;
      }
    }

    if (alloc || wcs->ctype == 0x0) {
      if (wcs->m_ctype) {
        /* In case the caller fiddled with it. */
        wcs->ctype = wcs->m_ctype;

      } else {
        if ((wcs->ctype = calloc(naxis, sizeof(char [72]))) == 0x0) {
          wcsfree(wcs);
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_ctype = wcs->ctype;
      }
    }

    if (alloc || wcs->pv == 0x0) {
      if (wcs->m_pv) {
        /* In case the caller fiddled with it. */
        wcs->pv = wcs->m_pv;

      } else {
        if (NPVMAX) {
          if ((wcs->pv = calloc(NPVMAX, sizeof(struct pvcard))) == 0x0) {
            wcsfree(wcs);
            return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
          }
        } else {
          wcs->pv = 0x0;
        }

        wcs->npvmax  = NPVMAX;

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_pv    = wcs->pv;
      }
    }

    if (alloc || wcs->ps == 0x0) {
      if (wcs->m_ps) {
        /* In case the caller fiddled with it. */
        wcs->ps = wcs->m_ps;

      } else {
        if (NPSMAX) {
          if ((wcs->ps = calloc(NPSMAX, sizeof(struct pscard))) == 0x0) {
            wcsfree(wcs);
            return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
          }
        } else {
          wcs->ps = 0x0;
        }

        wcs->npsmax  = NPSMAX;

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_ps    = wcs->ps;
      }
    }

    if (alloc || wcs->cd == 0x0) {
      if (wcs->m_cd) {
        /* In case the caller fiddled with it. */
        wcs->cd = wcs->m_cd;

      } else {
        if ((wcs->cd = calloc(naxis*naxis, sizeof(double))) == 0x0) {
          wcsfree(wcs);
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_cd    = wcs->cd;
      }
    }

    if (alloc || wcs->crota == 0x0) {
      if (wcs->m_crota) {
        /* In case the caller fiddled with it. */
        wcs->crota = wcs->m_crota;

      } else {
        if ((wcs->crota = calloc(naxis, sizeof(double))) == 0x0) {
          wcsfree(wcs);
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_crota = wcs->crota;
      }
    }

    if (alloc || wcs->colax == 0x0) {
      if (wcs->m_colax) {
        /* In case the caller fiddled with it. */
        wcs->colax = wcs->m_colax;

      } else {
        if ((wcs->colax = calloc(naxis, sizeof(int))) == 0x0) {
          wcsfree(wcs);
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_colax = wcs->colax;
      }
    }

    if (alloc || wcs->cname == 0x0) {
      if (wcs->m_cname) {
        /* In case the caller fiddled with it. */
        wcs->cname = wcs->m_cname;

      } else {
        if ((wcs->cname = calloc(naxis, sizeof(char [72]))) == 0x0) {
          wcsfree(wcs);
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_cname = wcs->cname;
      }
    }

    if (alloc || wcs->crder == 0x0) {
      if (wcs->m_crder) {
        /* In case the caller fiddled with it. */
        wcs->crder = wcs->m_crder;

      } else {
        if ((wcs->crder = calloc(naxis, sizeof(double))) == 0x0) {
          wcsfree(wcs);
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_crder = wcs->crder;
      }
    }

    if (alloc || wcs->csyer == 0x0) {
      if (wcs->m_csyer) {
        /* In case the caller fiddled with it. */
        wcs->csyer = wcs->m_csyer;

      } else {
        if ((wcs->csyer = calloc(naxis, sizeof(double))) == 0x0) {
          wcsfree(wcs);
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_csyer = wcs->csyer;
      }
    }
  }


  wcs->flag  = 0;
  wcs->naxis = naxis;


  /* Set defaults for the linear transformation. */
  wcs->lin.crpix  = wcs->crpix;
  wcs->lin.pc     = wcs->pc;
  wcs->lin.cdelt  = wcs->cdelt;
  if ((status = linini(0, naxis, &(wcs->lin)))) {
    return wcserr_set(WCS_ERRMSG(wcs_linerr[status]));
  }


  /* CRVALia defaults to 0.0. */
  for (i = 0; i < naxis; i++) {
    wcs->crval[i] = 0.0;
  }


  /* CUNITia and CTYPEia are blank by default. */
  for (i = 0; i < naxis; i++) {
    memset(wcs->cunit[i], 0, 72);
    memset(wcs->ctype[i], 0, 72);
  }


  /* Set defaults for the celestial transformation parameters. */
  wcs->lonpole = UNDEFINED;
  wcs->latpole = +90.0;

  /* Set defaults for the spectral transformation parameters. */
  wcs->restfrq = 0.0;
  wcs->restwav = 0.0;

  /* Default parameter values. */
  wcs->npv = 0;
  for (k = 0; k < wcs->npvmax; k++) {
    wcs->pv[k].i = 0;
    wcs->pv[k].m = 0;
    wcs->pv[k].value = 0.0;
  }

  wcs->nps = 0;
  for (k = 0; k < wcs->npsmax; k++) {
    wcs->ps[k].i = 0;
    wcs->ps[k].m = 0;
    memset(wcs->ps[k].value, 0, 72);
  }

  /* Defaults for alternate linear transformations. */
  cd = wcs->cd;
  for (i = 0; i < naxis; i++) {
    for (j = 0; j < naxis; j++) {
      *(cd++) = 0.0;
    }
  }
  for (i = 0; i < naxis; i++) {
    wcs->crota[i] = 0.0;
  }
  wcs->altlin = 0;
  wcs->velref = 0;

  /* Defaults for auxiliary coordinate system information. */
  memset(wcs->alt, 0, 4);
  wcs->alt[0] = ' ';
  wcs->colnum = 0;

  memset(wcs->wcsname, 0, 72);
  for (i = 0; i < naxis; i++) {
    wcs->colax[i] = 0;
    memset(wcs->cname[i], 0, 72);
    wcs->crder[i] = UNDEFINED;
    wcs->csyer[i] = UNDEFINED;
  }
  memset(wcs->radesys, 0, 72);
  wcs->equinox    = UNDEFINED;
  memset(wcs->specsys, 0, 72);
  memset(wcs->ssysobs, 0, 72);
  wcs->velosys    = UNDEFINED;
  memset(wcs->ssyssrc, 0, 72);
  wcs->zsource    = UNDEFINED;
  wcs->obsgeo[0]  = UNDEFINED;
  wcs->obsgeo[1]  = UNDEFINED;
  wcs->obsgeo[2]  = UNDEFINED;
  memset(wcs->dateobs, 0, 72);
  memset(wcs->dateavg, 0, 72);
  wcs->mjdobs     = UNDEFINED;
  wcs->mjdavg     = UNDEFINED;
  wcs->velangl    = UNDEFINED;

  wcs->ntab = 0;
  wcs->tab  = 0x0;
  wcs->nwtb = 0;
  wcs->wtb  = 0x0;

  /* Reset derived values. */
  strcpy(wcs->lngtyp, "    ");
  strcpy(wcs->lattyp, "    ");
  wcs->lng  = -1;
  wcs->lat  = -1;
  wcs->spec = -1;
  wcs->cubeface = -1;

  celini(&(wcs->cel));
  spcini(&(wcs->spc));

  return 0;
}

/*--------------------------------------------------------------------------*/

int wcssub(
  int alloc,
  const struct wcsprm *wcssrc,
  int *nsub,
  int axes[],
  struct wcsprm *wcsdst)

{
  static const char *function = "wcssub";

  char *c, ctypei[16];
  int  axis, cubeface, dealloc, dummy, i, itab, *itmp = 0x0, j, k, latitude,
       longitude, m, *map, msub, naxis, npv, nps, ntmp, other, spectral,
       status, stokes;
  const double *srcp;
  double *dstp;
  struct tabprm *tabp;
  struct wcserr **err;

  if (wcssrc == 0x0) return WCSERR_NULL_POINTER;
  if (wcsdst == 0x0) return WCSERR_NULL_POINTER;
  err = &(wcsdst->err);

  if ((naxis = wcssrc->naxis) <= 0) {
    return wcserr_set(WCSERR_SET(WCSERR_MEMORY),
      "naxis must be positive (got %d)", naxis);
  }

  if (nsub == 0x0) {
    nsub = &dummy;
    *nsub = naxis;
  } else if (*nsub == 0) {
    *nsub = naxis;
  }

  /* Allocate enough temporary storage to hold either axes[] or map[].*/
  ntmp = (*nsub <= naxis) ? naxis : *nsub;
  if ((itmp = calloc(ntmp, sizeof(int))) == 0x0) {
    return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
  }

  if ((dealloc = (axes == 0x0))) {
    /* Construct an index array. */
    if ((axes = calloc(naxis, sizeof(int))) == 0x0) {
      free(itmp);
      return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
    }

    for (i = 0; i < naxis; i++) {
      axes[i] = i+1;
    }
  }

  /* So that we don't try to free an uninitialized pointer on cleanup. */
  wcsdst->m_tab = 0x0;


  msub = 0;
  for (j = 0; j < *nsub; j++) {
    axis = axes[j];

    if (abs(axis) > 0x1000) {
      /* Subimage extraction by type. */
      k = abs(axis) & 0xFF;

      longitude = k & WCSSUB_LONGITUDE;
      latitude  = k & WCSSUB_LATITUDE;
      cubeface  = k & WCSSUB_CUBEFACE;
      spectral  = k & WCSSUB_SPECTRAL;
      stokes    = k & WCSSUB_STOKES;

      if ((other = (axis < 0))) {
        longitude = !longitude;
        latitude  = !latitude;
        cubeface  = !cubeface;
        spectral  = !spectral;
        stokes    = !stokes;
      }

      for (i = 0; i < naxis; i++) {
        strncpy (ctypei, (char *)(wcssrc->ctype + i), 8);
        ctypei[8] = '\0';

        /* Find the last non-blank character. */
        c = ctypei + 8;
        while (c-- > ctypei) {
          if (*c == ' ') *c = '\0';
          if (*c != '\0') break;
        }

        if (
          strcmp(ctypei,   "RA")  == 0 ||
          strcmp(ctypei+1, "LON") == 0 ||
          strcmp(ctypei+2, "LN")  == 0 ||
          strncmp(ctypei,   "RA---", 5) == 0 ||
          strncmp(ctypei+1, "LON-", 4) == 0 ||
          strncmp(ctypei+2, "LN-", 3) == 0) {
          if (!longitude) {
            continue;
          }

        } else if (
          strcmp(ctypei,   "DEC") == 0 ||
          strcmp(ctypei+1, "LAT") == 0 ||
          strcmp(ctypei+2, "LT")  == 0 ||
          strncmp(ctypei,   "DEC--", 5) == 0 ||
          strncmp(ctypei+1, "LAT-", 4) == 0 ||
          strncmp(ctypei+2, "LT-", 3) == 0) {
          if (!latitude) {
            continue;
          }

        } else if (strcmp(ctypei, "CUBEFACE") == 0) {
          if (!cubeface) {
            continue;
          }

        } else if ((
          strncmp(ctypei, "FREQ", 4) == 0 ||
          strncmp(ctypei, "ENER", 4) == 0 ||
          strncmp(ctypei, "WAVN", 4) == 0 ||
          strncmp(ctypei, "VRAD", 4) == 0 ||
          strncmp(ctypei, "WAVE", 4) == 0 ||
          strncmp(ctypei, "VOPT", 4) == 0 ||
          strncmp(ctypei, "ZOPT", 4) == 0 ||
          strncmp(ctypei, "AWAV", 4) == 0 ||
          strncmp(ctypei, "VELO", 4) == 0 ||
          strncmp(ctypei, "BETA", 4) == 0) &&
          (ctypei[4] == '\0' || ctypei[4] == '-')) {
          if (!spectral) {
            continue;
          }

        } else if (strcmp(ctypei, "STOKES") == 0) {
          if (!stokes) {
            continue;
          }

        } else if (!other) {
          continue;
        }

        /* This axis is wanted, but has it already been added? */
        for (k = 0; k < msub; k++) {
          if (itmp[k] == i+1) {
            break;
          }
        }
        if (k == msub) itmp[msub++] = i+1;
      }

    } else if (0 < axis && axis <= naxis) {
      /* Check that the requested axis has not already been added. */
      for (k = 0; k < msub; k++) {
        if (itmp[k] == axis) {
          break;
        }
      }
      if (k == msub) itmp[msub++] = axis;

    } else if (axis == 0) {
      /* Graft on a new axis. */
      itmp[msub++] = 0;

    } else {
      status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_SUBIMAGE));
      goto cleanup;
    }
  }

  if ((*nsub = msub) == 0) {
    status = wcsini(alloc, 0, wcsdst);
    goto cleanup;
  }

  for (i = 0; i < *nsub; i++) {
    axes[i] = itmp[i];
  }


  /* Construct the inverse axis map:
     axes[i] == j means that output axis i+1 comes from input axis j,
     axes[i] == 0 means to create a new axis,
      map[i] == j means that input axis i+1 goes to output axis j,
      map[i] == 0 means that input axis i+1 is not used. */
  map = itmp;
  for (i = 0; i < naxis; i++) {
    map[i] = 0;
  }

  for (i = 0; i < *nsub; i++) {
    if (axes[i] > 0) {
      map[axes[i]-1] = i+1;
    }
  }

  /* Check that the subimage coordinate system is separable. */
  srcp = wcssrc->pc;
  for (i = 0; i < naxis; i++) {
    for (j = 0; j < naxis; j++) {
      if (*(srcp++) == 0.0 || j == i) continue;

      if ((map[i] == 0) != (map[j] == 0)) {
        status = wcserr_set(WCS_ERRMSG(WCSERR_NON_SEPARABLE));
        goto cleanup;
      }
    }
  }


  /* Initialize the destination. */
  npv = NPVMAX;
  nps = NPSMAX;

  NPVMAX = 0;
  for (k = 0; k < wcssrc->npv; k++) {
    i = wcssrc->pv[k].i;
    if (i == 0 || (i > 0 && map[i-1])) {
      NPVMAX++;
    }
  }

  NPSMAX = 0;
  for (k = 0; k < wcssrc->nps; k++) {
    i = wcssrc->ps[k].i;
    if (i > 0 && map[i-1]) {
      NPSMAX++;
    }
  }

  status = wcsini(alloc, *nsub, wcsdst);

  NPVMAX = npv;
  NPSMAX = nps;

  if (status) {
    goto cleanup;
  }


  /* Linear transformation. */
  srcp = wcssrc->crpix;
  dstp = wcsdst->crpix;
  for (j = 0; j < *nsub; j++, dstp++) {
    if (axes[j] > 0) {
      k = axes[j] - 1;
      *dstp = *(srcp+k);
    }
  }

  srcp = wcssrc->pc;
  dstp = wcsdst->pc;
  for (i = 0; i < *nsub; i++) {
    if (axes[i] > 0) {
      for (j = 0; j < *nsub; j++, dstp++) {
        if (axes[j] > 0) {
          k = (axes[i]-1)*naxis + (axes[j]-1);
          *dstp = *(srcp+k);
        }
      }
    }
  }

  srcp = wcssrc->cdelt;
  dstp = wcsdst->cdelt;
  for (i = 0; i < *nsub; i++, dstp++) {
    if (axes[i] > 0) {
      k = axes[i] - 1;
      *dstp = *(srcp+k);
    }
  }

  /* Coordinate reference value. */
  srcp = wcssrc->crval;
  dstp = wcsdst->crval;
  for (i = 0; i < *nsub; i++, dstp++) {
    if (axes[i] > 0) {
      k = axes[i] - 1;
      *dstp = *(srcp+k);
    }
  }

  /* Coordinate units and type. */
  for (i = 0; i < *nsub; i++) {
    if (axes[i] > 0) {
      k = axes[i] - 1;
      strncpy(wcsdst->cunit[i], wcssrc->cunit[k], 72);
      strncpy(wcsdst->ctype[i], wcssrc->ctype[k], 72);
    }
  }

  /* Celestial and spectral transformation parameters. */
  wcsdst->lonpole = wcssrc->lonpole;
  wcsdst->latpole = wcssrc->latpole;
  wcsdst->restfrq = wcssrc->restfrq;
  wcsdst->restwav = wcssrc->restwav;

  /* Parameter values. */
  npv = 0;
  for (k = 0; k < wcssrc->npv; k++) {
    i = wcssrc->pv[k].i;
    if (i == 0) {
      /* i == 0 is a special code that means "the latitude axis". */
      wcsdst->pv[npv] = wcssrc->pv[k];
      wcsdst->pv[npv].i = 0;
      npv++;
    } else if (i > 0 && map[i-1]) {
      wcsdst->pv[npv] = wcssrc->pv[k];
      wcsdst->pv[npv].i = map[i-1];
      npv++;
    }
  }
  wcsdst->npv = npv;

  nps = 0;
  for (k = 0; k < wcssrc->nps; k++) {
    i = wcssrc->ps[k].i;
    if (i > 0 && map[i-1]) {
      wcsdst->ps[nps] = wcssrc->ps[k];
      wcsdst->ps[nps].i = map[i-1];
      nps++;
    }
  }
  wcsdst->nps = nps;

  /* Alternate linear transformations. */
  srcp = wcssrc->cd;
  dstp = wcsdst->cd;
  for (i = 0; i < *nsub; i++) {
    if (axes[i] > 0) {
      for (j = 0; j < *nsub; j++, dstp++) {
        if (axes[j] > 0) {
          k = (axes[i]-1)*naxis + (axes[j]-1);
          *dstp = *(srcp+k);
        }
      }
    }
  }

  srcp = wcssrc->crota;
  dstp = wcsdst->crota;
  for (i = 0; i < *nsub; i++, dstp++) {
    if (axes[i] > 0) {
      k = axes[i] - 1;
      *dstp = *(srcp+k);
    }
  }

  wcsdst->altlin = wcssrc->altlin;
  wcsdst->velref = wcssrc->velref;

  /* Auxiliary coordinate system information. */
  strncpy(wcsdst->alt, wcssrc->alt, 4);
  wcsdst->colnum = wcssrc->colnum;

  strncpy(wcsdst->wcsname, wcssrc->wcsname, 72);
  for (i = 0; i < *nsub; i++) {
    if (axes[i] > 0) {
      k = axes[i] - 1;
      wcsdst->colax[i] = wcssrc->colax[k];
      strncpy(wcsdst->cname[i], wcssrc->cname[k], 72);
      wcsdst->crder[i] = wcssrc->crder[k];
      wcsdst->csyer[i] = wcssrc->csyer[k];
    }
  }

  strncpy(wcsdst->radesys, wcssrc->radesys, 72);
  wcsdst->equinox = wcssrc->equinox;

  strncpy(wcsdst->specsys, wcssrc->specsys, 72);
  strncpy(wcsdst->ssysobs, wcssrc->ssysobs, 72);
  wcsdst->velosys = wcssrc->velosys;
  strncpy(wcsdst->ssyssrc, wcssrc->ssyssrc, 72);
  wcsdst->zsource = wcssrc->zsource;

  wcsdst->obsgeo[0] = wcssrc->obsgeo[0];
  wcsdst->obsgeo[1] = wcssrc->obsgeo[1];
  wcsdst->obsgeo[2] = wcssrc->obsgeo[2];

  strncpy(wcsdst->dateobs, wcssrc->dateobs, 72);
  strncpy(wcsdst->dateavg, wcssrc->dateavg, 72);
  wcsdst->mjdobs = wcssrc->mjdobs;
  wcsdst->mjdavg = wcssrc->mjdavg;


  /* Coordinate lookup tables; only copy what's needed. */
  wcsdst->ntab = 0;
  for (itab = 0; itab < wcssrc->ntab; itab++) {
    /* Is this table wanted? */
    for (m = 0; m < wcssrc->tab[itab].M; m++) {
      i = wcssrc->tab[itab].map[m];

      if (map[i-1]) {
        wcsdst->ntab++;
        break;
      }
    }
  }

  if (wcsdst->ntab) {
    /* Allocate memory for tabprm structs. */
    if ((wcsdst->tab = calloc(wcsdst->ntab, sizeof(struct tabprm))) == 0x0) {
      wcsdst->ntab = 0;

      status = wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
      goto cleanup;
    }

    wcsdst->m_tab = wcsdst->tab;
  }

  tabp = wcsdst->tab;
  for (itab = 0; itab < wcssrc->ntab; itab++) {
    for (m = 0; m < wcssrc->tab[itab].M; m++) {
      i = wcssrc->tab[itab].map[m];

      if (map[i-1]) {
        if ((status = tabcpy(1, wcssrc->tab + itab, tabp))) {
          wcserr_set(WCS_ERRMSG(wcs_taberr[status]));
          goto cleanup;
        }

        tabp++;
        break;
      }
    }
  }

  if (*nsub == naxis) {
      lincpy(1, &(wcssrc->lin), &(wcsdst->lin));
  }

cleanup:
  if (itmp) free(itmp);
  if (dealloc) {
    free(axes);
  }

  if (status && wcsdst->m_tab) {
    free(wcsdst->m_tab);
    wcsdst->tab   = 0x0;
    wcsdst->m_tab = 0x0;
  }

  return status;
}

/*--------------------------------------------------------------------------*/

int wcscompare(
  int cmp,
  double tol,
  const struct wcsprm *wcs1,
  const struct wcsprm *wcs2,
  int *equal)

{
  int i, j, naxis, naxis2;
  double diff;
  int tab_equal;
  int status;

  if (wcs1  == 0x0) return WCSERR_NULL_POINTER;
  if (wcs2  == 0x0) return WCSERR_NULL_POINTER;
  if (equal == 0x0) return WCSERR_NULL_POINTER;

  *equal = 0;

  if (wcs1->naxis != wcs2->naxis) {
    return 0;
  }

  naxis = wcs1->naxis;
  naxis2 = wcs1->naxis*wcs1->naxis;

  if (cmp & WCSCOMPARE_CRPIX) {
    /* Don't compare crpix. */
  } else if (cmp & WCSCOMPARE_TILING) {
    for (i = 0; i < naxis; ++i) {
      diff = wcs1->crpix[i] - wcs2->crpix[i];
      if ((double)(int)(diff) != diff) {
        return 0;
      }
    }
  } else {
    if (!wcsutil_Eq(naxis, tol, wcs1->crpix, wcs2->crpix)) {
      return 0;
    }
  }

  if (!wcsutil_Eq(naxis2, tol, wcs1->pc, wcs2->pc) ||
      !wcsutil_Eq(naxis, tol, wcs1->cdelt, wcs2->cdelt) ||
      !wcsutil_Eq(naxis, tol, wcs1->crval, wcs2->crval) ||
      !wcsutil_strEq(naxis, wcs1->cunit, wcs2->cunit) ||
      !wcsutil_strEq(naxis, wcs1->ctype, wcs2->ctype) ||
      !wcsutil_Eq(1, tol, &wcs1->lonpole, &wcs2->lonpole) ||
      !wcsutil_Eq(1, tol, &wcs1->latpole, &wcs2->latpole) ||
      !wcsutil_Eq(1, tol, &wcs1->restfrq, &wcs2->restfrq) ||
      !wcsutil_Eq(1, tol, &wcs1->restwav, &wcs2->restwav) ||
      wcs1->npv != wcs2->npv ||
      wcs1->nps != wcs2->nps) {
    return 0;
  }

  /* Compare pv cards, which may not be in the same order */
  for (i = 0; i < wcs1->npv; ++i) {
    for (j = 0; j < wcs2->npv; ++j) {
      if (wcs1->pv[i].i == wcs2->pv[j].i &&
          wcs1->pv[i].m == wcs2->pv[j].m) {
        if (!wcsutil_Eq(1, tol, &wcs1->pv[i].value, &wcs2->pv[j].value)) {
          return 0;
        }
        break;
      }
    }
    /* We didn't find a match, so they are not equal */
    if (j == wcs2->npv) {
      return 0;
    }
  }

  /* Compare ps cards, which may not be in the same order */
  for (i = 0; i < wcs1->nps; ++i) {
    for (j = 0; j < wcs2->nps; ++j) {
      if (wcs1->ps[i].i == wcs2->ps[j].i &&
          wcs1->ps[i].m == wcs2->ps[j].m) {
        if (strncmp(wcs1->ps[i].value, wcs2->ps[j].value, 72)) {
          return 0;
        }
        break;
      }
    }
    /* We didn't find a match, so they are not equal */
    if (j == wcs2->nps) {
      return 0;
    }
  }

  if (wcs1->flag != WCSSET || wcs2->flag != WCSSET) {
    if (!wcsutil_Eq(naxis2, tol, wcs1->cd, wcs2->cd) ||
        !wcsutil_Eq(naxis, tol, wcs1->crota, wcs2->crota) ||
        wcs1->altlin != wcs2->altlin ||
        wcs1->velref != wcs2->velref) {
      return 0;
    }
  }

  if (!(cmp & WCSCOMPARE_ANCILLARY)) {
    if (strncmp(wcs1->alt, wcs2->alt, 4) ||
        wcs1->colnum != wcs2->colnum ||
        !wcsutil_intEq(naxis, wcs1->colax, wcs2->colax) ||
        !wcsutil_strEq(naxis, wcs1->cname, wcs2->cname) ||
        !wcsutil_Eq(naxis, tol, wcs1->crder, wcs2->crder) ||
        !wcsutil_Eq(naxis, tol, wcs1->csyer, wcs2->csyer) ||
        strncmp(wcs1->dateavg, wcs2->dateavg, 72) ||
        strncmp(wcs1->dateobs, wcs2->dateobs, 72) ||
        !wcsutil_Eq(1, tol, &wcs1->equinox, &wcs2->equinox) ||
        !wcsutil_Eq(1, tol, &wcs1->mjdavg, &wcs2->mjdavg) ||
        !wcsutil_Eq(1, tol, &wcs1->mjdobs, &wcs2->mjdobs) ||
        !wcsutil_Eq(3, tol, wcs1->obsgeo, wcs2->obsgeo) ||
        strncmp(wcs1->radesys, wcs2->radesys, 72) ||
        strncmp(wcs1->specsys, wcs2->specsys, 72) ||
        strncmp(wcs1->ssysobs, wcs2->ssysobs, 72) ||
        !wcsutil_Eq(1, tol, &wcs1->velosys, &wcs2->velosys) ||
        !wcsutil_Eq(1, tol, &wcs1->zsource, &wcs2->zsource) ||
        strncmp(wcs1->ssyssrc, wcs2->ssyssrc, 72) ||
        !wcsutil_Eq(1, tol, &wcs1->velangl, &wcs2->velangl) ||
        strncmp(wcs1->wcsname, wcs2->wcsname, 72)) {
      return 0;
    }
  }

  /* Compare tabular parameters */
  if (wcs1->ntab != wcs2->ntab) {
    return 0;
  }

  for (i = 0; i < wcs1->ntab; ++i) {
    if ((status = tabcmp(0, tol, &wcs1->tab[i], &wcs2->tab[i], &tab_equal))) {
      return status;
    }
    if (!tab_equal) {
      return 0;
    }
  }

  *equal = 1;
  return 0;
}

/*--------------------------------------------------------------------------*/

int wcsfree(struct wcsprm *wcs)

{
  int j;

  if (wcs == 0x0) return WCSERR_NULL_POINTER;

  if (wcs->flag == -1) {
    wcs->lin.flag = -1;

  } else {
    /* Optionally allocated by wcsini() for given parameters. */
    if (wcs->m_flag == WCSSET) {
      if (wcs->crpix == wcs->m_crpix) wcs->crpix = 0x0;
      if (wcs->pc    == wcs->m_pc)    wcs->pc    = 0x0;
      if (wcs->cdelt == wcs->m_cdelt) wcs->cdelt = 0x0;
      if (wcs->crval == wcs->m_crval) wcs->crval = 0x0;
      if (wcs->cunit == wcs->m_cunit) wcs->cunit = 0x0;
      if (wcs->ctype == wcs->m_ctype) wcs->ctype = 0x0;
      if (wcs->pv    == wcs->m_pv)    wcs->pv    = 0x0;
      if (wcs->ps    == wcs->m_ps)    wcs->ps    = 0x0;
      if (wcs->cd    == wcs->m_cd)    wcs->cd    = 0x0;
      if (wcs->crota == wcs->m_crota) wcs->crota = 0x0;
      if (wcs->colax == wcs->m_colax) wcs->colax = 0x0;
      if (wcs->cname == wcs->m_cname) wcs->cname = 0x0;
      if (wcs->crder == wcs->m_crder) wcs->crder = 0x0;
      if (wcs->csyer == wcs->m_csyer) wcs->csyer = 0x0;
      if (wcs->tab   == wcs->m_tab)   wcs->tab   = 0x0;
      if (wcs->wtb   == wcs->m_wtb)   wcs->wtb   = 0x0;

      if (wcs->m_crpix)  free(wcs->m_crpix);
      if (wcs->m_pc)     free(wcs->m_pc);
      if (wcs->m_cdelt)  free(wcs->m_cdelt);
      if (wcs->m_crval)  free(wcs->m_crval);
      if (wcs->m_cunit)  free(wcs->m_cunit);
      if (wcs->m_ctype)  free(wcs->m_ctype);
      if (wcs->m_pv)     free(wcs->m_pv);
      if (wcs->m_ps)     free(wcs->m_ps);
      if (wcs->m_cd)     free(wcs->m_cd);
      if (wcs->m_crota)  free(wcs->m_crota);
      if (wcs->m_colax)  free(wcs->m_colax);
      if (wcs->m_cname)  free(wcs->m_cname);
      if (wcs->m_crder)  free(wcs->m_crder);
      if (wcs->m_csyer)  free(wcs->m_csyer);

      /* Allocated unconditionally by wcstab(). */
      if (wcs->m_tab) {
        for (j = 0; j < wcs->ntab; j++) {
          tabfree(wcs->m_tab + j);
        }

        free(wcs->m_tab);
      }
      if (wcs->m_wtb) free(wcs->m_wtb);
    }

    if (wcs->err) free(wcs->err);

    /* Allocated unconditionally by wcsset(). */
    if (wcs->types) free(wcs->types);

    if (wcs->lin.crpix == wcs->m_crpix) wcs->lin.crpix = 0x0;
    if (wcs->lin.pc    == wcs->m_pc)    wcs->lin.pc    = 0x0;
    if (wcs->lin.cdelt == wcs->m_cdelt) wcs->lin.cdelt = 0x0;
  }

  wcs->m_flag   = 0;
  wcs->m_naxis  = 0x0;
  wcs->m_crpix  = 0x0;
  wcs->m_pc     = 0x0;
  wcs->m_cdelt  = 0x0;
  wcs->m_crval  = 0x0;
  wcs->m_cunit  = 0x0;
  wcs->m_ctype  = 0x0;
  wcs->m_pv     = 0x0;
  wcs->m_ps     = 0x0;
  wcs->m_cd     = 0x0;
  wcs->m_crota  = 0x0;
  wcs->m_colax  = 0x0;
  wcs->m_cname  = 0x0;
  wcs->m_crder  = 0x0;
  wcs->m_csyer  = 0x0;

  wcs->ntab  = 0;
  wcs->m_tab = 0x0;
  wcs->nwtb  = 0;
  wcs->m_wtb = 0x0;

  wcs->types = 0x0;

  wcs->err  = 0x0;

  wcs->flag = 0;

  linfree(&(wcs->lin));
  celfree(&(wcs->cel));
  spcfree(&(wcs->spc));

  return 0;
}

/*--------------------------------------------------------------------------*/

int wcsprt(const struct wcsprm *wcs)

{
  int i, j, k;
  struct wtbarr *wtbp;

  if (wcs == 0x0) return WCSERR_NULL_POINTER;

  if (wcs->flag != WCSSET) {
    wcsprintf("The wcsprm struct is UNINITIALIZED.\n");
    return 0;
  }

  wcsprintf("       flag: %d\n", wcs->flag);
  wcsprintf("      naxis: %d\n", wcs->naxis);
  WCSPRINTF_PTR("      crpix: ", wcs->crpix, "\n");
  wcsprintf("            ");
  for (i = 0; i < wcs->naxis; i++) {
    wcsprintf("  %#- 11.5g", wcs->crpix[i]);
  }
  wcsprintf("\n");

  /* Linear transformation. */
  k = 0;
  WCSPRINTF_PTR("         pc: ", wcs->pc, "\n");
  for (i = 0; i < wcs->naxis; i++) {
    wcsprintf("    pc[%d][]:", i);
    for (j = 0; j < wcs->naxis; j++) {
      wcsprintf("  %#- 11.5g", wcs->pc[k++]);
    }
    wcsprintf("\n");
  }

  /* Coordinate increment at reference point. */
  WCSPRINTF_PTR("      cdelt: ", wcs->cdelt, "\n");
  wcsprintf("            ");
  for (i = 0; i < wcs->naxis; i++) {
    wcsprintf("  %#- 11.5g", wcs->cdelt[i]);
  }
  wcsprintf("\n");

  /* Coordinate value at reference point. */
  WCSPRINTF_PTR("      crval: ", wcs->crval, "\n");
  wcsprintf("            ");
  for (i = 0; i < wcs->naxis; i++) {
    wcsprintf("  %#- 11.5g", wcs->crval[i]);
  }
  wcsprintf("\n");

  /* Coordinate units and type. */
  WCSPRINTF_PTR("      cunit: ", wcs->cunit, "\n");
  for (i = 0; i < wcs->naxis; i++) {
    wcsprintf("             \"%s\"\n", wcs->cunit[i]);
  }

  WCSPRINTF_PTR("      ctype: ", wcs->ctype, "\n");
  for (i = 0; i < wcs->naxis; i++) {
    wcsprintf("             \"%s\"\n", wcs->ctype[i]);
  }

  /* Celestial and spectral transformation parameters. */
  if (undefined(wcs->lonpole)) {
    wcsprintf("    lonpole: UNDEFINED\n");
  } else {
    wcsprintf("    lonpole: %9f\n", wcs->lonpole);
  }
  wcsprintf("    latpole: %9f\n", wcs->latpole);
  wcsprintf("    restfrq: %f\n", wcs->restfrq);
  wcsprintf("    restwav: %f\n", wcs->restwav);

  /* Parameter values. */
  wcsprintf("        npv: %d\n", wcs->npv);
  wcsprintf("     npvmax: %d\n", wcs->npvmax);
  WCSPRINTF_PTR("         pv: ", wcs->pv, "\n");
  for (i = 0; i < wcs->npv; i++) {
    wcsprintf("             %3d%4d  %#- 11.5g\n", (wcs->pv[i]).i,
      (wcs->pv[i]).m, (wcs->pv[i]).value);
  }
  wcsprintf("        nps: %d\n", wcs->nps);
  wcsprintf("     npsmax: %d\n", wcs->npsmax);
  WCSPRINTF_PTR("         ps: ", wcs->ps, "\n");
  for (i = 0; i < wcs->nps; i++) {
    wcsprintf("             %3d%4d  %s\n", (wcs->ps[i]).i,
      (wcs->ps[i]).m, (wcs->ps[i]).value);
  }

  /* Alternate linear transformations. */
  k = 0;
  WCSPRINTF_PTR("         cd: ", wcs->cd, "\n");
  if (wcs->cd) {
    for (i = 0; i < wcs->naxis; i++) {
      wcsprintf("    cd[%d][]:", i);
      for (j = 0; j < wcs->naxis; j++) {
        wcsprintf("  %#- 11.5g", wcs->cd[k++]);
      }
      wcsprintf("\n");
    }
  }

  WCSPRINTF_PTR("      crota: ", wcs->crota, "\n");
  if (wcs->crota) {
    wcsprintf("            ");
    for (i = 0; i < wcs->naxis; i++) {
      wcsprintf("  %#- 11.5g", wcs->crota[i]);
    }
    wcsprintf("\n");
  }

  wcsprintf("     altlin: %d\n", wcs->altlin);
  wcsprintf("     velref: %d\n", wcs->velref);



  /* Auxiliary coordinate system information. */
  wcsprintf("        alt: '%c'\n", wcs->alt[0]);
  wcsprintf("     colnum: %d\n", wcs->colnum);

  WCSPRINTF_PTR("      colax: ", wcs->colax, "\n");
  if (wcs->colax) {
    wcsprintf("           ");
    for (i = 0; i < wcs->naxis; i++) {
      wcsprintf("  %5d", wcs->colax[i]);
    }
    wcsprintf("\n");
  }

  if (wcs->wcsname[0] == '\0') {
    wcsprintf("    wcsname: UNDEFINED\n");
  } else {
    wcsprintf("    wcsname: \"%s\"\n", wcs->wcsname);
  }

  WCSPRINTF_PTR("      cname: ", wcs->cname, "\n");
  if (wcs->cname) {
    for (i = 0; i < wcs->naxis; i++) {
      if (wcs->cname[i][0] == '\0') {
        wcsprintf("             UNDEFINED\n");
      } else {
        wcsprintf("             \"%s\"\n", wcs->cname[i]);
      }
    }
  }

  WCSPRINTF_PTR("      crder: ", wcs->crder, "\n");
  if (wcs->crder) {
    wcsprintf("           ");
    for (i = 0; i < wcs->naxis; i++) {
      if (undefined(wcs->crder[i])) {
        wcsprintf("  UNDEFINED   ");
      } else {
        wcsprintf("  %#- 11.5g", wcs->crder[i]);
      }
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("      csyer: ", wcs->csyer, "\n");
  if (wcs->csyer) {
    wcsprintf("           ");
    for (i = 0; i < wcs->naxis; i++) {
      if (undefined(wcs->csyer[i])) {
        wcsprintf("  UNDEFINED   ");
      } else {
        wcsprintf("  %#- 11.5g", wcs->csyer[i]);
      }
    }
    wcsprintf("\n");
  }

  if (wcs->radesys[0] == '\0') {
    wcsprintf("    radesys: UNDEFINED\n");
  } else {
    wcsprintf("    radesys: \"%s\"\n", wcs->radesys);
  }

  if (undefined(wcs->equinox)) {
    wcsprintf("    equinox: UNDEFINED\n");
  } else {
    wcsprintf("    equinox: %9f\n", wcs->equinox);
  }

  if (wcs->specsys[0] == '\0') {
    wcsprintf("    specsys: UNDEFINED\n");
  } else {
    wcsprintf("    specsys: \"%s\"\n", wcs->specsys);
  }

  if (wcs->ssysobs[0] == '\0') {
    wcsprintf("    ssysobs: UNDEFINED\n");
  } else {
    wcsprintf("    ssysobs: \"%s\"\n", wcs->ssysobs);
  }

  if (undefined(wcs->velosys)) {
    wcsprintf("    velosys: UNDEFINED\n");
  } else {
    wcsprintf("    velosys: %9f\n", wcs->velosys);
  }

  if (wcs->ssyssrc[0] == '\0') {
    wcsprintf("    ssyssrc: UNDEFINED\n");
  } else {
    wcsprintf("    ssyssrc: \"%s\"\n", wcs->ssyssrc);
  }

  if (undefined(wcs->zsource)) {
    wcsprintf("    zsource: UNDEFINED\n");
  } else {
    wcsprintf("    zsource: %9f\n", wcs->zsource);
  }

  wcsprintf("     obsgeo: ");
  for (i = 0; i < 3; i++) {
    if (undefined(wcs->obsgeo[i])) {
      wcsprintf("UNDEFINED     ");
    } else {
      wcsprintf("  %#- 11.5g", wcs->obsgeo[i]);
    }
  }
  wcsprintf("\n");

  if (wcs->dateobs[0] == '\0') {
    wcsprintf("    dateobs: UNDEFINED\n");
  } else {
    wcsprintf("    dateobs: \"%s\"\n", wcs->dateobs);
  }

  if (wcs->dateavg[0] == '\0') {
    wcsprintf("    dateavg: UNDEFINED\n");
  } else {
    wcsprintf("    dateavg: \"%s\"\n", wcs->dateavg);
  }

  if (undefined(wcs->mjdobs)) {
    wcsprintf("     mjdobs: UNDEFINED\n");
  } else {
    wcsprintf("     mjdobs: %9f\n", wcs->mjdobs);
  }

  if (undefined(wcs->mjdavg)) {
    wcsprintf("     mjdavg: UNDEFINED\n");
  } else {
    wcsprintf("     mjdavg: %9f\n", wcs->mjdavg);
  }

  wcsprintf("       ntab: %d\n", wcs->ntab);
  WCSPRINTF_PTR("        tab: ", wcs->tab, "");
  if (wcs->tab != 0x0) wcsprintf("  (see below)");
  wcsprintf("\n");
  wcsprintf("       nwtb: %d\n", wcs->nwtb);
  WCSPRINTF_PTR("        wtb: ", wcs->wtb, "");
  if (wcs->wtb != 0x0) wcsprintf("  (see below)");
  wcsprintf("\n");

  /* Derived values. */
  WCSPRINTF_PTR("      types: ", wcs->types, "\n           ");
  for (i = 0; i < wcs->naxis; i++) {
    wcsprintf("%5d", wcs->types[i]);
  }
  wcsprintf("\n");

  wcsprintf("     lngtyp: \"%s\"\n", wcs->lngtyp);
  wcsprintf("     lattyp: \"%s\"\n", wcs->lattyp);
  wcsprintf("        lng: %d\n", wcs->lng);
  wcsprintf("        lat: %d\n", wcs->lat);
  wcsprintf("       spec: %d\n", wcs->spec);
  wcsprintf("   cubeface: %d\n", wcs->cubeface);

  WCSPRINTF_PTR("        err: ", wcs->err, "\n");
  if (wcs->err) {
    wcserr_prt(wcs->err, "             ");
  }

  wcsprintf("        lin: (see below)\n");
  wcsprintf("        cel: (see below)\n");
  wcsprintf("        spc: (see below)\n");

  /* Memory management. */
  wcsprintf("     m_flag: %d\n", wcs->m_flag);
  wcsprintf("    m_naxis: %d\n", wcs->m_naxis);
  WCSPRINTF_PTR("    m_crpix: ", wcs->m_crpix, "");
  if (wcs->m_crpix == wcs->crpix) wcsprintf("  (= crpix)");
  wcsprintf("\n");
  WCSPRINTF_PTR("       m_pc: ", wcs->m_pc, "");
  if (wcs->m_pc == wcs->pc) wcsprintf("  (= pc)");
  wcsprintf("\n");
  WCSPRINTF_PTR("    m_cdelt: ", wcs->m_cdelt, "");
  if (wcs->m_cdelt == wcs->cdelt) wcsprintf("  (= cdelt)");
  wcsprintf("\n");
  WCSPRINTF_PTR("    m_crval: ", wcs->m_crval, "");
  if (wcs->m_crval == wcs->crval) wcsprintf("  (= crval)");
  wcsprintf("\n");
  WCSPRINTF_PTR("    m_cunit: ", wcs->m_cunit, "");
  if (wcs->m_cunit == wcs->cunit) wcsprintf("  (= cunit)");
  wcsprintf("\n");
  WCSPRINTF_PTR("    m_ctype: ", wcs->m_ctype, "");
  if (wcs->m_ctype == wcs->ctype) wcsprintf("  (= ctype)");
  wcsprintf("\n");
  WCSPRINTF_PTR("       m_pv: ", wcs->m_pv, "");
  if (wcs->m_pv == wcs->pv) wcsprintf("  (= pv)");
  wcsprintf("\n");
  WCSPRINTF_PTR("       m_ps: ", wcs->m_ps, "");
  if (wcs->m_ps == wcs->ps) wcsprintf("  (= ps)");
  wcsprintf("\n");
  WCSPRINTF_PTR("       m_cd: ", wcs->m_cd, "");
  if (wcs->m_cd == wcs->cd) wcsprintf("  (= cd)");
  wcsprintf("\n");
  WCSPRINTF_PTR("    m_crota: ", wcs->m_crota, "");
  if (wcs->m_crota == wcs->crota) wcsprintf("  (= crota)");
  wcsprintf("\n");
  wcsprintf("\n");
  WCSPRINTF_PTR("    m_colax: ", wcs->m_colax, "");
  if (wcs->m_colax == wcs->colax) wcsprintf("  (= colax)");
  wcsprintf("\n");
  WCSPRINTF_PTR("    m_cname: ", wcs->m_cname, "");
  if (wcs->m_cname == wcs->cname) wcsprintf("  (= cname)");
  wcsprintf("\n");
  WCSPRINTF_PTR("    m_crder: ", wcs->m_crder, "");
  if (wcs->m_crder == wcs->crder) wcsprintf("  (= crder)");
  wcsprintf("\n");
  WCSPRINTF_PTR("    m_csyer: ", wcs->m_csyer, "");
  if (wcs->m_csyer == wcs->csyer) wcsprintf("  (= csyer)");
  wcsprintf("\n");
  WCSPRINTF_PTR("      m_tab: ", wcs->m_tab, "");
  if (wcs->m_tab == wcs->tab) wcsprintf("  (= tab)");
  wcsprintf("\n");
  WCSPRINTF_PTR("      m_wtb: ", wcs->m_wtb, "");
  if (wcs->m_wtb == wcs->wtb) wcsprintf("  (= wtb)");
  wcsprintf("\n");

  /* Tabular transformation parameters. */
  if ((wtbp = wcs->wtb)) {
    for (j = 0; j < wcs->nwtb; j++, wtbp++) {
      wcsprintf("\n");
      wcsprintf("wtb[%d].*\n", j);
      wcsprintf("          i: %d\n", wtbp->i);
      wcsprintf("          m: %d\n", wtbp->m);
      wcsprintf("       kind: %c\n", wtbp->kind);
      wcsprintf("     extnam: %s\n", wtbp->extnam);
      wcsprintf("     extver: %d\n", wtbp->extver);
      wcsprintf("     extlev: %d\n", wtbp->extlev);
      wcsprintf("      ttype: %s\n", wtbp->ttype);
      wcsprintf("        row: %ld\n", wtbp->row);
      wcsprintf("       ndim: %d\n", wtbp->ndim);
      WCSPRINTF_PTR("     dimlen: ", wtbp->dimlen, "\n");
      WCSPRINTF_PTR("     arrayp: ", wtbp->arrayp, " -> ");
      WCSPRINTF_PTR("", *(wtbp->arrayp), "\n");
    }
  }

  if (wcs->tab) {
    for (j = 0; j < wcs->ntab; j++) {
      wcsprintf("\n");
      wcsprintf("tab[%d].*\n", j);
      tabprt(wcs->tab + j);
    }
  }

  /* Linear transformation parameters. */
  wcsprintf("\n");
  wcsprintf("   lin.*\n");
  linprt(&(wcs->lin));

  /* Celestial transformation parameters. */
  wcsprintf("\n");
  wcsprintf("   cel.*\n");
  celprt(&(wcs->cel));

  /* Spectral transformation parameters. */
  wcsprintf("\n");
  wcsprintf("   spc.*\n");
  spcprt(&(wcs->spc));

  return 0;
}

/*--------------------------------------------------------------------------*/

int wcsperr(const struct wcsprm *wcs, const char *prefix)

{
  int j;

  if (wcs == 0x0) return WCSERR_NULL_POINTER;

  if (wcs->err && wcserr_prt(wcs->err, prefix) == 0) {
    linperr(&(wcs->lin), prefix);
    celperr(&(wcs->cel), prefix);
    wcserr_prt(wcs->spc.err, prefix);
    if (wcs->tab) {
      for (j = 0; j < wcs->ntab; j++) {
        wcserr_prt((wcs->tab + j)->err, prefix);
      }
    }
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int wcsbchk(struct wcsprm *wcs, int bounds)

{
  int status;

  if (wcs == 0x0) return WCSERR_NULL_POINTER;

  if (wcs->flag != WCSSET) {
    if ((status = wcsset(wcs))) return status;
  }

  wcs->cel.prj.bounds = bounds;

  return 0;
}

/*--------------------------------------------------------------------------*/

int wcsset(struct wcsprm *wcs)

{
  static const char *function = "wcsset";

  char   dpq[8], scode[4], stype[5];
  int i, j, k, m, n, naxis, status;
  double lambda, rho;
  double *cd, *pc;
  struct disprm *dis;
  struct dpkey  *keyp;
  struct linprm *wcslin = &(wcs->lin);
  struct celprm *wcscel = &(wcs->cel);
  struct prjprm *wcsprj = &(wcscel->prj);
  struct spcprm *wcsspc = &(wcs->spc);
  struct wcserr **err;


  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  err = &(wcs->err);

  /* Determine axis types from CTYPEia. */
  if ((status = wcs_types(wcs))) {
    return status;
  }

  /* Convert to canonical units. */
  if ((status = wcs_units(wcs))) {
    return status;
  }

  naxis = wcs->naxis;


  /* Non-linear celestial axes present? */
  if (wcs->lng >= 0 && wcs->types[wcs->lng] == 2200) {
    celini(wcscel);

    /* CRVALia, LONPOLEa, and LATPOLEa keyvalues. */
    wcscel->ref[0] = wcs->crval[wcs->lng];
    wcscel->ref[1] = wcs->crval[wcs->lat];
    wcscel->ref[2] = wcs->lonpole;
    wcscel->ref[3] = wcs->latpole;

    /* Do alias translation for TPU/TPV before dealing with PVi_ma. */
    strncpy(wcsprj->code, wcs->ctype[wcs->lng]+5, 3);
    wcsprj->code[3] = '\0';
    if (strncmp(wcsprj->code, "TPU", 3) == 0 ||
        strncmp(wcsprj->code, "TPV", 3) == 0) {
      /* Translate the PV parameters. */
      if ((dis = calloc(1, sizeof(struct disprm))) == 0x0) {
        return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
      }

      disndp(6+wcs->npv);

      /* Attach it to linprm.  Also inits it. */
      dis->flag = -1;
      if (strncmp(wcsprj->code, "TPU", 3) == 0) {
        /* Prior distortion. */
        lindis(1, wcslin, dis);
        strcpy(dpq, "DP");
      } else {
        /* Sequent distortion. */
        lindis(2, wcslin, dis);
        strcpy(dpq, "DQ");
      }

      /* Yes, the distortion type is "TPV" even for TPU. */
      strcpy(dis->dtype[wcs->lng], "TPV");
      strcpy(dis->dtype[wcs->lat], "TPV");

      /* Keep the keywords in axis-order to aid debugging. */
      keyp = dis->dp;

      sprintf(dpq+2, "%d", wcs->lng+1);
      dpfill(keyp++, dpq, "NAXES",  0, 0, 2, 0.0);
      dpfill(keyp++, dpq, "AXIS.1", 0, 0, 1, 0.0);
      dpfill(keyp++, dpq, "AXIS.2", 0, 0, 2, 0.0);

      /* Copy distortion parameters for the longitude axis. */
      for (k = 0; k < wcs->npv; k++) {
        if (wcs->pv[k].i != wcs->lng+1) continue;
        sprintf(keyp->field, "%s.TPV.%d", dpq, wcs->pv[k].m);
        dpfill(keyp++, 0x0, 0x0, 0, 1, 0, wcs->pv[k].value);
      }

      /* Now the latitude axis. */
      sprintf(dpq+2, "%d", wcs->lat+1);
      dpfill(keyp++, dpq, "NAXES",  0, 0, 2, 0.0);
      dpfill(keyp++, dpq, "AXIS.1", 0, 0, 2, 0.0);
      dpfill(keyp++, dpq, "AXIS.2", 0, 0, 1, 0.0);

      for (k = 0; k < wcs->npv; k++) {
        if (wcs->pv[k].i != wcs->lat+1) continue;
        sprintf(keyp->field, "%s.TPV.%d", dpq, wcs->pv[k].m);
        dpfill(keyp++, 0x0, 0x0, 0, 1, 0, wcs->pv[k].value);
      }

      dis->ndp = keyp - dis->dp;

      /* Erase PVi_ma associated with the celestial axes. */
      n = 0;
      for (k = 0; k < wcs->npv; k++) {
        i = wcs->pv[k].i - 1;
        if (i == wcs->lng || i == wcs->lat) continue;

        wcs->pv[n].i = wcs->pv[k].i;
        wcs->pv[n].m = wcs->pv[k].m;
        wcs->pv[n].value = wcs->pv[k].value;

        n++;
      }

      wcs->npv = n;
      strcpy(wcsprj->code, "TAN");
    }

    /* PVi_ma keyvalues. */
    for (k = 0; k < wcs->npv; k++) {
      i = wcs->pv[k].i - 1;
      m = wcs->pv[k].m;

      if (i == -1) {
        /* From a PROJPn keyword. */
        i = wcs->lat;
      }

      if (i == wcs->lat) {
        /* PVi_ma associated with latitude axis. */
        if (m < 30) {
          wcsprj->pv[m] = wcs->pv[k].value;
        }

      } else if (i == wcs->lng) {
        /* PVi_ma associated with longitude axis. */
        switch (m) {
        case 0:
          wcscel->offset = (wcs->pv[k].value != 0.0);
          break;
        case 1:
          wcscel->phi0   = wcs->pv[k].value;
          break;
        case 2:
          wcscel->theta0 = wcs->pv[k].value;
          break;
        case 3:
          /* If present, overrides LONPOLEa. */
          wcscel->ref[2] = wcs->pv[k].value;
          break;
        case 4:
          /* If present, overrides LATPOLEa. */
          wcscel->ref[3] = wcs->pv[k].value;
          break;
        default:
          return wcserr_set(WCSERR_SET(WCSERR_BAD_COORD_TRANS),
            "PV%i_%i%s: Unrecognized coordinate transformation parameter",
            i+1, m, wcs->alt);
          break;
        }
      }
    }

    /* Do simple alias translations. */
    if (strncmp(wcs->ctype[wcs->lng]+5, "GLS", 3) == 0) {
      wcscel->offset = 1;
      wcscel->phi0   = 0.0;
      wcscel->theta0 = wcs->crval[wcs->lat];
      strcpy(wcsprj->code, "SFL");

    } else if (strncmp(wcs->ctype[wcs->lng]+5, "NCP", 3) == 0) {
      /* Convert NCP to SIN. */
      if (wcscel->ref[1] == 0.0) {
        return wcserr_set(WCSERR_SET(WCSERR_BAD_PARAM),
          "Invalid projection: NCP blows up on the equator");
      }

      strcpy(wcsprj->code, "SIN");
      wcsprj->pv[1] = 0.0;
      wcsprj->pv[2] = cosd(wcscel->ref[1])/sind(wcscel->ref[1]);
    }

    /* Initialize the celestial transformation routines. */
    wcsprj->r0 = 0.0;
    if ((status = celset(wcscel))) {
      return wcserr_set(WCS_ERRMSG(wcs_celerr[status]));
    }

    /* Update LONPOLE, LATPOLE, and PVi_ma keyvalues. */
    wcs->lonpole = wcscel->ref[2];
    wcs->latpole = wcscel->ref[3];

    for (k = 0; k < wcs->npv; k++) {
      i = wcs->pv[k].i - 1;
      m = wcs->pv[k].m;

      if (i == wcs->lng) {
        switch (m) {
        case 1:
          wcs->pv[k].value = wcscel->phi0;
          break;
        case 2:
          wcs->pv[k].value = wcscel->theta0;
          break;
        case 3:
          wcs->pv[k].value = wcscel->ref[2];
          break;
        case 4:
          wcs->pv[k].value = wcscel->ref[3];
          break;
        }
      }
    }
  }


  /* Non-linear spectral axis present? */
  if (wcs->spec >= 0 && wcs->types[wcs->spec] == 3300) {
    spcini(wcsspc);
    if ((status = spctype(wcs->ctype[wcs->spec], stype, scode, 0x0, 0x0, 0x0,
                          0x0, 0x0, err))) {
      return status;
    }
    strcpy(wcsspc->type, stype);
    strcpy(wcsspc->code, scode);

    /* CRVALia, RESTFRQa, and RESTWAVa keyvalues. */
    wcsspc->crval = wcs->crval[wcs->spec];
    wcsspc->restfrq = wcs->restfrq;
    wcsspc->restwav = wcs->restwav;

    /* PVi_ma keyvalues. */
    for (k = 0; k < wcs->npv; k++) {
      i = wcs->pv[k].i - 1;
      m = wcs->pv[k].m;

      if (i == wcs->spec) {
        /* PVi_ma associated with grism axis. */
        if (m < 7) {
          wcsspc->pv[m] = wcs->pv[k].value;
        }
      }
    }

    /* Initialize the spectral transformation routines. */
    if ((status = spcset(wcsspc))) {
      return wcserr_set(WCS_ERRMSG(wcs_spcerr[status]));
    }
  }


  /* Tabular axes present? */
  for (j = 0; j < wcs->ntab; j++) {
    if ((status = tabset(wcs->tab + j))) {
      return wcserr_set(WCS_ERRMSG(wcs_taberr[status]));
    }
  }


  /* Initialize the linear transformation. */
  wcs->altlin &= 7;
  if (wcs->altlin > 1 && !(wcs->altlin & 1)) {
    pc = wcs->pc;

    if (wcs->altlin & 2) {
      /* Copy CDi_ja to PCi_ja and reset CDELTia. */
      cd = wcs->cd;
      for (i = 0; i < naxis; i++) {
        for (j = 0; j < naxis; j++) {
          *(pc++) = *(cd++);
        }
        wcs->cdelt[i] = 1.0;
      }

    } else if (wcs->altlin & 4) {
      /* Construct PCi_ja from CROTAia. */
      if ((i = wcs->lng) >= 0 && (j = wcs->lat) >= 0) {
        rho = wcs->crota[j];

        if (wcs->cdelt[i] == 0.0) {
          return wcserr_set(WCSERR_SET(WCSERR_SINGULAR_MTX),
            "Singular transformation matrix, CDELT%d is zero", i+1);
        }
        lambda = wcs->cdelt[j]/wcs->cdelt[i];

        *(pc + i*naxis + i) = *(pc + j*naxis + j) = cosd(rho);
        *(pc + i*naxis + j) = *(pc + j*naxis + i) = sind(rho);
        *(pc + i*naxis + j) *= -lambda;
        *(pc + j*naxis + i) /=  lambda;
      }
    }
  }

  wcs->lin.crpix  = wcs->crpix;
  wcs->lin.pc     = wcs->pc;
  wcs->lin.cdelt  = wcs->cdelt;
  if ((status = linset(&(wcs->lin)))) {
    return wcserr_set(WCS_ERRMSG(wcs_linerr[status]));
  }


  /* Set defaults for radesys and equinox for equatorial or ecliptic. */
  if (strcmp(wcs->lngtyp, "RA")   == 0 ||
      strcmp(wcs->lngtyp, "ELON") == 0 ||
      strcmp(wcs->lngtyp, "HLON") == 0) {
    if (wcs->radesys[0] == '\0') {
      if (undefined(wcs->equinox)) {
        strcpy(wcs->radesys, "ICRS");
      } else if (wcs->equinox < 1984.0) {
        strcpy(wcs->radesys, "FK4");
      } else {
        strcpy(wcs->radesys, "FK5");
      }

    } else if (strcmp(wcs->radesys, "ICRS")  == 0 ||
               strcmp(wcs->radesys, "GAPPT") == 0) {
      /* Equinox is not applicable for these coordinate systems. */
      wcs->equinox = UNDEFINED;

    } else if (undefined(wcs->equinox)) {
      if (strcmp(wcs->radesys, "FK5") == 0) {
        wcs->equinox = 2000.0;
      } else if (strcmp(wcs->radesys, "FK4") == 0 ||
                 strcmp(wcs->radesys, "FK4-NO-E") == 0) {
        wcs->equinox = 1950.0;
      }
    }

  } else {
    /* No celestial axes, ensure that radesys and equinox are unset. */
    memset(wcs->radesys, 0, 72);
    wcs->equinox = UNDEFINED;
  }


  /* Strip off trailing blanks and null-fill auxiliary string members. */
  wcsutil_null_fill(4, wcs->alt);
  wcsutil_null_fill(72, wcs->wcsname);
  for (i = 0; i < naxis; i++) {
    wcsutil_null_fill(72, wcs->cname[i]);
  }
  wcsutil_null_fill(72, wcs->radesys);
  wcsutil_null_fill(72, wcs->specsys);
  wcsutil_null_fill(72, wcs->ssysobs);
  wcsutil_null_fill(72, wcs->ssyssrc);
  wcsutil_null_fill(72, wcs->dateobs);
  wcsutil_null_fill(72, wcs->dateavg);

  wcs->flag = WCSSET;

  return 0;
}

/* : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :  */

int wcs_types(struct wcsprm *wcs)

{
  static const char *function = "wcs_types";

  const int  nalias = 4;
  const char aliases [4][4] = {"NCP", "GLS", "TPU", "TPV"};

  const char *alt = "";
  char ctypei[16], pcode[4], requir[9], scode[4], specsys[9];
  int i, j, m, naxis, *ndx = 0x0, type;
  struct wcserr **err;

  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  err = &(wcs->err);

  /* Parse the CTYPEia keyvalues. */
  pcode[0]  = '\0';
  requir[0] = '\0';
  wcs->lng  = -1;
  wcs->lat  = -1;
  wcs->spec = -1;
  wcs->cubeface = -1;

  if (*(wcs->alt) != ' ') alt = wcs->alt;


  naxis = wcs->naxis;
  if (wcs->types) free(wcs->types);
  if ((wcs->types = calloc(naxis, sizeof(int))) == 0x0) {
    return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
  }

  for (i = 0; i < naxis; i++) {
    /* Null fill. */
    wcsutil_null_fill(72, wcs->ctype[i]);

    strncpy(ctypei, wcs->ctype[i], 15);
    ctypei[15] = '\0';

    /* Check for early Paper IV syntax (e.g. '-SIP' used by Spitzer). */
    if (strlen(ctypei) == 12 && ctypei[8] == '-') {
      /* Excise the "4-3-3" or "8-3"-form distortion code. */
      ctypei[8] = '\0';

      /* Remove trailing dashes from "8-3"-form codes. */
      for (j = 7; j > 0; j--) {
        if (ctypei[j] != '-') break;
        ctypei[j] = '\0';
      }
    }

    /* Logarithmic or tabular axis? */
    wcs->types[i] = 0;
    if (strcmp(ctypei+4, "-LOG") == 0) {
      /* Logarithmic axis. */
      wcs->types[i] = 400;

    } else if (strcmp(ctypei+4, "-TAB") == 0) {
      /* Tabular axis. */
      wcs->types[i] = 500;
    }

    if (wcs->types[i]) {
      /* Could have -LOG or -TAB with celestial or spectral types. */
      ctypei[4] = '\0';

      /* Take care of things like 'FREQ-LOG' or 'RA---TAB'. */
      for (j = 3; j >= 0; j--) {
        if (ctypei[j] != '-') break;
        ctypei[j] = '\0';
      }
    }

    /* Translate AIPS spectral types for spctyp(). */
    if (spcaips(ctypei, wcs->velref, ctypei, specsys) == 0) {
      strcpy(wcs->ctype[i], ctypei);
      if (wcs->specsys[0] == '\0') strcpy(wcs->specsys, specsys);
    }

    /* Process linear axes. */
    if (!(strlen(ctypei) == 8 && ctypei[4] == '-')) {
      /* Identify Stokes, celestial and spectral types. */
      if (strcmp(ctypei, "STOKES") == 0) {
        /* STOKES axis. */
        wcs->types[i] = 1100;

      } else if (strcmp(ctypei, "RA")  == 0 ||
        strcmp(ctypei+1, "LON") == 0 ||
        strcmp(ctypei+2, "LN")  == 0) {
        /* Longitude axis. */
        wcs->types[i] += 2000;
        if (wcs->lng < 0) {
          wcs->lng = i;
          strcpy(wcs->lngtyp, ctypei);
        }

      } else if (strcmp(ctypei,   "DEC") == 0 ||
                 strcmp(ctypei+1, "LAT") == 0 ||
                 strcmp(ctypei+2, "LT")  == 0) {
        /* Latitude axis. */
        wcs->types[i] += 2001;
        if (wcs->lat < 0) {
          wcs->lat = i;
          strcpy(wcs->lattyp, ctypei);
        }

      } else if (strcmp(ctypei, "CUBEFACE") == 0) {
        /* CUBEFACE axis. */
        if (wcs->cubeface == -1) {
          wcs->types[i] = 2102;
          wcs->cubeface = i;
        } else {
          /* Multiple CUBEFACE axes! */
          return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
            "Multiple CUBEFACE axes (in CTYPE%d%.1s and CTYPE%d%.1s)",
            wcs->cubeface+1, alt, i+1, alt);
        }

      } else if (spctyp(ctypei, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0) == 0) {
        /* Spectral axis. */
        if (wcs->spec < 0) wcs->spec = i;
        wcs->types[i] += 3000;
      }

      continue;
    }


    /* CTYPEia is in "4-3" form; is it a recognized spectral type? */
    if (spctyp(ctypei, 0x0, scode, 0x0, 0x0, 0x0, 0x0, 0x0) == 0) {
      /* Non-linear spectral axis found. */
      wcs->types[i] = 3300;

      /* Check uniqueness. */
      if (wcs->spec >= 0) {
        return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
          "Multiple spectral axes (in CTYPE%d%.1s and CTYPE%d%.1s)",
          wcs->spec+1, alt, i+1, alt);
      }

      wcs->spec = i;

      continue;
    }


    /* Is it a recognized celestial projection? */
    for (j = 0; j < prj_ncode; j++) {
      if (strncmp(ctypei+5, prj_codes[j], 3) == 0) break;
    }

    if (j == prj_ncode) {
      /* Not a standard projection code, maybe it's an alias. */
      for (j = 0; j < nalias; j++) {
        if (strncmp(ctypei+5, aliases[j], 3) == 0) break;
      }

      if (j == nalias) {
        /* Not a recognized algorithm code of any type. */
        wcs->types[i] = -1;
        return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
          "Unrecognized projection code (%s in CTYPE%d%.1s)",
          ctypei+5, i+1, alt);
      }
    }

    /* Parse the celestial axis type. */
    wcs->types[i] = 2200;
    if (*pcode == '\0') {
      /* The first of the two celestial axes. */
      sprintf(pcode, "%.3s", ctypei+5);

      if (strncmp(ctypei, "RA--", 4) == 0) {
        wcs->lng = i;
        strcpy(wcs->lngtyp, "RA");
        strcpy(wcs->lattyp, "DEC");
        ndx = &wcs->lat;
        sprintf(requir, "DEC--%s", pcode);
      } else if (strncmp(ctypei, "DEC-", 4) == 0) {
        wcs->lat = i;
        strcpy(wcs->lngtyp, "RA");
        strcpy(wcs->lattyp, "DEC");
        ndx = &wcs->lng;
        sprintf(requir, "RA---%s", pcode);
      } else if (strncmp(ctypei+1, "LON", 3) == 0) {
        wcs->lng = i;
        sprintf(wcs->lngtyp, "%cLON", ctypei[0]);
        sprintf(wcs->lattyp, "%cLAT", ctypei[0]);
        ndx = &wcs->lat;
        sprintf(requir, "%s-%s", wcs->lattyp, pcode);
      } else if (strncmp(ctypei+1, "LAT", 3) == 0) {
        wcs->lat = i;
        sprintf(wcs->lngtyp, "%cLON", ctypei[0]);
        sprintf(wcs->lattyp, "%cLAT", ctypei[0]);
        ndx = &wcs->lng;
        sprintf(requir, "%s-%s", wcs->lngtyp, pcode);
      } else if (strncmp(ctypei+2, "LN", 2) == 0) {
        wcs->lng = i;
        sprintf(wcs->lngtyp, "%c%cLN", ctypei[0], ctypei[1]);
        sprintf(wcs->lattyp, "%c%cLT", ctypei[0], ctypei[1]);
        ndx = &wcs->lat;
        sprintf(requir, "%s-%s", wcs->lattyp, pcode);
      } else if (strncmp(ctypei+2, "LT", 2) == 0) {
        wcs->lat = i;
        sprintf(wcs->lngtyp, "%c%cLN", ctypei[0], ctypei[1]);
        sprintf(wcs->lattyp, "%c%cLT", ctypei[0], ctypei[1]);
        ndx = &wcs->lng;
        sprintf(requir, "%s-%s", wcs->lngtyp, pcode);
      } else {
        /* Unrecognized celestial type. */
        wcs->types[i] = -1;

        wcs->lng = -1;
        wcs->lat = -1;
        return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
          "Unrecognized celestial type (%5s in CTYPE%d%.1s)",
          ctypei, i+1, alt);
      }

      if (wcs->lat >= 0) wcs->types[i]++;

    } else {
      /* Looking for the complementary celestial axis. */
      if (wcs->lat < 0) wcs->types[i]++;

      if (strncmp(ctypei, requir, 8) != 0) {
        /* Inconsistent projection types. */
        wcs->lng = -1;
        wcs->lat = -1;
        return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE), "Inconsistent "
          "projection types (expected %s, got %s in CTYPE%d%.1s)", requir,
          ctypei, i+1, alt);
      }

      *ndx = i;
      requir[0] = '\0';
    }
  }

  /* Do we have a complementary pair of celestial axes? */
  if (strcmp(requir, "")) {
    /* Unmatched celestial axis. */
    wcs->lng = -1;
    wcs->lat = -1;
    return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
      "Unmatched celestial axes");
  }

  /* Table group numbers. */
  for (j = 0; j < wcs->ntab; j++) {
    for (m = 0; m < wcs->tab[j].M; m++) {
      /* Get image axis number. */
      i = wcs->tab[j].map[m];

      type = (wcs->types[i] / 100) % 10;
      if (type != 5) {
        return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
          "Table parameters set for non-table axis type");
      }
      wcs->types[i] += 10 * j;
    }
  }

  return 0;
}

/* : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :  */

int wcs_units(struct wcsprm *wcs)

{
  static const char *function = "wcs_units";

  char ctype[9], units[16];
  int  i, j, naxis;
  double scale, offset, power;
  struct wcserr *uniterr = 0x0, **err;

  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  err = &(wcs->err);

  naxis = wcs->naxis;
  for (i = 0; i < naxis; i++) {
    /* Use types set by wcs_types(). */
    switch (wcs->types[i]/1000) {
    case 2:
      /* Celestial axis. */
      strcpy(units, "deg");
      break;

    case 3:
      /* Spectral axis. */
      strncpy(ctype, wcs->ctype[i], 8);
      ctype[8] = '\0';
      spctyp(ctype, 0x0, 0x0, 0x0, units, 0x0, 0x0, 0x0);
      break;

    default:
      continue;
    }

    /* Tabular axis, CDELTia and CRVALia relate to indices. */
    if ((wcs->types[i]/100)%10 == 5) {
      continue;
    }

    wcsutil_null_fill(72, wcs->cunit[i]);
    if (wcs->cunit[i][0]) {
      if (wcsunitse(wcs->cunit[i], units, &scale, &offset, &power,
                    &uniterr)) {
        wcserr_set(WCSERR_SET(WCSERR_BAD_COORD_TRANS),
          "In CUNIT%d%.1s: %s", i+1, (*wcs->alt)?wcs->alt:"", uniterr->msg);
        free(uniterr);
        return WCSERR_BAD_COORD_TRANS;
      }

      if (scale != 1.0) {
        wcs->cdelt[i] *= scale;
        wcs->crval[i] *= scale;

        for (j = 0; j < naxis; j++) {
          *(wcs->cd + i*naxis + j) *= scale;
        }

        strcpy(wcs->cunit[i], units);
      }

    } else {
      strcpy(wcs->cunit[i], units);
    }
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int wcsp2s(
  struct wcsprm *wcs,
  int ncoord,
  int nelem,
  const double pixcrd[],
  double imgcrd[],
  double phi[],
  double theta[],
  double world[],
  int stat[])

{
  static const char *function = "wcsp2s";

  int    bits, face, i, iso_x, iso_y, istat, *istatp, itab, k, m, nx, ny,
        *statp, status, type;
  double crvali, offset;
  register double *img, *wrl;
  struct celprm *wcscel = &(wcs->cel);
  struct prjprm *wcsprj = &(wcscel->prj);
  struct wcserr **err;

  /* Initialize if required. */
  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  err = &(wcs->err);

  if (wcs->flag != WCSSET) {
    if ((status = wcsset(wcs))) return status;
  }

  /* Sanity check. */
  if (ncoord < 1 || (ncoord > 1 && nelem < wcs->naxis)) {
    return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
      "ncoord and/or nelem inconsistent with the wcsprm");
  }


  /* Apply pixel-to-world linear transformation. */
  if ((status = linp2x(&(wcs->lin), ncoord, nelem, pixcrd, imgcrd))) {
    return wcserr_set(WCS_ERRMSG(wcs_linerr[status]));
  }

  /* Initialize status vectors. */
  if ((istatp = calloc(ncoord, sizeof(int))) == 0x0) {
    return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
  }

  stat[0] = 0;
  wcsutil_setAli(ncoord, 1, stat);


  /* Convert intermediate world coordinates to world coordinates. */
  for (i = 0; i < wcs->naxis; i++) {
    /* Extract the second digit of the axis type code. */
    type = (wcs->types[i] / 100) % 10;

    if (type <= 1) {
      /* Linear or quantized coordinate axis. */
      img = imgcrd + i;
      wrl = world  + i;
      crvali = wcs->crval[i];
      for (k = 0; k < ncoord; k++) {
        *wrl = *img + crvali;
        img += nelem;
        wrl += nelem;
      }

    } else if (wcs->types[i] == 2200) {
      /* Convert celestial coordinates; do we have a CUBEFACE axis? */
      if (wcs->cubeface != -1) {
        /* Separation between faces. */
        if (wcsprj->r0 == 0.0) {
          offset = 90.0;
        } else {
          offset = wcsprj->r0*PI/2.0;
        }

        /* Lay out faces in a plane. */
        img = imgcrd;
        statp = stat;
        bits = (1 << i) | (1 << wcs->lat);
        for (k = 0; k < ncoord; k++, statp++) {
          face = (int)(*(img+wcs->cubeface) + 0.5);
          if (fabs(*(img+wcs->cubeface) - face) > 1e-10) {
            *statp |= bits;
            status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_PIX));

          } else {
            *statp = 0;

            switch (face) {
            case 0:
              *(img+wcs->lat) += offset;
              break;
            case 1:
              break;
            case 2:
              *(img+i) += offset;
              break;
            case 3:
              *(img+i) += offset*2;
              break;
            case 4:
              *(img+i) += offset*3;
              break;
            case 5:
              *(img+wcs->lat) -= offset;
              break;
            default:
              *statp |= bits;
              status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_PIX));
            }
          }

          img += nelem;
        }
      }

      /* Check for constant x and/or y. */
      nx = ncoord;
      ny = 0;

      if ((iso_x = wcsutil_allEq(ncoord, nelem, imgcrd+i))) {
        nx = 1;
        ny = ncoord;
      }
      if ((iso_y = wcsutil_allEq(ncoord, nelem, imgcrd+wcs->lat))) {
        ny = 1;
      }

      /* Transform projection plane coordinates to celestial coordinates. */
      if ((istat = celx2s(wcscel, nx, ny, nelem, nelem, imgcrd+i,
                          imgcrd+wcs->lat, phi, theta, world+i,
                          world+wcs->lat, istatp))) {
        if (istat) {
          status = wcserr_set(WCS_ERRMSG(wcs_celerr[istat]));
          if (status != WCSERR_BAD_PIX) {
            goto cleanup;
          }
        }
      }

      /* If x and y were both constant, replicate values. */
      if (iso_x && iso_y) {
        wcsutil_setAll(ncoord, nelem, world+i);
        wcsutil_setAll(ncoord, nelem, world+wcs->lat);
        wcsutil_setAll(ncoord, 1, phi);
        wcsutil_setAll(ncoord, 1, theta);
        wcsutil_setAli(ncoord, 1, istatp);
      }

      if (istat == 5) {
        bits = (1 << i) | (1 << wcs->lat);
        wcsutil_setBit(ncoord, istatp, bits, stat);
      }

    } else if (type == 3 || type == 4) {
      /* Check for constant x. */
      nx = ncoord;
      if ((iso_x = wcsutil_allEq(ncoord, nelem, imgcrd+i))) {
        nx = 1;
      }

      istat = 0;
      if (wcs->types[i] == 3300) {
        /* Spectral coordinates. */
        istat = spcx2s(&(wcs->spc), nx, nelem, nelem, imgcrd+i, world+i,
                       istatp);
        if (istat) {
          status = wcserr_set(WCS_ERRMSG(wcs_spcerr[istat]));
          if (status != WCSERR_BAD_PIX) {
            goto cleanup;
          }
        }
      } else if (type == 4) {
        /* Logarithmic coordinates. */
        istat = logx2s(wcs->crval[i], nx, nelem, nelem, imgcrd+i, world+i,
                       istatp);
        if (istat) {
          status = wcserr_set(WCS_ERRMSG(wcs_logerr[istat]));
          if (status != WCSERR_BAD_PIX) {
            goto cleanup;
          }
        }
      }

      /* If x was constant, replicate values. */
      if (iso_x) {
        wcsutil_setAll(ncoord, nelem, world+i);
        wcsutil_setAli(ncoord, 1, istatp);
      }

      if (istat == 3) {
        wcsutil_setBit(ncoord, istatp, 1 << i, stat);
      }
    }
  }


  /* Do tabular coordinates. */
  for (itab = 0; itab < wcs->ntab; itab++) {
    istat = tabx2s(wcs->tab + itab, ncoord, nelem, imgcrd, world, istatp);

    if (istat) {
      status = wcserr_set(WCS_ERRMSG(wcs_taberr[istat]));

      if (status != WCSERR_BAD_PIX) {
        goto cleanup;

      } else {
        bits = 0;
        for (m = 0; m < wcs->tab[itab].M; m++) {
          bits |= 1 << wcs->tab[itab].map[m];
        }
        wcsutil_setBit(ncoord, istatp, bits, stat);
      }
    }
  }


  /* Zero the unused world coordinate elements. */
  for (i = wcs->naxis; i < nelem; i++) {
    world[i] = 0.0;
    wcsutil_setAll(ncoord, nelem, world+i);
  }

cleanup:
  free(istatp);
  return status;
}

/*--------------------------------------------------------------------------*/

int wcss2p(
  struct wcsprm* wcs,
  int ncoord,
  int nelem,
  const double world[],
  double phi[],
  double theta[],
  double imgcrd[],
  double pixcrd[],
  int stat[])

{
  static const char *function = "wcss2p";

  int    bits, i, isolat, isolng, isospec, istat, *istatp, itab, k, m, nlat,
         nlng, nwrld, status, type;
  double crvali, offset;
  register const double *wrl;
  register double *img;
  struct celprm *wcscel = &(wcs->cel);
  struct prjprm *wcsprj = &(wcscel->prj);
  struct wcserr **err;


  /* Initialize if required. */
  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  err = &(wcs->err);

  if (wcs->flag != WCSSET) {
    if ((status = wcsset(wcs))) return status;
  }

  /* Sanity check. */
  if (ncoord < 1 || (ncoord > 1 && nelem < wcs->naxis)) {
    return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
      "ncoord and/or nelem inconsistent with the wcsprm");
  }

  /* Initialize status vectors. */
  if ((istatp = calloc(ncoord, sizeof(int))) == 0x0) {
    return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
  }

  status = 0;
  stat[0] = 0;
  wcsutil_setAli(ncoord, 1, stat);


  /* Convert world coordinates to intermediate world coordinates. */
  for (i = 0; i < wcs->naxis; i++) {
    /* Extract the second digit of the axis type code. */
    type = (wcs->types[i] / 100) % 10;

    if (type <= 1) {
      /* Linear or quantized coordinate axis. */
      wrl = world  + i;
      img = imgcrd + i;
      crvali = wcs->crval[i];
      for (k = 0; k < ncoord; k++) {
        *img = *wrl - crvali;
        wrl += nelem;
        img += nelem;
      }

    } else if (wcs->types[i] == 2200) {
      /* Celestial coordinates; check for constant lng and/or lat. */
      nlng = ncoord;
      nlat = 0;

      if ((isolng = wcsutil_allEq(ncoord, nelem, world+i))) {
        nlng = 1;
        nlat = ncoord;
      }
      if ((isolat = wcsutil_allEq(ncoord, nelem, world+wcs->lat))) {
        nlat = 1;
      }

      /* Transform celestial coordinates to projection plane coordinates. */
      if ((istat = cels2x(wcscel, nlng, nlat, nelem, nelem, world+i,
                          world+wcs->lat, phi, theta, imgcrd+i,
                          imgcrd+wcs->lat, istatp))) {
        if (istat) {
          status = wcserr_set(WCS_ERRMSG(wcs_celerr[istat]));
          if (status != WCSERR_BAD_WORLD) {
            goto cleanup;
          }
        }
      }

      /* If lng and lat were both constant, replicate values. */
      if (isolng && isolat) {
        wcsutil_setAll(ncoord, nelem, imgcrd+i);
        wcsutil_setAll(ncoord, nelem, imgcrd+wcs->lat);
        wcsutil_setAll(ncoord, 1, phi);
        wcsutil_setAll(ncoord, 1, theta);
        wcsutil_setAli(ncoord, 1, istatp);
      }

      if (istat == CELERR_BAD_WORLD) {
        bits = (1 << i) | (1 << wcs->lat);
        wcsutil_setBit(ncoord, istatp, bits, stat);
      }

      /* Do we have a CUBEFACE axis? */
      if (wcs->cubeface != -1) {
        /* Separation between faces. */
        if (wcsprj->r0 == 0.0) {
          offset = 90.0;
        } else {
          offset = wcsprj->r0*PI/2.0;
        }

        /* Stack faces in a cube. */
        img = imgcrd;
        for (k = 0; k < ncoord; k++) {
          if (*(img+wcs->lat) < -0.5*offset) {
            *(img+wcs->lat) += offset;
            *(img+wcs->cubeface) = 5.0;
          } else if (*(img+wcs->lat) > 0.5*offset) {
            *(img+wcs->lat) -= offset;
            *(img+wcs->cubeface) = 0.0;
          } else if (*(img+i) > 2.5*offset) {
            *(img+i) -= 3.0*offset;
            *(img+wcs->cubeface) = 4.0;
          } else if (*(img+i) > 1.5*offset) {
            *(img+i) -= 2.0*offset;
            *(img+wcs->cubeface) = 3.0;
          } else if (*(img+i) > 0.5*offset) {
            *(img+i) -= offset;
            *(img+wcs->cubeface) = 2.0;
          } else {
            *(img+wcs->cubeface) = 1.0;
          }

          img += nelem;
        }
      }

    } else if (type == 3 || type == 4) {
      /* Check for constancy. */
      nwrld = ncoord;
      if ((isospec = wcsutil_allEq(ncoord, nelem, world+i))) {
        nwrld = 1;
      }

      istat = 0;
      if (wcs->types[i] == 3300) {
        /* Spectral coordinates. */
        istat = spcs2x(&(wcs->spc), nwrld, nelem, nelem, world+i,
                       imgcrd+i, istatp);
        if (istat) {
          status = wcserr_set(WCS_ERRMSG(wcs_spcerr[istat]));
          if (status != WCSERR_BAD_WORLD) {
            goto cleanup;
          }
        }
      } else if (type == 4) {
        /* Logarithmic coordinates. */
        istat = logs2x(wcs->crval[i], nwrld, nelem, nelem, world+i,
                       imgcrd+i, istatp);
        if (istat) {
          status = wcserr_set(WCS_ERRMSG(wcs_logerr[istat]));
          if (status != WCSERR_BAD_WORLD) {
            goto cleanup;
          }
        }
      }

      /* If constant, replicate values. */
      if (isospec) {
        wcsutil_setAll(ncoord, nelem, imgcrd+i);
        wcsutil_setAli(ncoord, 1, istatp);
      }

      if (istat == 4) {
        wcsutil_setBit(ncoord, istatp, 1 << i, stat);
      }
    }
  }


  /* Do tabular coordinates. */
  for (itab = 0; itab < wcs->ntab; itab++) {
    istat = tabs2x(wcs->tab + itab, ncoord, nelem, world, imgcrd, istatp);

    if (istat) {
      status = wcserr_set(WCS_ERRMSG(wcs_taberr[istat]));

      if (status == WCSERR_BAD_WORLD) {
        bits = 0;
        for (m = 0; m < wcs->tab[itab].M; m++) {
          bits |= 1 << wcs->tab[itab].map[m];
        }
        wcsutil_setBit(ncoord, istatp, bits, stat);

      } else {
        goto cleanup;
      }
    }
  }


  /* Zero the unused intermediate world coordinate elements. */
  for (i = wcs->naxis; i < nelem; i++) {
    imgcrd[i] = 0.0;
    wcsutil_setAll(ncoord, nelem, imgcrd+i);
  }


  /* Apply world-to-pixel linear transformation. */
  if ((istat = linx2p(&(wcs->lin), ncoord, nelem, imgcrd, pixcrd))) {
    status = wcserr_set(WCS_ERRMSG(wcs_linerr[istat]));
    goto cleanup;
  }

cleanup:
  free(istatp);
  return status;
}

/*--------------------------------------------------------------------------*/

int wcsmix(
  struct wcsprm *wcs,
  int mixpix,
  int mixcel,
  const double vspan[2],
  double vstep,
  int viter,
  double world[],
  double phi[],
  double theta[],
  double imgcrd[],
  double pixcrd[])

{
  static const char *function = "wcsmix";

  const int niter = 60;
  int    crossed, istep, iter, j, k, nstep, retry, stat[1], status;
  const double tol  = 1.0e-10;
  const double tol2 = 100.0*tol;
  double *worldlat, *worldlng;
  double lambda, span[2], step;
  double pixmix;
  double dlng, lng, lng0, lng0m, lng1, lng1m;
  double dlat, lat, lat0, lat0m, lat1, lat1m;
  double d, d0, d0m, d1, d1m, dx = 0.0;
  double dabs, dmin, lmin;
  double dphi, phi0, phi1;
  struct celprm *wcscel = &(wcs->cel);
  struct wcsprm wcs0;
  struct wcserr **err;

  /* Initialize if required. */
  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  err = &(wcs->err);

  if (wcs->flag != WCSSET) {
    if ((status = wcsset(wcs))) return status;
  }

  worldlng = world + wcs->lng;
  worldlat = world + wcs->lat;


  /* Check vspan. */
  if (vspan[0] <= vspan[1]) {
    span[0] = vspan[0];
    span[1] = vspan[1];
  } else {
    /* Swap them. */
    span[0] = vspan[1];
    span[1] = vspan[0];
  }

  /* Check vstep. */
  step = fabs(vstep);
  if (step == 0.0) {
    step = (span[1] - span[0])/10.0;
    if (step > 1.0 || step == 0.0) step = 1.0;
  }

  /* Check viter. */
  nstep = viter;
  if (nstep < 5) {
    nstep = 5;
  } else if (nstep > 10) {
    nstep = 10;
  }

  /* Given pixel element. */
  pixmix = pixcrd[mixpix];

  /* Iterate on the step size. */
  for (istep = 0; istep <= nstep; istep++) {
    if (istep) step /= 2.0;

    /* Iterate on the sky coordinate between the specified range. */
    if (mixcel == 1) {
      /* Celestial longitude is given. */

      /* Check whether the solution interval is a crossing interval. */
      lat0 = span[0];
      *worldlat = lat0;
      if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                           stat))) {
        if (status == WCSERR_BAD_WORLD) {
          status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
        }
        return status;
      }
      d0 = pixcrd[mixpix] - pixmix;

      dabs = fabs(d0);
      if (dabs < tol) return 0;

      lat1 = span[1];
      *worldlat = lat1;
      if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                           stat))) {
        if (status == WCSERR_BAD_WORLD) {
          status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
        }
        return status;
      }
      d1 = pixcrd[mixpix] - pixmix;

      dabs = fabs(d1);
      if (dabs < tol) return 0;

      lmin = lat1;
      dmin = dabs;

      /* Check for a crossing point. */
      if (signbit(d0) != signbit(d1)) {
        crossed = 1;
        dx = d1;
      } else {
        crossed = 0;
        lat0 = span[1];
      }

      for (retry = 0; retry < 4; retry++) {
        /* Refine the solution interval. */
        while (lat0 > span[0]) {
          lat0 -= step;
          if (lat0 < span[0]) lat0 = span[0];
          *worldlat = lat0;
          if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                               stat))) {
            if (status == WCSERR_BAD_WORLD) {
              status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
            }
            return status;
          }
          d0 = pixcrd[mixpix] - pixmix;

          /* Check for a solution. */
          dabs = fabs(d0);
          if (dabs < tol) return 0;

          /* Record the point of closest approach. */
          if (dabs < dmin) {
            lmin = lat0;
            dmin = dabs;
          }

          /* Check for a crossing point. */
          if (signbit(d0) != signbit(d1)) {
            crossed = 2;
            dx = d0;
            break;
          }

          /* Advance to the next subinterval. */
          lat1 = lat0;
          d1 = d0;
        }

        if (crossed) {
          /* A crossing point was found. */
          for (iter = 0; iter < niter; iter++) {
            /* Use regula falsi division of the interval. */
            lambda = d0/(d0-d1);
            if (lambda < 0.1) {
              lambda = 0.1;
            } else if (lambda > 0.9) {
              lambda = 0.9;
            }

            dlat = lat1 - lat0;
            lat = lat0 + lambda*dlat;
            *worldlat = lat;
            if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                                 stat))) {
              if (status == WCSERR_BAD_WORLD) {
                status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
              }
              return status;
            }

            /* Check for a solution. */
            d = pixcrd[mixpix] - pixmix;
            dabs = fabs(d);
            if (dabs < tol) return 0;

            if (dlat < tol) {
              /* An artifact of numerical imprecision. */
              if (dabs < tol2) return 0;

              /* Must be a discontinuity. */
              break;
            }

            /* Record the point of closest approach. */
            if (dabs < dmin) {
              lmin = lat;
              dmin = dabs;
            }

            if (signbit(d0) == signbit(d)) {
              lat0 = lat;
              d0 = d;
            } else {
              lat1 = lat;
              d1 = d;
            }
          }

          /* No convergence, must have been a discontinuity. */
          if (crossed == 1) lat0 = span[1];
          lat1 = lat0;
          d1 = dx;
          crossed = 0;

        } else {
          /* No crossing point; look for a tangent point. */
          if (lmin == span[0]) break;
          if (lmin == span[1]) break;

          lat = lmin;
          lat0 = lat - step;
          if (lat0 < span[0]) lat0 = span[0];
          lat1 = lat + step;
          if (lat1 > span[1]) lat1 = span[1];

          *worldlat = lat0;
          if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                               stat))) {
            if (status == WCSERR_BAD_WORLD) {
              status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
            }
            return status;
          }
          d0 = fabs(pixcrd[mixpix] - pixmix);

          d  = dmin;

          *worldlat = lat1;
          if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                               stat))) {
            if (status == WCSERR_BAD_WORLD) {
              status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
            }
            return status;
          }
          d1 = fabs(pixcrd[mixpix] - pixmix);

          for (iter = 0; iter < niter; iter++) {
            lat0m = (lat0 + lat)/2.0;
            *worldlat = lat0m;
            if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                                 stat))) {
              if (status == WCSERR_BAD_WORLD) {
                status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
              }
              return status;
            }
            d0m = fabs(pixcrd[mixpix] - pixmix);

            if (d0m < tol) return 0;

            lat1m = (lat1 + lat)/2.0;
            *worldlat = lat1m;
            if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                                 stat))) {
              if (status == WCSERR_BAD_WORLD) {
                status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
              }
              return status;
            }
            d1m = fabs(pixcrd[mixpix] - pixmix);

            if (d1m < tol) return 0;

            if (d0m < d && d0m <= d1m) {
              lat1 = lat;
              d1   = d;
              lat  = lat0m;
              d    = d0m;
            } else if (d1m < d) {
              lat0 = lat;
              d0   = d;
              lat  = lat1m;
              d    = d1m;
            } else {
              lat0 = lat0m;
              d0   = d0m;
              lat1 = lat1m;
              d1   = d1m;
            }
          }
        }
      }

    } else {
      /* Celestial latitude is given. */

      /* Check whether the solution interval is a crossing interval. */
      lng0 = span[0];
      *worldlng = lng0;
      if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                           stat))) {
        if (status == WCSERR_BAD_WORLD) {
          status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
        }
        return status;
      }
      d0 = pixcrd[mixpix] - pixmix;

      dabs = fabs(d0);
      if (dabs < tol) return 0;

      lng1 = span[1];
      *worldlng = lng1;
      if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                           stat))) {
        if (status == WCSERR_BAD_WORLD) {
          status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
        }
        return status;
      }
      d1 = pixcrd[mixpix] - pixmix;

      dabs = fabs(d1);
      if (dabs < tol) return 0;
      lmin = lng1;
      dmin = dabs;

      /* Check for a crossing point. */
      if (signbit(d0) != signbit(d1)) {
        crossed = 1;
        dx = d1;
      } else {
        crossed = 0;
        lng0 = span[1];
      }

      for (retry = 0; retry < 4; retry++) {
        /* Refine the solution interval. */
        while (lng0 > span[0]) {
          lng0 -= step;
          if (lng0 < span[0]) lng0 = span[0];
          *worldlng = lng0;
          if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                               stat))) {
            if (status == WCSERR_BAD_WORLD) {
              status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
            }
            return status;
          }
          d0 = pixcrd[mixpix] - pixmix;

          /* Check for a solution. */
          dabs = fabs(d0);
          if (dabs < tol) return 0;

          /* Record the point of closest approach. */
          if (dabs < dmin) {
            lmin = lng0;
            dmin = dabs;
          }

          /* Check for a crossing point. */
          if (signbit(d0) != signbit(d1)) {
            crossed = 2;
            dx = d0;
            break;
          }

          /* Advance to the next subinterval. */
          lng1 = lng0;
          d1 = d0;
        }

        if (crossed) {
          /* A crossing point was found. */
          for (iter = 0; iter < niter; iter++) {
            /* Use regula falsi division of the interval. */
            lambda = d0/(d0-d1);
            if (lambda < 0.1) {
              lambda = 0.1;
            } else if (lambda > 0.9) {
              lambda = 0.9;
            }

            dlng = lng1 - lng0;
            lng = lng0 + lambda*dlng;
            *worldlng = lng;
            if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                                 stat))) {
              if (status == WCSERR_BAD_WORLD) {
                status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
              }
              return status;
            }

            /* Check for a solution. */
            d = pixcrd[mixpix] - pixmix;
            dabs = fabs(d);
            if (dabs < tol) return 0;

            if (dlng < tol) {
              /* An artifact of numerical imprecision. */
              if (dabs < tol2) return 0;

              /* Must be a discontinuity. */
              break;
            }

            /* Record the point of closest approach. */
            if (dabs < dmin) {
              lmin = lng;
              dmin = dabs;
            }

            if (signbit(d0) == signbit(d)) {
              lng0 = lng;
              d0 = d;
            } else {
              lng1 = lng;
              d1 = d;
            }
          }

          /* No convergence, must have been a discontinuity. */
          if (crossed == 1) lng0 = span[1];
          lng1 = lng0;
          d1 = dx;
          crossed = 0;

        } else {
          /* No crossing point; look for a tangent point. */
          if (lmin == span[0]) break;
          if (lmin == span[1]) break;

          lng = lmin;
          lng0 = lng - step;
          if (lng0 < span[0]) lng0 = span[0];
          lng1 = lng + step;
          if (lng1 > span[1]) lng1 = span[1];

          *worldlng = lng0;
          if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                               stat))) {
            if (status == WCSERR_BAD_WORLD) {
              status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
            }
            return status;
          }
          d0 = fabs(pixcrd[mixpix] - pixmix);

          d  = dmin;

          *worldlng = lng1;
          if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                               stat))) {
            if (status == WCSERR_BAD_WORLD) {
              status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
            }
            return status;
          }
          d1 = fabs(pixcrd[mixpix] - pixmix);

          for (iter = 0; iter < niter; iter++) {
            lng0m = (lng0 + lng)/2.0;
            *worldlng = lng0m;
            if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                                 stat))) {
              if (status == WCSERR_BAD_WORLD) {
                status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
              }
              return status;
            }
            d0m = fabs(pixcrd[mixpix] - pixmix);

            if (d0m < tol) return 0;

            lng1m = (lng1 + lng)/2.0;
            *worldlng = lng1m;
            if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                                 stat))) {
              if (status == WCSERR_BAD_WORLD) {
                status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
              }
              return status;
            }
            d1m = fabs(pixcrd[mixpix] - pixmix);

            if (d1m < tol) return 0;

            if (d0m < d && d0m <= d1m) {
              lng1 = lng;
              d1   = d;
              lng  = lng0m;
              d    = d0m;
            } else if (d1m < d) {
              lng0 = lng;
              d0   = d;
              lng  = lng1m;
              d    = d1m;
            } else {
              lng0 = lng0m;
              d0   = d0m;
              lng1 = lng1m;
              d1   = d1m;
            }
          }
        }
      }
    }
  }


  /* Set cel0 to the unity transformation. */
  wcs0 = *wcs;
  wcs0.cel.euler[0] = -90.0;
  wcs0.cel.euler[1] =   0.0;
  wcs0.cel.euler[2] =  90.0;
  wcs0.cel.euler[3] =   1.0;
  wcs0.cel.euler[4] =   0.0;

  /* No convergence, check for aberrant behaviour at a native pole. */
  *theta = -90.0;
  for (j = 1; j <= 2; j++) {
    /* Could the celestial coordinate element map to a native pole? */
    *phi = 0.0;
    *theta = -*theta;
    sphx2s(wcscel->euler, 1, 1, 1, 1, phi, theta, &lng, &lat);

    if (mixcel == 1) {
      if (fabs(fmod(*worldlng-lng, 360.0)) > tol) continue;
      if (lat < span[0]) continue;
      if (lat > span[1]) continue;
      *worldlat = lat;
    } else {
      if (fabs(*worldlat-lat) > tol) continue;
      if (lng < span[0]) lng += 360.0;
      if (lng > span[1]) lng -= 360.0;
      if (lng < span[0]) continue;
      if (lng > span[1]) continue;
      *worldlng = lng;
    }

    /* Is there a solution for the given pixel coordinate element? */
    lng = *worldlng;
    lat = *worldlat;

    /* Feed native coordinates to wcss2p() with cel0 set to unity. */
    *worldlng = -180.0;
    *worldlat = *theta;
    if ((status = wcss2p(&wcs0, 1, 0, world, phi, theta, imgcrd, pixcrd,
                         stat))) {
      if (wcs->err) free(wcs->err);
      wcs->err = wcs0.err;
      if (status == WCSERR_BAD_WORLD) {
        status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
      }
      return status;
    }
    d0 = pixcrd[mixpix] - pixmix;

    /* Check for a solution. */
    if (fabs(d0) < tol) {
      /* Recall saved world coordinates. */
      *worldlng = lng;
      *worldlat = lat;
      return 0;
    }

    /* Search for a crossing interval. */
    phi0 = -180.0;
    for (k = -179; k <= 180; k++) {
      phi1 = (double) k;
      *worldlng = phi1;
      if ((status = wcss2p(&wcs0, 1, 0, world, phi, theta, imgcrd, pixcrd,
                           stat))) {
        if (wcs->err) free(wcs->err);
        wcs->err = wcs0.err;
        if (status == WCSERR_BAD_WORLD) {
          status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
        }
        return status;
      }
      d1 = pixcrd[mixpix] - pixmix;

      /* Check for a solution. */
      dabs = fabs(d1);
      if (dabs < tol) {
        /* Recall saved world coordinates. */
        *worldlng = lng;
        *worldlat = lat;
        return 0;
      }

      /* Is it a crossing interval? */
      if (signbit(d0) != signbit(d1)) break;

      phi0 = phi1;
      d0 = d1;
    }

    for (iter = 1; iter <= niter; iter++) {
      /* Use regula falsi division of the interval. */
      lambda = d0/(d0-d1);
      if (lambda < 0.1) {
        lambda = 0.1;
      } else if (lambda > 0.9) {
        lambda = 0.9;
      }

      dphi = phi1 - phi0;
      *worldlng = phi0 + lambda*dphi;
      if ((status = wcss2p(&wcs0, 1, 0, world, phi, theta, imgcrd, pixcrd,
                           stat))) {
        if (wcs->err) free(wcs->err);
        wcs->err = wcs0.err;
        if (status == WCSERR_BAD_WORLD) {
          status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
        }
        return status;
      }

      /* Check for a solution. */
      d = pixcrd[mixpix] - pixmix;
      dabs = fabs(d);
      if (dabs < tol || (dphi < tol && dabs < tol2)) {
        /* Recall saved world coordinates. */
        *worldlng = lng;
        *worldlat = lat;
        return 0;
      }

      if (signbit(d0) == signbit(d)) {
        phi0 = *worldlng;
        d0 = d;
      } else {
        phi1 = *worldlng;
        d1 = d;
      }
    }
  }


  /* No solution. */
  return wcserr_set(WCS_ERRMSG(WCSERR_NO_SOLUTION));
}

/*--------------------------------------------------------------------------*/

int wcssptr(
  struct wcsprm *wcs,
  int  *i,
  char ctype[9])

{
  static const char *function = "wcssptr";

  int    j, status;
  double cdelt, crval;
  struct wcserr **err;

  /* Initialize if required. */
  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  err = &(wcs->err);

  if (wcs->flag != WCSSET) {
    if ((status = wcsset(wcs))) return status;
  }

  if ((j = *i) < 0) {
    if ((j = wcs->spec) < 0) {
      /* Look for a linear spectral axis. */
      for (j = 0; j < wcs->naxis; j++) {
        if (wcs->types[j]/100 == 30) {
          break;
        }
      }

      if (j >= wcs->naxis) {
        /* No spectral axis. */
        return wcserr_set(WCSERR_SET(WCSERR_BAD_SUBIMAGE),
          "No spectral axis found.");
      }
    }

    *i = j;
  }

  /* Translate the spectral axis. */
  if ((status = spctrne(wcs->ctype[j], wcs->crval[j], wcs->cdelt[j],
                        wcs->restfrq, wcs->restwav, ctype, &crval, &cdelt,
                        &(wcs->spc.err)))) {
    return wcserr_set(WCS_ERRMSG(wcs_spcerr[status]));
  }


  /* Translate keyvalues. */
  wcs->flag = 0;
  wcs->cdelt[j] = cdelt;
  wcs->crval[j] = crval;
  spctyp(ctype, 0x0, 0x0, 0x0, wcs->cunit[j], 0x0, 0x0, 0x0);
  strcpy(wcs->ctype[j], ctype);

  /* This keeps things tidy if the spectral axis is linear. */
  spcfree(&(wcs->spc));
  spcini(&(wcs->spc));

  return 0;
}

/*--------------------------------------------------------------------------*/

#define STRINGIZE(s) STRINGIFY(s)
#define STRINGIFY(s) #s

const char *wcslib_version(
  int  vers[3])

{
  static const char *wcsver = STRINGIZE(WCSLIB_VERSION);

  if (vers != 0x0) {
    vers[2] = 0;
    sscanf(wcsver, "%d.%d.%d", vers, vers+1, vers+2);
  }

  return wcsver;
}
