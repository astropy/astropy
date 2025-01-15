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
  $Id: wcs.c,v 7.12 2022/09/09 04:57:58 mcalabre Exp $
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
#include "wtbarr.h"
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

// Maximum number of PVi_ma and PSi_ma keywords.
int NPVMAX = 64;
int NPSMAX =  8;

// Map status return value to message.
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
  "Non-separable subimage coordinate system",
  "wcsprm struct is unset, use wcsset()"};

// Map error returns for lower-level routines.
const int wcs_linerr[] = {
  WCSERR_SUCCESS,		//  0: LINERR_SUCCESS
  WCSERR_NULL_POINTER,		//  1: LINERR_NULL_POINTER
  WCSERR_MEMORY,		//  2: LINERR_MEMORY
  WCSERR_SINGULAR_MTX,		//  3: LINERR_SINGULAR_MTX
  WCSERR_BAD_PARAM,		//  4: LINERR_DISTORT_INIT
  WCSERR_BAD_PIX,		//  5: LINERR_DISTORT
  WCSERR_BAD_WORLD		//  6: LINERR_DEDISTORT
};

const int wcs_logerr[] = {
  WCSERR_SUCCESS,		//  0: LOGERR_SUCCESS
  WCSERR_NULL_POINTER,		//  1: LOGERR_NULL_POINTER
  WCSERR_BAD_PARAM,		//  2: LOGERR_BAD_LOG_REF_VAL
  WCSERR_BAD_PIX,		//  3: LOGERR_BAD_X
  WCSERR_BAD_WORLD		//  4: LOGERR_BAD_WORLD
};

const int wcs_spcerr[] = {
				// -1: SPCERR_NO_CHANGE
  WCSERR_SUCCESS,		//  0: SPCERR_SUCCESS
  WCSERR_NULL_POINTER,		//  1: SPCERR_NULL_POINTER
  WCSERR_BAD_PARAM,		//  2: SPCERR_BAD_SPEC_PARAMS
  WCSERR_BAD_PIX,		//  3: SPCERR_BAD_X
  WCSERR_BAD_WORLD		//  4: SPCERR_BAD_SPEC
};

const int wcs_celerr[] = {
  WCSERR_SUCCESS,		//  0: CELERR_SUCCESS
  WCSERR_NULL_POINTER,		//  1: CELERR_NULL_POINTER
  WCSERR_BAD_PARAM,		//  2: CELERR_BAD_PARAM
  WCSERR_BAD_COORD_TRANS,	//  3: CELERR_BAD_COORD_TRANS
  WCSERR_ILL_COORD_TRANS,	//  4: CELERR_ILL_COORD_TRANS
  WCSERR_BAD_PIX,		//  5: CELERR_BAD_PIX
  WCSERR_BAD_WORLD		//  6: CELERR_BAD_WORLD
};

const int wcs_taberr[] = {
  WCSERR_SUCCESS,		//  0: TABERR_SUCCESS
  WCSERR_NULL_POINTER,		//  1: TABERR_NULL_POINTER
  WCSERR_MEMORY,		//  2: TABERR_MEMORY
  WCSERR_BAD_PARAM,		//  3: TABERR_BAD_PARAMS
  WCSERR_BAD_PIX,		//  4: TABERR_BAD_X
  WCSERR_BAD_WORLD		//  5: TABERR_BAD_WORLD
};

// Convenience macro for invoking wcserr_set().
#define WCS_ERRMSG(status) WCSERR_SET(status), wcs_errmsg[status]

#ifndef signbit
#define signbit(X) ((X) < 0.0 ? 1 : 0)
#endif

// Internal helper functions, not for general use.
static int wcs_types(struct wcsprm *);
static int time_type(const char *);
static int time_code(const char *ctype, int nc);
static int wcs_units(struct wcsprm *);

//----------------------------------------------------------------------------

int wcsnpv(int npvmax) { if (npvmax >= 0) NPVMAX = npvmax; return NPVMAX; }
int wcsnps(int npsmax) { if (npsmax >= 0) NPSMAX = npsmax; return NPSMAX; }

//----------------------------------------------------------------------------

int wcsini(int alloc, int naxis, struct wcsprm *wcs)

{
  return wcsinit(alloc, naxis, wcs, -1, -1, -1);
}

//----------------------------------------------------------------------------

int wcsinit(
  int alloc,
  int naxis,
  struct wcsprm *wcs,
  int npvmax,
  int npsmax,
  int ndpmax)

{
  static const char *function = "wcsinit";

  int status;

  // Check inputs.
  if (wcs == 0x0) return WCSERR_NULL_POINTER;

  if (npvmax < 0) npvmax = wcsnpv(-1);
  if (npsmax < 0) npsmax = wcsnps(-1);


  // Initialize error message handling...
  if (wcs->flag == -1) {
    wcs->err = 0x0;
  }
  struct wcserr **err = &(wcs->err);
  wcserr_clear(err);

  // ...and also in the contained structs in case we have to return due to
  // an error before they can be initialized by their specialized routines,
  // since wcsperr() assumes their wcserr pointers are valid.
  if (wcs->flag == -1) {
    wcs->lin.err = 0x0;
    wcs->cel.err = 0x0;
    wcs->spc.err = 0x0;
  }
  wcserr_clear(&(wcs->lin.err));
  wcserr_clear(&(wcs->cel.err));
  wcserr_clear(&(wcs->spc.err));


  // Initialize pointers.
  if (wcs->flag == -1 || wcs->m_flag != WCSSET) {
    if (wcs->flag == -1) {
      wcs->tab   = 0x0;
      wcs->types = 0x0;
      wcs->lin.flag = -1;
    }

    // Initialize memory management.
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
    wcs->m_czphs = 0x0;
    wcs->m_cperi = 0x0;
    wcs->m_aux   = 0x0;
    wcs->m_tab   = 0x0;
    wcs->m_wtb   = 0x0;
  }

  if (naxis < 0) {
    return wcserr_set(WCSERR_SET(WCSERR_MEMORY),
      "naxis must not be negative (got %d)", naxis);
  }


  // Allocate memory for arrays if required.
  if (alloc ||
     wcs->crpix == 0x0 ||
     wcs->pc    == 0x0 ||
     wcs->cdelt == 0x0 ||
     wcs->crval == 0x0 ||
     wcs->cunit == 0x0 ||
     wcs->ctype == 0x0 ||
     (npvmax && wcs->pv == 0x0) ||
     (npsmax && wcs->ps == 0x0) ||
     wcs->cd    == 0x0 ||
     wcs->crota == 0x0 ||
     wcs->colax == 0x0 ||
     wcs->cname == 0x0 ||
     wcs->crder == 0x0 ||
     wcs->csyer == 0x0 ||
     wcs->czphs == 0x0 ||
     wcs->cperi == 0x0) {

    // Was sufficient allocated previously?
    if (wcs->m_flag == WCSSET &&
       (wcs->m_naxis < naxis  ||
        wcs->npvmax  < npvmax ||
        wcs->npsmax  < npsmax)) {
      // No, free it.
      wcsfree(wcs);
    }

    if (alloc || wcs->crpix == 0x0) {
      if (wcs->m_crpix) {
        // In case the caller fiddled with it.
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
        // In case the caller fiddled with it.
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
        // In case the caller fiddled with it.
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
        // In case the caller fiddled with it.
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
        // In case the caller fiddled with it.
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
        // In case the caller fiddled with it.
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
        // In case the caller fiddled with it.
        wcs->pv = wcs->m_pv;

      } else {
        if (npvmax) {
          if ((wcs->pv = calloc(npvmax, sizeof(struct pvcard))) == 0x0) {
            wcsfree(wcs);
            return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
          }
        } else {
          wcs->pv = 0x0;
        }

        wcs->npvmax  = npvmax;

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_pv    = wcs->pv;
      }
    }

    if (alloc || wcs->ps == 0x0) {
      if (wcs->m_ps) {
        // In case the caller fiddled with it.
        wcs->ps = wcs->m_ps;

      } else {
        if (npsmax) {
          if ((wcs->ps = calloc(npsmax, sizeof(struct pscard))) == 0x0) {
            wcsfree(wcs);
            return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
          }
        } else {
          wcs->ps = 0x0;
        }

        wcs->npsmax  = npsmax;

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_ps    = wcs->ps;
      }
    }

    if (alloc || wcs->cd == 0x0) {
      if (wcs->m_cd) {
        // In case the caller fiddled with it.
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
        // In case the caller fiddled with it.
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
        // In case the caller fiddled with it.
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
        // In case the caller fiddled with it.
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
        // In case the caller fiddled with it.
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
        // In case the caller fiddled with it.
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

    if (alloc || wcs->czphs == 0x0) {
      if (wcs->m_czphs) {
        // In case the caller fiddled with it.
        wcs->czphs = wcs->m_czphs;

      } else {
        if ((wcs->czphs = calloc(naxis, sizeof(double))) == 0x0) {
          wcsfree(wcs);
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_czphs = wcs->czphs;
      }
    }

    if (alloc || wcs->cperi == 0x0) {
      if (wcs->m_cperi) {
        // In case the caller fiddled with it.
        wcs->cperi = wcs->m_cperi;

      } else {
        if ((wcs->cperi = calloc(naxis, sizeof(double))) == 0x0) {
          wcsfree(wcs);
          return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
        }

        wcs->m_flag  = WCSSET;
        wcs->m_naxis = naxis;
        wcs->m_cperi = wcs->cperi;
      }
    }
  }


  wcs->flag  = 0;
  wcs->naxis = naxis;


  // Set defaults for the linear transformation.
  wcs->lin.crpix  = wcs->crpix;
  wcs->lin.pc     = wcs->pc;
  wcs->lin.cdelt  = wcs->cdelt;
  if ((status = lininit(0, naxis, &(wcs->lin), ndpmax))) {
    return wcserr_set(WCS_ERRMSG(wcs_linerr[status]));
  }


  // CRVALia defaults to 0.0.
  for (int i = 0; i < naxis; i++) {
    wcs->crval[i] = 0.0;
  }


  // CUNITia and CTYPEia are blank by default.
  for (int i = 0; i < naxis; i++) {
    memset(wcs->cunit[i], 0, 72);
    memset(wcs->ctype[i], 0, 72);
  }


  // Set defaults for the celestial transformation parameters.
  wcs->lonpole = UNDEFINED;
  wcs->latpole = +90.0;

  // Set defaults for the spectral transformation parameters.
  wcs->restfrq = 0.0;
  wcs->restwav = 0.0;

  // Default parameter values.
  wcs->npv = 0;
  for (int k = 0; k < wcs->npvmax; k++) {
    wcs->pv[k].i = 0;
    wcs->pv[k].m = 0;
    wcs->pv[k].value = 0.0;
  }

  wcs->nps = 0;
  for (int k = 0; k < wcs->npsmax; k++) {
    wcs->ps[k].i = 0;
    wcs->ps[k].m = 0;
    memset(wcs->ps[k].value, 0, 72);
  }

  // Defaults for alternate linear transformations.
  double *cd = wcs->cd;
  for (int i = 0; i < naxis; i++) {
    for (int j = 0; j < naxis; j++) {
      *(cd++) = 0.0;
    }
  }
  for (int i = 0; i < naxis; i++) {
    wcs->crota[i] = 0.0;
  }
  wcs->altlin = 0;
  wcs->velref = 0;

  // Defaults for auxiliary coordinate system information.
  memset(wcs->alt, 0, 4);
  wcs->alt[0] = ' ';
  wcs->colnum = 0;

  for (int i = 0; i < naxis; i++) {
    wcs->colax[i] = 0;
    memset(wcs->cname[i], 0, 72);
    wcs->crder[i] = UNDEFINED;
    wcs->csyer[i] = UNDEFINED;
    wcs->czphs[i] = UNDEFINED;
    wcs->cperi[i] = UNDEFINED;
  }

  memset(wcs->wcsname, 0, 72);

  memset(wcs->timesys,  0, 72);
  memset(wcs->trefpos,  0, 72);
  memset(wcs->trefdir,  0, 72);
  memset(wcs->plephem,  0, 72);

  memset(wcs->timeunit, 0, 72);
  memset(wcs->dateref,  0, 72);
  wcs->mjdref[0]  = UNDEFINED;
  wcs->mjdref[1]  = UNDEFINED;
  wcs->timeoffs   = UNDEFINED;

  memset(wcs->dateobs, 0, 72);
  memset(wcs->datebeg, 0, 72);
  memset(wcs->dateavg, 0, 72);
  memset(wcs->dateend, 0, 72);
  wcs->mjdobs     = UNDEFINED;
  wcs->mjdbeg     = UNDEFINED;
  wcs->mjdavg     = UNDEFINED;
  wcs->mjdend     = UNDEFINED;
  wcs->jepoch     = UNDEFINED;
  wcs->bepoch     = UNDEFINED;
  wcs->tstart     = UNDEFINED;
  wcs->tstop      = UNDEFINED;
  wcs->xposure    = UNDEFINED;
  wcs->telapse    = UNDEFINED;

  wcs->timsyer    = UNDEFINED;
  wcs->timrder    = UNDEFINED;
  wcs->timedel    = UNDEFINED;
  wcs->timepixr   = UNDEFINED;

  wcs->obsgeo[0]  = UNDEFINED;
  wcs->obsgeo[1]  = UNDEFINED;
  wcs->obsgeo[2]  = UNDEFINED;
  wcs->obsgeo[3]  = UNDEFINED;
  wcs->obsgeo[4]  = UNDEFINED;
  wcs->obsgeo[5]  = UNDEFINED;
  memset(wcs->obsorbit, 0, 72);
  memset(wcs->radesys,  0, 72);
  wcs->equinox    = UNDEFINED;
  memset(wcs->specsys,  0, 72);
  memset(wcs->ssysobs,  0, 72);
  wcs->velosys    = UNDEFINED;
  wcs->zsource    = UNDEFINED;
  memset(wcs->ssyssrc,  0, 72);
  wcs->velangl    = UNDEFINED;

  // No additional auxiliary coordinate system information.
  wcs->aux  = 0x0;

  // Tabular parameters.
  wcs->ntab = 0;
  wcs->tab  = 0x0;
  wcs->nwtb = 0;
  wcs->wtb  = 0x0;

  // Reset derived values.
  strcpy(wcs->lngtyp, "    ");
  strcpy(wcs->lattyp, "    ");
  wcs->lng  = -1;
  wcs->lat  = -1;
  wcs->spec = -1;
  wcs->cubeface = -1;

  celini(&(wcs->cel));
  spcini(&(wcs->spc));

  return WCSERR_SUCCESS;
}

//----------------------------------------------------------------------------

int wcsauxi(
  int alloc,
  struct wcsprm *wcs)

{
  static const char *function = "wcsauxi";

  // Check inputs.
  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  struct wcserr **err = &(wcs->err);

  // Allocate memory if required.
  if (alloc || wcs->aux == 0x0) {
    if (wcs->m_aux) {
      // In case the caller fiddled with it.
      wcs->aux = wcs->m_aux;

    } else {
      if ((wcs->aux = malloc(sizeof(struct auxprm))) == 0x0) {
        return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
      }

      wcs->m_aux = wcs->aux;
    }
  }

  struct auxprm *aux = wcs->aux;
  aux->rsun_ref = UNDEFINED;
  aux->dsun_obs = UNDEFINED;
  aux->crln_obs = UNDEFINED;
  aux->hgln_obs = UNDEFINED;
  aux->hglt_obs = UNDEFINED;

  return WCSERR_SUCCESS;
}

//----------------------------------------------------------------------------

int wcssub(
  int alloc,
  const struct wcsprm *wcssrc,
  int *nsub,
  int axes[],
  struct wcsprm *wcsdst)

{
  static const char *function = "wcssub";

  const char *pq = "PQ";
  int  status;

  if (wcssrc == 0x0) return WCSERR_NULL_POINTER;
  if (wcsdst == 0x0) return WCSERR_NULL_POINTER;
  struct wcserr **err = &(wcsdst->err);

  // N.B. we do not rely on the wcsprm struct having been set up.
  int naxis;
  if ((naxis = wcssrc->naxis) <= 0) {
    return wcserr_set(WCSERR_SET(WCSERR_MEMORY),
      "naxis must be positive (got %d)", naxis);
  }

  int dummy;
  if (nsub == 0x0) {
    nsub = &dummy;
    *nsub = naxis;
  } else if (*nsub == 0) {
    *nsub = naxis;
  }

  // Allocate enough temporary storage to hold either axes[] xor map[].
  int *itmp;
  int  ntmp = (*nsub <= naxis) ? naxis : *nsub;
  if ((itmp = calloc(ntmp, sizeof(int))) == 0x0) {
    return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
  }

  int dealloc;
  if ((dealloc = (axes == 0x0))) {
    // Construct an index array.
    if ((axes = calloc(naxis, sizeof(int))) == 0x0) {
      free(itmp);
      return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
    }

    for (int i = 0; i < naxis; i++) {
      axes[i] = i+1;
    }
  }

  // So that we don't try to free uninitialized pointers on cleanup.
  wcsdst->m_aux = 0x0;
  wcsdst->m_tab = 0x0;


  int msub = 0;
  for (int j = 0; j < *nsub; j++) {
    int axis = axes[j];

    if (abs(axis) > 0x1000) {
      // Subimage extraction by type.
      int k = abs(axis) & 0xFF;

      int longitude = k & WCSSUB_LONGITUDE;
      int latitude  = k & WCSSUB_LATITUDE;
      int cubeface  = k & WCSSUB_CUBEFACE;
      int spectral  = k & WCSSUB_SPECTRAL;
      int stokes    = k & WCSSUB_STOKES;
      int time      = k & WCSSUB_TIME;

      int other;
      if ((other = (axis < 0))) {
        longitude = !longitude;
        latitude  = !latitude;
        cubeface  = !cubeface;
        spectral  = !spectral;
        stokes    = !stokes;
        time      = !time;
      }

      for (int i = 0; i < naxis; i++) {
        char ctypei[16];
        strncpy (ctypei, (char *)(wcssrc->ctype + i), 8);
        ctypei[8] = '\0';

        // Find the last non-blank character.
        char *c = ctypei + 8;
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

        } else if (time_type(ctypei)) {
          if (!time) {
            continue;
          }

        } else if (strcmp(ctypei, "STOKES") == 0) {
          if (!stokes) {
            continue;
          }

        } else if (strcmp(ctypei, "CUBEFACE") == 0) {
          if (!cubeface) {
            continue;
          }

        } else if (!other) {
          continue;
        }

        // This axis is wanted, but has it already been added?
        int k;
        for (k = 0; k < msub; k++) {
          if (itmp[k] == i+1) {
            break;
          }
        }
        if (k == msub) itmp[msub++] = i+1;
      }

    } else if (0 < axis && axis <= naxis) {
      // Check that the requested axis has not already been added.
      int k;
      for (k = 0; k < msub; k++) {
        if (itmp[k] == axis) {
          break;
        }
      }
      if (k == msub) itmp[msub++] = axis;

    } else if (axis == 0) {
      // Graft on a new axis.
      itmp[msub++] = 0;

    } else {
      status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_SUBIMAGE));
      goto cleanup;
    }
  }

  if ((*nsub = msub) == 0) {
    // Zero out this struct.
    status = wcsinit(alloc, 0, wcsdst, 0, 0, 0);
    goto cleanup;
  }

  for (int i = 0; i < *nsub; i++) {
    axes[i] = itmp[i];
  }


  // Construct the inverse axis map (i is 0-relative, j is 1-relative):
  // axes[i] == j means that output axis i+1 comes from input axis j,
  // axes[i] == 0 means to create a new axis,
  //  map[i] == j means that input axis i+1 goes to output axis j,
  //  map[i] == 0 means that input axis i+1 is not used.
  int *map = itmp;
  for (int i = 0; i < naxis; i++) {
    map[i] = 0;
  }

  for (int i = 0; i < *nsub; i++) {
    if (axes[i] > 0) {
      map[axes[i]-1] = i+1;
    }
  }


  // Check that the subimage coordinate system is separable.  First check
  // non-zero, off-diagonal elements of the linear transformation matrix.
  double *dstp;
  const double *srcp = wcssrc->pc;
  for (int i = 0; i < naxis; i++) {
    for (int j = 0; j < naxis; j++) {
      if (*(srcp++) == 0.0 || j == i) continue;

      if ((map[i] == 0) != (map[j] == 0)) {
        status = wcserr_set(WCSERR_SET(WCSERR_NON_SEPARABLE),
          "Non-zero off-diagonal matrix elements excluded from the subimage");
        goto cleanup;
      }
    }
  }

  // Tabular coordinates, if any, will be checked below.

  // Now check for distortions that depend on other axes.  As the disprm
  // struct may not have been initialized, we must parse the dpkey entries.
  int ndpmax = 0;
  for (int m = 0; m < 2; m++) {
    struct disprm *dissrc;
    if (m == 0) {
      dissrc = wcssrc->lin.dispre;
    } else {
      dissrc = wcssrc->lin.disseq;
    }

    int ndp = 0;
    if (dissrc != 0x0) {
      for (int j = 0; j < naxis; j++) {
        if (map[j] == 0) continue;

        // Axis numbers in axmap[] are 0-relative.
        int axmap[32];
        for (int jhat = 0; jhat < 32; jhat++) {
          axmap[jhat] = -1;
        }

        int Nhat = 0;
        struct dpkey *dpsrc = dissrc->dp;
        for (int idp = 0; idp < dissrc->ndp; idp++, dpsrc++) {
          // Thorough error checking will be done later by disset().
          if (dpsrc->j != j+1) continue;
          if (dpsrc->field[1] != pq[m]) continue;
          char *fp;
          if ((fp = strchr(dpsrc->field, '.')) == 0x0) continue;
          fp++;

          ndp++;

          if (strncmp(fp, "NAXES", 6) == 0) {
            Nhat = dpkeyi(dpsrc);
          } else if (strncmp(fp, "AXIS.", 5) == 0) {
            int jhat;
            sscanf(fp+5, "%d", &jhat);
            axmap[jhat-1] = dpkeyi(dpsrc) - 1;
          }
        }

        if (Nhat < 0 || (Nhat == 0 && 1 < ndp) || naxis < Nhat || 32 < Nhat) {
          status = wcserr_set(WCSERR_SET(WCSERR_BAD_PARAM),
            "NAXES was not set (or bad) for %s distortion on axis %d",
            dissrc->dtype[j], j+1);
          goto cleanup;
        }

        for (int jhat = 0; jhat < Nhat; jhat++) {
          if (axmap[jhat] < 0) {
            axmap[jhat] = jhat;

            // Make room for an additional DPja.AXIS.j record.
            ndp++;
          }

          if (map[axmap[jhat]] == 0) {
            // Distortion depends on an axis excluded from the subimage.
            status = wcserr_set(WCSERR_SET(WCSERR_NON_SEPARABLE),
              "Distortion depends on an axis excluded from the subimage.");
            goto cleanup;
          }
        }
      }
    }

    if (ndpmax < ndp) ndpmax = ndp;
  }


  // Number of PVi_ma records in the subimage.
  int npvmax = 0;
  for (int m = 0; m < wcssrc->npv; m++) {
    int i = wcssrc->pv[m].i;
    if (i == 0 || (i > 0 && map[i-1])) {
      npvmax++;
    }
  }

  // Number of PSi_ma records in the subimage.
  int npsmax = 0;
  for (int m = 0; m < wcssrc->nps; m++) {
    int i = wcssrc->ps[m].i;
    if (i > 0 && map[i-1]) {
      npsmax++;
    }
  }

  // Initialize the destination.
  status = wcsinit(alloc, *nsub, wcsdst, npvmax, npsmax, ndpmax);

  for (int m = 0; m < 2; m++) {
    struct disprm *dissrc, *disdst;
    if (m == 0) {
      dissrc = wcssrc->lin.dispre;
      disdst = wcsdst->lin.dispre;
    } else {
      dissrc = wcssrc->lin.disseq;
      disdst = wcsdst->lin.disseq;
    }

    if (dissrc && !disdst) {
      if ((disdst = calloc(1, sizeof(struct disprm))) == 0x0) {
        return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
      }

      // Also inits disdst.
      disdst->flag = -1;
      lindist(m+1, &(wcsdst->lin), disdst, ndpmax);
    }
  }

  if (status) {
    goto cleanup;
  }


  // Linear transformation.
  srcp = wcssrc->crpix;
  dstp = wcsdst->crpix;
  for (int j = 0; j < *nsub; j++, dstp++) {
    if (axes[j] > 0) {
      int k = axes[j] - 1;
      *dstp = *(srcp+k);
    }
  }

  srcp = wcssrc->pc;
  dstp = wcsdst->pc;
  for (int i = 0; i < *nsub; i++) {
    for (int j = 0; j < *nsub; j++, dstp++) {
      if (axes[i] > 0 && axes[j] > 0) {
        int k = (axes[i]-1)*naxis + (axes[j]-1);
        *dstp = *(srcp+k);
      }
    }
  }

  srcp = wcssrc->cdelt;
  dstp = wcsdst->cdelt;
  for (int i = 0; i < *nsub; i++, dstp++) {
    if (axes[i] > 0) {
      int k = axes[i] - 1;
      *dstp = *(srcp+k);
    }
  }

  // Coordinate reference value.
  srcp = wcssrc->crval;
  dstp = wcsdst->crval;
  for (int i = 0; i < *nsub; i++, dstp++) {
    if (axes[i] > 0) {
      int k = axes[i] - 1;
      *dstp = *(srcp+k);
    }
  }

  // Coordinate units and type.
  for (int i = 0; i < *nsub; i++) {
    if (axes[i] > 0) {
      int k = axes[i] - 1;
      strncpy(wcsdst->cunit[i], wcssrc->cunit[k], 72);
      strncpy(wcsdst->ctype[i], wcssrc->ctype[k], 72);
    }
  }

  // Celestial and spectral transformation parameters.
  wcsdst->lonpole = wcssrc->lonpole;
  wcsdst->latpole = wcssrc->latpole;
  wcsdst->restfrq = wcssrc->restfrq;
  wcsdst->restwav = wcssrc->restwav;

  // Parameter values.
  int npv = 0;
  for (int m = 0; m < wcssrc->npv; m++) {
    int i = wcssrc->pv[m].i;
    if (i == 0) {
      // i == 0 is a special code that means "the latitude axis".
      wcsdst->pv[npv] = wcssrc->pv[m];
      wcsdst->pv[npv].i = 0;
      npv++;
    } else if (i > 0 && map[i-1]) {
      wcsdst->pv[npv] = wcssrc->pv[m];
      wcsdst->pv[npv].i = map[i-1];
      npv++;
    }
  }
  wcsdst->npv = npv;

  int nps = 0;
  for (int m = 0; m < wcssrc->nps; m++) {
    int i = wcssrc->ps[m].i;
    if (i > 0 && map[i-1]) {
      wcsdst->ps[nps] = wcssrc->ps[m];
      wcsdst->ps[nps].i = map[i-1];
      nps++;
    }
  }
  wcsdst->nps = nps;

  // Alternate linear transformations.
  if (wcssrc->cd) {
    srcp = wcssrc->cd;
    dstp = wcsdst->cd;
    for (int i = 0; i < *nsub; i++) {
      for (int j = 0; j < *nsub; j++, dstp++) {
        if (axes[i] > 0 && axes[j] > 0) {
          int k = (axes[i]-1)*naxis + (axes[j]-1);
          *dstp = *(srcp+k);
        } else if (i == j && wcssrc->altlin & 2) {
          // A new axis is being created where CDi_ja was present in the input
          // header, so override the default value of 0 set by wcsinit().
          *dstp = 1.0;
        }
      }
    }
  }

  if (wcssrc->crota) {
    srcp = wcssrc->crota;
    dstp = wcsdst->crota;
    for (int i = 0; i < *nsub; i++, dstp++) {
      if (axes[i] > 0) {
        int k = axes[i] - 1;
        *dstp = *(srcp+k);
      }
    }
  }

  wcsdst->altlin = wcssrc->altlin;
  wcsdst->velref = wcssrc->velref;

  // Auxiliary coordinate system information.
  strncpy(wcsdst->alt, wcssrc->alt, 4);
  wcsdst->colnum = wcssrc->colnum;

  for (int i = 0; i < *nsub; i++) {
    if (axes[i] > 0) {
      int k = axes[i] - 1;
      if (wcssrc->colax) wcsdst->colax[i] = wcssrc->colax[k];
      if (wcssrc->cname) strncpy(wcsdst->cname[i], wcssrc->cname[k], 72);
      if (wcssrc->crder) wcsdst->crder[i] = wcssrc->crder[k];
      if (wcssrc->csyer) wcsdst->csyer[i] = wcssrc->csyer[k];
      if (wcssrc->czphs) wcsdst->czphs[i] = wcssrc->czphs[k];
      if (wcssrc->cperi) wcsdst->cperi[i] = wcssrc->cperi[k];
    }
  }

  strncpy(wcsdst->wcsname, wcssrc->wcsname, 72);

  strncpy(wcsdst->timesys, wcssrc->timesys, 72);
  strncpy(wcsdst->trefpos, wcssrc->trefpos, 72);
  strncpy(wcsdst->trefdir, wcssrc->trefdir, 72);
  strncpy(wcsdst->plephem, wcssrc->plephem, 72);

  strncpy(wcsdst->timeunit, wcssrc->timeunit, 72);
  strncpy(wcsdst->dateref,  wcssrc->dateref, 72);
  wcsdst->mjdref[0] = wcssrc->mjdref[0];
  wcsdst->mjdref[1] = wcssrc->mjdref[1];
  wcsdst->timeoffs  = wcssrc->timeoffs;

  strncpy(wcsdst->dateobs, wcssrc->dateobs, 72);
  strncpy(wcsdst->datebeg, wcssrc->datebeg, 72);
  strncpy(wcsdst->dateavg, wcssrc->dateavg, 72);
  strncpy(wcsdst->dateend, wcssrc->dateend, 72);

  wcsdst->mjdobs  = wcssrc->mjdobs;
  wcsdst->mjdbeg  = wcssrc->mjdbeg;
  wcsdst->mjdavg  = wcssrc->mjdavg;
  wcsdst->mjdend  = wcssrc->mjdend;
  wcsdst->jepoch  = wcssrc->jepoch;
  wcsdst->bepoch  = wcssrc->bepoch;
  wcsdst->tstart  = wcssrc->tstart;
  wcsdst->tstop   = wcssrc->tstop;
  wcsdst->xposure = wcssrc->xposure;
  wcsdst->telapse = wcssrc->telapse;

  wcsdst->timsyer  = wcssrc->timsyer;
  wcsdst->timrder  = wcssrc->timrder;
  wcsdst->timedel  = wcssrc->timedel;
  wcsdst->timepixr = wcssrc->timepixr;

  wcsdst->obsgeo[0] = wcssrc->obsgeo[0];
  wcsdst->obsgeo[1] = wcssrc->obsgeo[1];
  wcsdst->obsgeo[2] = wcssrc->obsgeo[2];
  wcsdst->obsgeo[3] = wcssrc->obsgeo[3];
  wcsdst->obsgeo[4] = wcssrc->obsgeo[4];
  wcsdst->obsgeo[5] = wcssrc->obsgeo[5];

  strncpy(wcsdst->obsorbit, wcssrc->obsorbit, 72);
  strncpy(wcsdst->radesys,  wcssrc->radesys, 72);
  wcsdst->equinox = wcssrc->equinox;
  strncpy(wcsdst->specsys,  wcssrc->specsys, 72);
  strncpy(wcsdst->ssysobs,  wcssrc->ssysobs, 72);
  wcsdst->velosys = wcssrc->velosys;
  wcsdst->zsource = wcssrc->zsource;
  strncpy(wcsdst->ssyssrc,  wcssrc->ssyssrc, 72);
  wcsdst->velangl = wcssrc->velangl;


  // Additional auxiliary coordinate system information.
  if (wcssrc->aux && !wcsdst->aux) {
    if ((wcsdst->aux = calloc(1, sizeof(struct auxprm))) == 0x0) {
      status = wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
      goto cleanup;
    }

    wcsdst->m_aux = wcsdst->aux;

    wcsdst->aux->rsun_ref = wcssrc->aux->rsun_ref;
    wcsdst->aux->dsun_obs = wcssrc->aux->dsun_obs;
    wcsdst->aux->crln_obs = wcssrc->aux->crln_obs;
    wcsdst->aux->hgln_obs = wcssrc->aux->hgln_obs;
    wcsdst->aux->hglt_obs = wcssrc->aux->hglt_obs;
  }


  // Coordinate lookup tables; only copy what's needed.
  wcsdst->ntab = 0;
  for (int itab = 0; itab < wcssrc->ntab; itab++) {
    // Is this table wanted?
    for (int m = 0; m < wcssrc->tab[itab].M; m++) {
      int i = wcssrc->tab[itab].map[m];

      if (map[i]) {
        wcsdst->ntab++;
        break;
      }
    }
  }

  if (wcsdst->ntab) {
    // Allocate memory for tabprm structs.
    if ((wcsdst->tab = calloc(wcsdst->ntab, sizeof(struct tabprm))) == 0x0) {
      wcsdst->ntab = 0;

      status = wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
      goto cleanup;
    }

    wcsdst->m_tab = wcsdst->tab;
  }

  struct tabprm *tab = wcsdst->tab;
  for (int itab = 0; itab < wcssrc->ntab; itab++) {
    for (int m = 0; m < wcssrc->tab[itab].M; m++) {
      int i = wcssrc->tab[itab].map[m];

      if (map[i]) {
        tab->flag = -1;
        if ((status = tabcpy(1, wcssrc->tab + itab, tab))) {
          wcserr_set(WCS_ERRMSG(wcs_taberr[status]));
          goto cleanup;
        }

        // Translate axis numbers in tab->map[].
        for (int m = 0; m < wcssrc->tab[itab].M; m++) {
          // Table axis mapping, followed by...
          int i = wcssrc->tab[itab].map[m];

          // ...subimaging axis mapping.
          int j = map[i] - 1;
          if (j < 0) {
            // In general, multi-dimensional tables (i.e. with M > 1) are not
            // separable, so if one axis is selected then all must be.
            status = wcserr_set(WCSERR_SET(WCSERR_NON_SEPARABLE),
              "Table with M>1 depends on axis excluded from the subimage");
            goto cleanup;
          }

          tab->map[m] = j;
        }

        tab++;
        break;
      }
    }
  }


  // Distortion parameters (in linprm).
  for (int m = 0; m < 2; m++) {
    struct disprm *dissrc, *disdst;
    if (m == 0) {
      dissrc = wcssrc->lin.dispre;
      disdst = wcsdst->lin.dispre;
    } else {
      dissrc = wcssrc->lin.disseq;
      disdst = wcsdst->lin.disseq;
    }

    if (dissrc) {
      disdst->naxis = *nsub;

      // Distortion type and maximum distortion (but not total distortion).
      for (int j = 0; j < *nsub; j++) {
        if (axes[j] > 0) {
          int k = axes[j] - 1;
          strncpy(disdst->dtype[j], dissrc->dtype[k], 72);
          disdst->maxdis[j] = dissrc->maxdis[k];
        }
      }

      // DPja or DQia keyvalues.
      int ndp = 0;
      struct dpkey *dpdst = disdst->dp;
      for (int j = 0; j < *nsub; j++) {
        if (axes[j] == 0) continue;

        // Determine the axis mapping.
        int axmap[32];
        for (int jhat = 0; jhat < 32; jhat++) {
          axmap[jhat] = -1;
        }

        int Nhat = 0;
        struct dpkey *dpsrc = dissrc->dp;
        for (int idp = 0; idp < dissrc->ndp; idp++, dpsrc++) {
          if (dpsrc->j != axes[j]) continue;
          if (dpsrc->field[1] != pq[m]) continue;
          char *fp;
          if ((fp = strchr(dpsrc->field, '.')) == 0x0) continue;
          fp++;

          if (strncmp(fp, "NAXES", 6) == 0) {
            Nhat = dpkeyi(dpsrc);
          } else if (strncmp(fp, "AXIS.", 5) == 0) {
            int jhat;
            sscanf(fp+5, "%d", &jhat);
            axmap[jhat-1] = dpkeyi(dpsrc) - 1;
          }
        }

        for (int jhat = 0; jhat < Nhat; jhat++) {
          if (axmap[jhat] < 0) {
            axmap[jhat] = jhat;
          }
        }

        // Copy the DPja or DQia keyvalues.
        dpsrc = dissrc->dp;
        for (int idp = 0; idp < dissrc->ndp; idp++, dpsrc++) {
          if (dpsrc->j != axes[j]) continue;
          if (dpsrc->field[1] != pq[m]) continue;
          char *fp;
          if ((fp = strchr(dpsrc->field, '.')) == 0x0) continue;
          fp++;

          if (strncmp(fp, "AXIS.", 5) == 0) {
            // Skip it, we will create our own later.
            continue;
          }

          *dpdst = *dpsrc;
          char ctmp[16];
          sprintf(ctmp, "%d", j+1);
          dpdst->field[2] = ctmp[0];
          dpdst->j = j+1;

          ndp++;
          dpdst++;

          if (strncmp(fp, "NAXES", 6) == 0) {
            for (int jhat = 0; jhat < Nhat; jhat++) {
              strcpy(dpdst->field, dpsrc->field);
              dpdst->field[2] = ctmp[0];
              fp = strchr(dpdst->field, '.') + 1;
              sprintf(fp, "AXIS.%d", jhat+1);
              dpdst->j = j+1;
              dpdst->type = 0;
              dpdst->value.i = map[axmap[jhat]];

              ndp++;
              dpdst++;
            }
          }
        }
      }

      disdst->ndp = ndp;
    }
  }


cleanup:
  if (itmp) free(itmp);
  if (dealloc) {
    free(axes);
  }

  if (status && wcsdst->m_aux) {
    free(wcsdst->m_aux);
    wcsdst->aux   = 0x0;
    wcsdst->m_aux = 0x0;
  }

  if (status && wcsdst->m_tab) {
    tabfree(wcsdst->m_tab);
  }

  return status;
}

//----------------------------------------------------------------------------

int wcscompare(
  int cmp,
  double tol,
  const struct wcsprm *wcs1,
  const struct wcsprm *wcs2,
  int *equal)

{
  int status;

  if (wcs1  == 0x0) return WCSERR_NULL_POINTER;
  if (wcs2  == 0x0) return WCSERR_NULL_POINTER;
  if (equal == 0x0) return WCSERR_NULL_POINTER;

  *equal = 0;

  if (wcs1->naxis != wcs2->naxis) {
    return 0;
  }

  int naxis = wcs1->naxis;
  int naxis2 = wcs1->naxis*wcs1->naxis;

  if (cmp & WCSCOMPARE_CRPIX) {
    // Don't compare crpix.
  } else if (cmp & WCSCOMPARE_TILING) {
    for (int i = 0; i < naxis; ++i) {
      double diff = wcs1->crpix[i] - wcs2->crpix[i];
      if ((double)(int)(diff) != diff) {
        return 0;
      }
    }
  } else {
    if (!wcsutil_dblEq(naxis, tol, wcs1->crpix, wcs2->crpix)) {
      return 0;
    }
  }

  if (!wcsutil_dblEq(naxis2, tol, wcs1->pc, wcs2->pc) ||
      !wcsutil_dblEq(naxis, tol, wcs1->cdelt, wcs2->cdelt) ||
      !wcsutil_dblEq(naxis, tol, wcs1->crval, wcs2->crval) ||
      !wcsutil_strEq(naxis, wcs1->cunit, wcs2->cunit) ||
      !wcsutil_strEq(naxis, wcs1->ctype, wcs2->ctype) ||
      !wcsutil_dblEq(1, tol, &wcs1->lonpole, &wcs2->lonpole) ||
      !wcsutil_dblEq(1, tol, &wcs1->latpole, &wcs2->latpole) ||
      !wcsutil_dblEq(1, tol, &wcs1->restfrq, &wcs2->restfrq) ||
      !wcsutil_dblEq(1, tol, &wcs1->restwav, &wcs2->restwav) ||
      wcs1->npv != wcs2->npv ||
      wcs1->nps != wcs2->nps) {
    return 0;
  }

  // Compare pv cards, which may not be in the same order
  for (int i = 0; i < wcs1->npv; ++i) {
    int j;
    for (j = 0; j < wcs2->npv; ++j) {
      if (wcs1->pv[i].i == wcs2->pv[j].i &&
          wcs1->pv[i].m == wcs2->pv[j].m) {
        if (!wcsutil_dblEq(1, tol, &wcs1->pv[i].value, &wcs2->pv[j].value)) {
          return 0;
        }
        break;
      }
    }
    // We didn't find a match, so they are not equal
    if (j == wcs2->npv) {
      return 0;
    }
  }

  // Compare ps cards, which may not be in the same order
  for (int i = 0; i < wcs1->nps; ++i) {
    int j;
    for (j = 0; j < wcs2->nps; ++j) {
      if (wcs1->ps[i].i == wcs2->ps[j].i &&
          wcs1->ps[i].m == wcs2->ps[j].m) {
        if (strncmp(wcs1->ps[i].value, wcs2->ps[j].value, 72)) {
          return 0;
        }
        break;
      }
    }
    // We didn't find a match, so they are not equal
    if (j == wcs2->nps) {
      return 0;
    }
  }

  if (wcs1->flag != WCSSET || wcs2->flag != WCSSET) {
    if (!wcsutil_dblEq(naxis2, tol, wcs1->cd, wcs2->cd) ||
        !wcsutil_dblEq(naxis, tol, wcs1->crota, wcs2->crota) ||
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
        !wcsutil_dblEq(naxis, tol, wcs1->crder, wcs2->crder) ||
        !wcsutil_dblEq(naxis, tol, wcs1->csyer, wcs2->csyer) ||
        !wcsutil_dblEq(naxis, tol, wcs1->czphs, wcs2->czphs) ||
        !wcsutil_dblEq(naxis, tol, wcs1->cperi, wcs2->cperi) ||
        strncmp(wcs1->wcsname,  wcs2->wcsname,  72) ||
        strncmp(wcs1->timesys,  wcs2->timesys,  72) ||
        strncmp(wcs1->trefpos,  wcs2->trefpos,  72) ||
        strncmp(wcs1->trefdir,  wcs2->trefdir,  72) ||
        strncmp(wcs1->plephem,  wcs2->plephem,  72) ||
        strncmp(wcs1->timeunit, wcs2->timeunit, 72) ||
        strncmp(wcs1->dateref,  wcs2->dateref,  72) ||
        !wcsutil_dblEq(2, tol,  wcs1->mjdref,    wcs2->mjdref)   ||
        !wcsutil_dblEq(1, tol, &wcs1->timeoffs, &wcs2->timeoffs) ||
        strncmp(wcs1->dateobs,  wcs2->dateobs, 72) ||
        strncmp(wcs1->datebeg,  wcs2->datebeg, 72) ||
        strncmp(wcs1->dateavg,  wcs2->dateavg, 72) ||
        strncmp(wcs1->dateend,  wcs2->dateend, 72) ||
        !wcsutil_dblEq(1, tol, &wcs1->mjdobs,   &wcs2->mjdobs)   ||
        !wcsutil_dblEq(1, tol, &wcs1->mjdbeg,   &wcs2->mjdbeg)   ||
        !wcsutil_dblEq(1, tol, &wcs1->mjdavg,   &wcs2->mjdavg)   ||
        !wcsutil_dblEq(1, tol, &wcs1->mjdend,   &wcs2->mjdend)   ||
        !wcsutil_dblEq(1, tol, &wcs1->jepoch,   &wcs2->jepoch)   ||
        !wcsutil_dblEq(1, tol, &wcs1->bepoch,   &wcs2->bepoch)   ||
        !wcsutil_dblEq(1, tol, &wcs1->tstart,   &wcs2->tstart)   ||
        !wcsutil_dblEq(1, tol, &wcs1->tstop,    &wcs2->tstop)    ||
        !wcsutil_dblEq(1, tol, &wcs1->xposure,  &wcs2->xposure)  ||
        !wcsutil_dblEq(1, tol, &wcs1->telapse,  &wcs2->telapse)  ||
        !wcsutil_dblEq(1, tol, &wcs1->timsyer,  &wcs2->timsyer)  ||
        !wcsutil_dblEq(1, tol, &wcs1->timrder,  &wcs2->timrder)  ||
        !wcsutil_dblEq(1, tol, &wcs1->timedel,  &wcs2->timedel)  ||
        !wcsutil_dblEq(1, tol, &wcs1->timepixr, &wcs2->timepixr) ||
        !wcsutil_dblEq(6, tol,  wcs1->obsgeo,    wcs2->obsgeo)   ||
        strncmp(wcs1->obsorbit, wcs2->obsorbit, 72) ||
        strncmp(wcs1->radesys,  wcs2->radesys,  72) ||
        !wcsutil_dblEq(1, tol, &wcs1->equinox,  &wcs2->equinox)  ||
        strncmp(wcs1->specsys,  wcs2->specsys,  72) ||
        strncmp(wcs1->ssysobs,  wcs2->ssysobs,  72) ||
        !wcsutil_dblEq(1, tol, &wcs1->velosys,  &wcs2->velosys)  ||
        !wcsutil_dblEq(1, tol, &wcs1->zsource,  &wcs2->zsource)  ||
        strncmp(wcs1->ssyssrc,  wcs2->ssyssrc,  72) ||
        !wcsutil_dblEq(1, tol, &wcs1->velangl,  &wcs2->velangl)) {
      return 0;
    }

    // Compare additional auxiliary parameters.
    if (wcs1->aux && wcs2->aux) {
      if (!wcsutil_dblEq(1, tol, &wcs1->aux->rsun_ref, &wcs2->aux->rsun_ref) ||
          !wcsutil_dblEq(1, tol, &wcs1->aux->dsun_obs, &wcs2->aux->dsun_obs) ||
          !wcsutil_dblEq(1, tol, &wcs1->aux->crln_obs, &wcs2->aux->crln_obs) ||
          !wcsutil_dblEq(1, tol, &wcs1->aux->hgln_obs, &wcs2->aux->hgln_obs) ||
          !wcsutil_dblEq(1, tol, &wcs1->aux->hglt_obs, &wcs2->aux->hglt_obs)) {
        return 0;
      }
    } else if (wcs1->aux || wcs2->aux) {
      return 0;
    }
  }

  // Compare tabular parameters
  if (wcs1->ntab != wcs2->ntab) {
    return 0;
  }

  for (int i = 0; i < wcs1->ntab; ++i) {
    int tab_equal;
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

//----------------------------------------------------------------------------

int wcsfree(struct wcsprm *wcs)

{
  if (wcs == 0x0) return WCSERR_NULL_POINTER;

  if (wcs->flag == -1) {
    wcs->lin.flag = -1;

  } else {
    // Optionally allocated by wcsinit() for given parameters.
    if (wcs->m_flag == WCSSET) {
      // Start by cleaning the slate.
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
      if (wcs->czphs == wcs->m_czphs) wcs->czphs = 0x0;
      if (wcs->cperi == wcs->m_cperi) wcs->cperi = 0x0;

      if (wcs->aux   == wcs->m_aux)   wcs->aux   = 0x0;
      if (wcs->tab   == wcs->m_tab)   wcs->tab   = 0x0;
      if (wcs->wtb   == wcs->m_wtb)   wcs->wtb   = 0x0;

      // Now release the memory.
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
      if (wcs->m_czphs)  free(wcs->m_czphs);
      if (wcs->m_cperi)  free(wcs->m_cperi);

      // May have been allocated by wcspih() or wcssub().
      if (wcs->m_aux) free(wcs->m_aux);

      // Allocated unconditionally by wcstab().
      if (wcs->m_tab) {
        for (int itab = 0; itab < wcs->ntab; itab++) {
          tabfree(wcs->m_tab + itab);
        }

        free(wcs->m_tab);
      }
      if (wcs->m_wtb) free(wcs->m_wtb);
    }

    // Allocated unconditionally by wcsset().
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
  wcs->m_czphs  = 0x0;
  wcs->m_cperi  = 0x0;

  wcs->m_aux    = 0x0;

  wcs->ntab  = 0;
  wcs->m_tab = 0x0;
  wcs->nwtb  = 0;
  wcs->m_wtb = 0x0;

  wcs->types = 0x0;

  wcserr_clear(&(wcs->err));

  wcs->flag = 0;

  linfree(&(wcs->lin));
  celfree(&(wcs->cel));
  spcfree(&(wcs->spc));

  return WCSERR_SUCCESS;
}

//----------------------------------------------------------------------------

int wcstrim(struct wcsprm *wcs)

{
  if (wcs == 0x0) return WCSERR_NULL_POINTER;

  if (wcs->m_flag != WCSSET) {
    // Nothing to do.
    return WCSERR_SUCCESS;
  }

  if (wcs->flag != WCSSET) {
    return WCSERR_UNSET;
  }

  if (wcs->npv < wcs->npvmax) {
    if (wcs->m_pv) {
      if (wcs->npv == 0) {
        free(wcs->m_pv);
        wcs->pv = wcs->m_pv = 0x0;
      } else {
        size_t size = wcs->npv * sizeof(struct pvcard);
        // No error if realloc() fails, it will leave the array untouched.
        if ((wcs->pv = wcs->m_pv = realloc(wcs->m_pv, size))) {
          wcs->npvmax = wcs->npv;
        }
      }
    }
  }

  if (wcs->nps < wcs->npsmax) {
    if (wcs->m_ps) {
      if (wcs->nps == 0) {
        free(wcs->m_ps);
        wcs->ps = wcs->m_ps = 0x0;
      } else {
        size_t size = wcs->nps * sizeof(struct pscard);
        // No error if realloc() fails, it will leave the array untouched.
        if ((wcs->ps = wcs->m_ps = realloc(wcs->m_ps, size))) {
          wcs->npsmax = wcs->nps;
        }
      }
    }
  }

  if (!(wcs->altlin & 2)) {
    if (wcs->m_cd) {
      free(wcs->m_cd);
      wcs->cd = wcs->m_cd = 0x0;
    }
  }

  if (!(wcs->altlin & 4)) {
    if (wcs->m_crota) {
      free(wcs->m_crota);
      wcs->crota = wcs->m_crota = 0x0;
    }
  }

  if (wcs->colax) {
    if (wcsutil_all_ival(wcs->naxis, 0, wcs->colax)) {
      free(wcs->m_colax);
      wcs->colax = wcs->m_colax = 0x0;
    }
  }

  if (wcs->cname) {
    if (wcsutil_all_sval(wcs->naxis, "", (const char (*)[72])wcs->cname)) {
      free(wcs->m_cname);
      wcs->cname = wcs->m_cname = 0x0;
    }
  }

  if (wcs->crder) {
    if (wcsutil_all_dval(wcs->naxis, UNDEFINED, wcs->crder)) {
      free(wcs->m_crder);
      wcs->crder = wcs->m_crder = 0x0;
    }
  }

  if (wcs->csyer) {
    if (wcsutil_all_dval(wcs->naxis, UNDEFINED, wcs->csyer)) {
      free(wcs->m_csyer);
      wcs->csyer = wcs->m_csyer = 0x0;
    }
  }

  if (wcs->czphs) {
    if (wcsutil_all_dval(wcs->naxis, UNDEFINED, wcs->czphs)) {
      free(wcs->m_czphs);
      wcs->czphs = wcs->m_czphs = 0x0;
    }
  }

  if (wcs->cperi) {
    if (wcsutil_all_dval(wcs->naxis, UNDEFINED, wcs->cperi)) {
      free(wcs->m_cperi);
      wcs->cperi = wcs->m_cperi = 0x0;
    }
  }

  return WCSERR_SUCCESS;
}

//----------------------------------------------------------------------------

int wcssize(const struct wcsprm *wcs, int sizes[2])

{
  if (wcs == 0x0) {
    sizes[0] = sizes[1] = 0;
    return WCSERR_SUCCESS;
  }

  // Base size, in bytes.
  sizes[0] = sizeof(struct wcsprm);

  // Total size of allocated memory, in bytes.
  sizes[1] = 0;

  int exsizes[2];
  int naxis = wcs->naxis;

  // wcsprm::crpix[].
  sizes[1] += naxis * sizeof(double);

  // wcsprm::pc[].
  sizes[1] += naxis*naxis * sizeof(double);

  // wcsprm::cdelt[].
  sizes[1] += naxis * sizeof(double);

  // wcsprm::crval[].
  sizes[1] += naxis * sizeof(double);

  // wcsprm::cunit[].
  if (wcs->cunit) {
    sizes[1] += naxis * sizeof(char [72]);
  }

  // wcsprm::ctype[].
  sizes[1] += naxis * sizeof(char [72]);

  // wcsprm::pv[].
  if (wcs->pv) {
    sizes[1] += wcs->npvmax * sizeof(struct pvcard);
  }

  // wcsprm::ps[].
  if (wcs->ps) {
    sizes[1] += wcs->npsmax * sizeof(struct pscard);
  }

  // wcsprm::cd[].
  if (wcs->cd) {
    sizes[1] += naxis*naxis * sizeof(double);
  }

  // wcsprm::crota[].
  if (wcs->crota) {
    sizes[1] += naxis * sizeof(double);
  }

  // wcsprm::colax[].
  if (wcs->colax) {
    sizes[1] += naxis * sizeof(int);
  }

  // wcsprm::cname[].
  if (wcs->cname) {
    sizes[1] += naxis * sizeof(char [72]);
  }

  // wcsprm::crder[].
  if (wcs->crder) {
    sizes[1] += naxis * sizeof(double);
  }

  // wcsprm::csyer[].
  if (wcs->csyer) {
    sizes[1] += naxis * sizeof(double);
  }

  // wcsprm::czphs[].
  if (wcs->czphs) {
    sizes[1] += naxis * sizeof(double);
  }

  // wcsprm::cperi[].
  if (wcs->cperi) {
    sizes[1] += naxis * sizeof(double);
  }

  // wcsprm::aux.
  if (wcs->aux) {
    sizes[1] += sizeof(struct auxprm);
  }

  // wcsprm::tab.
  for (int itab = 0; itab < wcs->ntab; itab++) {
    tabsize(wcs->tab + itab, exsizes);
    sizes[1] += exsizes[0] + exsizes[1];
  }

  // wcsprm::wtb.
  if (wcs->wtb) {
    sizes[1] += wcs->nwtb * sizeof(struct wtbarr);
  }

  // wcsprm::lin.
  linsize(&(wcs->lin), exsizes);
  sizes[1] += exsizes[1];

  // wcsprm::err.
  wcserr_size(wcs->err, exsizes);
  sizes[1] += exsizes[0] + exsizes[1];

  return WCSERR_SUCCESS;
}

//----------------------------------------------------------------------------

int auxsize(const struct auxprm *aux, int sizes[2])

{
  if (aux == 0x0) {
    sizes[0] = sizes[1] = 0;
    return WCSERR_SUCCESS;
  }

  // Base size, in bytes.
  sizes[0] = sizeof(struct auxprm);

  // Total size of allocated memory, in bytes.
  sizes[1] = 0;

  return WCSERR_SUCCESS;
}


//----------------------------------------------------------------------------

static void wcsprt_auxc(const char *name, const char *value)
{
  if (value[0] == '\0') {
    wcsprintf("   %s: UNDEFINED\n", name);
  } else {
    wcsprintf("   %s: \"%s\"\n", name, value);
  }
}

static void wcsprt_auxd(const char *name, double value)
{
  if (undefined(value)) {
    wcsprintf("   %s: UNDEFINED\n", name);
  } else {
    wcsprintf("   %s:  %15.9f\n", name, value);
  }
}

int wcsprt(const struct wcsprm *wcs)

{
  if (wcs == 0x0) return WCSERR_NULL_POINTER;

  if (wcs->flag != WCSSET) {
    wcsprintf("The wcsprm struct is UNINITIALIZED.\n");
    return WCSERR_SUCCESS;
  }

  wcsprintf("       flag: %d\n", wcs->flag);
  wcsprintf("      naxis: %d\n", wcs->naxis);
  WCSPRINTF_PTR("      crpix: ", wcs->crpix, "\n");
  wcsprintf("            ");
  for (int i = 0; i < wcs->naxis; i++) {
    wcsprintf("  %#- 11.5g", wcs->crpix[i]);
  }
  wcsprintf("\n");

  // Linear transformation.
  int k = 0;
  WCSPRINTF_PTR("         pc: ", wcs->pc, "\n");
  for (int i = 0; i < wcs->naxis; i++) {
    wcsprintf("    pc[%d][]:", i);
    for (int j = 0; j < wcs->naxis; j++) {
      wcsprintf("  %#- 11.5g", wcs->pc[k++]);
    }
    wcsprintf("\n");
  }

  // Coordinate increment at reference point.
  WCSPRINTF_PTR("      cdelt: ", wcs->cdelt, "\n");
  wcsprintf("            ");
  for (int i = 0; i < wcs->naxis; i++) {
    wcsprintf("  %#- 11.5g", wcs->cdelt[i]);
  }
  wcsprintf("\n");

  // Coordinate value at reference point.
  WCSPRINTF_PTR("      crval: ", wcs->crval, "\n");
  wcsprintf("            ");
  for (int i = 0; i < wcs->naxis; i++) {
    wcsprintf("  %#- 11.5g", wcs->crval[i]);
  }
  wcsprintf("\n");

  // Coordinate units and type.
  WCSPRINTF_PTR("      cunit: ", wcs->cunit, "\n");
  for (int i = 0; i < wcs->naxis; i++) {
    wcsprintf("             \"%s\"\n", wcs->cunit[i]);
  }

  WCSPRINTF_PTR("      ctype: ", wcs->ctype, "\n");
  for (int i = 0; i < wcs->naxis; i++) {
    wcsprintf("             \"%s\"\n", wcs->ctype[i]);
  }

  // Celestial and spectral transformation parameters.
  if (undefined(wcs->lonpole)) {
    wcsprintf("    lonpole: UNDEFINED\n");
  } else {
    wcsprintf("    lonpole: %9f\n", wcs->lonpole);
  }
  wcsprintf("    latpole: %9f\n", wcs->latpole);
  wcsprintf("    restfrq: %f\n", wcs->restfrq);
  wcsprintf("    restwav: %f\n", wcs->restwav);

  // Parameter values.
  wcsprintf("        npv: %d\n", wcs->npv);
  wcsprintf("     npvmax: %d\n", wcs->npvmax);
  WCSPRINTF_PTR("         pv: ", wcs->pv, "\n");
  for (int k = 0; k < wcs->npv; k++) {
    wcsprintf("             %3d%4d  %#- 11.5g\n", (wcs->pv[k]).i,
      (wcs->pv[k]).m, (wcs->pv[k]).value);
  }
  wcsprintf("        nps: %d\n", wcs->nps);
  wcsprintf("     npsmax: %d\n", wcs->npsmax);
  WCSPRINTF_PTR("         ps: ", wcs->ps, "\n");
  for (int k = 0; k < wcs->nps; k++) {
    wcsprintf("             %3d%4d  %s\n", (wcs->ps[k]).i,
      (wcs->ps[k]).m, (wcs->ps[k]).value);
  }

  // Alternate linear transformations.
  k = 0;
  WCSPRINTF_PTR("         cd: ", wcs->cd, "\n");
  if (wcs->cd) {
    for (int i = 0; i < wcs->naxis; i++) {
      wcsprintf("    cd[%d][]:", i);
      for (int j = 0; j < wcs->naxis; j++) {
        wcsprintf("  %#- 11.5g", wcs->cd[k++]);
      }
      wcsprintf("\n");
    }
  }

  WCSPRINTF_PTR("      crota: ", wcs->crota, "\n");
  if (wcs->crota) {
    wcsprintf("            ");
    for (int i = 0; i < wcs->naxis; i++) {
      wcsprintf("  %#- 11.5g", wcs->crota[i]);
    }
    wcsprintf("\n");
  }

  wcsprintf("     altlin: %d\n", wcs->altlin);
  wcsprintf("     velref: %d\n", wcs->velref);



  // Auxiliary coordinate system information.
  wcsprintf("        alt: '%c'\n", wcs->alt[0]);
  wcsprintf("     colnum: %d\n", wcs->colnum);

  WCSPRINTF_PTR("      colax: ", wcs->colax, "\n");
  if (wcs->colax) {
    wcsprintf("           ");
    for (int i = 0; i < wcs->naxis; i++) {
      wcsprintf("  %5d", wcs->colax[i]);
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("      cname: ", wcs->cname, "\n");
  if (wcs->cname) {
    for (int i = 0; i < wcs->naxis; i++) {
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
    for (int i = 0; i < wcs->naxis; i++) {
      if (undefined(wcs->crder[i])) {
        wcsprintf("    UNDEFINED");
      } else {
        wcsprintf("  %#- 11.5g", wcs->crder[i]);
      }
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("      csyer: ", wcs->csyer, "\n");
  if (wcs->csyer) {
    wcsprintf("           ");
    for (int i = 0; i < wcs->naxis; i++) {
      if (undefined(wcs->csyer[i])) {
        wcsprintf("    UNDEFINED");
      } else {
        wcsprintf("  %#- 11.5g", wcs->csyer[i]);
      }
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("      czphs: ", wcs->czphs, "\n");
  if (wcs->czphs) {
    wcsprintf("           ");
    for (int i = 0; i < wcs->naxis; i++) {
      if (undefined(wcs->czphs[i])) {
        wcsprintf("    UNDEFINED");
      } else {
        wcsprintf("  %#- 11.5g", wcs->czphs[i]);
      }
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("      cperi: ", wcs->cperi, "\n");
  if (wcs->cperi) {
    wcsprintf("           ");
    for (int i = 0; i < wcs->naxis; i++) {
      if (undefined(wcs->cperi[i])) {
        wcsprintf("    UNDEFINED");
      } else {
        wcsprintf("  %#- 11.5g", wcs->cperi[i]);
      }
    }
    wcsprintf("\n");
  }

  wcsprt_auxc(" wcsname", wcs->wcsname);

  wcsprt_auxc(" timesys", wcs->timesys);
  wcsprt_auxc(" trefpos", wcs->trefpos);
  wcsprt_auxc(" trefdir", wcs->trefdir);
  wcsprt_auxc(" plephem", wcs->plephem);
  wcsprt_auxc("timeunit", wcs->timeunit);
  wcsprt_auxc(" dateref", wcs->dateref);
  wcsprintf("     mjdref: ");
  for (int k = 0; k < 2; k++) {
    if (undefined(wcs->mjdref[k])) {
      wcsprintf("       UNDEFINED");
    } else {
      wcsprintf(" %15.9f", wcs->mjdref[k]);
    }
  }
  wcsprintf("\n");
  wcsprt_auxd("timeoffs", wcs->timeoffs);

  wcsprt_auxc(" dateobs", wcs->dateobs);
  wcsprt_auxc(" datebeg", wcs->datebeg);
  wcsprt_auxc(" dateavg", wcs->dateavg);
  wcsprt_auxc(" dateend", wcs->dateend);
  wcsprt_auxd("  mjdobs", wcs->mjdobs);
  wcsprt_auxd("  mjdbeg", wcs->mjdbeg);
  wcsprt_auxd("  mjdavg", wcs->mjdavg);
  wcsprt_auxd("  mjdend", wcs->mjdend);
  wcsprt_auxd("  jepoch", wcs->jepoch);
  wcsprt_auxd("  bepoch", wcs->bepoch);
  wcsprt_auxd("  tstart", wcs->tstart);
  wcsprt_auxd("   tstop", wcs->tstop);
  wcsprt_auxd(" xposure", wcs->xposure);
  wcsprt_auxd(" telapse", wcs->telapse);


  wcsprt_auxd(" timsyer", wcs->timsyer);
  wcsprt_auxd(" timrder", wcs->timrder);
  wcsprt_auxd(" timedel", wcs->timedel);
  wcsprt_auxd("timepixr", wcs->timepixr);

  wcsprintf("     obsgeo: ");
  for (int k = 0; k < 3; k++) {
    if (undefined(wcs->obsgeo[k])) {
      wcsprintf("       UNDEFINED");
    } else {
      wcsprintf(" %15.6f", wcs->obsgeo[k]);
    }
  }
  wcsprintf("\n             ");
  for (int k = 3; k < 6; k++) {
    if (undefined(wcs->obsgeo[k])) {
      wcsprintf("       UNDEFINED");
    } else {
      wcsprintf(" %15.6f", wcs->obsgeo[k]);
    }
  }
  wcsprintf("\n");

  wcsprt_auxc("obsorbit", wcs->obsorbit);
  wcsprt_auxc(" radesys", wcs->radesys);
  wcsprt_auxd(" equinox", wcs->equinox);
  wcsprt_auxc(" specsys", wcs->specsys);
  wcsprt_auxc(" ssysobs", wcs->ssysobs);
  wcsprt_auxd(" velosys", wcs->velosys);
  wcsprt_auxd(" zsource", wcs->zsource);
  wcsprt_auxc(" ssyssrc", wcs->ssyssrc);
  wcsprt_auxd(" velangl", wcs->velangl);

  // Additional auxiliary coordinate system information.
  WCSPRINTF_PTR("        aux: ", wcs->aux, "\n");
  if (wcs->aux) {
    wcsprt_auxd("rsun_ref", wcs->aux->rsun_ref);
    wcsprt_auxd("dsun_obs", wcs->aux->dsun_obs);
    wcsprt_auxd("crln_obs", wcs->aux->crln_obs);
    wcsprt_auxd("hgln_obs", wcs->aux->hgln_obs);
    wcsprt_auxd("hglt_obs", wcs->aux->hglt_obs);
  }

  wcsprintf("       ntab: %d\n", wcs->ntab);
  WCSPRINTF_PTR("        tab: ", wcs->tab, "");
  if (wcs->tab != 0x0) wcsprintf("  (see below)");
  wcsprintf("\n");
  wcsprintf("       nwtb: %d\n", wcs->nwtb);
  WCSPRINTF_PTR("        wtb: ", wcs->wtb, "");
  if (wcs->wtb != 0x0) wcsprintf("  (see below)");
  wcsprintf("\n");

  // Derived values.
  WCSPRINTF_PTR("      types: ", wcs->types, "\n           ");
  for (int i = 0; i < wcs->naxis; i++) {
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

  // Memory management.
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
  WCSPRINTF_PTR("    m_czphs: ", wcs->m_czphs, "");
  if (wcs->m_czphs == wcs->czphs) wcsprintf("  (= czphs)");
  wcsprintf("\n");
  WCSPRINTF_PTR("    m_cperi: ", wcs->m_cperi, "");
  if (wcs->m_cperi == wcs->cperi) wcsprintf("  (= cperi)");
  wcsprintf("\n");
  WCSPRINTF_PTR("      m_aux: ", wcs->m_aux, "");
  if (wcs->m_aux == wcs->aux) wcsprintf("  (= aux)");
  wcsprintf("\n");
  WCSPRINTF_PTR("      m_tab: ", wcs->m_tab, "");
  if (wcs->m_tab == wcs->tab) wcsprintf("  (= tab)");
  wcsprintf("\n");
  WCSPRINTF_PTR("      m_wtb: ", wcs->m_wtb, "");
  if (wcs->m_wtb == wcs->wtb) wcsprintf("  (= wtb)");
  wcsprintf("\n");

  // Tabular transformation parameters.
  struct wtbarr *wtbp = wcs->wtb;
  if (wtbp) {
    for (int iwtb = 0; iwtb < wcs->nwtb; iwtb++, wtbp++) {
      wcsprintf("\n");
      wcsprintf("wtb[%d].*\n", iwtb);
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
    for (int itab = 0; itab < wcs->ntab; itab++) {
      wcsprintf("\n");
      wcsprintf("tab[%d].*\n", itab);
      tabprt(wcs->tab + itab);
    }
  }

  // Linear transformation parameters.
  wcsprintf("\n");
  wcsprintf("   lin.*\n");
  linprt(&(wcs->lin));

  // Celestial transformation parameters.
  wcsprintf("\n");
  wcsprintf("   cel.*\n");
  celprt(&(wcs->cel));

  // Spectral transformation parameters.
  wcsprintf("\n");
  wcsprintf("   spc.*\n");
  spcprt(&(wcs->spc));

  return WCSERR_SUCCESS;
}

//----------------------------------------------------------------------------

int wcsperr(const struct wcsprm *wcs, const char *prefix)

{
  if (wcs == 0x0) return WCSERR_NULL_POINTER;

  if (wcs->err && wcserr_prt(wcs->err, prefix) == 0) {
    linperr(&(wcs->lin), prefix);
    celperr(&(wcs->cel), prefix);
    wcserr_prt(wcs->spc.err, prefix);
    if (wcs->tab) {
      for (int itab = 0; itab < wcs->ntab; itab++) {
        wcserr_prt((wcs->tab + itab)->err, prefix);
      }
    }
  }

  return WCSERR_SUCCESS;
}

//----------------------------------------------------------------------------

int wcsbchk(struct wcsprm *wcs, int bounds)

{
  if (wcs == 0x0) return WCSERR_NULL_POINTER;

  if (wcs->flag != WCSSET) {
    int status;
    if ((status = wcsset(wcs))) return status;
  }

  wcs->cel.prj.bounds = bounds;

  return WCSERR_SUCCESS;
}

//----------------------------------------------------------------------------

int wcsset(struct wcsprm *wcs)

{
  static const char *function = "wcsset";

  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  struct wcserr **err = &(wcs->err);

  // Determine axis types from CTYPEia.
  int status;
  if ((status = wcs_types(wcs))) {
    return status;
  }

  // Convert to canonical units.
  if ((status = wcs_units(wcs))) {
    return status;
  }

  int naxis = wcs->naxis;
  if (32 < naxis) {
    return wcserr_set(WCSERR_SET(WCSERR_BAD_PARAM),
      "naxis must not exceed 32 (got %d)", naxis);
  }


  // Non-linear celestial axes present?
  if (wcs->lng >= 0 && wcs->types[wcs->lng] == 2200) {
    struct celprm *wcscel = &(wcs->cel);
    celini(wcscel);

    // CRVALia, LONPOLEa, and LATPOLEa keyvalues.
    wcscel->ref[0] = wcs->crval[wcs->lng];
    wcscel->ref[1] = wcs->crval[wcs->lat];
    wcscel->ref[2] = wcs->lonpole;
    wcscel->ref[3] = wcs->latpole;

    // Do alias translation for TPU/TPV before dealing with PVi_ma.
    struct prjprm *wcsprj = &(wcscel->prj);
    strncpy(wcsprj->code, wcs->ctype[wcs->lng]+5, 3);
    wcsprj->code[3] = '\0';
    if (strncmp(wcsprj->code, "TPU", 3) == 0 ||
        strncmp(wcsprj->code, "TPV", 3) == 0) {
      // Translate the PV parameters.
      struct disprm *dis;
      if ((dis = calloc(1, sizeof(struct disprm))) == 0x0) {
        return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
      }

      int ndpmax = 6 + wcs->npv;

      // Attach it to linprm.  Also inits it.
      char dpq[16];
      struct linprm *wcslin = &(wcs->lin);
      dis->flag = -1;
      if (strncmp(wcsprj->code, "TPU", 3) == 0) {
        // Prior distortion.
        lindist(1, wcslin, dis, ndpmax);
        strcpy(dpq, "DP");
      } else {
        // Sequent distortion.
        lindist(2, wcslin, dis, ndpmax);
        strcpy(dpq, "DQ");
      }

      // Yes, the distortion type is "TPV" even for TPU.
      strcpy(dis->dtype[wcs->lng], "TPV");
      strcpy(dis->dtype[wcs->lat], "TPV");

      // Keep the keywords in axis-order to aid debugging.
      struct dpkey *keyp = dis->dp;
      dis->ndp = 0;

      sprintf(dpq+2, "%d", wcs->lng+1);
      dpfill(keyp++, dpq, "NAXES",  0, 0, 2, 0.0);
      dpfill(keyp++, dpq, "AXIS.1", 0, 0, 1, 0.0);
      dpfill(keyp++, dpq, "AXIS.2", 0, 0, 2, 0.0);
      dis->ndp += 3;

      // Copy distortion parameters for the longitude axis.
      for (int k = 0; k < wcs->npv; k++) {
        if (wcs->pv[k].i != wcs->lng+1) continue;
        sprintf(keyp->field, "%s.TPV.%d", dpq, wcs->pv[k].m);
        dpfill(keyp++, 0x0, 0x0, 0, 1, 0, wcs->pv[k].value);
        dis->ndp++;
      }

      // Now the latitude axis.
      sprintf(dpq+2, "%d", wcs->lat+1);
      dpfill(keyp++, dpq, "NAXES",  0, 0, 2, 0.0);
      dpfill(keyp++, dpq, "AXIS.1", 0, 0, 2, 0.0);
      dpfill(keyp++, dpq, "AXIS.2", 0, 0, 1, 0.0);
      dis->ndp += 3;

      for (int k = 0; k < wcs->npv; k++) {
        if (wcs->pv[k].i != wcs->lat+1) continue;
        sprintf(keyp->field, "%s.TPV.%d", dpq, wcs->pv[k].m);
        dpfill(keyp++, 0x0, 0x0, 0, 1, 0, wcs->pv[k].value);
        dis->ndp++;
      }

      // Erase PVi_ma associated with the celestial axes.
      int n = 0;
      for (int k = 0; k < wcs->npv; k++) {
        int i = wcs->pv[k].i - 1;
        if (i == wcs->lng || i == wcs->lat) continue;

        wcs->pv[n].i = wcs->pv[k].i;
        wcs->pv[n].m = wcs->pv[k].m;
        wcs->pv[n].value = wcs->pv[k].value;

        n++;
      }

      wcs->npv = n;
      strcpy(wcsprj->code, "TAN");

      // As the PVi_ma have now been erased, ctype must be reset to prevent
      // this translation from re-occurring if wcsset() is called again.
      strcpy(wcs->ctype[wcs->lng]+5, "TAN");
      strcpy(wcs->ctype[wcs->lat]+5, "TAN");

    } else if (strncmp(wcsprj->code, "TNX", 3) == 0) {
      // The WAT distortion should already have been encoded in disseq.
      strcpy(wcsprj->code, "TAN");
      strcpy(wcs->ctype[wcs->lng]+5, "TAN");
      strcpy(wcs->ctype[wcs->lat]+5, "TAN");

    } else if (strncmp(wcsprj->code, "ZPX", 3) == 0) {
      // The WAT distortion should already have been encoded in disseq.
      strcpy(wcsprj->code, "ZPN");
      strcpy(wcs->ctype[wcs->lng]+5, "ZPN");
      strcpy(wcs->ctype[wcs->lat]+5, "ZPN");
    }

    // PVi_ma keyvalues.
    for (int k = 0; k < wcs->npv; k++) {
      if (wcs->pv[k].i == 0) {
        // From a PROJPn keyword.
        wcs->pv[k].i = wcs->lat + 1;
      }

      int i = wcs->pv[k].i - 1;
      int m = wcs->pv[k].m;

      if (i == wcs->lat) {
        // PVi_ma associated with latitude axis.
        if (m < 30) {
          wcsprj->pv[m] = wcs->pv[k].value;
        }

      } else if (i == wcs->lng) {
        // PVi_ma associated with longitude axis.
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
          // If present, overrides LONPOLEa.
          wcscel->ref[2] = wcs->pv[k].value;
          break;
        case 4:
          // If present, overrides LATPOLEa.
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

    // Do simple alias translations.
    if (strncmp(wcs->ctype[wcs->lng]+5, "GLS", 3) == 0) {
      wcscel->offset = 1;
      wcscel->phi0   = 0.0;
      wcscel->theta0 = wcs->crval[wcs->lat];
      strcpy(wcsprj->code, "SFL");

    } else if (strncmp(wcs->ctype[wcs->lng]+5, "NCP", 3) == 0) {
      // Convert NCP to SIN.
      if (wcscel->ref[1] == 0.0) {
        return wcserr_set(WCSERR_SET(WCSERR_BAD_PARAM),
          "Invalid projection: NCP blows up on the equator");
      }

      strcpy(wcsprj->code, "SIN");
      wcsprj->pv[1] = 0.0;
      wcsprj->pv[2] = cosd(wcscel->ref[1])/sind(wcscel->ref[1]);
    }

    // Initialize the celestial transformation routines.
    wcsprj->r0 = 0.0;
    if ((status = celset(wcscel))) {
      return wcserr_set(WCS_ERRMSG(wcs_celerr[status]));
    }

    // Update LONPOLE, LATPOLE, and PVi_ma keyvalues.
    wcs->lonpole = wcscel->ref[2];
    wcs->latpole = wcscel->ref[3];

    for (int k = 0; k < wcs->npv; k++) {
      int i = wcs->pv[k].i - 1;
      int m = wcs->pv[k].m;

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


  // Non-linear spectral axis present?
  if (wcs->spec >= 0 && wcs->types[wcs->spec] == 3300) {
    char scode[4], stype[5];
    struct spcprm *wcsspc = &(wcs->spc);
    spcini(wcsspc);
    if ((status = spctype(wcs->ctype[wcs->spec], stype, scode, 0x0, 0x0, 0x0,
                          0x0, 0x0, err))) {
      return status;
    }
    strcpy(wcsspc->type, stype);
    strcpy(wcsspc->code, scode);

    // CRVALia, RESTFRQa, and RESTWAVa keyvalues.
    wcsspc->crval = wcs->crval[wcs->spec];
    wcsspc->restfrq = wcs->restfrq;
    wcsspc->restwav = wcs->restwav;

    // PVi_ma keyvalues.
    for (int k = 0; k < wcs->npv; k++) {
      int i = wcs->pv[k].i - 1;
      int m = wcs->pv[k].m;

      if (i == wcs->spec) {
        // PVi_ma associated with grism axis.
        if (m < 7) {
          wcsspc->pv[m] = wcs->pv[k].value;
        }
      }
    }

    // Initialize the spectral transformation routines.
    if ((status = spcset(wcsspc))) {
      return wcserr_set(WCS_ERRMSG(wcs_spcerr[status]));
    }
  }


  // Tabular axes present?
  for (int itab = 0; itab < wcs->ntab; itab++) {
    if ((status = tabset(wcs->tab + itab))) {
      return wcserr_set(WCS_ERRMSG(wcs_taberr[status]));
    }
  }


  // Initialize the linear transformation.
  wcs->altlin &= 15;
  if (wcs->altlin > 1 && !(wcs->altlin & 1)) {
    double *pc = wcs->pc;

    if ((wcs->altlin & 2) && !(wcs->altlin & 8)) {
      // Copy CDi_ja to PCi_ja and reset CDELTia.
      double *cd = wcs->cd;
      for (int i = 0; i < naxis; i++) {
        for (int j = 0; j < naxis; j++) {
          *(pc++) = *(cd++);
        }
        wcs->cdelt[i] = 1.0;
      }

    } else if (wcs->altlin & 4) {
      // Construct PCi_ja from CROTAia.
      int i, j;
      if ((i = wcs->lng) >= 0 && (j = wcs->lat) >= 0) {
        double rho = wcs->crota[j];

        if (wcs->cdelt[i] == 0.0) {
          return wcserr_set(WCSERR_SET(WCSERR_SINGULAR_MTX),
            "Singular transformation matrix, CDELT%d is zero", i+1);
        }
        double lambda = wcs->cdelt[j]/wcs->cdelt[i];

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


  // Set defaults for radesys and equinox for equatorial or ecliptic.
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
      // Equinox is not applicable for these coordinate systems.
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
    // No celestial axes, ensure that radesys and equinox are unset.
    memset(wcs->radesys, 0, 72);
    wcs->equinox = UNDEFINED;
  }


  // Strip off trailing blanks and null-fill auxiliary string members.
  if (wcs->alt[0] == '\0') wcs->alt[0] = ' ';
  memset(wcs->alt+1, '\0', 3);

  for (int i = 0; i < naxis; i++) {
    wcsutil_null_fill(72, wcs->cname[i]);
  }
  wcsutil_null_fill(72, wcs->wcsname);
  wcsutil_null_fill(72, wcs->timesys);
  wcsutil_null_fill(72, wcs->trefpos);
  wcsutil_null_fill(72, wcs->trefdir);
  wcsutil_null_fill(72, wcs->plephem);
  wcsutil_null_fill(72, wcs->timeunit);
  wcsutil_null_fill(72, wcs->dateref);
  wcsutil_null_fill(72, wcs->dateobs);
  wcsutil_null_fill(72, wcs->datebeg);
  wcsutil_null_fill(72, wcs->dateavg);
  wcsutil_null_fill(72, wcs->dateend);
  wcsutil_null_fill(72, wcs->obsorbit);
  wcsutil_null_fill(72, wcs->radesys);
  wcsutil_null_fill(72, wcs->specsys);
  wcsutil_null_fill(72, wcs->ssysobs);
  wcsutil_null_fill(72, wcs->ssyssrc);

  // MJDREF defaults to zero if no reference date keywords were defined.
  if (wcs->dateref[0] == '\0') {
    if (undefined(wcs->mjdref[0])) {
      wcs->mjdref[0] = 0.0;
    }
    if (undefined(wcs->mjdref[1])) {
      wcs->mjdref[1] = 0.0;
    }
  }

  wcs->flag = WCSSET;

  return WCSERR_SUCCESS;
}

// : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :

int wcs_types(struct wcsprm *wcs)

{
  static const char *function = "wcs_types";

  const int  nalias = 6;
  const char aliases [6][4] = {"NCP", "GLS", "TPU", "TPV", "TNX", "ZPX"};

  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  struct wcserr **err = &(wcs->err);

  // Parse the CTYPEia keyvalues.
  char pcode[4], requir[16];
  pcode[0]  = '\0';
  requir[0] = '\0';
  wcs->lng  = -1;
  wcs->lat  = -1;
  wcs->spec = -1;
  wcs->cubeface = -1;

  const char *alt = "";
  if (*(wcs->alt) != ' ') alt = wcs->alt;


  int naxis = wcs->naxis;
  if (wcs->types) free(wcs->types);
  if ((wcs->types = calloc(naxis, sizeof(int))) == 0x0) {
    return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
  }

  int *ndx = 0x0;
  for (int i = 0; i < naxis; i++) {
    // Null fill.
    wcsutil_null_fill(72, wcs->ctype[i]);

    char ctypei[16];
    strncpy(ctypei, wcs->ctype[i], 15);
    ctypei[15] = '\0';

    // Check for early Paper IV syntax (e.g. '-SIP' used by Spitzer).
    if (strlen(ctypei) == 12 && ctypei[8] == '-') {
      // Excise the "4-3-3" or "8-3"-form distortion code.
      ctypei[8] = '\0';

      // Remove trailing dashes from "8-3"-form codes.
      for (int j = 7; j > 0; j--) {
        if (ctypei[j] != '-') break;
        ctypei[j] = '\0';
      }
    }

    // Logarithmic or tabular axis?
    wcs->types[i] = 0;
    if (strcmp(ctypei+4, "-LOG") == 0) {
      // Logarithmic axis.
      wcs->types[i] = 400;

    } else if (strcmp(ctypei+4, "-TAB") == 0) {
      // Tabular axis.
      wcs->types[i] = 500;
    }

    if (wcs->types[i]) {
      // Could have -LOG or -TAB with celestial or spectral types.
      ctypei[4] = '\0';

      // Take care of things like 'FREQ-LOG' or 'RA---TAB'.
      for (int j = 3; j >= 0; j--) {
        if (ctypei[j] != '-') break;
        ctypei[j] = '\0';
      }
    }

    // Translate AIPS spectral types for spctyp().
    char specsys[9];
    if (spcaips(ctypei, wcs->velref, ctypei, specsys) == 0) {
      strcpy(wcs->ctype[i], ctypei);
      if (wcs->specsys[0] == '\0') strcpy(wcs->specsys, specsys);
    }

    // Process linear axes.
    if (!(strlen(ctypei) == 8 && ctypei[4] == '-')) {
      // Identify Stokes, celestial, spectral, and time types.
      if (strcmp(ctypei, "STOKES") == 0) {
        // STOKES axis.
        wcs->types[i] = 1100;

      } else if (strcmp(ctypei, "RA")  == 0 ||
        strcmp(ctypei+1, "LON") == 0 ||
        strcmp(ctypei+2, "LN")  == 0) {
        // Longitude axis.
        wcs->types[i] += 2000;
        if (wcs->lng < 0) {
          wcs->lng = i;
          strcpy(wcs->lngtyp, ctypei);
        }

      } else if (strcmp(ctypei,   "DEC") == 0 ||
                 strcmp(ctypei+1, "LAT") == 0 ||
                 strcmp(ctypei+2, "LT")  == 0) {
        // Latitude axis.
        wcs->types[i] += 2001;
        if (wcs->lat < 0) {
          wcs->lat = i;
          strcpy(wcs->lattyp, ctypei);
        }

      } else if (strcmp(ctypei, "CUBEFACE") == 0) {
        // CUBEFACE axis.
        if (wcs->cubeface == -1) {
          wcs->types[i] = 2102;
          wcs->cubeface = i;
        } else {
          // Multiple CUBEFACE axes!
          return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
            "Multiple CUBEFACE axes (in CTYPE%d%.1s and CTYPE%d%.1s)",
            wcs->cubeface+1, alt, i+1, alt);
        }

      } else if (spctyp(ctypei, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0) == 0) {
        // Spectral axis.
        if (wcs->spec < 0) wcs->spec = i;
        wcs->types[i] += 3000;

      } else if (time_type(ctypei)) {
        // Time axis.
        wcs->types[i] += 4000;
      }

      continue;
    }


    // CTYPEia is in "4-3" form; is it a recognized spectral type?
    char scode[4];
    if (spctyp(ctypei, 0x0, scode, 0x0, 0x0, 0x0, 0x0, 0x0) == 0) {
      // Non-linear spectral axis found.
      wcs->types[i] = 3300;

      // Check uniqueness.
      if (wcs->spec >= 0) {
        return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
          "Multiple spectral axes (in CTYPE%d%.1s and CTYPE%d%.1s)",
          wcs->spec+1, alt, i+1, alt);
      }

      wcs->spec = i;

      continue;
    }


    // Is it a recognized celestial projection?
    int j;
    for (j = 0; j < prj_ncode; j++) {
      if (strncmp(ctypei+5, prj_codes[j], 3) == 0) break;
    }

    if (j == prj_ncode) {
      // Not a standard projection code, maybe it's an alias.
      for (j = 0; j < nalias; j++) {
        if (strncmp(ctypei+5, aliases[j], 3) == 0) break;
      }

      if (j == nalias) {
        // Not a recognized algorithm code of any type.
        wcs->types[i] = -1;
        return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
          "Unrecognized projection code (%s in CTYPE%d%.1s)",
          ctypei+5, i+1, alt);
      }
    }

    // Parse the celestial axis type.
    wcs->types[i] = 2200;
    if (*pcode == '\0') {
      // The first of the two celestial axes.
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
        // Unrecognized celestial type.
        wcs->types[i] = -1;

        wcs->lng = -1;
        wcs->lat = -1;
        return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
          "Unrecognized celestial type (%5s in CTYPE%d%.1s)",
          ctypei, i+1, alt);
      }

      if (wcs->lat >= 0) wcs->types[i]++;

    } else {
      // Looking for the complementary celestial axis.
      if (wcs->lat < 0) wcs->types[i]++;

      if (strncmp(ctypei, requir, 8) != 0) {
        // Inconsistent projection types.
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

  // Do we have a complementary pair of celestial axes?
  if (strcmp(requir, "")) {
    // Unmatched celestial axis.
    wcs->lng = -1;
    wcs->lat = -1;
    return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
      "Unmatched celestial axes");
  }

  // Table group numbers.
  for (int j = 0; j < wcs->ntab; j++) {
    for (int m = 0; m < wcs->tab[j].M; m++) {
      // Get image axis number.
      int i = wcs->tab[j].map[m];

      int type = (wcs->types[i] / 100) % 10;
      if (type != 5) {
        return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
          "Table parameters set for non-table axis type");
      }
      wcs->types[i] += 10 * j;
    }
  }

  return WCSERR_SUCCESS;
}

// : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :

int time_type(const char *ctype)

{
  // Is it a recognised time system as listed in Table 2 of WCS Paper VII?
  if (strncmp(ctype, "TIME", 4) == 0) return time_code(ctype, 4);
  if (strncmp(ctype, "UTC",  3) == 0) return time_code(ctype, 3);
  if (strncmp(ctype, "TAI",  3) == 0) return time_code(ctype, 3);
  if (strncmp(ctype, "IAT",  3) == 0) return time_code(ctype, 3);
  if (strncmp(ctype, "TT",   2) == 0) return time_code(ctype, 2);
  if (strncmp(ctype, "TDB",  3) == 0) return time_code(ctype, 3);
  if (strncmp(ctype, "TDT",  3) == 0) return time_code(ctype, 3);
  if (strncmp(ctype, "GPS",  3) == 0) return time_code(ctype, 3);
  if (strncmp(ctype, "TCB",  3) == 0) return time_code(ctype, 3);
  if (strncmp(ctype, "TCG",  3) == 0) return time_code(ctype, 3);
  if (strncmp(ctype, "GMT",  3) == 0) return time_code(ctype, 3);
  if (strncmp(ctype, "UT1",  3) == 0) return time_code(ctype, 3);
  if (strncmp(ctype, "UT",   2) == 0) return time_code(ctype, 2);
  if (strncmp(ctype, "ET",   2) == 0) return time_code(ctype, 2);
  if (strncmp(ctype, "LOCAL",5) == 0) return 1;

  return 0;
}

// : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :

static int time_code(const char *ctype, int nc)

{
  // If no algorithm code then we're finished.
  if (*(ctype+nc) == '\0') return 1;

  // Check the correct number of hyphens for things like "TT---TAB".
  while (nc < 4) {
    if (*(ctype+nc) != '-') return 0;
    nc++;
  }

  // Is it a code applicable to time-like axes?
  const char *code = ctype + 4;
  if (*code == '-') {
    if (strncmp(code, "-LOG", 5) == 0) return 1;
    if (strncmp(code, "-TAB", 5) == 0) return 1;
  }

  return 0;
}

// : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :

int wcs_units(struct wcsprm *wcs)

{
  static const char *function = "wcs_units";

  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  struct wcserr **err = &(wcs->err);

  int naxis = wcs->naxis;
  for (int i = 0; i < naxis; i++) {
    // Squeeze out trailing blanks.
    wcsutil_null_fill(72, wcs->cunit[i]);

    // Use types set by wcs_types().
    char ctype[9], units[16];
    switch (wcs->types[i]/1000) {
    case 2:
      // Celestial axis.
      strcpy(units, "deg");
      break;

    case 3:
      // Spectral axis.
      strncpy(ctype, wcs->ctype[i], 8);
      ctype[8] = '\0';
      spctyp(ctype, 0x0, 0x0, 0x0, units, 0x0, 0x0, 0x0);
      break;

    default:
      continue;
    }

    // Tabular axis, CDELTia and CRVALia relate to indices.
    if ((wcs->types[i]/100)%10 == 5) {
      continue;
    }

    if (wcs->cunit[i][0]) {
      double scale, offset, power;
      struct wcserr *uniterr;
      if (wcsunitse(wcs->cunit[i], units, &scale, &offset, &power,
                    &uniterr)) {
        if (uniterr) {
          // uniterr will not be set if wcserr is not enabled.
          wcserr_set(WCSERR_SET(WCSERR_BAD_COORD_TRANS),
            "In CUNIT%d%.1s: %s", i+1, (*wcs->alt)?wcs->alt:"", uniterr->msg);
          free(uniterr);
        }
        return WCSERR_BAD_COORD_TRANS;
      }

      if (scale != 1.0) {
        wcs->cdelt[i] *= scale;
        wcs->crval[i] *= scale;

        for (int j = 0; j < naxis; j++) {
          *(wcs->cd + i*naxis + j) *= scale;
        }

        strcpy(wcs->cunit[i], units);
      }

    } else {
      strcpy(wcs->cunit[i], units);
    }
  }

  return WCSERR_SUCCESS;
}

//----------------------------------------------------------------------------

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

  // Initialize if required.
  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  struct wcserr **err = &(wcs->err);

  int status = 0;
  if (wcs->flag != WCSSET) {
    if ((status = wcsset(wcs))) return status;
  }

  // Sanity check.
  if (ncoord < 1 || (ncoord > 1 && nelem < wcs->naxis)) {
    return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
      "ncoord and/or nelem inconsistent with the wcsprm");
  }


  // Initialize status vectors.
  int *istatp;
  if ((istatp = calloc(ncoord, sizeof(int))) == 0x0) {
    return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
  }

  stat[0] = 0;
  wcsutil_setAli(ncoord, 1, stat);


  // Apply pixel-to-world linear transformation.
  struct linprm *lin = &(wcs->lin);
  if (!(lin->dispre || lin->disseq)) {
    // No distortions present, do vector call.
    int istat = linp2x(lin, ncoord, nelem, pixcrd, imgcrd);
    if (istat) {
      // If one fails then all fail.
      status = wcserr_set(WCS_ERRMSG(wcs_linerr[istat]));
      goto cleanup;
    }

  } else {
    // Distortions present, get the status return for each coordinate.
    int disaxes = 0;

    register const double *pix = pixcrd;
    register double *img = imgcrd;
    register int *statp = stat;
    for (int k = 0 ; k < ncoord; k++, pix += nelem, img += nelem, statp++) {
      int istat = linp2x(lin, 1, nelem, pix, img);
      if (istat) {
        status = wcserr_set(WCS_ERRMSG(wcs_linerr[istat]));
        if (status != WCSERR_BAD_PIX) {
          goto cleanup;
        }

        if (disaxes == 0) {
          // Which axes have distortions?
          struct disprm *dispre = lin->dispre;
          struct disprm *disseq = lin->disseq;
          for (int i = 0; i < wcs->naxis; i++) {
            if (dispre && dispre->disp2x[i]) {
              disaxes |= (1 << i);
            } else if (disseq && disseq->disp2x[i]) {
              disaxes |= (1 << i);
            }
          }

          if (disaxes == 0) {
            // Shouldn't happen.
            disaxes = (2 << wcs->naxis) - 1;
          }
        }

        // WCSERR_BAD_PIX stat[] vector accounting.
        *statp = disaxes;
      }
    }
  }


  // Convert intermediate world coordinates to world coordinates.
  struct celprm *wcscel = &(wcs->cel);
  struct prjprm *wcsprj = &(wcscel->prj);
  for (int i = 0; i < wcs->naxis; i++) {
    // Extract the second digit of the axis type code.
    int type = (wcs->types[i] / 100) % 10;

    register double *img, *wrl;
    if (type <= 1) {
      // Linear or quantized coordinate axis.
      img = imgcrd + i;
      wrl = world  + i;
      double crvali = wcs->crval[i];
      for (int k = 0; k < ncoord; k++) {
        *wrl = *img + crvali;
        img += nelem;
        wrl += nelem;
      }

    } else if (wcs->types[i] == 2200) {
      // Convert celestial coordinates; do we have a CUBEFACE axis?
      if (wcs->cubeface != -1) {
        // Separation between faces.
        double offset;
        if (wcsprj->r0 == 0.0) {
          offset = 90.0;
        } else {
          offset = wcsprj->r0*PI/2.0;
        }

        // Lay out faces in a plane.
        img = imgcrd;
        int *statp = stat;
        int bits = (1 << i) | (1 << wcs->lat);
        for (int k = 0; k < ncoord; k++, statp++) {
          int face = (int)(*(img+wcs->cubeface) + 0.5);
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

      // Check for constant x and/or y.
      int iso_x = 0;
      int iso_y = 0;
      int nx = ncoord;
      int ny = 0;

      if (ncoord > 1) {
        if ((iso_x = wcsutil_allEq(ncoord, nelem, imgcrd+i))) {
          nx = 1;
          ny = ncoord;
        }
        if ((iso_y = wcsutil_allEq(ncoord, nelem, imgcrd+wcs->lat))) {
          ny = 1;
        }
      }

      // Transform projection plane coordinates to celestial coordinates.
      int istat = celx2s(wcscel, nx, ny, nelem, nelem, imgcrd+i,
                         imgcrd+wcs->lat, phi, theta, world+i,
                         world+wcs->lat, istatp);
      if (istat) {
        status = wcserr_set(WCS_ERRMSG(wcs_celerr[istat]));
        if (status != WCSERR_BAD_PIX) {
          goto cleanup;
        }
      }

      // If x and y were both constant, replicate values.
      if (iso_x && iso_y) {
        wcsutil_setAll(ncoord, nelem, world+i);
        wcsutil_setAll(ncoord, nelem, world+wcs->lat);
        wcsutil_setAll(ncoord, 1, phi);
        wcsutil_setAll(ncoord, 1, theta);
        wcsutil_setAli(ncoord, 1, istatp);
      }

      // WCSERR_BAD_PIX stat[] vector accounting.
      if (istat) {
        int bits = (1 << i) | (1 << wcs->lat);
        wcsutil_setBit(ncoord, istatp, bits, stat);
      }

    } else if (type == 3 || type == 4) {
      // Spectral and logarithmic coordinates; check for constant x.
      int iso_x;
      int nx = ncoord;

      if (ncoord > 1) {
        if ((iso_x = wcsutil_allEq(ncoord, nelem, imgcrd+i))) {
          nx = 1;
        }
      }

      int istat = 0;
      if (wcs->types[i] == 3300) {
        // Spectral coordinates.
        istat = spcx2s(&(wcs->spc), nx, nelem, nelem, imgcrd+i, world+i,
                       istatp);
        if (istat) {
          status = wcserr_set(WCS_ERRMSG(wcs_spcerr[istat]));
          if (status != WCSERR_BAD_PIX) {
            goto cleanup;
          }
        }
      } else if (type == 4) {
        // Logarithmic coordinates.
        istat = logx2s(wcs->crval[i], nx, nelem, nelem, imgcrd+i, world+i,
                       istatp);
        if (istat) {
          status = wcserr_set(WCS_ERRMSG(wcs_logerr[istat]));
          if (status != WCSERR_BAD_PIX) {
            goto cleanup;
          }
        }
      }

      // If x was constant, replicate values.
      if (iso_x) {
        wcsutil_setAll(ncoord, nelem, world+i);
        wcsutil_setAli(ncoord, 1, istatp);
      }

      // WCSERR_BAD_PIX stat[] vector accounting.
      if (istat) {
        wcsutil_setBit(ncoord, istatp, 1 << i, stat);
      }
    }
  }


  // Do tabular coordinates.
  for (int itab = 0; itab < wcs->ntab; itab++) {
    int istat = tabx2s(wcs->tab + itab, ncoord, nelem, imgcrd, world, istatp);

    if (istat) {
      status = wcserr_set(WCS_ERRMSG(wcs_taberr[istat]));
      if (status != WCSERR_BAD_PIX) {
        goto cleanup;
      }

      // WCSERR_BAD_PIX stat[] vector accounting.
      int bits = 0;
      for (int m = 0; m < wcs->tab[itab].M; m++) {
        bits |= 1 << wcs->tab[itab].map[m];
      }
      wcsutil_setBit(ncoord, istatp, bits, stat);
    }
  }


  // Zero the unused world coordinate elements.
  for (int i = wcs->naxis; i < nelem; i++) {
    world[i] = 0.0;
    wcsutil_setAll(ncoord, nelem, world+i);
  }

cleanup:
  free(istatp);
  return status;
}

//----------------------------------------------------------------------------

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

  // Initialize if required.
  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  struct wcserr **err = &(wcs->err);

  int status = 0;
  if (wcs->flag != WCSSET) {
    if ((status = wcsset(wcs))) return status;
  }

  // Sanity check.
  if (ncoord < 1 || (ncoord > 1 && nelem < wcs->naxis)) {
    return wcserr_set(WCSERR_SET(WCSERR_BAD_CTYPE),
      "ncoord and/or nelem inconsistent with the wcsprm");
  }

  // Initialize status vectors.
  int *istatp;
  if ((istatp = calloc(ncoord, sizeof(int))) == 0x0) {
    return wcserr_set(WCS_ERRMSG(WCSERR_MEMORY));
  }

  stat[0] = 0;
  wcsutil_setAli(ncoord, 1, stat);


  // Convert world coordinates to intermediate world coordinates.
  struct celprm *wcscel = &(wcs->cel);
  struct prjprm *wcsprj = &(wcscel->prj);
  for (int i = 0; i < wcs->naxis; i++) {
    // Extract the second digit of the axis type code.
    int type = (wcs->types[i] / 100) % 10;

    if (type <= 1) {
      // Linear or quantized coordinate axis.
      register const double *wrl = world  + i;
      register double *img = imgcrd + i;
      double crvali = wcs->crval[i];
      for (int k = 0; k < ncoord; k++) {
        *img = *wrl - crvali;
        wrl += nelem;
        img += nelem;
      }

    } else if (wcs->types[i] == 2200) {
      // Celestial coordinates; check for constant lng and/or lat.
      int isolng = 0;
      int isolat = 0;
      int nlng = ncoord;
      int nlat = 0;

      if (ncoord > 1) {
        if ((isolng = wcsutil_allEq(ncoord, nelem, world+i))) {
          nlng = 1;
          nlat = ncoord;
        }
        if ((isolat = wcsutil_allEq(ncoord, nelem, world+wcs->lat))) {
          nlat = 1;
        }
      }

      // Transform celestial coordinates to projection plane coordinates.
      int istat = cels2x(wcscel, nlng, nlat, nelem, nelem, world+i,
                         world+wcs->lat, phi, theta, imgcrd+i,
                         imgcrd+wcs->lat, istatp);
      if (istat) {
        status = wcserr_set(WCS_ERRMSG(wcs_celerr[istat]));
        if (status != WCSERR_BAD_WORLD) {
          goto cleanup;
        }
      }

      // If lng and lat were both constant, replicate values.
      if (isolng && isolat) {
        wcsutil_setAll(ncoord, nelem, imgcrd+i);
        wcsutil_setAll(ncoord, nelem, imgcrd+wcs->lat);
        wcsutil_setAll(ncoord, 1, phi);
        wcsutil_setAll(ncoord, 1, theta);
        wcsutil_setAli(ncoord, 1, istatp);
      }

      // WCSERR_BAD_WORLD stat[] vector accounting.
      if (istat) {
        int bits = (1 << i) | (1 << wcs->lat);
        wcsutil_setBit(ncoord, istatp, bits, stat);
      }

      // Do we have a CUBEFACE axis?
      if (wcs->cubeface != -1) {
        // Separation between faces.
        double offset;
        if (wcsprj->r0 == 0.0) {
          offset = 90.0;
        } else {
          offset = wcsprj->r0*PI/2.0;
        }

        // Stack faces in a cube.
        register double *img = imgcrd;
        for (int k = 0; k < ncoord; k++) {
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
      // Spectral and logarithmic coordinates; check for constancy.
      int isospec = 0;
      int nwrld = ncoord;

      if (ncoord > 1) {
        if ((isospec = wcsutil_allEq(ncoord, nelem, world+i))) {
          nwrld = 1;
        }
      }

      int istat = 0;
      if (wcs->types[i] == 3300) {
        // Spectral coordinates.
        istat = spcs2x(&(wcs->spc), nwrld, nelem, nelem, world+i,
                       imgcrd+i, istatp);
        if (istat) {
          status = wcserr_set(WCS_ERRMSG(wcs_spcerr[istat]));
          if (status != WCSERR_BAD_WORLD) {
            goto cleanup;
          }
        }
      } else if (type == 4) {
        // Logarithmic coordinates.
        istat = logs2x(wcs->crval[i], nwrld, nelem, nelem, world+i,
                       imgcrd+i, istatp);
        if (istat) {
          status = wcserr_set(WCS_ERRMSG(wcs_logerr[istat]));
          if (status != WCSERR_BAD_WORLD) {
            goto cleanup;
          }
        }
      }

      // If constant, replicate values.
      if (isospec) {
        wcsutil_setAll(ncoord, nelem, imgcrd+i);
        wcsutil_setAli(ncoord, 1, istatp);
      }

      // WCSERR_BAD_WORLD stat[] vector accounting.
      if (istat) {
        wcsutil_setBit(ncoord, istatp, 1 << i, stat);
      }
    }
  }


  // Do tabular coordinates.
  for (int itab = 0; itab < wcs->ntab; itab++) {
    int istat = tabs2x(wcs->tab + itab, ncoord, nelem, world, imgcrd, istatp);

    if (istat) {
      status = wcserr_set(WCS_ERRMSG(wcs_taberr[istat]));

      if (status != WCSERR_BAD_WORLD) {
        goto cleanup;
      }

      int bits = 0;
      for (int m = 0; m < wcs->tab[itab].M; m++) {
        bits |= 1 << wcs->tab[itab].map[m];
      }
      wcsutil_setBit(ncoord, istatp, bits, stat);
    }
  }


  // Zero the unused intermediate world coordinate elements.
  for (int i = wcs->naxis; i < nelem; i++) {
    imgcrd[i] = 0.0;
    wcsutil_setAll(ncoord, nelem, imgcrd+i);
  }


  // Apply world-to-pixel linear transformation.
  struct linprm *lin = &(wcs->lin);
  if (!(lin->dispre || lin->disseq)) {
    // No distortions present, do vector call.
    int istat = linx2p(lin, ncoord, nelem, imgcrd, pixcrd);
    if (istat) {
      status = wcserr_set(WCS_ERRMSG(wcs_linerr[istat]));
      goto cleanup;
    }

  } else {
    // Distortions present, get the status return for each coordinate.
    int disaxes = 0;

    register const double *img = imgcrd;
    register double *pix = pixcrd;
    register int *statp = stat;
    for (int k = 0 ; k < ncoord; k++, pix += nelem, img += nelem, statp++) {
      int istat = linx2p(lin, 1, nelem, img, pix);
      if (istat) {
        status = wcserr_set(WCS_ERRMSG(wcs_linerr[istat]));
        if (status != WCSERR_BAD_WORLD) {
          goto cleanup;
        }

        if (disaxes == 0) {
          // Which axes have distortions?
          struct disprm *dispre = lin->dispre;
          struct disprm *disseq = lin->disseq;
          for (int i = 0; i < wcs->naxis; i++) {
            if (dispre && dispre->disp2x[i]) {
              disaxes |= (1 << i);
            } else if (disseq && disseq->disp2x[i]) {
              disaxes |= (1 << i);
            }
          }

          if (disaxes == 0) {
            // Shouldn't happen.
            disaxes = (2 << wcs->naxis) - 1;
          }
        }

        // WCSERR_BAD_WORLD stat[] vector accounting.
        *statp = disaxes;
      }
    }
  }

cleanup:
  free(istatp);
  return status;
}

//----------------------------------------------------------------------------

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
  const double tol  = 1.0e-10;
  const double tol2 = 100.0*tol;

  // Initialize if required.
  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  struct wcserr **err = &(wcs->err);

  int status;
  if (wcs->flag != WCSSET) {
    if ((status = wcsset(wcs))) return status;
  }

  if (wcs->lng < 0 || wcs->lat < 0) {
    return wcserr_set(WCSERR_SET(WCSERR_BAD_SUBIMAGE),
      "Image does not have celestial axes");
  }

  double *worldlng = world + wcs->lng;
  double *worldlat = world + wcs->lat;


  // Check vspan.
  double span[2];
  if (vspan[0] <= vspan[1]) {
    span[0] = vspan[0];
    span[1] = vspan[1];
  } else {
    // Swap them.
    span[0] = vspan[1];
    span[1] = vspan[0];
  }

  // Check vstep.
  double step = fabs(vstep);
  if (step == 0.0) {
    step = (span[1] - span[0])/10.0;
    if (step > 1.0 || step == 0.0) step = 1.0;
  }

  // Check viter.
  int nstep = viter;
  if (nstep < 5) {
    nstep = 5;
  } else if (nstep > 10) {
    nstep = 10;
  }

  // Given pixel element.
  double pixmix = pixcrd[mixpix];

  // Iterate on the step size.
  for (int istep = 0; istep <= nstep; istep++) {
    if (istep) step /= 2.0;

    // Iterate on the sky coordinate between the specified range.
    if (mixcel == 1) {
      // Celestial longitude is given.

      // Check whether the solution interval is a crossing interval.
      double lat0 = span[0];
      *worldlat = lat0;
      int stat[1];
      if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                           stat))) {
        if (status == WCSERR_BAD_WORLD) {
          status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
        }
        return status;
      }
      double d0 = pixcrd[mixpix] - pixmix;

      double dabs = fabs(d0);
      if (dabs < tol) return WCSERR_SUCCESS;

      double lat1 = span[1];
      *worldlat = lat1;
      if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                           stat))) {
        if (status == WCSERR_BAD_WORLD) {
          status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
        }
        return status;
      }
      double d1 = pixcrd[mixpix] - pixmix;

      dabs = fabs(d1);
      if (dabs < tol) return WCSERR_SUCCESS;

      double lmin = lat1;
      double dmin = dabs;

      // Check for a crossing point.
      int crossed;
      double dx = 0.0;
      if (signbit(d0) != signbit(d1)) {
        crossed = 1;
        dx = d1;
      } else {
        crossed = 0;
        lat0 = span[1];
      }

      for (int retry = 0; retry < 4; retry++) {
        // Refine the solution interval.
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

          // Check for a solution.
          dabs = fabs(d0);
          if (dabs < tol) return WCSERR_SUCCESS;

          // Record the point of closest approach.
          if (dabs < dmin) {
            lmin = lat0;
            dmin = dabs;
          }

          // Check for a crossing point.
          if (signbit(d0) != signbit(d1)) {
            crossed = 2;
            dx = d0;
            break;
          }

          // Advance to the next subinterval.
          lat1 = lat0;
          d1 = d0;
        }

        if (crossed) {
          // A crossing point was found.
          for (int iter = 0; iter < niter; iter++) {
            // Use regula falsi division of the interval.
            double lambda = d0/(d0-d1);
            if (lambda < 0.1) {
              lambda = 0.1;
            } else if (lambda > 0.9) {
              lambda = 0.9;
            }

            double dlat = lat1 - lat0;
            double lat = lat0 + lambda*dlat;
            *worldlat = lat;
            if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                                 stat))) {
              if (status == WCSERR_BAD_WORLD) {
                status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
              }
              return status;
            }

            // Check for a solution.
            double d = pixcrd[mixpix] - pixmix;
            dabs = fabs(d);
            if (dabs < tol) return WCSERR_SUCCESS;

            if (dlat < tol) {
              // An artifact of numerical imprecision.
              if (dabs < tol2) return WCSERR_SUCCESS;

              // Must be a discontinuity.
              break;
            }

            // Record the point of closest approach.
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

          // No convergence, must have been a discontinuity.
          if (crossed == 1) lat0 = span[1];
          lat1 = lat0;
          d1 = dx;
          crossed = 0;

        } else {
          // No crossing point; look for a tangent point.
          if (lmin == span[0]) break;
          if (lmin == span[1]) break;

          double lat = lmin;
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

          double d  = dmin;

          *worldlat = lat1;
          if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                               stat))) {
            if (status == WCSERR_BAD_WORLD) {
              status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
            }
            return status;
          }
          d1 = fabs(pixcrd[mixpix] - pixmix);

          for (int iter = 0; iter < niter; iter++) {
            double lat0m = (lat0 + lat)/2.0;
            *worldlat = lat0m;
            if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                                 stat))) {
              if (status == WCSERR_BAD_WORLD) {
                status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
              }
              return status;
            }
            double d0m = fabs(pixcrd[mixpix] - pixmix);

            if (d0m < tol) return WCSERR_SUCCESS;

            double lat1m = (lat1 + lat)/2.0;
            *worldlat = lat1m;
            if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                                 stat))) {
              if (status == WCSERR_BAD_WORLD) {
                status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
              }
              return status;
            }
            double d1m = fabs(pixcrd[mixpix] - pixmix);

            if (d1m < tol) return WCSERR_SUCCESS;

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
      // Celestial latitude is given.

      // Check whether the solution interval is a crossing interval.
      double lng0 = span[0];
      *worldlng = lng0;
      int stat[1];
      if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                           stat))) {
        if (status == WCSERR_BAD_WORLD) {
          status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
        }
        return status;
      }
      double d0 = pixcrd[mixpix] - pixmix;

      double dabs = fabs(d0);
      if (dabs < tol) return WCSERR_SUCCESS;

      double lng1 = span[1];
      *worldlng = lng1;
      if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                           stat))) {
        if (status == WCSERR_BAD_WORLD) {
          status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
        }
        return status;
      }
      double d1 = pixcrd[mixpix] - pixmix;

      dabs = fabs(d1);
      if (dabs < tol) return WCSERR_SUCCESS;
      double lmin = lng1;
      double dmin = dabs;

      // Check for a crossing point.
      int crossed;
      double dx = 0.0;
      if (signbit(d0) != signbit(d1)) {
        crossed = 1;
        dx = d1;
      } else {
        crossed = 0;
        lng0 = span[1];
      }

      for (int retry = 0; retry < 4; retry++) {
        // Refine the solution interval.
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

          // Check for a solution.
          dabs = fabs(d0);
          if (dabs < tol) return WCSERR_SUCCESS;

          // Record the point of closest approach.
          if (dabs < dmin) {
            lmin = lng0;
            dmin = dabs;
          }

          // Check for a crossing point.
          if (signbit(d0) != signbit(d1)) {
            crossed = 2;
            dx = d0;
            break;
          }

          // Advance to the next subinterval.
          lng1 = lng0;
          d1 = d0;
        }

        if (crossed) {
          // A crossing point was found.
          for (int iter = 0; iter < niter; iter++) {
            // Use regula falsi division of the interval.
            double lambda = d0/(d0-d1);
            if (lambda < 0.1) {
              lambda = 0.1;
            } else if (lambda > 0.9) {
              lambda = 0.9;
            }

            double dlng = lng1 - lng0;
            double lng = lng0 + lambda*dlng;
            *worldlng = lng;
            if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                                 stat))) {
              if (status == WCSERR_BAD_WORLD) {
                status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
              }
              return status;
            }

            // Check for a solution.
            double d = pixcrd[mixpix] - pixmix;
            dabs = fabs(d);
            if (dabs < tol) return WCSERR_SUCCESS;

            if (dlng < tol) {
              // An artifact of numerical imprecision.
              if (dabs < tol2) return WCSERR_SUCCESS;

              // Must be a discontinuity.
              break;
            }

            // Record the point of closest approach.
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

          // No convergence, must have been a discontinuity.
          if (crossed == 1) lng0 = span[1];
          lng1 = lng0;
          d1 = dx;
          crossed = 0;

        } else {
          // No crossing point; look for a tangent point.
          if (lmin == span[0]) break;
          if (lmin == span[1]) break;

          double lng = lmin;
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

          double d  = dmin;

          *worldlng = lng1;
          if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                               stat))) {
            if (status == WCSERR_BAD_WORLD) {
              status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
            }
            return status;
          }
          d1 = fabs(pixcrd[mixpix] - pixmix);

          for (int iter = 0; iter < niter; iter++) {
            double lng0m = (lng0 + lng)/2.0;
            *worldlng = lng0m;
            if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                                 stat))) {
              if (status == WCSERR_BAD_WORLD) {
                status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
              }
              return status;
            }
            double d0m = fabs(pixcrd[mixpix] - pixmix);

            if (d0m < tol) return WCSERR_SUCCESS;

            double lng1m = (lng1 + lng)/2.0;
            *worldlng = lng1m;
            if ((status = wcss2p(wcs, 1, 0, world, phi, theta, imgcrd, pixcrd,
                                 stat))) {
              if (status == WCSERR_BAD_WORLD) {
                status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
              }
              return status;
            }
            double d1m = fabs(pixcrd[mixpix] - pixmix);

            if (d1m < tol) return WCSERR_SUCCESS;

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


  // Set cel0 to the unity transformation.
  struct wcsprm wcs0 = *wcs;
  wcs0.cel.euler[0] = -90.0;
  wcs0.cel.euler[1] =   0.0;
  wcs0.cel.euler[2] =  90.0;
  wcs0.cel.euler[3] =   1.0;
  wcs0.cel.euler[4] =   0.0;

  // No convergence, check for aberrant behaviour at a native pole.
  *theta = -90.0;
  for (int j = 1; j <= 2; j++) {
    // Could the celestial coordinate element map to a native pole?
    *phi = 0.0;
    *theta = -*theta;
    struct celprm *wcscel = &(wcs->cel);
    double lng, lat;
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

    // Is there a solution for the given pixel coordinate element?
    lng = *worldlng;
    lat = *worldlat;

    // Feed native coordinates to wcss2p() with cel0 set to unity.
    *worldlng = -180.0;
    *worldlat = *theta;
    int stat[1];
    if ((status = wcss2p(&wcs0, 1, 0, world, phi, theta, imgcrd, pixcrd,
                         stat))) {
      wcserr_clear(err);
      wcs->err = wcs0.err;
      if (status == WCSERR_BAD_WORLD) {
        status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
      }
      return status;
    }
    double d0 = pixcrd[mixpix] - pixmix;

    // Check for a solution.
    if (fabs(d0) < tol) {
      // Recall saved world coordinates.
      *worldlng = lng;
      *worldlat = lat;
      return WCSERR_SUCCESS;
    }

    // Search for a crossing interval.
    double phi0 = -180.0, phi1;
    double d1;
    for (int k = -179; k <= 180; k++) {
      phi1 = (double)k;
      *worldlng = phi1;
      if ((status = wcss2p(&wcs0, 1, 0, world, phi, theta, imgcrd, pixcrd,
                           stat))) {
        wcserr_clear(err);
        wcs->err = wcs0.err;
        if (status == WCSERR_BAD_WORLD) {
          status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
        }
        return status;
      }
      d1 = pixcrd[mixpix] - pixmix;

      // Check for a solution.
      double dabs = fabs(d1);
      if (dabs < tol) {
        // Recall saved world coordinates.
        *worldlng = lng;
        *worldlat = lat;
        return WCSERR_SUCCESS;
      }

      // Is it a crossing interval?
      if (signbit(d0) != signbit(d1)) break;

      phi0 = phi1;
      d0 = d1;
    }

    for (int iter = 1; iter <= niter; iter++) {
      // Use regula falsi division of the interval.
      double lambda = d0/(d0-d1);
      if (lambda < 0.1) {
        lambda = 0.1;
      } else if (lambda > 0.9) {
        lambda = 0.9;
      }

      double dphi = phi1 - phi0;
      *worldlng = phi0 + lambda*dphi;
      if ((status = wcss2p(&wcs0, 1, 0, world, phi, theta, imgcrd, pixcrd,
                           stat))) {
        wcserr_clear(err);
        wcs->err = wcs0.err;
        if (status == WCSERR_BAD_WORLD) {
          status = wcserr_set(WCS_ERRMSG(WCSERR_BAD_WORLD_COORD));
        }
        return status;
      }

      // Check for a solution.
      double d = pixcrd[mixpix] - pixmix;
      double dabs = fabs(d);
      if (dabs < tol || (dphi < tol && dabs < tol2)) {
        // Recall saved world coordinates.
        *worldlng = lng;
        *worldlat = lat;
        return WCSERR_SUCCESS;
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


  // No solution.
  return wcserr_set(WCS_ERRMSG(WCSERR_NO_SOLUTION));
}

//----------------------------------------------------------------------------

int wcsccs(
  struct wcsprm *wcs,
  double lng2P1,
  double lat2P1,
  double lng1P2,
  const char *clng,
  const char *clat,
  const char *radesys,
  double equinox,
  const char *alt)

{
  static const char *function = "wcsccs";

  int status;

  // Initialize if required.
  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  struct wcserr **err = &(wcs->err);

  if (wcs->flag != WCSSET) {
    if ((status = wcsset(wcs))) return status;
  }

  if (wcs->lng < 0 || wcs->lat < 0) {
    return wcserr_set(WCSERR_SET(WCSERR_BAD_SUBIMAGE),
      "Image does not have celestial axes");
  }

  // (lng1XX,lat1XX)  ...longitude and latitude of XX in the old system.
  // (lng2XX,lat2XX)  ...longitude and latitude of XX in the new system.
  // XX = NP          ...natuve pole,
  //      P1          ...pole of the old system,
  //      P2          ...pole of the new system,
  //      FP          ...fiducial point.

  // Set up the transformation from the old to the new system.
  double euler12[5];
  euler12[0] = lng2P1;
  euler12[1] = 90.0 - lat2P1;
  euler12[2] = lng1P2;
  euler12[3] = cosd(euler12[1]);
  euler12[4] = sind(euler12[1]);

  // Transform coordinates of the fiducial point (FP) to the new system.
  double lng1FP = wcs->crval[wcs->lng];
  double lat1FP = wcs->crval[wcs->lat];
  double lng2FP, lat2FP;
  (void)sphx2s(euler12, 1, 1, 1, 1, &lng1FP, &lat1FP, &lng2FP, &lat2FP);

  // Compute native coordinates of the new pole (noting lat1P2 == lat2P1).
  double phiP2, thetaP2;
  (void)sphs2x(wcs->cel.euler, 1, 1, 1, 1, &lng1P2, &lat2P1,
               &phiP2, &thetaP2);

  if (fabs(lat2FP) == 90.0 || fabs(thetaP2) == 90.0) {
    // If one of the poles of the new system is at the fiducial point, then
    // lng2FP is indeterminate, and if one of them is at the native pole, then
    // phiP2 is indeterminate.  We have to work harder to obtain these values.

    // Compute coordinates of the native pole (NP) in the old and new systems.
    double phiNP = 0.0, thetaNP = 90.0;
    double lng1NP, lat1NP;
    (void)sphx2s(wcs->cel.euler, 1, 1, 1, 1, &phiNP, &thetaNP,
                 &lng1NP, &lat1NP);

    double lng2NP, lat2NP;
    (void)sphx2s(euler12, 1, 1, 1, 1, &lng1NP, &lat1NP, &lng2NP, &lat2NP);

    // Native latitude and longitude of the fiducial point, (phi0,theta0).
    double phiFP   = wcs->cel.prj.phi0;
    double thetaFP = wcs->cel.prj.theta0;

    if (fabs(lat2NP) == 90.0) {
      // Following WCS Paper II equations (3) and (4), we are free to choose
      // phiP2 and set lng2NP accordingly.  So set phiP2 to its default value
      // for the projection.
      if (thetaFP < lat2FP) {
        phiP2 = 0.0;
      } else {
        phiP2 = 180.0;
      }

      // Compute coordinates in the old system of test point X.
      double phiX = 0.0, thetaX = 0.0;
      double lng1X, lat1X;
      (void)sphx2s(wcs->cel.euler, 1, 1, 1, 1, &phiX, &thetaX,
                   &lng1X, &lat1X);

      // Ensure that lng1X is not indeterminate.
      if (fabs(lat1X) == 90.0) {
        phiX = 90.0;
        (void)sphx2s(wcs->cel.euler, 1, 1, 1, 1, &phiX, &thetaX,
                     &lng1X, &lat1X);
      }

      // Compute coordinates in the new system of test point X.
      double lng2X, lat2X;
      (void)sphx2s(euler12, 1, 1, 1, 1, &lng1X, &lat1X, &lng2X, &lat2X);

      // Apply WCS Paper II equations (3) and (4).
      if (lat2NP == +90.0) {
        lng2NP = lng2X + (phiP2 - phiX) + 180.0;
      } else {
        lng2NP = lng2X - (phiP2 - phiX);
      }

    } else {
      // For (lng2NP + 90, 0), WCS Paper II equation (5) reduces to
      // phi = phiP2 - 90.
      double lng2X = lng2NP + 90.0;
      double lat2X = 0.0;
      double lng1X, lat1X;
      (void)sphs2x(euler12, 1, 1, 1, 1, &lng2X, &lat2X, &lng1X, &lat1X);

      double phiX, thetaX;
      (void)sphs2x(wcs->cel.euler, 1, 1, 1, 1, &lng1X, &lat1X,
                   &phiX, &thetaX);

      phiP2 = phiX + 90.0;
    }

    // Compute the longitude of the fiducial point in the new system.
    double eulerN2[5];
    eulerN2[0] = lng2NP;
    eulerN2[1] = 90.0 - lat2NP;
    eulerN2[2] = phiP2;
    eulerN2[3] = cosd(eulerN2[1]);
    eulerN2[4] = sind(eulerN2[1]);

    (void)sphx2s(eulerN2, 1, 1, 1, 1, &phiFP, &thetaFP, &lng2FP, &lat2FP);
  }

  // Update reference values in wcsprm.
  wcs->flag = 0;
  wcs->crval[wcs->lng] = lng2FP;
  wcs->crval[wcs->lat] = lat2FP;
  wcs->lonpole = phiP2;
  wcs->latpole = thetaP2;

  // Update wcsprm::ctype.
  if (clng) {
    strncpy(wcs->ctype[wcs->lng], clng, 4);
    for (int i = 0; i < 4; i++) {
      if (wcs->ctype[wcs->lng][i] == '\0') {
        wcs->ctype[wcs->lng][i] = '-';
      }
    }
  }

  if (clat) {
    strncpy(wcs->ctype[wcs->lat], clat, 4);
    for (int i = 0; i < 4; i++) {
      if (wcs->ctype[wcs->lat][i] == '\0') {
        wcs->ctype[wcs->lat][i] = '-';
      }
    }
  }

  // Update auxiliary values.
  if (strncmp(wcs->ctype[wcs->lng], "RA--", 4) == 0 &&
      strncmp(wcs->ctype[wcs->lat], "DEC-", 4) == 0) {
    // Transforming to equatorial coordinates.
    if (radesys) {
      strncpy(wcs->radesys, radesys, 71);
    }

    if (equinox != 0.0) {
      wcs->equinox = equinox;
    }
  } else {
    // Meaningless for other than equatorial coordinates.
    memset(wcs->radesys, 0, 72);
    wcs->equinox = UNDEFINED;
  }

  if (alt && *alt) {
    wcs->alt[0] = *alt;
  }

  // Reset the struct.
  if ((status = wcsset(wcs))) return status;

  return WCSERR_SUCCESS;
}

//----------------------------------------------------------------------------

int wcssptr(
  struct wcsprm *wcs,
  int  *i,
  char ctype[9])

{
  static const char *function = "wcssptr";

  int status;

  // Initialize if required.
  if (wcs == 0x0) return WCSERR_NULL_POINTER;
  struct wcserr **err = &(wcs->err);

  if (wcs->flag != WCSSET) {
    if ((status = wcsset(wcs))) return status;
  }

  int j;
  if ((j = *i) < 0) {
    if ((j = wcs->spec) < 0) {
      // Look for a linear spectral axis.
      for (j = 0; j < wcs->naxis; j++) {
        if (wcs->types[j]/100 == 30) {
          break;
        }
      }

      if (j >= wcs->naxis) {
        // No spectral axis.
        return wcserr_set(WCSERR_SET(WCSERR_BAD_SUBIMAGE),
          "No spectral axis found");
      }
    }

    *i = j;
  }

  // Translate the spectral axis.
  double cdelt, crval;
  if ((status = spctrne(wcs->ctype[j], wcs->crval[j], wcs->cdelt[j],
                        wcs->restfrq, wcs->restwav, ctype, &crval, &cdelt,
                        &(wcs->spc.err)))) {
    return wcserr_set(WCS_ERRMSG(wcs_spcerr[status]));
  }


  // Translate keyvalues.
  wcs->flag = 0;
  wcs->cdelt[j] = cdelt;
  wcs->crval[j] = crval;
  spctyp(ctype, 0x0, 0x0, 0x0, wcs->cunit[j], 0x0, 0x0, 0x0);
  strcpy(wcs->ctype[j], ctype);

  // This keeps things tidy if the spectral axis is linear.
  spcfree(&(wcs->spc));
  spcini(&(wcs->spc));

  // Reset the struct.
  if ((status = wcsset(wcs))) return status;

  return WCSERR_SUCCESS;
}

//----------------------------------------------------------------------------

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
