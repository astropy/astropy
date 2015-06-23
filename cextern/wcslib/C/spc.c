/*============================================================================

  WCSLIB 5.6 - an implementation of the FITS WCS standard.
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
  $Id: spc.c,v 5.6 2015/06/14 07:11:24 mcalabre Exp $
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
#include "spc.h"
#include "spx.h"

/* Spectral algorithm codes. */
#define F2S 100;		/* Axis linear in frequency.          */
#define W2S 200;		/* Axis linear in vacuum wavelengths. */
#define A2S 300;		/* Axis linear in air wavelengths.    */
#define V2S 400;		/* Axis linear in velocity.           */
#define GRI 500;		/* Grism in vacuum.                   */
#define GRA 600;		/* Grism in air.                      */

/* S-type spectral variables. */
#define FREQ  0;		/* Frequency-like.                    */
#define AFRQ  1;		/* Frequency-like.                    */
#define ENER  2;		/* Frequency-like.                    */
#define WAVN  3;		/* Frequency-like.                    */
#define VRAD  4;		/* Frequency-like.                    */
#define WAVE 10;		/* Vacuum wavelength-like.            */
#define VOPT 11;		/* Vacuum wavelength-like.            */
#define ZOPT 12;		/* Vacuum wavelength-like.            */
#define AWAV 20;		/* Air wavelength-like.               */
#define VELO 30;		/* Velocity-like.                     */
#define BETA 31;		/* Velocity-like.                     */


/* Map status return value to message. */
const char *spc_errmsg[] = {
  "Success",
  "Null spcprm pointer passed",
  "Invalid spectral parameters",
  "One or more of x coordinates were invalid",
  "One or more of the spec coordinates were invalid"};

/* Map error returns for lower-level routines.  SPXERR_BAD_INSPEC_COORD */
/* maps to either SPCERR_BAD_X or SPCERR_BAD_SPEC depending on context. */
const int spc_spxerr[] = {

  SPCERR_SUCCESS,		/*  0: SPXERR_SUCCESS          */
  SPCERR_NULL_POINTER,		/*  1: SPXERR_NULL_POINTER     */
  SPCERR_BAD_SPEC_PARAMS,	/*  2: SPXERR_BAD_SPEC_PARAMS  */
  SPCERR_BAD_SPEC_PARAMS	/*  3: SPXERR_BAD_SPEC_VAR     */
				/*  4: SPXERR_BAD_INSPEC_COORD */
};

/* Convenience macro for invoking wcserr_set(). */
#define SPC_ERRMSG(status) WCSERR_SET(status), spc_errmsg[status]


#define C 2.99792458e8

/*--------------------------------------------------------------------------*/

int spcini(struct spcprm *spc)

{
  register int k;

  if (spc == 0x0) return SPCERR_NULL_POINTER;

  spc->flag = 0;

  memset(spc->type, 0, 8);
  strcpy(spc->type, "    ");
  strcpy(spc->code, "   ");

  spc->crval = UNDEFINED;
  spc->restfrq =  0.0;
  spc->restwav =  0.0;

  for (k = 0; k < 7; k++) {
    spc->pv[k] = UNDEFINED;
  }

  for (k = 0; k < 6; k++) {
    spc->w[k] = 0.0;
  }

  spc->isGrism  = 0;
  spc->padding1 = 0;

  spc->err = 0x0;

  spc->padding2 = 0x0;
  spc->spxX2P = 0x0;
  spc->spxP2S = 0x0;
  spc->spxS2P = 0x0;
  spc->spxP2X = 0x0;

  return 0;
}

/*--------------------------------------------------------------------------*/

int spcfree(struct spcprm *spc)

{
  if (spc == 0x0) return SPCERR_NULL_POINTER;

  if (spc->err) {
    free(spc->err);
    spc->err = 0x0;
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int spcprt(const struct spcprm *spc)

{
  char hext[32];
  int  i;

  if (spc == 0x0) return SPCERR_NULL_POINTER;

  wcsprintf("       flag: %d\n", spc->flag);
  wcsprintf("       type: \"%s\"\n", spc->type);
  wcsprintf("       code: \"%s\"\n", spc->code);
  if (undefined(spc->crval)) {
    wcsprintf("      crval: UNDEFINED\n");
  } else {
    wcsprintf("      crval: %#- 11.5g\n", spc->crval);
  }
  wcsprintf("    restfrq: %f\n", spc->restfrq);
  wcsprintf("    restwav: %f\n", spc->restwav);

  wcsprintf("         pv:");
  if (spc->isGrism) {
    for (i = 0; i < 5; i++) {
      if (undefined(spc->pv[i])) {
        wcsprintf("  UNDEFINED   ");
      } else {
        wcsprintf("  %#- 11.5g", spc->pv[i]);
      }
    }
    wcsprintf("\n            ");
    for (i = 5; i < 7; i++) {
      if (undefined(spc->pv[i])) {
        wcsprintf("  UNDEFINED   ");
      } else {
        wcsprintf("  %#- 11.5g", spc->pv[i]);
      }
    }
    wcsprintf("\n");

  } else {
    wcsprintf(" (not used)\n");
  }

  wcsprintf("          w:");
  for (i = 0; i < 3; i++) {
    wcsprintf("  %#- 11.5g", spc->w[i]);
  }
  if (spc->isGrism) {
    wcsprintf("\n            ");
    for (i = 3; i < 6; i++) {
      wcsprintf("  %#- 11.5g", spc->w[i]);
    }
    wcsprintf("\n");
  } else {
    wcsprintf("  (remainder unused)\n");
  }

  wcsprintf("    isGrism: %d\n", spc->isGrism);

  WCSPRINTF_PTR("        err: ", spc->err, "\n");
  if (spc->err) {
    wcserr_prt(spc->err, "             ");
  }

  wcsprintf("     spxX2P: %s\n",
    wcsutil_fptr2str((int (*)(void))spc->spxX2P, hext));
  wcsprintf("     spxP2S: %s\n",
    wcsutil_fptr2str((int (*)(void))spc->spxP2S, hext));
  wcsprintf("     spxS2P: %s\n",
    wcsutil_fptr2str((int (*)(void))spc->spxS2P, hext));
  wcsprintf("     spxP2X: %s\n",
    wcsutil_fptr2str((int (*)(void))spc->spxP2X, hext));

  return 0;
}

/*--------------------------------------------------------------------------*/

int spcset(struct spcprm *spc)

{
  static const char *function = "spcset";

  char   ctype[9], ptype, xtype;
  int    restreq, status;
  double alpha, beta_r, crvalX, dn_r, dXdS, epsilon, G, m, lambda_r, n_r,
         t, restfrq, restwav, theta;
  struct wcserr **err;

  if (spc == 0x0) return SPCERR_NULL_POINTER;
  err = &(spc->err);

  if (undefined(spc->crval)) {
    return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
      "Spectral crval is undefined");
  }

  memset((spc->type)+4, 0, 4);
  spc->code[3] = '\0';
  wcsutil_blank_fill(4, spc->type);
  wcsutil_blank_fill(3, spc->code);
  spc->w[0] = 0.0;


  /* Analyse the spectral axis type. */
  memset(ctype, 0, 9);
  strncpy(ctype, spc->type, 4);
  if (*(spc->code) != ' ') {
    sprintf(ctype+4, "-%s", spc->code);
  }
  restfrq = spc->restfrq;
  restwav = spc->restwav;
  if ((status = spcspxe(ctype, spc->crval, restfrq, restwav, &ptype, &xtype,
                        &restreq, &crvalX, &dXdS, &(spc->err)))) {
    return status;
  }

  /* Satisfy rest frequency/wavelength requirements. */
  if (restreq) {
    if (restreq == 3 && restfrq == 0.0 && restwav == 0.0) {
      /* VRAD-V2F, VOPT-V2W, and ZOPT-V2W require the rest frequency or */
      /* wavelength for the S-P and P-X transformations but not for S-X */
      /* so supply a phoney value. */
      restwav = 1.0;
    }

    if (restfrq == 0.0) {
      restfrq = C/restwav;
    } else {
      restwav = C/restfrq;
    }

    if (ptype == 'F') {
      spc->w[0] = restfrq;
    } else if (ptype != 'V') {
      spc->w[0] = restwav;
    } else {
      if (xtype == 'F') {
        spc->w[0] = restfrq;
      } else {
        spc->w[0] = restwav;
      }
    }
  }

  spc->w[1] = crvalX;
  spc->w[2] = dXdS;


  /* Set pointers-to-functions for the linear part of the transformation. */
  if (ptype == 'F') {
    if (strcmp(spc->type, "FREQ") == 0) {
      /* Frequency. */
      spc->flag = FREQ;
      spc->spxP2S = 0x0;
      spc->spxS2P = 0x0;

    } else if (strcmp(spc->type, "AFRQ") == 0) {
      /* Angular frequency. */
      spc->flag = AFRQ;
      spc->spxP2S = freqafrq;
      spc->spxS2P = afrqfreq;

    } else if (strcmp(spc->type, "ENER") == 0) {
      /* Photon energy. */
      spc->flag = ENER;
      spc->spxP2S = freqener;
      spc->spxS2P = enerfreq;

    } else if (strcmp(spc->type, "WAVN") == 0) {
      /* Wave number. */
      spc->flag = WAVN;
      spc->spxP2S = freqwavn;
      spc->spxS2P = wavnfreq;

    } else if (strcmp(spc->type, "VRAD") == 0) {
      /* Radio velocity. */
      spc->flag = VRAD;
      spc->spxP2S = freqvrad;
      spc->spxS2P = vradfreq;
    }

  } else if (ptype == 'W') {
    if (strcmp(spc->type, "WAVE") == 0) {
      /* Vacuum wavelengths. */
      spc->flag = WAVE;
      spc->spxP2S = 0x0;
      spc->spxS2P = 0x0;

    } else if (strcmp(spc->type, "VOPT") == 0) {
      /* Optical velocity. */
      spc->flag = VOPT;
      spc->spxP2S = wavevopt;
      spc->spxS2P = voptwave;

    } else if (strcmp(spc->type, "ZOPT") == 0) {
      /* Redshift. */
      spc->flag = ZOPT;
      spc->spxP2S = wavezopt;
      spc->spxS2P = zoptwave;
    }

  } else if (ptype == 'A') {
    if (strcmp(spc->type, "AWAV") == 0) {
      /* Air wavelengths. */
      spc->flag = AWAV;
      spc->spxP2S = 0x0;
      spc->spxS2P = 0x0;
    }

  } else if (ptype == 'V') {
    if (strcmp(spc->type, "VELO") == 0) {
      /* Relativistic velocity. */
      spc->flag = VELO;
      spc->spxP2S = 0x0;
      spc->spxS2P = 0x0;

    } else if (strcmp(spc->type, "BETA") == 0) {
      /* Velocity ratio (v/c). */
      spc->flag = BETA;
      spc->spxP2S = velobeta;
      spc->spxS2P = betavelo;
    }
  }


  /* Set pointers-to-functions for the non-linear part of the spectral */
  /* transformation.                                                   */
  spc->isGrism = 0;
  if (xtype == 'F') {
    /* Axis is linear in frequency. */
    if (ptype == 'F') {
      spc->spxX2P = 0x0;
      spc->spxP2X = 0x0;

    } else if (ptype == 'W') {
      spc->spxX2P = freqwave;
      spc->spxP2X = wavefreq;

    } else if (ptype == 'A') {
      spc->spxX2P = freqawav;
      spc->spxP2X = awavfreq;

    } else if (ptype == 'V') {
      spc->spxX2P = freqvelo;
      spc->spxP2X = velofreq;
    }

    spc->flag += F2S;

  } else if (xtype == 'W' || xtype == 'w') {
    /* Axis is linear in vacuum wavelengths. */
    if (ptype == 'F') {
      spc->spxX2P = wavefreq;
      spc->spxP2X = freqwave;

    } else if (ptype == 'W') {
      spc->spxX2P = 0x0;
      spc->spxP2X = 0x0;

    } else if (ptype == 'A') {
      spc->spxX2P = waveawav;
      spc->spxP2X = awavwave;

    } else if (ptype == 'V') {
      spc->spxX2P = wavevelo;
      spc->spxP2X = velowave;
    }

    if (xtype == 'W') {
      spc->flag += W2S;
    } else {
      /* Grism in vacuum. */
      spc->isGrism = 1;
      spc->flag += GRI;
    }

  } else if (xtype == 'A' || xtype == 'a') {
    /* Axis is linear in air wavelengths. */
    if (ptype == 'F') {
      spc->spxX2P = awavfreq;
      spc->spxP2X = freqawav;

    } else if (ptype == 'W') {
      spc->spxX2P = awavwave;
      spc->spxP2X = waveawav;

    } else if (ptype == 'A') {
      spc->spxX2P = 0x0;
      spc->spxP2X = 0x0;

    } else if (ptype == 'V') {
      spc->spxX2P = awavvelo;
      spc->spxP2X = veloawav;
    }

    if (xtype == 'A') {
      spc->flag += A2S;
    } else {
      /* Grism in air. */
      spc->isGrism = 2;
      spc->flag += GRA;
    }

  } else if (xtype == 'V') {
    /* Axis is linear in relativistic velocity. */
    if (ptype == 'F') {
      spc->spxX2P = velofreq;
      spc->spxP2X = freqvelo;

    } else if (ptype == 'W') {
      spc->spxX2P = velowave;
      spc->spxP2X = wavevelo;

    } else if (ptype == 'A') {
      spc->spxX2P = veloawav;
      spc->spxP2X = awavvelo;

    } else if (ptype == 'V') {
      spc->spxX2P = 0x0;
      spc->spxP2X = 0x0;
    }

    spc->flag += V2S;
  }


  /* Check for grism axes. */
  if (spc->isGrism) {
    /* Axis is linear in "grism parameter"; work in wavelength. */
    lambda_r = crvalX;

    /* Set defaults. */
    if (undefined(spc->pv[0])) spc->pv[0] = 0.0;
    if (undefined(spc->pv[1])) spc->pv[1] = 0.0;
    if (undefined(spc->pv[2])) spc->pv[2] = 0.0;
    if (undefined(spc->pv[3])) spc->pv[3] = 1.0;
    if (undefined(spc->pv[4])) spc->pv[4] = 0.0;
    if (undefined(spc->pv[5])) spc->pv[5] = 0.0;
    if (undefined(spc->pv[6])) spc->pv[6] = 0.0;

    /* Compute intermediaries. */
    G       = spc->pv[0];
    m       = spc->pv[1];
    alpha   = spc->pv[2];
    n_r     = spc->pv[3];
    dn_r    = spc->pv[4];
    epsilon = spc->pv[5];
    theta   = spc->pv[6];

    t = G*m/cosd(epsilon);
    beta_r = asind(t*lambda_r - n_r*sind(alpha));

    t -= dn_r*sind(alpha);

    spc->w[1] = -tand(theta);
    spc->w[2] *= t / (cosd(beta_r)*cosd(theta)*cosd(theta));
    spc->w[3] = beta_r + theta;
    spc->w[4] = (n_r - dn_r*lambda_r)*sind(alpha);
    spc->w[5] = 1.0 / t;
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int spcx2s(
  struct spcprm *spc,
  int nx,
  int sx,
  int sspec,
  const double x[],
  double spec[],
  int stat[])

{
  static const char *function = "spcx2s";

  int statP2S, status = 0, statX2P;
  double beta;
  register int ix;
  register int *statp;
  register const double *xp;
  register double *specp;
  struct wcserr **err;

  /* Initialize. */
  if (spc == 0x0) return SPCERR_NULL_POINTER;
  err = &(spc->err);

  if (spc->flag == 0) {
    if ((status = spcset(spc))) return status;
  }

  /* Convert intermediate world coordinate x to X. */
  xp = x;
  specp = spec;
  statp = stat;
  for (ix = 0; ix < nx; ix++, xp += sx, specp += sspec) {
    *specp = spc->w[1] + (*xp)*spc->w[2];
    *(statp++) = 0;
  }

  /* If X is the grism parameter then convert it to wavelength. */
  if (spc->isGrism) {
    specp = spec;
    for (ix = 0; ix < nx; ix++, specp += sspec) {
      beta = atand(*specp) + spc->w[3];
      *specp = (sind(beta) + spc->w[4]) * spc->w[5];
    }
  }

  /* Apply the non-linear step of the algorithm chain to convert the    */
  /* X-type spectral variable to P-type intermediate spectral variable. */
  if (spc->spxX2P) {
    if ((statX2P = spc->spxX2P(spc->w[0], nx, sspec, sspec, spec, spec,
                               stat))) {
      if (statX2P == SPXERR_BAD_INSPEC_COORD) {
        status = SPCERR_BAD_X;
      } else if (statX2P == SPXERR_BAD_SPEC_PARAMS) {
        return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
          "Invalid spectral parameters: Frequency or wavelength is 0");
      } else {
        return wcserr_set(SPC_ERRMSG(spc_spxerr[statX2P]));
      }
    }
  }

  /* Apply the linear step of the algorithm chain to convert P-type  */
  /* intermediate spectral variable to the required S-type variable. */
  if (spc->spxP2S) {
    if ((statP2S = spc->spxP2S(spc->w[0], nx, sspec, sspec, spec, spec,
                               stat))) {
      if (statP2S == SPXERR_BAD_INSPEC_COORD) {
        status = SPCERR_BAD_X;
      } else if (statP2S == SPXERR_BAD_SPEC_PARAMS) {
        return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
          "Invalid spectral parameters: Frequency or wavelength is 0");
      } else {
        return wcserr_set(SPC_ERRMSG(spc_spxerr[statP2S]));
      }
    }
  }

  if (status) {
    wcserr_set(SPC_ERRMSG(status));
  }
  return status;
}

/*--------------------------------------------------------------------------*/

int spcs2x(
  struct spcprm *spc,
  int nspec,
  int sspec,
  int sx,
  const double spec[],
  double x[],
  int stat[])

{
  static const char *function = "spcs2x";

  int statP2X, status = 0, statS2P;
  double beta, s;
  register int ispec;
  register int *statp;
  register const double *specp;
  register double *xp;
  struct wcserr **err;

  /* Initialize. */
  if (spc == 0x0) return SPCERR_NULL_POINTER;
  err = &(spc->err);

  if (spc->flag == 0) {
    if ((status = spcset(spc))) return status;
  }

  /* Apply the linear step of the algorithm chain to convert the S-type */
  /* spectral variable to P-type intermediate spectral variable.        */
  if (spc->spxS2P) {
    if ((statS2P = spc->spxS2P(spc->w[0], nspec, sspec, sx, spec, x, stat))) {
      if (statS2P == SPXERR_BAD_INSPEC_COORD) {
        status = SPCERR_BAD_SPEC;
      } else if (statS2P == SPXERR_BAD_SPEC_PARAMS) {
        return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
          "Invalid spectral parameters: Frequency or wavelength is 0");
      } else {
        return wcserr_set(SPC_ERRMSG(spc_spxerr[statS2P]));
      }
    }

  } else {
    /* Just a copy. */
    xp = x;
    specp = spec;
    statp = stat;
    for (ispec = 0; ispec < nspec; ispec++, specp += sspec, xp += sx) {
      *xp = *specp;
      *(statp++) = 0;
    }
  }


  /* Apply the non-linear step of the algorithm chain to convert P-type */
  /* intermediate spectral variable to X-type spectral variable. */
  if (spc->spxP2X) {
    if ((statP2X = spc->spxP2X(spc->w[0], nspec, sx, sx, x, x, stat))) {
      if (statP2X == SPCERR_BAD_SPEC) {
        status = SPCERR_BAD_SPEC;
      } else if (statP2X == SPXERR_BAD_SPEC_PARAMS) {
        return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
          "Invalid spectral parameters: Frequency or wavelength is 0");
      } else {
        return wcserr_set(SPC_ERRMSG(spc_spxerr[statP2X]));
      }
    }
  }

  if (spc->isGrism) {
    /* Convert X-type spectral variable (wavelength) to grism parameter. */
    xp = x;
    statp = stat;
    for (ispec = 0; ispec < nspec; ispec++, xp += sx, statp++) {
      if (*statp) continue;

      s = *xp/spc->w[5] - spc->w[4];
      if (fabs(s) <= 1.0) {
        beta = asind(s);
        *xp = tand(beta - spc->w[3]);
      } else {
        *statp = 1;
      }
    }
  }


  /* Convert X-type spectral variable to intermediate world coordinate x. */
  xp = x;
  statp = stat;
  for (ispec = 0; ispec < nspec; ispec++, xp += sx) {
    if (*(statp++)) continue;

    *xp -= spc->w[1];
    *xp /= spc->w[2];
  }

  if (status) {
    wcserr_set(SPC_ERRMSG(status));
  }
  return status;
}

/*--------------------------------------------------------------------------*/

int spctyp(
  const char ctypei[9],
  char stype[],
  char scode[],
  char sname[],
  char units[],
  char *ptype,
  char *xtype,
  int  *restreq)

{
  return spctype(
    ctypei, stype, scode, sname, units, ptype, xtype, restreq, NULL);
}

/* : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :  */

int spctype(
  const char ctypei[9],
  char stype[],
  char scode[],
  char sname[],
  char units[],
  char *ptype,
  char *xtype,
  int  *restreq,
  struct wcserr **err)

{
  static const char *function = "spctype";

  char ctype[9], ptype_t, sname_t[32], units_t[8], xtype_t;
  int  restreq_t = 0;

  if (err) *err = 0x0;

  /* Copy with blank padding. */
  sprintf(ctype, "%-8.8s", ctypei);
  ctype[8] = '\0';

  /* Validate the S-type spectral variable. */
  if (strncmp(ctype, "FREQ", 4) == 0) {
    strcpy(sname_t, "Frequency");
    strcpy(units_t, "Hz");
    ptype_t = 'F';
  } else if (strncmp(ctype, "AFRQ", 4) == 0) {
    strcpy(sname_t, "Angular frequency");
    strcpy(units_t, "rad/s");
    ptype_t = 'F';
  } else if (strncmp(ctype, "ENER", 4) == 0) {
    strcpy(sname_t, "Photon energy");
    strcpy(units_t, "J");
    ptype_t = 'F';
  } else if (strncmp(ctype, "WAVN", 4) == 0) {
    strcpy(sname_t, "Wavenumber");
    strcpy(units_t, "/m");
    ptype_t = 'F';
  } else if (strncmp(ctype, "VRAD", 4) == 0) {
    strcpy(sname_t, "Radio velocity");
    strcpy(units_t, "m/s");
    ptype_t = 'F';
    restreq_t = 1;
  } else if (strncmp(ctype, "WAVE", 4) == 0) {
    strcpy(sname_t, "Vacuum wavelength");
    strcpy(units_t, "m");
    ptype_t = 'W';
  } else if (strncmp(ctype, "VOPT", 4) == 0) {
    strcpy(sname_t, "Optical velocity");
    strcpy(units_t, "m/s");
    ptype_t = 'W';
    restreq_t = 1;
  } else if (strncmp(ctype, "ZOPT", 4) == 0) {
    strcpy(sname_t, "Redshift");
    strcpy(units_t, "");
    ptype_t = 'W';
    restreq_t = 1;
  } else if (strncmp(ctype, "AWAV", 4) == 0) {
    strcpy(sname_t, "Air wavelength");
    strcpy(units_t, "m");
    ptype_t = 'A';
  } else if (strncmp(ctype, "VELO", 4) == 0) {
    strcpy(sname_t, "Relativistic velocity");
    strcpy(units_t, "m/s");
    ptype_t = 'V';
  } else if (strncmp(ctype, "BETA", 4) == 0) {
    strcpy(sname_t, "Velocity ratio (v/c)");
    strcpy(units_t, "");
    ptype_t = 'V';
  } else {
    return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
      "Unknown spectral type '%s'", ctype);
  }


  /* Determine X-type and validate the spectral algorithm code. */
  if ((xtype_t = ctype[5]) == ' ') {
    /* The algorithm code must be completely blank. */
    if (strcmp(ctype+4, "    ") != 0) {
      return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
        "Invalid spectral algorithm '%s'", ctype+4);
    }

    xtype_t = ptype_t;

  } else if (ctype[4] != '-') {
    return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
      "Invalid spectral type '%s'", ctype);

  } else if (strcmp(ctype+5, "LOG") == 0 || strcmp(ctype+5, "TAB") == 0) {
    /* Logarithmic or tabular axis, not linear in any spectral type. */

  } else if (xtype_t == 'G') {
    /* Validate the algorithm code. */
    if (ctype[6] != 'R') {
      return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
        "Invalid spectral algorithm '%s'", xtype_t);
    }

    /* Grism coordinates... */
    if (ctype[7] == 'I') {
      /* ...in vacuum. */
      xtype_t = 'w';
    } else if (ctype[7] == 'A') {
      /* ...in air. */
      xtype_t = 'a';
    } else {
      return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
        "Invalid spectral algorithm '%s'", xtype_t);
    }

  } else if (ctype[6] != '2') {
    /* Algorithm code has invalid syntax. */
    return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
      "Invalid spectral algorithm syntax '%s'", xtype_t);
  } else if (ctype[7] != ptype_t && ctype[7] != '?') {
    /* The P-, and S-type variables are inconsistent. */
    return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
      "In spectral type '%s', P- and S-type variables are inconsistent",
      ctype);

  } else if (ctype[7] == ctype[5]) {
    /* Degenerate algorithm code. */
    sprintf(ctype+4, "    ");
  }


  /* Rest freq/wavelength required for transformation between P and X? */
  if (strchr("FWAwa", (int)xtype_t)) {
    if (ptype_t == 'V') {
      restreq_t += 2;
    }
  } else if (xtype_t == 'V') {
    if (strchr("FWAwa", (int)ptype_t)) {
      restreq_t += 2;
    }
  } else if (strchr("LT", (int)xtype_t) == 0) {
    /* Invalid X-type variable code. */
    return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
      "In spectral type '%s', invalid X-type variable code", ctype);
  }


  /* Copy results. */
  if (stype) {
    strncpy(stype, ctype, 4);
    stype[4] = '\0';
  }
  if (scode) strcpy(scode, ctype+5);
  if (sname) strcpy(sname, sname_t);
  if (units) strcpy(units, units_t);
  if (ptype) *ptype = ptype_t;
  if (xtype) *xtype = xtype_t;
  if (restreq) *restreq = restreq_t;


  return 0;
}

/*--------------------------------------------------------------------------*/

int spcspx(
  const char ctypeS[9],
  double crvalS,
  double restfrq,
  double restwav,
  char *ptype,
  char *xtype,
  int *restreq,
  double *crvalX,
  double *dXdS)

{
  return spcspxe(ctypeS, crvalS, restfrq, restwav, ptype, xtype, restreq,
                 crvalX, dXdS, 0x0);
}

/* : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :  */

int spcspxe(
  const char ctypeS[9],
  double crvalS,
  double restfrq,
  double restwav,
  char *ptype,
  char *xtype,
  int *restreq,
  double *crvalX,
  double *dXdS,
  struct wcserr **err)

{
  static const char *function = "spcspxe";

  char scode[4], stype[5], type[8];
  int  status;
  double dPdS, dXdP;
  struct spxprm spx;


  /* Analyse the spectral axis code. */
  if ((status = spctype(ctypeS, stype, scode, 0x0, 0x0, ptype, xtype, restreq,
                        err))) {
    return status;
  }

  if (strchr("LT", (int)(*xtype))) {
    /* Can't handle logarithmic or tabular coordinates. */
    return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
      "Can't handle logarithmic or tabular coordinates");
  }

  /* Do we have rest frequency and/or wavelength as required? */
  if ((*restreq)%3 && restfrq == 0.0 && restwav == 0.0) {
    return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
      "Missing required rest frequency or wavelength");
  }

  /* Compute all spectral parameters and their derivatives. */
  strcpy(type, stype);
  spx.err = (err ? *err : 0x0);
  if ((status = specx(type, crvalS, restfrq, restwav, &spx))) {
    status = spc_spxerr[status];
    if (err) {
      if ((*err = spx.err)) {
        (*err)->status = status;
      }
    } else {
      free(spx.err);
    }
    return status;
  }


  /* Transform S-P (linear) and P-X (non-linear). */
  dPdS = 0.0;
  dXdP = 0.0;
  if (*ptype == 'F') {
    if (strcmp(stype, "FREQ") == 0) {
      dPdS = 1.0;
    } else if (strcmp(stype, "AFRQ") == 0) {
      dPdS = spx.dfreqafrq;
    } else if (strcmp(stype, "ENER") == 0) {
      dPdS = spx.dfreqener;
    } else if (strcmp(stype, "WAVN") == 0) {
      dPdS = spx.dfreqwavn;
    } else if (strcmp(stype, "VRAD") == 0) {
      dPdS = spx.dfreqvrad;
    }

    if (*xtype == 'F') {
      *crvalX = spx.freq;
      dXdP = 1.0;
    } else if (*xtype == 'W' || *xtype == 'w') {
      *crvalX = spx.wave;
      dXdP = spx.dwavefreq;
    } else if (*xtype == 'A' || *xtype == 'a') {
      *crvalX = spx.awav;
      dXdP = spx.dawavfreq;
    } else if (*xtype == 'V') {
      *crvalX = spx.velo;
      dXdP = spx.dvelofreq;
    }

  } else if (*ptype == 'W' || *ptype == 'w') {
    if (strcmp(stype, "WAVE") == 0) {
      dPdS = 1.0;
    } else if (strcmp(stype, "VOPT") == 0) {
      dPdS = spx.dwavevopt;
    } else if (strcmp(stype, "ZOPT") == 0) {
      dPdS = spx.dwavezopt;
    }

    if (*xtype == 'F') {
      *crvalX = spx.freq;
      dXdP = spx.dfreqwave;
    } else if (*xtype == 'W' || *xtype == 'w') {
      *crvalX = spx.wave;
      dXdP = 1.0;
    } else if (*xtype == 'A' || *xtype == 'a') {
      *crvalX = spx.awav;
      dXdP = spx.dawavwave;
    } else if (*xtype == 'V') {
      *crvalX = spx.velo;
      dXdP = spx.dvelowave;
    }

  } else if (*ptype == 'A' || *ptype == 'a') {
    if (strcmp(stype, "AWAV") == 0) {
      dPdS = 1.0;
    }

    if (*xtype == 'F') {
      *crvalX = spx.freq;
      dXdP = spx.dfreqawav;
    } else if (*xtype == 'W' || *xtype == 'w') {
      *crvalX = spx.wave;
      dXdP = spx.dwaveawav;
    } else if (*xtype == 'A' || *xtype == 'a') {
      *crvalX = spx.awav;
      dXdP = 1.0;
    } else if (*xtype == 'V') {
      *crvalX = spx.velo;
      dXdP = spx.dveloawav;
    }

  } else if (*ptype == 'V') {
    if (strcmp(stype, "VELO") == 0) {
      dPdS = 1.0;
    } else if (strcmp(stype, "BETA") == 0) {
      dPdS = spx.dvelobeta;
    }

    if (*xtype == 'F') {
      *crvalX = spx.freq;
      dXdP = spx.dfreqvelo;
    } else if (*xtype == 'W' || *xtype == 'w') {
      *crvalX = spx.wave;
      dXdP = spx.dwavevelo;
    } else if (*xtype == 'A' || *xtype == 'a') {
      *crvalX = spx.awav;
      dXdP = spx.dawavvelo;
    } else if (*xtype == 'V') {
      *crvalX = spx.velo;
      dXdP = 1.0;
    }
  }

  *dXdS = dXdP * dPdS;

  return 0;
}

/*--------------------------------------------------------------------------*/

int spcxps(
  const char ctypeS[9],
  double crvalX,
  double restfrq,
  double restwav,
  char *ptype,
  char *xtype,
  int *restreq,
  double *crvalS,
  double *dSdX)

{
  return spcxpse(ctypeS, crvalX, restfrq, restwav, ptype, xtype, restreq,
                 crvalS, dSdX, NULL);
}

/* : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :  */

int spcxpse(
  const char ctypeS[9],
  double crvalX,
  double restfrq,
  double restwav,
  char *ptype,
  char *xtype,
  int *restreq,
  double *crvalS,
  double *dSdX,
  struct wcserr **err)

{
  static const char *function = "spcxpse";

  char scode[4], stype[5], type[8];
  int  status;
  double dPdX, dSdP;
  struct spxprm spx;

  /* Analyse the spectral axis type. */
  if ((status = spctype(ctypeS, stype, scode, 0x0, 0x0, ptype, xtype, restreq,
                        err))) {
    return status;
  }

  if (strchr("LT", (int)(*xtype))) {
    /* Can't handle logarithmic or tabular coordinates. */
    return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
      "Can't handle logarithmic or tabular coordinates");
  }

  /* Do we have rest frequency and/or wavelength as required? */
  if ((*restreq)%3 && restfrq == 0.0 && restwav == 0.0) {
    return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
      "Missing required rest frequency or wavelength");
  }

  /* Compute all spectral parameters and their derivatives. */
  if (*xtype == 'F') {
    strcpy(type, "FREQ");
  } else if (*xtype == 'W' || *xtype == 'w') {
    strcpy(type, "WAVE");
  } else if (*xtype == 'A' || *xtype == 'a') {
    strcpy(type, "AWAV");
  } else if (*xtype == 'V') {
    strcpy(type, "VELO");
  }

  spx.err = (err ? *err : 0x0);
  if (specx(type, crvalX, restfrq, restwav, &spx)) {
    status = spc_spxerr[status];
    if (err) {
      if ((*err = spx.err)) {
        (*err)->status = status;
      }
    } else {
      free(spx.err);
    }
    return status;
  }


  /* Transform X-P (non-linear) and P-S (linear). */
  dPdX = 0.0;
  dSdP = 0.0;
  if (*ptype == 'F') {
    if (*xtype == 'F') {
      dPdX = 1.0;
    } else if (*xtype == 'W' || *xtype == 'w') {
      dPdX = spx.dfreqwave;
    } else if (*xtype == 'A' || *xtype == 'a') {
      dPdX = spx.dfreqawav;
    } else if (*xtype == 'V') {
      dPdX = spx.dfreqvelo;
    }

    if (strcmp(stype, "FREQ") == 0) {
      *crvalS = spx.freq;
      dSdP = 1.0;
    } else if (strcmp(stype, "AFRQ") == 0) {
      *crvalS = spx.afrq;
      dSdP = spx.dafrqfreq;
    } else if (strcmp(stype, "ENER") == 0) {
      *crvalS = spx.ener;
      dSdP = spx.denerfreq;
    } else if (strcmp(stype, "WAVN") == 0) {
      *crvalS = spx.wavn;
      dSdP = spx.dwavnfreq;
    } else if (strcmp(stype, "VRAD") == 0) {
      *crvalS = spx.vrad;
      dSdP = spx.dvradfreq;
    }

  } else if (*ptype == 'W') {
    if (*xtype == 'F') {
      dPdX = spx.dwavefreq;
    } else if (*xtype == 'W' || *xtype == 'w') {
      dPdX = 1.0;
    } else if (*xtype == 'A' || *xtype == 'a') {
      dPdX = spx.dwaveawav;
    } else if (*xtype == 'V') {
      dPdX = spx.dwavevelo;
    }

    if (strcmp(stype, "WAVE") == 0) {
      *crvalS = spx.wave;
      dSdP = 1.0;
    } else if (strcmp(stype, "VOPT") == 0) {
      *crvalS = spx.vopt;
      dSdP = spx.dvoptwave;
    } else if (strcmp(stype, "ZOPT") == 0) {
      *crvalS = spx.zopt;
      dSdP = spx.dzoptwave;
    }

  } else if (*ptype == 'A') {
    if (*xtype == 'F') {
      dPdX = spx.dawavfreq;
    } else if (*xtype == 'W' || *xtype == 'w') {
      dPdX = spx.dawavwave;
    } else if (*xtype == 'A' || *xtype == 'a') {
      dPdX = 1.0;
    } else if (*xtype == 'V') {
      dPdX = spx.dawavvelo;
    }

    if (strcmp(stype, "AWAV") == 0) {
      *crvalS = spx.awav;
      dSdP = 1.0;
    }

  } else if (*ptype == 'V') {
    if (*xtype == 'F') {
      dPdX = spx.dvelofreq;
    } else if (*xtype == 'W' || *xtype == 'w') {
      dPdX = spx.dvelowave;
    } else if (*xtype == 'A' || *xtype == 'a') {
      dPdX = spx.dveloawav;
    } else if (*xtype == 'V') {
      dPdX = 1.0;
    }

    if (strcmp(stype, "VELO") == 0) {
      *crvalS = spx.velo;
      dSdP = 1.0;
    } else if (strcmp(stype, "BETA") == 0) {
      *crvalS = spx.beta;
      dSdP = spx.dbetavelo;
    }
  }

  *dSdX = dSdP * dPdX;

  return 0;
}

/*--------------------------------------------------------------------------*/

int spctrn(
  const char ctypeS1[9],
  double crvalS1,
  double cdeltS1,
  double restfrq,
  double restwav,
  char   ctypeS2[9],
  double *crvalS2,
  double *cdeltS2)

{
  return spctrne(ctypeS1, crvalS1, cdeltS1, restfrq, restwav,
                 ctypeS2, crvalS2, cdeltS2, NULL);
}

/* : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :  */

int spctrne(
  const char ctypeS1[9],
  double crvalS1,
  double cdeltS1,
  double restfrq,
  double restwav,
  char   ctypeS2[9],
  double *crvalS2,
  double *cdeltS2,
  struct wcserr **err)

{
  static const char *function = "spctrne";

  char *cp, ptype1, ptype2, stype1[5], stype2[5], xtype1, xtype2;
  int  restreq, status;
  double crvalX, dS2dX, dXdS1;

  if (restfrq == 0.0 && restwav == 0.0) {
    /* If translating between two velocity-characteristic types, or between
       two wave-characteristic types, then we may need to set a dummy rest
       frequency or wavelength to perform the calculations. */
    strncpy(stype1, ctypeS1, 4);
    strncpy(stype2, ctypeS2, 4);
    stype1[4] = stype2[4] = '\0';
    if ((strstr("VRAD VOPT ZOPT VELO BETA", stype1) != 0x0) ==
        (strstr("VRAD VOPT ZOPT VELO BETA", stype2) != 0x0)) {
      restwav = 1.0;
    }
  }

  if ((status = spcspxe(ctypeS1, crvalS1, restfrq, restwav, &ptype1, &xtype1,
                        &restreq, &crvalX, &dXdS1, err))) {
    return status;
  }

  /* Pad with blanks. */
  ctypeS2[8] = '\0';
  for (cp = ctypeS2; *cp; cp++);
  while (cp < ctypeS2+8) *(cp++) = ' ';

  if (strncmp(ctypeS2+5, "???", 3) == 0) {
    /* Set the algorithm code if required. */
    if (xtype1 == 'w') {
      strcpy(ctypeS2+5, "GRI");
    } else if (xtype1 == 'a') {
      strcpy(ctypeS2+5, "GRA");
    } else {
      ctypeS2[5] = xtype1;
      ctypeS2[6] = '2';
    }
  }

  if ((status = spcxpse(ctypeS2, crvalX, restfrq, restwav, &ptype2, &xtype2,
                        &restreq, crvalS2, &dS2dX, err))) {
    return status;
  }

  /* Are the X-types compatible? */
  if (xtype2 != xtype1) {
    return wcserr_set(WCSERR_SET(SPCERR_BAD_SPEC_PARAMS),
      "Incompatible X-types '%c' and '%c'", xtype1, xtype2);
  }

  if (ctypeS2[7] == '?') {
    if (ptype2 == xtype2) {
      strcpy(ctypeS2+4, "    ");
    } else {
      ctypeS2[7] = ptype2;
    }
  }

  *cdeltS2 = dS2dX * dXdS1 * cdeltS1;

  return 0;
}

/*--------------------------------------------------------------------------*/

int spcaips(
  const char ctypeA[9],
  int  velref,
  char ctype[9],
  char specsys[9])

{
  const char *frames[] = {"LSRK", "BARYCENT", "TOPOCENT",
                          "LSRD", "GEOCENTR", "SOURCE", "GALACTOC"};
  char *fcode;
  int  ivf, status;

  /* Make a null-filled copy of ctypeA. */
  if (ctype != ctypeA) strncpy(ctype, ctypeA, 8);
  ctype[8] = '\0';
  wcsutil_null_fill(9, ctype);
  *specsys = '\0';

  /* Is it a recognized AIPS-convention type? */
  status = SPCERR_NO_CHANGE;
  if (strncmp(ctype, "FREQ", 4) == 0 ||
      strncmp(ctype, "VELO", 4) == 0 ||
      strncmp(ctype, "FELO", 4) == 0) {
    /* Look for the Doppler frame. */
    if (*(fcode = ctype+4)) {
      if (strcmp(fcode, "-LSR") == 0) {
        strcpy(specsys, "LSRK");
      } else if (strcmp(fcode, "-HEL") == 0) {
        strcpy(specsys, "BARYCENT");
      } else if (strcmp(fcode, "-OBS") == 0) {
        strcpy(specsys, "TOPOCENT");
      } else {
        /* Not a recognized AIPS spectral type. */
        return SPCERR_NO_CHANGE;
      }

      *fcode = '\0';
      status = 0;
    }

    /* VELREF takes precedence if present. */
    ivf = velref%256;
    if (0 < ivf && ivf <= 7) {
      strcpy(specsys, frames[ivf-1]);
      status = 0;
    } else if (ivf) {
      status = SPCERR_BAD_SPEC_PARAMS;
    }

    if (strcmp(ctype, "VELO") == 0) {
      /* Check that we found an AIPS-convention Doppler frame. */
      if (*specsys) {
        /* 'VELO' in AIPS means radio or optical depending on VELREF. */
        ivf = velref/256;
        if (ivf == 0) {
          strcpy(ctype, "VOPT");
        } else if (ivf == 1) {
          strcpy(ctype, "VRAD");
        } else {
          status = SPCERR_BAD_SPEC_PARAMS;
        }
      }
    } else if (strcmp(ctype, "FELO") == 0) {
      /* Uniform in frequency but expressed as an optical velocity (strictly
         we should also have found an AIPS-convention Doppler frame). */
      strcpy(ctype, "VOPT-F2W");
      if (status < 0) status = 0;
    }
  }

  return status;
}
