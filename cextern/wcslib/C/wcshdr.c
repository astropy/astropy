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
  $Id: wcshdr.c,v 5.11 2015/10/18 09:13:05 mcalabre Exp $
*===========================================================================*/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcserr.h"
#include "wcsmath.h"
#include "wcsutil.h"
#include "wcshdr.h"
#include "tab.h"
#include "dis.h"
#include "wcs.h"

extern const int WCSSET;

extern const int DIS_DOTPD;

/* Map status return value to message. */
const char *wcshdr_errmsg[] = {
  "Success",
  "Null wcsprm pointer passed",
  "Memory allocation failed",
  "Invalid column selection",
  "Fatal error returned by Flex parser",
  "Invalid tabular parameters"};

/* Map error returns for lower-level routines. */
const int wcshdr_taberr[] = {
  WCSHDRERR_SUCCESS,		/*  0: TABERR_SUCCESS         */
  WCSHDRERR_NULL_POINTER,	/*  1: TABERR_NULL_POINTER    */
  WCSHDRERR_MEMORY,		/*  2: TABERR_MEMORY          */
  WCSHDRERR_BAD_TABULAR_PARAMS	/*  3: TABERR_BAD_PARAMS      */
				/*  4: TABERR_BAD_X           */
				/*  5: TABERR_BAD_WORLD       */
};

/* Convenience macro for invoking wcserr_set(). */
#define WCSHDR_ERRMSG(status) WCSERR_SET(status), wcshdr_errmsg[status]

/* Internal helper functions, not for general use. */
static void wcshdo_format(int, int, const double [], char *);
static void wcshdo_tpdterm(int, int, char *);
static void wcshdo_util(int, const char [], const char [], int, const char [],
  int, int, int, char, int, int [], char [], const char [], int *, char **,
  int *);

/*--------------------------------------------------------------------------*/

int wcstab(struct wcsprm *wcs)

{
  static const char *function = "wcstab";

  char (*PSi_0a)[72] = 0x0, (*PSi_1a)[72] = 0x0, (*PSi_2a)[72] = 0x0;
  int  *PVi_1a = 0x0, *PVi_2a = 0x0, *PVi_3a = 0x0, *tabax, *tabidx = 0x0;
  int   getcrd, i, ip, itab, itabax, j, jtabax, m, naxis, ntabax, status;
  struct wtbarr *wtbp;
  struct tabprm *tabp;
  struct wcserr **err;

  if (wcs == 0x0) return WCSHDRERR_NULL_POINTER;
  err = &(wcs->err);

  /* Free memory previously allocated by wcstab(). */
  if (wcs->flag != -1 && wcs->m_flag == WCSSET) {
    if (wcs->wtb == wcs->m_wtb) wcs->wtb = 0x0;
    if (wcs->tab == wcs->m_tab) wcs->tab = 0x0;

    if (wcs->m_wtb) free(wcs->m_wtb);
    if (wcs->m_tab) {
      for (j = 0; j < wcs->ntab; j++) {
        tabfree(wcs->m_tab + j);
      }

      free(wcs->m_tab);
    }
  }

  wcs->ntab = 0;
  wcs->nwtb = 0;
  wcs->wtb  = 0x0;
  wcs->tab  = 0x0;


  /* Determine the number of -TAB axes. */
  naxis = wcs->naxis;
  if (!(tabax = calloc(naxis, sizeof(int)))) {
    return wcserr_set(WCSHDR_ERRMSG(WCSHDRERR_MEMORY));
  }

  ntabax = 0;
  for (i = 0; i < naxis; i++) {
    /* Null fill. */
    wcsutil_null_fill(72, wcs->ctype[i]);

    if (!strcmp(wcs->ctype[i]+4, "-TAB")) {
      tabax[i] = ntabax++;
    } else {
      tabax[i] = -1;
    }
  }

  if (ntabax == 0) {
    /* No lookup tables. */
    status = 0;
    goto cleanup;
  }


  /* Collect information from the PSi_ma and PVi_ma keyvalues. */
  if (!((PSi_0a = calloc(ntabax, sizeof(char[72]))) &&
        (PVi_1a = calloc(ntabax, sizeof(int)))      &&
        (PVi_2a = calloc(ntabax, sizeof(int)))      &&
        (PSi_1a = calloc(ntabax, sizeof(char[72]))) &&
        (PSi_2a = calloc(ntabax, sizeof(char[72]))) &&
        (PVi_3a = calloc(ntabax, sizeof(int)))      &&
        (tabidx = calloc(ntabax, sizeof(int))))) {
    status = wcserr_set(WCSHDR_ERRMSG(WCSHDRERR_MEMORY));
    goto cleanup;
  }

  for (itabax = 0; itabax < ntabax; itabax++) {
    /* Remember that calloc() zeroes allocated memory. */
    PVi_1a[itabax] = 1;
    PVi_2a[itabax] = 1;
    PVi_3a[itabax] = 1;
  }

  for (ip = 0; ip < wcs->nps; ip++) {
    itabax = tabax[wcs->ps[ip].i - 1];
    if (itabax >= 0) {
      switch (wcs->ps[ip].m) {
      case 0:
        /* EXTNAME. */
        strcpy(PSi_0a[itabax], wcs->ps[ip].value);
        wcsutil_null_fill(72, PSi_0a[itabax]);
        break;
      case 1:
        /* TTYPEn for coordinate array. */
        strcpy(PSi_1a[itabax], wcs->ps[ip].value);
        wcsutil_null_fill(72, PSi_1a[itabax]);
        break;
      case 2:
        /* TTYPEn for index vector. */
        strcpy(PSi_2a[itabax], wcs->ps[ip].value);
        wcsutil_null_fill(72, PSi_2a[itabax]);
        break;
      }
    }
  }

  for (ip = 0; ip < wcs->npv; ip++) {
    itabax = tabax[wcs->pv[ip].i - 1];
    if (itabax >= 0) {
      switch (wcs->pv[ip].m) {
      case 1:
        /* EXTVER. */
        PVi_1a[itabax] = (int)(wcs->pv[ip].value + 0.5);
        break;
      case 2:
        /* EXTLEVEL. */
        PVi_2a[itabax] = (int)(wcs->pv[ip].value + 0.5);
        break;
      case 3:
        /* Table axis number. */
        PVi_3a[itabax] = (int)(wcs->pv[ip].value + 0.5);
        break;
      }
    }
  }


  /* Determine the number of independent tables. */
  for (itabax = 0; itabax < ntabax; itabax++) {
    /* These have no defaults. */
    if (!PSi_0a[itabax][0] || !PSi_1a[itabax][0]) {
      status = wcserr_set(WCSERR_SET(WCSHDRERR_BAD_TABULAR_PARAMS),
        "Invalid tabular parameters: PSi_0a and PSi_1a must be specified");
      goto cleanup;
    }

    tabidx[itabax] = -1;
    for (jtabax = 0; jtabax < i; jtabax++) {
      /* EXTNAME, EXTVER, EXTLEVEL, and TTYPEn for the coordinate array */
      /* must match for each axis of a multi-dimensional lookup table.  */
      if (strcmp(PSi_0a[itabax], PSi_0a[jtabax]) == 0 &&
          strcmp(PSi_1a[itabax], PSi_1a[jtabax]) == 0 &&
          PVi_1a[itabax] == PVi_1a[jtabax] &&
          PVi_2a[itabax] == PVi_2a[jtabax]) {
        tabidx[itabax] = tabidx[jtabax];
        break;
      }
    }

    if (jtabax == itabax) {
      tabidx[itabax] = wcs->ntab;
      wcs->ntab++;
    }
  }

  if (!(wcs->tab = calloc(wcs->ntab, sizeof(struct tabprm)))) {
    status = wcserr_set(WCSHDR_ERRMSG(WCSHDRERR_MEMORY));
    goto cleanup;
  }
  wcs->m_tab = wcs->tab;

  /* Table dimensionality; find the largest axis number. */
  for (itabax = 0; itabax < ntabax; itabax++) {
    tabp = wcs->tab + tabidx[itabax];

    /* PVi_3a records the 1-relative table axis number. */
    if (PVi_3a[itabax] > tabp->M) {
      tabp->M = PVi_3a[itabax];
    }
  }

  for (itab = 0; itab < wcs->ntab; itab++) {
    if ((status = tabini(1, wcs->tab[itab].M, 0, wcs->tab + itab))) {
      status = wcserr_set(WCSHDR_ERRMSG(wcshdr_taberr[status]));
      goto cleanup;
    }
  }


  /* Copy parameters into the tabprm structs. */
  for (i = 0; i < naxis; i++) {
    if ((itabax = tabax[i]) < 0) {
      /* Not a -TAB axis. */
      continue;
    }

    /* PVi_3a records the 1-relative table axis number. */
    m = PVi_3a[itabax] - 1;

    tabp = wcs->tab + tabidx[itabax];
    tabp->map[m] = i;
    tabp->crval[m] = wcs->crval[i];
  }

  /* Check for completeness. */
  for (itab = 0; itab < wcs->ntab; itab++) {
    for (m = 0; m < wcs->tab[itab].M; m++) {
      if (wcs->tab[itab].map[m] < 0) {
        status = wcserr_set(WCSERR_SET(WCSHDRERR_BAD_TABULAR_PARAMS),
          "Invalid tabular parameters: the axis mapping is undefined");
        goto cleanup;
      }
    }
  }


  /* Set up for reading the arrays; how many arrays are there? */
  for (itabax = 0; itabax < ntabax; itabax++) {
    /* Does this -TAB axis have a non-degenerate index array? */
    if (PSi_2a[itabax][0]) {
      wcs->nwtb++;
    }
  }

  /* Add one coordinate array for each table. */
  wcs->nwtb += wcs->ntab;

  /* Allocate memory for structs to be returned. */
  if (!(wcs->wtb = calloc(wcs->nwtb, sizeof(struct wtbarr)))) {
    wcs->nwtb = 0;

    status = wcserr_set(WCSHDR_ERRMSG(WCSHDRERR_MEMORY));
    goto cleanup;
  }
  wcs->m_wtb = wcs->wtb;

  /* Set pointers for the index and coordinate arrays. */
  wtbp = wcs->wtb;
  for (itab = 0; itab < wcs->ntab; itab++) {
    getcrd = 1;
    for (itabax = 0; itabax < ntabax; itabax++) {
      if (tabidx[itabax] != itab) continue;

      if (getcrd) {
        /* Coordinate array. */
        wtbp->i = itabax + 1;
        wtbp->m = PVi_3a[itabax];
        wtbp->kind = 'c';

        strcpy(wtbp->extnam, PSi_0a[itabax]);
        wtbp->extver = PVi_1a[itabax];
        wtbp->extlev = PVi_2a[itabax];
        strcpy(wtbp->ttype, PSi_1a[itabax]);
        wtbp->row    = 1L;
        wtbp->ndim   = wcs->tab[itab].M + 1;
        wtbp->dimlen = wcs->tab[itab].K;
        wtbp->arrayp = &(wcs->tab[itab].coord);

        /* Signal for tabset() to take this memory. */
        wcs->tab[itab].m_coord = (double *)0x1;

        wtbp++;
        getcrd = 0;
      }

      if (PSi_2a[itabax][0]) {
        /* Index array. */
        wtbp->i = itabax + 1;
        wtbp->m = PVi_3a[itabax];
        wtbp->kind = 'i';

        m = wtbp->m - 1;
        strcpy(wtbp->extnam, PSi_0a[itabax]);
        wtbp->extver = PVi_1a[itabax];
        wtbp->extlev = PVi_2a[itabax];
        strcpy(wtbp->ttype, PSi_2a[itabax]);
        wtbp->row    = 1L;
        wtbp->ndim   = 1;
        wtbp->dimlen = wcs->tab[itab].K + m;
        wtbp->arrayp = wcs->tab[itab].index + m;

        /* Signal for tabset() to take this memory. */
        wcs->tab[itab].m_indxs[m] = (double *)0x1;

        wtbp++;
      }
    }
  }

  status = 0;

cleanup:
  if (tabax)  free(tabax);
  if (tabidx) free(tabidx);
  if (PSi_0a) free(PSi_0a);
  if (PVi_1a) free(PVi_1a);
  if (PVi_2a) free(PVi_2a);
  if (PSi_1a) free(PSi_1a);
  if (PSi_2a) free(PSi_2a);
  if (PVi_3a) free(PVi_3a);

  if (status) {
    if (wcs->tab) free(wcs->tab);
    if (wcs->wtb) free(wcs->wtb);
  }

  return status;
}

/*--------------------------------------------------------------------------*/

int wcsidx(int nwcs, struct wcsprm **wcs, int alts[27])

{
  int a, iwcs;
  struct wcsprm *wcsp;

  for (a = 0; a < 27; a++) {
    alts[a] = -1;
  }

  if (wcs == 0x0) {
    return WCSHDRERR_NULL_POINTER;
  }

  wcsp = *wcs;
  for (iwcs = 0; iwcs < nwcs; iwcs++, wcsp++) {
    if (wcsp->colnum || wcsp->colax[0]) continue;

    if (wcsp->alt[0] == ' ') {
      a = 0;
    } else {
      a = wcsp->alt[0] - 'A' + 1;
    }

    alts[a] = iwcs;
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int wcsbdx(int nwcs, struct wcsprm **wcs, int type, short alts[1000][28])

{
  short  *ip;
  int    a, i, icol, iwcs;
  struct wcsprm *wcsp;

  for (ip = alts[0]; ip < alts[0] + 28*1000; ip++) {
    *ip = -1;
  }

  for (icol = 0; icol < 1000; icol++) {
    alts[icol][27] = 0;
  }

  if (wcs == 0x0) {
    return WCSHDRERR_NULL_POINTER;
  }

  wcsp = *wcs;
  for (iwcs = 0; iwcs < nwcs; iwcs++, wcsp++) {
    if (wcsp->alt[0] == ' ') {
      a = 0;
    } else {
      a = wcsp->alt[0] - 'A' + 1;
    }

    if (type) {
      /* Pixel list. */
      if (wcsp->colax[0]) {
        for (i = 0; i < wcsp->naxis; i++) {
          alts[wcsp->colax[i]][a]  = iwcs;
          alts[wcsp->colax[i]][27]++;
        }
      } else if (!wcsp->colnum) {
        alts[0][a]  = iwcs;
        alts[0][27]++;
      }

    } else {
      /* Binary table image array. */
      if (wcsp->colnum) {
        alts[wcsp->colnum][a] = iwcs;
        alts[wcsp->colnum][27]++;
      } else if (!wcsp->colax[0]) {
        alts[0][a]  = iwcs;
        alts[0][27]++;
      }
    }
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int wcsvfree(int *nwcs, struct wcsprm **wcs)

{
  int a, status = 0;
  struct wcsprm *wcsp;

  if (wcs == 0x0) {
    return WCSHDRERR_NULL_POINTER;
  }

  wcsp = *wcs;
  for (a = 0; a < *nwcs; a++, wcsp++) {
    status |= wcsfree(wcsp);
  }

  free(*wcs);

  *nwcs = 0;
  *wcs = 0x0;

  return status;
}

/*--------------------------------------------------------------------------*/
#define I_DTYPE   0	/* Distortion type code.                            */
#define I_NIPARM  1	/* Full (allocated) length of iparm[].              */
#define I_NDPARM  2	/* No. of parameters in dparm[], excl. work space.  */
#define I_DOCORR  3	/* True if distortion func computes a correction.   */
#define I_TPDNCO  4	/* No. of TPD coefficients, forward...              */
#define I_TPDINV  5     /* ...and inverse.                                  */
#define I_TPDAUX  6	/* True if auxiliary variables are used.            */
#define I_TPDRAD  7	/* True if the radial variable is used.             */

int wcshdo(int relax, struct wcsprm *wcs, int *nkeyrec, char **header)

/* ::: CUBEFACE and STOKES handling? */

{
  static const char *function = "wcshdo";

  const char axid[] = "xyxuvu", *cp;
  const int  nTPD[] = {1, 4, 7, 12, 17, 24, 31, 40, 49, 60};

  char alt, comment[72], ctemp[32], *ctypei, format[8], fmt01[8],
       keyvalue[72], keyword[16], *kp, *kp0, obsg[8] = "OBSG?",
       obsgeo[8] = "OBSGEO-?", pq, ptype, xtype, term[16], tpdsrc[24],
       xyz[] = "XYZ";
  int  *axmap, bintab, col0, *colax, colnum, degree, direct, doaux, dosip,
       dotpd, dotpv, i, idis, idp, *iparm, j, jhat, k, m, naxis, ncoeff, Nhat,
       p, pixlist, precision, primage, q, status = 0;
  double *dparm, keyval;
  struct disprm *dis;
  struct dpkey  *keyp;
  struct wcserr **err;

  *nkeyrec = 0;
  *header  = 0x0;

  if (wcs == 0x0) return WCSHDRERR_NULL_POINTER;
  err = &(wcs->err);

  if (wcs->flag != WCSSET) {
    if ((status = wcsset(wcs))) return status;
  }

  if ((naxis = wcs->naxis) == 0) {
    return 0;
  }


  /* These are mainly for convenience. */
  alt = wcs->alt[0];
  if (alt == ' ') alt = '\0';
  colnum = wcs->colnum;
  colax  = wcs->colax;

  primage = 0;
  bintab  = 0;
  pixlist = 0;
  if (colnum) {
    bintab  = 1;
    col0 = colnum;
  } else if (colax[0]) {
    pixlist = 1;
    col0 = colax[0];
  } else {
    primage = 1;
  }


  /* WCS dimension. */
  if (!pixlist) {
    sprintf(keyvalue, "%20d", naxis);
    wcshdo_util(relax, "WCSAXES", "WCAX", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "Number of coordinate axes", nkeyrec, header, &status);
  }

  /* Reference pixel coordinates. */
  wcshdo_format('G', naxis, wcs->crpix, format);
  for (j = 0; j < naxis; j++) {
    wcsutil_double2str(keyvalue, format, wcs->crpix[j]);
    wcshdo_util(relax, "CRPIX", "CRP", WCSHDO_CRPXna, "CRPX", 0, j+1, 0, alt,
      colnum, colax, keyvalue, "Pixel coordinate of reference point", nkeyrec,
      header, &status);
  }

  /* Linear transformation matrix. */
  wcshdo_format('G', naxis*naxis, wcs->pc, format);
  k = 0;
  for (i = 0; i < naxis; i++) {
    for (j = 0; j < naxis; j++, k++) {
      if (i == j) {
        if (wcs->pc[k] == 1.0) continue;
      } else {
        if (wcs->pc[k] == 0.0) continue;
      }

      wcsutil_double2str(keyvalue, format, wcs->pc[k]);
      wcshdo_util(relax, "PC", bintab ? "PC" : "P", WCSHDO_TPCn_ka,
        bintab ? 0x0 : "PC", i+1, j+1, 0, alt, colnum, colax,
        keyvalue, "Coordinate transformation matrix element",
        nkeyrec, header, &status);
    }
  }

  /* Coordinate increment at reference point. */
  wcshdo_format('G', naxis, wcs->cdelt, format);
  for (i = 0; i < naxis; i++) {
    wcsutil_double2str(keyvalue, format, wcs->cdelt[i]);
    comment[0] = '\0';
    if (wcs->cunit[i][0]) sprintf(comment, "[%s] ", wcs->cunit[i]);
    strcat(comment, "Coordinate increment at reference point");
    wcshdo_util(relax, "CDELT", "CDE", WCSHDO_CRPXna, "CDLT", i+1, 0, 0, alt,
      colnum, colax, keyvalue, comment, nkeyrec, header, &status);
  }

  /* Units of coordinate increment and reference value. */
  for (i = 0; i < naxis; i++) {
    if (wcs->cunit[i][0] == '\0') continue;

    sprintf(keyvalue, "'%s'", wcs->cunit[i]);
    wcshdo_util(relax, "CUNIT", "CUN", WCSHDO_CRPXna, "CUNI", i+1, 0, 0, alt,
      colnum, colax, keyvalue, "Units of coordinate increment and value",
      nkeyrec, header, &status);
  }

  /* May need to alter ctype for particular distortions so do basic checks */
  /* now.  Note that SIP, TPV, DSS, TNX, and ZPX are restricted to exactly */
  /* two axes and cannot coexist with other distortion types.              */
  dosip = 0;
  dotpv = 0;
  dotpd = 0;

  if ((dis = wcs->lin.dispre)) {
    for (i = 0; i < naxis; i++) {
      if (strcmp(dis->dtype[i], "SIP") == 0) {
        /* Simple Imaging Polynomial (SIP).  Write it in its native form  */
        /* if possible, unless specifically requested to write it as TPD. */
        dotpd = (dis->iparm[i][I_DTYPE] & DIS_DOTPD);

        if (!dotpd) {;
          if (alt ||
              dis->Nhat[0]      != 2 ||
              dis->Nhat[1]      != 2 ||
              dis->axmap[0][0]  != 0 ||
              dis->axmap[0][1]  != 1 ||
              dis->axmap[1][0]  != 0 ||
              dis->axmap[1][1]  != 1 ||
              dis->offset[0][0] != wcs->crpix[0] ||
              dis->offset[0][1] != wcs->crpix[1] ||
              dis->offset[1][0] != wcs->crpix[0] ||
              dis->offset[1][1] != wcs->crpix[1] ||
              dis->scale[0][0]  != 1.0 ||
              dis->scale[0][1]  != 1.0 ||
              dis->scale[1][0]  != 1.0 ||
              dis->scale[1][1]  != 1.0) {
            /* Must have been read as a 'SIP' distortion, CPDISja = 'SIP'. */
            /* Cannot be written as native SIP so write it as TPD.         */
            dotpd = DIS_DOTPD;
          } else if (strncmp(wcs->ctype[0], "RA---TAN", 8) ||
                     strncmp(wcs->ctype[1], "DEC--TAN", 8)) {
            /* Must have been permuted by wcssub(). */
            /* Native SIP doesn't have axis mapping so write it as TPD. */
            dotpd = DIS_DOTPD;
          }

          if (dotpd) {
            strcpy(tpdsrc, "SIP coordinates");
          } else {
            dosip = 1;
          }
        }

        break;
      }
    }
  }

  if ((dis = wcs->lin.disseq)) {
    for (i = 0; i < naxis; i++) {
      if (strcmp(dis->dtype[i], "TPV") == 0) {
        /* TPV "projection".  Write it in its native form if possible, */
        /* unless specifically requested to write it as TPD.           */
        dotpd = (dis->iparm[i][I_DTYPE] & DIS_DOTPD);

        if (!dotpd) {;
          if (dis->axmap[wcs->lng][0] != wcs->lng ||
              dis->axmap[wcs->lng][1] != wcs->lat ||
              dis->axmap[wcs->lat][0] != wcs->lat ||
              dis->axmap[wcs->lat][1] != wcs->lng ||
              dis->offset[wcs->lng][wcs->lng] != 0.0 ||
              dis->offset[wcs->lng][wcs->lat] != 0.0 ||
              dis->offset[wcs->lat][wcs->lng] != 0.0 ||
              dis->offset[wcs->lat][wcs->lat] != 0.0 ||
              dis->scale[wcs->lng][wcs->lng]  != 1.0 ||
              dis->scale[wcs->lng][wcs->lat]  != 1.0 ||
              dis->scale[wcs->lat][wcs->lng]  != 1.0 ||
              dis->scale[wcs->lat][wcs->lat]  != 1.0) {
            /* Must have been read as a 'TPV' distortion, CPDISja = 'TPV'. */
            /* Cannot be written as native TPV so write it as TPD.         */
            dotpd = DIS_DOTPD;
          }

          if (dotpd) {
            strcpy(tpdsrc, "TPV \"projection\"");
          } else {
            dotpv = 1;
          }
        }

        break;

      } else if (strcmp(dis->dtype[i], "DSS") == 0) {
        /* Always written as TPD. */
        dotpd = DIS_DOTPD;
        strcpy(tpdsrc, dis->dtype[i]);

      } else if (strncmp(dis->dtype[i], "WAT", 3) == 0) {
        /* Always written as TPD. */
        dotpd = DIS_DOTPD;
        strcpy(tpdsrc, dis->dtype[i]+4);

        if (strcmp(dis->dtype[i], "DSS") == 0) {
          strcpy(tpdsrc, wcs->wcsname);
        } else {
          strcat(tpdsrc, " \"projection\"");
        }

        break;
      }
    }
  }

  /* Coordinate type. */
  for (i = 0; i < naxis; i++) {
    if (wcs->ctype[i][0] == '\0') continue;

    sprintf(keyvalue, "'%s'", wcs->ctype[i]);
    strcpy(comment, "Coordinate type code");

    ctypei = keyvalue + 1;
    if (i == wcs->lng || i == wcs->lat) {
      /* Alter ctype for particular distortions. */
      if (dosip) {
        /* It could have come in as CPDISja = 'SIP'. */
        strcpy(ctypei+8, "-SIP'");
      } else if (dotpv) {
        /* Reinstate projection code edited by wcsset(). */
        strcpy(ctypei+4, "-TPV'");
      }

      if (strncmp(ctypei+8, "-SIP", 4) == 0) {
        strcpy(comment, "TAN (gnomonic) projection + SIP distortions");

      } else if (strncmp(ctypei+4, "-TPV", 4) == 0) {
        strcpy(comment, "TAN (gnomonic) projection + distortions");

      } else {
        if (strncmp(ctypei, "RA--", 4) == 0) {
          strcpy(comment, "Right ascension, ");

        } else if (strncmp(ctypei, "DEC-", 4) == 0) {
          strcpy(comment, "Declination, ");

        } else if (strncmp(ctypei+1, "LON", 3) == 0 ||
                   strncmp(ctypei+1, "LAT", 3) == 0) {
          ctypei[0] = toupper(ctypei[0]);

          switch (ctypei[0]) {
          case 'G':
            strcpy(comment, "galactic ");
            break;
          case 'E':
            strcpy(comment, "ecliptic ");
            break;
          case 'H':
            strcpy(comment, "helioecliptic ");
            break;
          case 'S':
            strcpy(comment, "supergalactic ");
            break;
          }

          if (i == wcs->lng) {
            strcat(comment, "longitude, ");
          } else {
            strcat(comment, "latitude, ");
          }
        }

        strcat(comment, wcs->cel.prj.name);
        strcat(comment, " projection");
      }

    } else if (i == wcs->spec) {
      spctyp(wcs->ctype[i], 0x0, 0x0, comment, 0x0, &ptype, &xtype, 0x0);
      if (ptype == xtype) {
        strcat(comment, " (linear)");
      } else {
        switch (xtype) {
        case 'F':
          strcat(comment, " (linear in frequency)");
          break;
        case 'V':
          strcat(comment, " (linear in velocity)");
          break;
        case 'W':
          strcat(comment, " (linear in wavelength)");
          break;
        }
      }
    }

    wcshdo_util(relax, "CTYPE", "CTY", WCSHDO_CRPXna, "CTYP", i+1, 0, 0, alt,
      colnum, colax, keyvalue, comment, nkeyrec, header, &status);
  }

  /* Coordinate value at reference point. */
  wcshdo_format('G', naxis, wcs->crval, format);
  for (i = 0; i < naxis; i++) {
    wcsutil_double2str(keyvalue, format, wcs->crval[i]);
    comment[0] = '\0';
    if (wcs->cunit[i][0]) sprintf(comment, "[%s] ", wcs->cunit[i]);
    strcat(comment, "Coordinate value at reference point");
    wcshdo_util(relax, "CRVAL", "CRV", WCSHDO_CRPXna, "CRVL", i+1, 0, 0, alt,
      colnum, colax, keyvalue, comment, nkeyrec, header, &status);
  }

  /* Parameter values. */
  for (k = 0; k < wcs->npv; k++) {
    wcsutil_double2str(keyvalue, "%20.12G", (wcs->pv[k]).value);
    if ((wcs->pv[k]).i == (wcs->lng + 1)) {
      switch ((wcs->pv[k]).m) {
      case 1:
        strcpy(comment, "[deg] Native longitude of the reference point");
        break;
      case 2:
        strcpy(comment, "[deg] Native latitude  of the reference point");
        break;
      case 3:
        if (primage) {
          sprintf(keyword, "LONPOLE%c", alt);
        } else if (bintab) {
          sprintf(keyword, "LONP%d%c", colnum, alt);
        } else {
          sprintf(keyword, "LONP%d%c", colax[(wcs->pv[k]).i - 1], alt);
        }
        sprintf(comment, "[deg] alias for %s (has precedence)", keyword);
        break;
      case 4:
        if (primage) {
          sprintf(keyword, "LATPOLE%c", alt);
        } else if (bintab) {
          sprintf(keyword, "LATP%d%c", colnum, alt);
        } else {
          sprintf(keyword, "LATP%d%c", colax[(wcs->pv[k]).i - 1], alt);
        }
        sprintf(comment, "[deg] alias for %s (has precedence)", keyword);
        break;
      }
    } else if ((wcs->pv[k]).i == (wcs->lat + 1)) {
      sprintf(comment, "%s projection parameter", wcs->cel.prj.code);
    } else {
      strcpy(comment, "Coordinate transformation parameter");
    }

    wcshdo_util(relax, "PV", "V", WCSHDO_PVn_ma, "PV", wcs->pv[k].i, -1,
      wcs->pv[k].m, alt, colnum, colax, keyvalue, comment,
      nkeyrec, header, &status);
  }

  for (k = 0; k < wcs->nps; k++) {
    sprintf(keyvalue, "'%s'", (wcs->ps[k]).value);
    wcshdo_util(relax, "PS", "S", WCSHDO_PVn_ma, "PS", wcs->ps[k].i, -1,
      wcs->ps[k].m, alt, colnum, colax, keyvalue,
      "Coordinate transformation parameter",
      nkeyrec, header, &status);
  }

  /* Celestial and spectral transformation parameters. */
  if (!undefined(wcs->lonpole)) {
    wcsutil_double2str(keyvalue, "%20.12G", wcs->lonpole);
    wcshdo_util(relax, "LONPOLE", "LONP", 0, 0x0, 0, 0, 0, alt,
      colnum, colax, keyvalue, "[deg] Native longitude of celestial pole",
      nkeyrec, header, &status);
  }

  if (!undefined(wcs->latpole)) {
    wcsutil_double2str(keyvalue, "%20.12G", wcs->latpole);
    wcshdo_util(relax, "LATPOLE", "LATP", 0, 0x0, 0, 0, 0, alt,
      colnum, colax, keyvalue, "[deg] Native latitude of celestial pole",
      nkeyrec, header, &status);
  }

  if (wcs->restfrq != 0.0) {
    wcsutil_double2str(keyvalue, "%20.12G", wcs->restfrq);
    wcshdo_util(relax, "RESTFRQ", "RFRQ", 0, 0x0, 0, 0, 0, alt,
      colnum, colax, keyvalue, "[Hz] Line rest frequency",
      nkeyrec, header, &status);
  }

  if (wcs->restwav != 0.0) {
    wcsutil_double2str(keyvalue, "%20.12G", wcs->restwav);
    wcshdo_util(relax, "RESTWAV", "RWAV", 0, 0x0, 0, 0, 0, alt,
      colnum, colax, keyvalue, "[Hz] Line rest wavelength",
      nkeyrec, header, &status);
  }

  /* Coordinate system title. */
  if (wcs->wcsname[0]) {
    sprintf(keyvalue, "'%s'", wcs->wcsname);
    if (bintab) {
      wcshdo_util(relax, "WCSNAME", "WCSN", 0, 0x0, 0, 0, 0, alt,
        colnum, colax, keyvalue, "Coordinate system title",
        nkeyrec, header, &status);
    } else {
      /* TWCS was a mistake. */
      wcshdo_util(relax, "WCSNAME", "TWCS", WCSHDO_WCSNna, "WCSN", 0, 0, 0,
        alt, colnum, colax, keyvalue, "Coordinate system title",
        nkeyrec, header, &status);
    }
  }

  /* Coordinate axis title. */
  if (wcs->cname) {
    for (i = 0; i < naxis; i++) {
      if (wcs->cname[i][0] == '\0') continue;

      sprintf(keyvalue, "'%s'", wcs->cname[i]);
      wcshdo_util(relax, "CNAME", "CNA", WCSHDO_CNAMna, "CNAM", i+1, 0, 0,
        alt, colnum, colax, keyvalue, "Axis name for labelling purposes",
        nkeyrec, header, &status);
    }
  }

  /* Random error in coordinate. */
  if (wcs->crder) {
    for (i = 0; i < naxis; i++) {
      if (undefined(wcs->crder[i])) continue;

      wcsutil_double2str(keyvalue, "%20.12G", wcs->crder[i]);
      comment[0] = '\0';
      if (wcs->cunit[i][0]) sprintf(comment, "[%s] ", wcs->cunit[i]);
      strcat(comment, "Random error in coordinate");
      wcshdo_util(relax, "CRDER", "CRD", WCSHDO_CNAMna, "CRDE", i+1, 0, 0,
        alt, colnum, colax, keyvalue, comment, nkeyrec, header, &status);
    }
  }

  /* Systematic error in coordinate. */
  if (wcs->csyer) {
    for (i = 0; i < naxis; i++) {
      if (undefined(wcs->csyer[i])) continue;

      wcsutil_double2str(keyvalue, "%20.12G", wcs->csyer[i]);
      comment[0] = '\0';
      if (wcs->cunit[i][0]) sprintf(comment, "[%s] ", wcs->cunit[i]);
      strcat(comment, "Systematic error in coordinate");
      wcshdo_util(relax, "CSYER", "CSY", WCSHDO_CNAMna, "CSYE", i+1, 0, 0,
        alt, colnum, colax, keyvalue, comment, nkeyrec, header, &status);
    }
  }

  /* Equatorial coordinate system type. */
  if (wcs->radesys[0]) {
    sprintf(keyvalue, "'%s'", wcs->radesys);
    wcshdo_util(relax, "RADESYS", "RADE", 0, 0x0, 0, 0, 0, alt,
      colnum, colax, keyvalue, "Equatorial coordinate system",
      nkeyrec, header, &status);
  }

  /* Equinox of equatorial coordinate system. */
  if (!undefined(wcs->equinox)) {
    wcsutil_double2str(keyvalue, "%20.12G", wcs->equinox);
    wcshdo_util(relax, "EQUINOX", "EQUI", 0, 0x0, 0, 0, 0, alt,
      colnum, colax, keyvalue, "[yr] Equinox of equatorial coordinates",
      nkeyrec, header, &status);
  }

  /* Reference frame of spectral coordinates. */
  if (wcs->specsys[0]) {
    sprintf(keyvalue, "'%s'", wcs->specsys);
    wcshdo_util(relax, "SPECSYS", "SPEC", 0, 0x0, 0, 0, 0, alt,
      colnum, colax, keyvalue, "Reference frame of spectral coordinates",
      nkeyrec, header, &status);
  }

  /* Reference frame of spectral observation. */
  if (wcs->ssysobs[0]) {
    sprintf(keyvalue, "'%s'", wcs->ssysobs);
    wcshdo_util(relax, "SSYSOBS", "SOBS", 0, 0x0, 0, 0, 0, alt,
      colnum, colax, keyvalue, "Reference frame of spectral observation",
      nkeyrec, header, &status);
  }

  /* Observer's velocity towards source. */
  if (!undefined(wcs->velosys)) {
    wcsutil_double2str(keyvalue, "%20.12G", wcs->velosys);
    wcshdo_util(relax, "VELOSYS", "VSYS", 0, 0x0, 0, 0, 0, alt,
      colnum, colax, keyvalue, "[m/s] Velocity towards source",
      nkeyrec, header, &status);
  }

  /* Reference frame of source redshift. */
  if (wcs->ssyssrc[0]) {
    sprintf(keyvalue, "'%s'", wcs->ssyssrc);
    wcshdo_util(relax, "SSYSSRC", "SSRC", 0, 0x0, 0, 0, 0, alt,
      colnum, colax, keyvalue, "Reference frame of source redshift",
      nkeyrec, header, &status);
  }

  /* Redshift of the source. */
  if (!undefined(wcs->zsource)) {
    wcsutil_double2str(keyvalue, "%20.12G", wcs->zsource);
    wcshdo_util(relax, "ZSOURCE", "ZSOU", 0, 0x0, 0, 0, 0, alt,
      colnum, colax, keyvalue, "Redshift of the source",
      nkeyrec, header, &status);
  }

  /* Observatory coordinates. */
  for (k = 0; k < 3; k++) {
    if (undefined(wcs->obsgeo[k])) continue;

    wcsutil_double2str(keyvalue, "%20.12G", wcs->obsgeo[k]);
    sprintf(comment, "[m] ITRF observatory %c-coordinate", xyz[k]);
    obsgeo[7] = xyz[k];
    obsg[4]   = xyz[k];
    wcshdo_util(relax, obsgeo, obsg, 0, 0x0, 0, 0, 0, ' ',
      colnum, colax, keyvalue, comment, nkeyrec, header, &status);
  }

  /* MJD of observation. */
  if (!undefined(wcs->mjdobs)) {
    wcsutil_double2str(keyvalue, "%20.12G", wcs->mjdobs);

    strcpy(comment, "[d] MJD of observation");
    if (wcs->dateobs[0]) {
      if (primage || (relax & 1) == 0) {
        sprintf(comment+22, " matching DATE-OBS");
      } else {
        sprintf(comment+22, " matching DOBS%d", col0);
      }
    }

    wcshdo_util(relax, "MJD-OBS", "MJDOB", 0, 0x0, 0, 0, 0, ' ',
      colnum, colax, keyvalue, comment, nkeyrec, header, &status);
  }

  /* MJD mid-observation time. */
  if (!undefined(wcs->mjdavg)) {
    wcsutil_double2str(keyvalue, "%20.12G", wcs->mjdavg);

    strcpy(comment, "[d] MJD mid-observation");
    if (wcs->dateavg[0]) {
      if (primage) {
        sprintf(comment+23, " matching DATE-AVG");
      } else {
        sprintf(comment+23, " matching DAVG%d", col0);
      }
    }

    wcshdo_util(relax, "MJD-AVG", "MJDA", 0, 0x0, 0, 0, 0, ' ',
      colnum, colax, keyvalue, comment, nkeyrec, header, &status);
  }

  /* ISO-8601 date corresponding to MJD-OBS. */
  if (wcs->dateobs[0]) {
    sprintf(keyvalue, "'%s'", wcs->dateobs);

    strcpy(comment, "ISO-8601 observation date");
    if (!undefined(wcs->mjdobs)) {
      if (primage) {
        sprintf(comment+25, " matching MJD-OBS");
      } else {
        sprintf(comment+25, " matching MJDOB%d", col0);
      }
    }

    if (relax & 1) {
      /* Allow DOBSn. */
      wcshdo_util(relax, "DATE-OBS", "DOBS", WCSHDO_DOBSn, 0x0, 0, 0, 0,
        ' ', colnum, colax, keyvalue, comment, nkeyrec, header, &status);
    } else {
      /* Force DATE-OBS. */
      wcshdo_util(relax, "DATE-OBS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0,
        0x0, keyvalue, comment, nkeyrec, header, &status);
    }
  }

  /* ISO-8601 date corresponding to MJD-OBS. */
  if (wcs->dateavg[0]) {
    sprintf(keyvalue, "'%s'", wcs->dateavg);

    strcpy(comment, "ISO-8601 mid-observation date");
    if (!undefined(wcs->mjdavg)) {
      if (primage) {
        sprintf(comment+29, " matching MJD-AVG");
      } else {
        sprintf(comment+29, " matching MJDA%d", col0);
      }
    }

    wcshdo_util(relax, "DATE-AVG", "DAVG", 0, 0x0, 0, 0, 0, ' ',
      colnum, colax, keyvalue, comment, nkeyrec, header, &status);
  }


  /* Write distortion function parameters. */
  if (dosip) {
    /* Simple Imaging Polynomial (SIP) is handled by translating its dpkey */
    /* records.  Determine a suitable numerical precision for the          */
    /* polynomial coefficients to avoid trailing zeroes common to all of   */
    /* them.                                                               */
    dis  = wcs->lin.dispre;
    keyp = dis->dp;
    kp0  = keyvalue + 2;
    for (idp = 0; idp < dis->ndp; idp++, keyp++) {
      cp = strchr(keyp->field, '.') + 1;
      if (strncmp(cp, "SIP.", 4) != 0) continue;
      wcsutil_double2str(keyvalue, "%20.13E", wcsutil_dpkey_double(keyp));

      kp = keyvalue + 15;
      while (kp0 < kp && *kp == '0') kp--;
      kp0 = kp;
    }

    precision = kp - (keyvalue + 2);
    if (precision < 1)  precision = 1;
    if (13 < precision) precision = 13;
    sprintf(format, "%%20.%dE", precision);

    /* Ensure the coefficients are written in a human-readable sequence. */
    for (j = 0; j <= 1; j++) {
      /* Distortion function polynomial coefficients. */
      wcshdo_util(relax, "", "", 0, 0x0, 0, 0, 0, ' ', 0, 0, "", "",
        nkeyrec, header, &status);

      if (j == 0) {
        strcpy(keyword, "A_");
      } else {
        strcpy(keyword, "B_");
      }

      ncoeff = dis->iparm[j][I_TPDNCO];
      for (degree = 0; degree <= 9; degree++) {
        if (ncoeff <= nTPD[degree]) break;
      }

      strcpy(keyword+2, "ORDER");
      sprintf(keyvalue, "%20d", degree);
      sprintf(comment, "SIP polynomial degree, axis %d, pixel-to-sky", j+1);
      wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, ' ', 0, 0,
        keyvalue, comment, nkeyrec, header, &status);

      keyp = dis->dp;
      for (idp = 0; idp < dis->ndp; idp++, keyp++) {
        if (keyp->j != j+1) continue;
        if ((keyval = wcsutil_dpkey_double(keyp)) == 0.0) continue;

        cp = strchr(keyp->field, '.') + 1;
        if (strncmp(cp, "SIP.FWD.", 8) != 0) continue;
        cp += 8;
        strcpy(keyword+2, cp);
        sscanf(cp, "%d_%d", &p, &q);
        strncpy(term, "xxxxxxxxx", p);
        strncpy(term+p, "yyyyyyyyy", q);
        term[p+q] = '\0';

        wcsutil_double2str(keyvalue, format, keyval);
        sprintf(comment, "SIP distortion coefficient: %s", term);
        wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, ' ', 0, 0,
          keyvalue, comment, nkeyrec, header, &status);
      }

      if (dis->maxdis[j] != 0.0) {
        strcpy(keyword+2, "DMAX");
        wcsutil_double2str(keyvalue, "%20.3f", dis->maxdis[j]);
        wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, ' ', 0, 0,
          keyvalue, "Maximum value of distortion function", nkeyrec,
          header, &status);
      }

      /* Inverse distortion function polynomial coefficients. */
      if (dis->disx2p == 0x0) continue;

      wcshdo_util(relax, "", "", 0, 0x0, 0, 0, 0, ' ', 0, 0, "", "",
        nkeyrec, header, &status);

      if (j == 0) {
        strcpy(keyword, "AP_");
      } else {
        strcpy(keyword, "BP_");
      }

      ncoeff = dis->iparm[j][I_NDPARM] - dis->iparm[j][I_TPDNCO];
      for (degree = 0; degree <= 9; degree++) {
        if (ncoeff <= nTPD[degree]) break;
      }

      strcpy(keyword+3, "ORDER");
      sprintf(keyvalue, "%20d", degree);
      sprintf(comment, "SIP polynomial degree, axis %d, sky-to-pixel", j+1);
      wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, ' ', 0, 0,
        keyvalue, comment, nkeyrec, header, &status);

      keyp = dis->dp;
      for (idp = 0; idp < dis->ndp; idp++, keyp++) {
        if (keyp->j != j+1) continue;
        if ((keyval = wcsutil_dpkey_double(keyp)) == 0.0) continue;

        cp = strchr(keyp->field, '.') + 1;
        if (strncmp(cp, "SIP.REV.", 8) != 0) continue;
        cp += 8;
        strcpy(keyword+3, cp);
        sscanf(cp, "%d_%d", &p, &q);
        strncpy(term, "xxxxxxxxx", p);
        strncpy(term+p, "yyyyyyyyy", q);
        term[p+q] = '\0';

        wcsutil_double2str(keyvalue, format, keyval);
        sprintf(comment, "SIP inverse coefficient: %s", term);
        wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, ' ', 0, 0,
          keyvalue, comment, nkeyrec, header, &status);
      }
    }
  }

  for (idis = 0; idis < 2; idis++) {
    if (idis == 0 && (dis = wcs->lin.dispre) == 0x0) continue;
    if (idis == 1 && (dis = wcs->lin.disseq) == 0x0) continue;

    for (j = 0; j < naxis; j++) {
      if (dis->disp2x[j] == 0x0) continue;

      iparm = dis->iparm[j];
      dparm = dis->dparm[j];

      /* Identify the distortion type. */
      if (dotpv) {
        /* TPV "projection" is handled by translating its dpkey records, */
        /* which were originally translated from PVi_ma by wcsset(), or  */
        /* possibly input directly as a CQDISia = 'TPV' distortion type. */
        /* Determine a suitable numerical precision for the polynomial   */
        /* coefficients to avoid trailing zeroes common to all of them.  */
        wcshdo_format('E', iparm[I_NDPARM], dparm, format);
        sprintf(fmt01, "%.3ss", format);

        wcshdo_util(relax, "", "", 0, 0x0, 0, 0, 0, ' ', 0, 0, "", "",
          nkeyrec, header, &status);

        /* Distortion function polynomial coefficients. */
        sprintf(keyword, "PV%d_", j+1);
        kp = keyword + strlen(keyword);

        keyp = dis->dp;
        for (idp = 0; idp < dis->ndp; idp++, keyp++) {
          if (keyp->j != j+1) continue;
          if ((keyval = wcsutil_dpkey_double(keyp)) == 0.0) continue;

          cp = strchr(keyp->field, '.') + 1;
          if (strncmp(cp, "TPV.", 4) != 0) continue;
          strcpy(kp, cp+4);

          /* Identify the term of the TPV polynomial for human readers. */
          sscanf(cp+4, "%d", &m);
          wcshdo_tpdterm(m, j == wcs->lng, term);
          sprintf(comment, "TPV coefficient: %s", term);

          if (keyval == 1.0) {
            sprintf(keyvalue, fmt01, "1.0");
          } else {
            wcsutil_double2str(keyvalue, format, keyval);
          }
          wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
            keyvalue, comment, nkeyrec, header, &status);
        }

      } else if (strcmp(dis->dtype[j], "TPD") == 0 || dotpd ||
                 strcmp(dis->dtype[j], "Polynomial")  == 0 ||
                 strcmp(dis->dtype[j], "Polynomial*") == 0) {
        /* One of the Paper IV type polynomial distortions. */
        wcshdo_util(relax, "", "", 0, 0x0, 0, 0, 0, ' ', 0, 0, "", "",
          nkeyrec, header, &status);

        if (strcmp(dis->dtype[j], "TPD") == 0) {
          /* Pure TPD. */
          dotpd = 1;
        } else if (strncmp(dis->dtype[j], "Polynomial", 10) == 0) {
          /* Polynomial distortion.  Write it as TPD by request? */
          dotpd = (iparm[I_DTYPE] & DIS_DOTPD);
          strcpy(tpdsrc, "Polynomial distortion");
        }

        pq = idis ? 'Q' : 'P';
        Nhat = dis->Nhat[j];

        /* CPDISja/CQDISia */
        sprintf(keyword, "C%cDIS%d", pq, j+1);
        if (idis == 0) {
          strcpy(comment, "P = prior, ");
        } else {
          strcpy(comment, "Q = sequent, ");
        }

        if (dotpd) {
          strcpy(keyvalue, "'TPD'");
          strcat(comment, "Template Polynomial Distortion");

          /* For identifying terms of the TPD polynomial. */
          axmap  = dis->axmap[j];
          direct = 1;
          doaux  = iparm[I_TPDAUX];
          if (Nhat == 2) {
            /* Associate x with longitude, y with latitude. */
            if (axmap[0] == wcs->lng && axmap[1] == wcs->lat) {
              direct = 1;
            } else if (axmap[0] == wcs->lat && axmap[1] == wcs->lng) {
              direct = 0;
            } else {
              /* Non-celestial. */
              direct = (axmap[0] < axmap[1]);
            }
          }
        } else {
          strcpy(keyvalue, "'Polynomial'");
          strcat(comment, "general Polynomial distortion");
        }

        wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
          keyvalue, comment, nkeyrec, header, &status);

        /* NAXES. */
        sprintf(keyword,  "D%c%d", pq, j+1);
        sprintf(keyvalue, "'NAXES:  %d'", Nhat);
        if (Nhat == 1) {
          strcpy(comment,  "One independent variable");
        } else if (Nhat == 2) {
          strcpy(comment,  "Two independent variables");
        } else {
          strcpy(comment,  "Number of independent variables");
        }
        wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
          keyvalue, comment, nkeyrec, header, &status);

        /* AXIS.jhat */
        for (jhat = 0; jhat < Nhat; jhat++) {
          axmap = dis->axmap[j];
          sprintf(keyvalue, "'AXIS.%d: %d'", jhat+1, axmap[jhat]+1);
          if (jhat == 0) {
            strcpy(comment, "1st");
          } else if (jhat == 1) {
            strcpy(comment, "2nd");
          } else if (jhat == 2) {
            strcpy(comment, "3rd");
          } else {
            sprintf(comment, "%dth", jhat+1);
          }

          sprintf(comment+strlen(comment), " independent variable: axis %d",
            axmap[jhat]+1);
          if (dotpd) {
            /* axid is "xyxuvu". */
            cp = axid;
            if (!direct) cp++;
            if (doaux) cp += 3;
            sprintf(comment+strlen(comment), " (= %c)", cp[jhat]);
          }

          wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
            keyvalue, comment, nkeyrec, header, &status);
        }

        /* OFFSET.jhat */
        wcshdo_format('f', Nhat, dis->offset[j], format);
        for (jhat = 0; jhat < Nhat; jhat++) {
          if (dis->offset[j][jhat] == 0.0) continue;

          wcsutil_double2str(ctemp, format, dis->offset[j][jhat]);
          sprintf(keyvalue, "'OFFSET.%d: %s'", jhat+1, ctemp);
          sprintf(comment, "Variable %d renormalization offset", jhat+1);

          wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
            keyvalue, comment, nkeyrec, header, &status);
        }

        /* SCALE.jhat */
        wcshdo_format('f', Nhat, dis->scale[j], format);
        for (jhat = 0; jhat < Nhat; jhat++) {
          if (dis->scale[j][jhat] == 1.0) continue;

          wcsutil_double2str(ctemp, format, dis->scale[j][jhat]);
          sprintf(keyvalue, "'SCALE.%d: %s'", jhat+1, ctemp);
          sprintf(comment, "Variable %d renormalization scale", jhat+1);

          wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
            keyvalue, comment, nkeyrec, header, &status);
        }

        /* Does the distortion function compute a correction? */
        if (iparm[I_DOCORR]) {
          wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
            "'DOCORR: 1'", "Distortion function computes a correction",
            nkeyrec, header, &status);
        }

        if (dotpd) {
          /* Template Polynomial Distortion (TPD).  As it may have been */
          /* translated from SIP, TPV, DSS, TNX, ZPX, or perhaps        */
          /* Polynomial, the dpkey records may not relate to TPD.       */
          /* Output is therefore handled via dparm.                     */
          wcshdo_format('E', iparm[I_NDPARM], dparm, format);
          sprintf(fmt01, "%.3ss", format);

          /* AUX.jhat.COEFF.m */
          if (doaux) {
            for (idp = 0; idp < 6; idp++) {
              if (dparm[idp] == 0.0) {
                sprintf(ctemp, fmt01, "0.0");
              } else if (dparm[idp] == 1.0) {
                sprintf(ctemp, fmt01, "1.0");
              } else {
                wcsutil_double2str(ctemp, format, dparm[idp]);
              }

              if (idp < 3) {
                sprintf(keyvalue, "'AUX.1.COEFF.%d: %s'", idp%3, ctemp);
                sprintf(comment, "TPD: x = c0 + c1*u + c2*v");
              } else {
                sprintf(keyvalue, "'AUX.2.COEFF.%d: %s'", idp%3, ctemp);
                sprintf(comment, "TPD: y = d0 + d1*u + d2*v");
              }

              wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
                keyvalue, comment, nkeyrec, header, &status);

            }

            dparm += 6;
          }

          /* TPD.FWD.m */
          for (idp = 0; idp < iparm[I_TPDNCO]; idp++) {
            if (dparm[idp] == 0.0) continue;

            if (dparm[idp] == 1.0) {
              sprintf(ctemp, fmt01, "1.0");
            } else {
              wcsutil_double2str(ctemp, format, dparm[idp]);
            }

            m = idp;
            sprintf(keyvalue, "'TPD.FWD.%d:%s %s'", m, (m<10)?" ":"", ctemp);
            wcshdo_tpdterm(m, direct, term);
            sprintf(comment, "TPD coefficient: %s", term);

            wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
              keyvalue, comment, nkeyrec, header, &status);
          }

          /* CPERRja/CQERRia */
          if (dis->maxdis[j] != 0.0) {
            sprintf(keyword,  "C%cERR%d", pq, j+1);
            sprintf(keyvalue, "%20.2f", dis->maxdis[j]);
            sprintf(comment, "%sMaximum absolute value of distortion",
              idis?"":"[pix] ");
            wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
              keyvalue, comment, nkeyrec, header, &status);
          }

          /* Inverse distortion function polynomial coefficients. */
          if (dis->disx2p[j] == 0x0) continue;

          wcshdo_util(relax, "", "", 0, 0x0, 0, 0, 0, ' ', 0, 0, "", "",
            nkeyrec, header, &status);

          /* TPD.REV.m */
          sprintf(keyword,  "D%c%d", pq, j+1);
          for (idp = iparm[I_TPDNCO]; idp < iparm[I_NDPARM]; idp++) {
            if (dparm[idp] == 0.0) continue;

            wcsutil_double2str(ctemp, format, dparm[idp]);
            m = idp - iparm[I_TPDNCO];
            sprintf(keyvalue, "'TPD.REV.%d:%s %s'", m, (m<10)?" ":"", ctemp);
            wcshdo_tpdterm(m, direct, term);
            sprintf(comment, "TPD coefficient: %s", term);

            wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
              keyvalue, comment, nkeyrec, header, &status);
          }

        } else {
          /* General polynomial distortion, handled via its dpkey records */
          /* since iparm and dparm may hold a translation to TPD.         */

          /* Do auxiliary variables first. */
          keyp = dis->dp;
          for (idp = 0; idp < dis->ndp; idp++, keyp++) {
            if (keyp->j != j+1) continue;

            cp = strchr(keyp->field, '.') + 1;
            if (strncmp(cp, "NAUX", 4) != 0) continue;

            sprintf(keyvalue, "'%s: %d'", cp, wcsutil_dpkey_int(keyp));
            wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
              keyvalue, "Number of auxiliary variables", nkeyrec, header,
              &status);

            keyp = dis->dp;
            for (idp = 0; idp < dis->ndp; idp++, keyp++) {
              if (keyp->j != j+1) continue;

              keyval = wcsutil_dpkey_double(keyp);

              cp = strchr(keyp->field, '.') + 1;
              if (strncmp(cp, "AUX.", 4) != 0) continue;

              sscanf(cp+4, "%d", &m);
              sprintf(keyvalue, "'%s:", cp);

              cp = strchr(cp+4, '.') + 1;
              kp = keyvalue + strlen(keyvalue);

              if ((double)((int)keyval) == keyval) {
                sprintf(kp, "%4d'", (int)keyval);
              } else if (keyval == 0.5) {
                strcat(kp, " 0.5'");
              } else {
                wcsutil_double2str(kp, "%21.13E", keyval);
                strcat(keyvalue, "'");
              }

              sscanf(cp+6, "%d", &p);
              if (strncmp(cp, "POWER.", 4) == 0) {
                if (p) {
                  sprintf(comment, "Aux %d: var %d power", m, p);
                } else {
                  sprintf(comment, "Aux %d: power of sum of terms", m);
                }
              } else {
                if (p) {
                  sprintf(comment, "Aux %d: var %d coefficient", m, p);
                } else {
                  sprintf(comment, "Aux %d: offset term", m);
                }
              }

              wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
                keyvalue, comment, nkeyrec, header, &status);
            }

            break;
          }

          /* Do polynomial terms. */
          keyp = dis->dp;
          for (idp = 0; idp < dis->ndp; idp++, keyp++) {
            if (keyp->j != j+1) continue;

            cp = strchr(keyp->field, '.') + 1;
            if (strncmp(cp, "NTERMS", 6) != 0) continue;

            sprintf(keyvalue, "'%s: %d'", cp, wcsutil_dpkey_int(keyp));
            wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
              keyvalue, "Number of terms in the polynomial", nkeyrec, header,
              &status);
          }

          keyp = dis->dp;
          for (idp = 0; idp < dis->ndp; idp++, keyp++) {
            if (keyp->j != j+1) continue;

            if ((keyval = wcsutil_dpkey_double(keyp)) == 0.0) continue;

            cp = strchr(keyp->field, '.') + 1;
            if (strncmp(cp, "TERM.", 5) != 0) continue;

            sscanf(cp+5, "%d", &m);
            sprintf(keyvalue, "'%s:%s ", cp, (m<10)?" ":"");

            cp = strchr(cp+5, '.') + 1;
            kp = keyvalue + strlen(keyvalue);
            if (strncmp(cp, "VAR.", 4) == 0) {
              if ((double)((int)keyval) == keyval) {
                sprintf(kp, "%20d", (int)keyval);
              } else {
                wcsutil_double2str(kp, "%20.13f", keyval);
              }

              sscanf(cp+4, "%d", &p);
              if (p <= Nhat) {
                sprintf(comment, "Poly term %d: var %d power", m, p);
              } else {
                sprintf(comment, "Poly term %d: aux %d power", m, p-Nhat);
              }

            } else {
              wcsutil_double2str(kp, "%20.13E", keyval);
              sprintf(comment, "Poly term %d: coefficient", m);
            }
            strcat(keyvalue, "'");

            wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
              keyvalue, comment, nkeyrec, header, &status);
          }

          /* CPERRja/CQERRia */
          if (dis->maxdis[j] != 0.0) {
            sprintf(keyword,  "C%cERR%d", pq, j+1);
            sprintf(keyvalue, "%20.2f", dis->maxdis[j]);
            sprintf(comment, "%sMaximum absolute value of distortion",
              idis?"":"[pix] ");
            wcshdo_util(relax, keyword, "", 0, 0x0, 0, 0, 0, alt, 0, 0,
              keyvalue, comment, nkeyrec, header, &status);
          }
        }
      }
    }

    /* DVERRa */
    if (dis->totdis != 0.0) {
      sprintf(keyvalue, "%20.2f", dis->totdis);
      sprintf(comment, "Maximum combined distortion");
      wcshdo_util(relax, "DVERR", "", 0, 0x0, 0, 0, 0, alt, 0, 0,
        keyvalue, comment, nkeyrec, header, &status);
    }
  }


  /* Add identification. */
  wcshdo_util(relax, "", "", 0, 0x0, 0, 0, 0, ' ', 0, 0, "", "",
    nkeyrec, header, &status);

  if (dotpd == DIS_DOTPD) {
    /* TPD by translation. */
    sprintf(comment, "Translated from %s to TPD by WCSLIB %s", tpdsrc,
      wcslib_version(0x0));
  } else {
    sprintf(comment, "WCS header keyrecords produced by WCSLIB %s",
      wcslib_version(0x0));
  }

  wcshdo_util(relax, "COMMENT", "", 0, 0x0, 0, 0, 0, ' ', 0, 0,
    "", comment, nkeyrec, header, &status);


  if (status == WCSHDRERR_MEMORY) {
    wcserr_set(WCSHDR_ERRMSG(status));
  }
  return status;
}

/*--------------------------------------------------------------------------*/
/* Determine a suitable floating point format for a set of parameters.      */

void wcshdo_format(
  int fmt,
  int nval,
  const double val[],
  char *format)

{
  char *cp, *cp0, cval[24];
  int  i, expmax, expon, nsig, precision;

  if (fmt == 'G') {
    fmt = 'f';
    for (i = 0; i < nval; i++) {
      if (fabs(val[i]) < 1e-4 || 1e12 < val[i]) {
        fmt = 'E';
        break;
      }
    }
  }

  cp0 = cval + 2;
  expmax = -999;
  for (i = 0; i < nval; i++) {
    /* Double precision has at least 15 significant digits, and up to 17:  */
    /* http://en.wikipedia.org/wiki/Double-precision_floating-point_format */
    wcsutil_double2str(cval, "%21.14E", val[i]);

    cp = cval + 16;
    while (cp0 < cp && *cp == '0') cp--;
    cp0 = cp;

    sscanf(cval+18, "%d", &expon);
    if (expmax < expon) expmax = expon;
  }

  nsig = cp0 - (cval + 2) + 1;


  if (fmt == 'f') {
    precision = nsig - (expmax + 1);
    if (precision < 1)  precision = 1;
    if (17 < precision) precision = 17;
    sprintf(format, "%%20.%df", precision);

  } else {
    precision = nsig - 1;
    if (precision < 1)  precision = 1;
    if (14 < precision) precision = 14;
    if (precision < 14) {
      sprintf(format, "%%20.%dE", precision);
    } else {
      sprintf(format, "%%21.%dE", precision);
    }
  }
}

/*--------------------------------------------------------------------------*/
/* Construct a string that identifies the term of a TPD or TPV polynomial.  */

void wcshdo_tpdterm(
  int m,
  int direct,
  char *term)

{
  const int nTPD[] = {1, 4, 7, 12, 17, 24, 31, 40, 49, 60};

  int degree, k;

  for (degree = 0; degree <= 9; degree++) {
    if (m < nTPD[degree]) break;
  }

  if (degree == 0) {
    strcpy(term, "1");

  } else {
    k = degree - (m - nTPD[degree-1]);

    if (k < 0) {
      strncpy(term, "rrrrrrrrr", degree);
    } else if (direct) {
      strncpy(term, "xxxxxxxxx", k);
      strncpy(term+k, "yyyyyyyyy", degree-k);
    } else {
      strncpy(term, "yyyyyyyyy", k);
      strncpy(term+k, "xxxxxxxxx", degree-k);
    }

    term[degree] = '\0';
  }
}

/*--------------------------------------------------------------------------*/
/* Construct a keyrecord from the components given.                         */

void wcshdo_util(
  int relax,
  const char pikey[],
  const char tbkey[],
  int level,
  const char tlkey[],
  int i,
  int j,
  int m,
  char alt,
  int  btcol,
  int  plcol[],
  char keyvalue[],
  const char keycomment[],
  int  *nkeyrec,
  char **header,
  int  *status)

{
  char ch0, ch1, *hptr, keyword[16], *kptr;
  int  nbyte, nc = 47, nv;

  if (*status) return;

  /* Reallocate memory in blocks of 2880 bytes. */
  if ((*nkeyrec)%32 == 0) {
    nbyte = ((*nkeyrec)/32 + 1) * 2880;
    if (!(hptr = realloc(*header, nbyte))) {
      *status = WCSHDRERR_MEMORY;
      return;
    }

    *header = hptr;
  }

  /* Construct the keyword. */
  if (alt == ' ') alt = '\0';
  if (btcol) {
    /* Binary table image array. */
    if (i > 0 && j) {
      if (j > 0) {
        sprintf(keyword, "%d%d%s%d%c", i, j, tbkey, btcol, alt);
      } else {
        sprintf(keyword, "%d%s%d_%d%c", i, tbkey, btcol, m, alt);
      }
    } else if (i > 0) {
      sprintf(keyword, "%d%s%d%c", i, tbkey, btcol, alt);
    } else if (j > 0) {
      sprintf(keyword, "%d%s%d%c", j, tbkey, btcol, alt);
    } else {
      sprintf(keyword, "%s%d%c", tbkey, btcol, alt);
    }

    if ((strlen(keyword) < 8) && tlkey && (relax & level)) {
      /* Use the long form. */
      if (i > 0 && j) {
        if (j > 0) {
          sprintf(keyword, "%d%d%s%d%c", i, j, tlkey, btcol, alt);
        } else {
          sprintf(keyword, "%d%s%d_%d%c", i, tlkey, btcol, m, alt);
        }
      } else if (i > 0) {
        sprintf(keyword, "%d%s%d%c", i, tlkey, btcol, alt);
      } else if (j > 0) {
        sprintf(keyword, "%d%s%d%c", j, tlkey, btcol, alt);
      } else {
        sprintf(keyword, "%s%d%c", tlkey, btcol, alt);
      }
    }

  } else if (plcol && plcol[0]) {
    /* Pixel list. */
    if (i > 0 && j) {
      if (j > 0) {
        sprintf(keyword, "T%s%d_%d%c", tbkey, plcol[i-1], plcol[j-1], alt);
      } else {
        sprintf(keyword, "T%s%d_%d%c", tbkey, plcol[i-1], m, alt);
      }
    } else if (i > 0) {
      sprintf(keyword, "T%s%d%c", tbkey, plcol[i-1], alt);
    } else if (j > 0) {
      sprintf(keyword, "T%s%d%c", tbkey, plcol[j-1], alt);
    } else {
      sprintf(keyword, "%s%d%c", tbkey, plcol[0], alt);
    }

    if ((strlen(keyword) < 8) && tlkey && (relax & level)) {
      /* Use the long form. */
      if (i > 0 && j) {
        if (j > 0) {
          sprintf(keyword, "T%s%d_%d%c", tlkey, plcol[i-1], plcol[j-1], alt);
        } else {
          sprintf(keyword, "T%s%d_%d%c", tlkey, plcol[i-1], m, alt);
        }
      } else if (i > 0) {
        sprintf(keyword, "T%s%d%c", tlkey, plcol[i-1], alt);
      } else if (j > 0) {
        sprintf(keyword, "T%s%d%c", tlkey, plcol[j-1], alt);
      } else {
        sprintf(keyword, "%s%d%c", tlkey, plcol[0], alt);
      }
    }
  } else {
    if (i > 0 && j) {
      if (j > 0) {
        sprintf(keyword, "%s%d_%d%c", pikey, i, j, alt);
      } else {
        sprintf(keyword, "%s%d_%d%c", pikey, i, m, alt);
      }
    } else if (i > 0) {
      sprintf(keyword, "%s%d%c", pikey, i, alt);
    } else if (j > 0) {
      sprintf(keyword, "%s%d%c", pikey, j, alt);
    } else {
      sprintf(keyword, "%s%c", pikey, alt);
    }
  }

  /* Double-up single-quotes in string keyvalues. */
  if (*keyvalue == '\'') {
    hptr = keyvalue + 1;
    while (*hptr) {
      if (*hptr == '\'') {
        kptr = hptr++;
        if (*hptr) {
          ch0 = *kptr;
          while (*kptr) {
            ch1 = *(++kptr);
            *kptr = ch0;
            ch0 = ch1;
          }
        } else {
          break;
        }
      }

      hptr++;
    }

    /* Check length. */
    if (strlen(keyvalue) > 70) {
      /* Truncate. */
      keyvalue[69] = '\'';
      keyvalue[70] = '\0';
    }

  } else {
    /* Check length. */
    if (strlen(keyvalue) > 70) {
      /* Truncate. */
      keyvalue[70] = '\0';
    }
  }

  if ((nv = strlen(keyvalue) > 20)) {
    /* Rob the keycomment to make space for the keyvalue. */
    nc -= (nv - 20);
  }

  hptr = *header + (80 * ((*nkeyrec)++));
  if (*keyword == '\0') {
    sprintf(hptr, "%80.80s", " ");
  } else if (strcmp(keyword, "COMMENT") == 0) {
    sprintf(hptr, "%-8.8s %-71.71s", keyword, keycomment);
  } else {
    sprintf(hptr, "%-8.8s= %-20s / %-*.*s", keyword, keyvalue, nc, nc,
      keycomment);
  }
}
