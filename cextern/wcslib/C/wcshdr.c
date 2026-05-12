/*============================================================================
  WCSLIB 8.7 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2026, Mark Calabretta

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
  http://www.atnf.csiro.au/computing/software/wcs
  $Id: wcshdr.c,v 8.7 2026/05/11 12:01:10 mcalabre Exp $
*===========================================================================*/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "wcserr.h"
#include "wcsmath.h"
#include "wcsutil.h"
#include "wcshdr.h"
#include "wtbarr.h"
#include "tab.h"
#include "dis.h"
#include "wcs.h"

// Map status return value to message.
const char *wcshdr_errmsg[] = {
  "Success",
  "Null wcsprm pointer passed",
  "Memory allocation failed",
  "Invalid column selection",
  "Fatal error returned by Flex parser",
  "Invalid tabular parameters"};

// Map error returns for lower-level routines.
const int wcshdr_taberr[] = {
  WCSHDRERR_SUCCESS,		//  0: TABERR_SUCCESS
  WCSHDRERR_NULL_POINTER,	//  1: TABERR_NULL_POINTER
  WCSHDRERR_MEMORY,		//  2: TABERR_MEMORY
  WCSHDRERR_BAD_TABULAR_PARAMS	//  3: TABERR_BAD_PARAMS
				//  4: TABERR_BAD_X
				//  5: TABERR_BAD_WORLD
};

static const int WCSSET    =  137;	// Matching wcs.c

static const int DIS_DOTPD = 1024;	// Matching dis.c

// Internal helper functions, not for general use.
static void wcshdo_format(int, int, const double [], char [16]);
static void wcshdo_tpdterm(int, int, char [16]);
static void wcshdo_util(int, const char [], const char [], int, const char [],
  int, int, int, char, int, int [], char [], const char [], int *, char **,
  int *);

// Convenience macro for invoking wcserr_set().
#define WCSHDR_ERRMSG(status) WCSERR_SET(status), wcshdr_errmsg[status]

//----------------------------------------------------------------------------

int wcstab(struct wcsprm *wcs)

{
  static const char *function = "wcstab";

  // Pointers to allocated memory.
  char (*PSi_0a)[72] = 0x0, (*PSi_1a)[72] = 0x0, (*PSi_2a)[72] = 0x0;
  int  *PVi_1a = 0x0, *PVi_2a = 0x0, *PVi_3a = 0x0;
  int  *tabidx = 0x0;

  int status = 0;

  if (wcs == 0x0) return WCSHDRERR_NULL_POINTER;
  struct wcserr **err = &(wcs->err);

  // Free memory previously allocated by wcstab().
  if (wcs->flag != -1 && wcs->m_flag == WCSSET) {
    if (wcs->wtb == wcs->m_wtb) wcs->wtb = 0x0;
    if (wcs->tab == wcs->m_tab) wcs->tab = 0x0;

    if (wcs->m_wtb) free(wcs->m_wtb);
    if (wcs->m_tab) {
      for (int j = 0; j < wcs->ntab; j++) {
        tabfree(wcs->m_tab + j);
      }

      free(wcs->m_tab);
    }
  }

  wcs->ntab = 0;
  wcs->nwtb = 0;
  wcs->wtb  = 0x0;
  wcs->tab  = 0x0;


  // Determine the number of -TAB axes.
  int naxis = wcs->naxis;
  int *tabax;
  if (!(tabax = calloc(naxis, sizeof(int)))) {
    return wcserr_set(WCSHDR_ERRMSG(WCSHDRERR_MEMORY));
  }

  int ntabax = 0;
  for (int i = 0; i < naxis; i++) {
    // Null fill.
    wcsutil_null_fill(72, wcs->ctype[i]);

    if (!strcmp(wcs->ctype[i]+4, "-TAB")) {
      tabax[i] = ntabax++;
    } else {
      tabax[i] = -1;
    }
  }

  if (ntabax == 0) {
    // No lookup tables.
    goto cleanup;
  }


  // Collect information from the PSi_ma and PVi_ma keyvalues.
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

  for (int itabax = 0; itabax < ntabax; itabax++) {
    // Remember that calloc() zeroes allocated memory.
    PVi_1a[itabax] = 1;
    PVi_2a[itabax] = 1;
    PVi_3a[itabax] = 1;
  }

  for (int ip = 0; ip < wcs->nps; ip++) {
    int itabax = tabax[wcs->ps[ip].i - 1];
    if (itabax >= 0) {
      switch (wcs->ps[ip].m) {
      case 0:
        // EXTNAME.
        strncpy(PSi_0a[itabax], wcs->ps[ip].value, 72);
        wcsutil_null_fill(72, PSi_0a[itabax]);
        break;
      case 1:
        // TTYPEn for coordinate array.
        strncpy(PSi_1a[itabax], wcs->ps[ip].value, 72);
        wcsutil_null_fill(72, PSi_1a[itabax]);
        break;
      case 2:
        // TTYPEn for index vector.
        strncpy(PSi_2a[itabax], wcs->ps[ip].value, 72);
        wcsutil_null_fill(72, PSi_2a[itabax]);
        break;
      default:
        // Undefined (ignored).
        break;
      }
    }
  }

  for (int ip = 0; ip < wcs->npv; ip++) {
    int itabax = tabax[wcs->pv[ip].i - 1];
    if (itabax >= 0) {
      switch (wcs->pv[ip].m) {
      case 1:
        // EXTVER.
        PVi_1a[itabax] = (int)(wcs->pv[ip].value + 0.5);
        break;
      case 2:
        // EXTLEVEL.
        PVi_2a[itabax] = (int)(wcs->pv[ip].value + 0.5);
        break;
      case 3:
        // Table axis number.
        PVi_3a[itabax] = (int)(wcs->pv[ip].value + 0.5);
        break;
      default:
        // Undefined (ignored).
        break;
      }
    }
  }


  // Determine the number of independent tables.
  for (int itabax = 0; itabax < ntabax; itabax++) {
    // These have no defaults.
    if (!PSi_0a[itabax][0] || !PSi_1a[itabax][0]) {
      status = wcserr_set(WCSERR_SET(WCSHDRERR_BAD_TABULAR_PARAMS),
        "Invalid tabular parameters: PSi_0a and PSi_1a must be specified");
      goto cleanup;
    }

    tabidx[itabax] = -1;
    int jtabax;
    for (jtabax = 0; jtabax < itabax; jtabax++) {
      // EXTNAME, EXTVER, EXTLEVEL, and TTYPEn for the coordinate array
      // must match for each axis of a multi-dimensional lookup table.
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

  // Table dimensionality; find the largest axis number.
  for (int itabax = 0; itabax < ntabax; itabax++) {
    struct tabprm *tabp = wcs->tab + tabidx[itabax];

    // PVi_3a records the 1-relative table axis number.
    if (PVi_3a[itabax] > tabp->M) {
      tabp->M = PVi_3a[itabax];
    }
  }

  for (int itab = 0; itab < wcs->ntab; itab++) {
    if ((status = tabini(1, wcs->tab[itab].M, 0, wcs->tab + itab))) {
      status = wcserr_set(WCSHDR_ERRMSG(wcshdr_taberr[status]));
      goto cleanup;
    }
  }


  // Copy parameters into the tabprm structs.
  for (int i = 0; i < naxis; i++) {
    int itabax;
    if ((itabax = tabax[i]) < 0) {
      // Not a -TAB axis.
      continue;
    }

    // PVi_3a records the 1-relative table axis number.
    int m = PVi_3a[itabax] - 1;

    struct tabprm *tabp;
    tabp = wcs->tab + tabidx[itabax];
    tabp->map[m] = i;
    tabp->crval[m] = wcs->crval[i];
  }

  // Check for completeness.
  for (int itab = 0; itab < wcs->ntab; itab++) {
    for (int m = 0; m < wcs->tab[itab].M; m++) {
      if (wcs->tab[itab].map[m] < 0) {
        status = wcserr_set(WCSERR_SET(WCSHDRERR_BAD_TABULAR_PARAMS),
          "Invalid tabular parameters: the axis mapping is undefined");
        goto cleanup;
      }
    }
  }


  // Set up for reading the arrays; how many arrays are there?
  for (int itabax = 0; itabax < ntabax; itabax++) {
    // Does this -TAB axis have a non-degenerate index array?
    if (PSi_2a[itabax][0]) {
      wcs->nwtb++;
    }
  }

  // Add one coordinate array for each table.
  wcs->nwtb += wcs->ntab;

  // Allocate memory for structs to be returned.
  if (!(wcs->wtb = calloc(wcs->nwtb, sizeof(struct wtbarr)))) {
    wcs->nwtb = 0;

    status = wcserr_set(WCSHDR_ERRMSG(WCSHDRERR_MEMORY));
    goto cleanup;
  }
  wcs->m_wtb = wcs->wtb;

  // Set pointers for the index and coordinate arrays.
  struct wtbarr *wtbp = wcs->wtb;
  for (int itab = 0; itab < wcs->ntab; itab++) {
    int getcrd = 1;
    for (int itabax = 0; itabax < ntabax; itabax++) {
      if (tabidx[itabax] != itab) continue;

      if (getcrd) {
        // Coordinate array.
        wtbp->i = itabax + 1;
        wtbp->m = PVi_3a[itabax];
        wtbp->kind = 'c';

        strncpy(wtbp->extnam, PSi_0a[itabax], 72);
        wtbp->extver = PVi_1a[itabax];
        wtbp->extlev = PVi_2a[itabax];
        strncpy(wtbp->ttype, PSi_1a[itabax], 72);
        wtbp->row    = 1L;
        wtbp->ndim   = wcs->tab[itab].M + 1;
        wtbp->dimlen = wcs->tab[itab].K;
        wtbp->arrayp = &(wcs->tab[itab].coord);

        // Signal for tabset() to take this memory.
        wcs->tab[itab].m_coord = (double *)0x1;

        wtbp++;
        getcrd = 0;
      }

      if (PSi_2a[itabax][0]) {
        // Index array.
        wtbp->i = itabax + 1;
        wtbp->m = PVi_3a[itabax];
        wtbp->kind = 'i';

        int m = wtbp->m - 1;
        strncpy(wtbp->extnam, PSi_0a[itabax], 72);
        wtbp->extver = PVi_1a[itabax];
        wtbp->extlev = PVi_2a[itabax];
        strncpy(wtbp->ttype, PSi_2a[itabax], 72);
        wtbp->row    = 1L;
        wtbp->ndim   = 1;
        wtbp->dimlen = wcs->tab[itab].K + m;
        wtbp->arrayp = wcs->tab[itab].index + m;

        // Signal for tabset() to take this memory.
        wcs->tab[itab].m_indxs[m] = (double *)0x1;

        wtbp++;
      }
    }
  }

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

//----------------------------------------------------------------------------

int wcsidx(int nwcs, struct wcsprm **wcs, int alts[27])

{
  for (int a = 0; a < 27; a++) {
    alts[a] = -1;
  }

  if (wcs == 0x0) {
    return WCSHDRERR_NULL_POINTER;
  }

  struct wcsprm *wcsp = *wcs;
  for (int iwcs = 0; iwcs < nwcs; iwcs++, wcsp++) {
    if (wcsp->colnum || wcsp->colax[0]) continue;

    int a;
    if (wcsp->alt[0] == ' ') {
      a = 0;
    } else {
      a = wcsp->alt[0] - 'A' + 1;
    }

    alts[a] = iwcs;
  }

  return 0;
}

//----------------------------------------------------------------------------

int wcsbdx(int nwcs, struct wcsprm **wcs, int type, short alts[1000][28])

{
  // End of the array.
  short *ipend = alts[0] + (ptrdiff_t)(28*1000);
  for (short *ip = alts[0]; ip < ipend; ip++) {
    *ip = -1;
  }

  for (int icol = 0; icol < 1000; icol++) {
    alts[icol][27] = 0;
  }

  if (wcs == 0x0) {
    return WCSHDRERR_NULL_POINTER;
  }

  struct wcsprm *wcsp = *wcs;
  for (short iwcs = 0; iwcs < (short)nwcs; iwcs++, wcsp++) {
    int a;
    if (wcsp->alt[0] == ' ') {
      a = 0;
    } else {
      a = wcsp->alt[0] - 'A' + 1;
    }

    if (type) {
      // Pixel list.
      if (wcsp->colax[0]) {
        for (int i = 0; i < wcsp->naxis; i++) {
          alts[wcsp->colax[i]][a] = iwcs;
          alts[wcsp->colax[i]][27]++;
        }
      } else if (!wcsp->colnum) {
        alts[0][a] = iwcs;
        alts[0][27]++;
      }

    } else {
      // Binary table image array.
      if (wcsp->colnum) {
        alts[wcsp->colnum][a] = iwcs;
        alts[wcsp->colnum][27]++;
      } else if (!wcsp->colax[0]) {
        alts[0][a] = iwcs;
        alts[0][27]++;
      }
    }
  }

  return 0;
}

//----------------------------------------------------------------------------

int wcsvfree(int *nwcs, struct wcsprm **wcs)

{
  int status = 0;

  if (wcs == 0x0) {
    return WCSHDRERR_NULL_POINTER;
  }

  struct wcsprm *wcsp = *wcs;
  for (int a = 0; a < *nwcs; a++, wcsp++) {
    status |= wcsfree(wcsp);
  }

  free(*wcs);

  *nwcs = 0;
  *wcs = 0x0;

  return status;
}

//----------------------------------------------------------------------------
// Matching the definitions in dis.c.
#define I_DTYPE   0	// Distortion type code.
#define I_NIPARM  1	// Full (allocated) length of iparm[].
#define I_NDPARM  2	// No. of parameters in dparm[], excl. work space.
#define I_TPDNCO  3	// No. of TPD coefficients, forward...
#define I_TPDINV  4	// ...and inverse.
#define I_TPDAUX  5	// True if auxiliary variables are used.
#define I_TPDRAD  6	// True if the radial variable is used.

int wcshdo(int ctrl, struct wcsprm *wcs, int *nkeyrec, char **header)

// ::: CUBEFACE and STOKES handling?

{
  static const char *function = "wcshdo";

  const char axid[] = "xyxuvu", *cp;
  const int  nTPD[] = {1, 4, 7, 12, 17, 24, 31, 40, 49, 60};

  *nkeyrec = 0;
  *header  = 0x0;

  if (wcs == 0x0) return WCSHDRERR_NULL_POINTER;
  struct wcserr **err = &(wcs->err);

  int status = 0;
  if (abs(wcs->flag) != WCSSET) {
    if ((status = wcsset(wcs))) return status;
  }

  int naxis;
  if ((naxis = wcs->naxis) == 0) {
    return 0;
  }


  // These are mainly for convenience.
  char alt = wcs->alt[0];
  if (alt == ' ') alt = '\0';
  int colnum = wcs->colnum;
  int *colax = wcs->colax;

  int primage = 0;
  int bintab  = 0;
  int pixlist = 0;
  if (colnum) {
    bintab  = 1;
  } else if (colax[0]) {
    pixlist = 1;
  } else {
    primage = 1;
  }


  // Initialize floating point format control.
  char format[16];
  *format = '\0';
  if (ctrl & WCSHDO_P17) {
    strncpy(format, "% 20.17G", 16);
  } else if (ctrl & WCSHDO_P16) {
    strncpy(format, "% 20.16G", 16);
  } else if (ctrl & WCSHDO_P15) {
    strncpy(format, "% 20.15G", 16);
  } else if (ctrl & WCSHDO_P14) {
    strncpy(format, "% 20.14G", 16);
  } else if (ctrl & WCSHDO_P13) {
    strncpy(format, "% 20.13G", 16);
  } else if (ctrl & WCSHDO_P12) {
    strncpy(format, "%20.12G", 16);
  }

  if (*format && (ctrl & WCSHDO_EFMT)) {
    if (format[6] == 'G') {
      format[6] = 'E';
    } else {
      format[7] = 'E';
    }
  }

  int dofmt = (*format == '\0');

  char keyword[16],  *kwe = keyword  + 16;
  char keyvalue[96], *kve = keyvalue + 96;
  char comment[128], *cme = comment + 128;


  // WCS dimension.
  if (!pixlist) {
    snprintf(keyvalue, 96, "%20d", naxis);
    wcshdo_util(ctrl, "WCSAXES", "WCAX", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "Number of coordinate axes", nkeyrec, header, &status);
  }

  // Reference pixel coordinates.
  if (dofmt) wcshdo_format('G', naxis, wcs->crpix, format);
  for (int j = 0; j < naxis; j++) {
    wcsutil_double2str(keyvalue, 96, format, wcs->crpix[j]);
    wcshdo_util(ctrl, "CRPIX", "CRP", WCSHDO_CRPXna, "CRPX", 0, j+1, 0, alt,
      colnum, colax, keyvalue, "Pixel coordinate of reference point", nkeyrec,
      header, &status);
  }

  // Linear transformation matrix.
  if (dofmt) wcshdo_format('G', naxis*naxis, wcs->pc, format);
  int k = 0;
  for (int i = 0; i < naxis; i++) {
    for (int j = 0; j < naxis; j++, k++) {
      if (i == j) {
        if (wcs->pc[k] == 1.0) continue;
      } else {
        if (wcs->pc[k] == 0.0) continue;
      }

      wcsutil_double2str(keyvalue, 96, format, wcs->pc[k]);
      wcshdo_util(ctrl, "PC", bintab ? "PC" : "P", WCSHDO_TPCn_ka,
        bintab ? 0x0 : "PC", i+1, j+1, 0, alt, colnum, colax,
        keyvalue, "Coordinate transformation matrix element",
        nkeyrec, header, &status);
    }
  }

  // Coordinate increment at reference point.
  if (dofmt) wcshdo_format('G', naxis, wcs->cdelt, format);
  for (int i = 0; i < naxis; i++) {
    wcsutil_double2str(keyvalue, 96, format, wcs->cdelt[i]);
    comment[0] = '\0';
    if (wcs->cunit[i][0]) snprintf(comment, 128, "[%s] ", wcs->cunit[i]);
    char *cm = comment + strlen(comment);
    size_t csize = cme - cm;
    strncpy(cm, "Coordinate increment at reference point", csize);
    wcshdo_util(ctrl, "CDELT", "CDE", WCSHDO_CRPXna, "CDLT", i+1, 0, 0, alt,
      colnum, colax, keyvalue, comment, nkeyrec, header, &status);
  }

  // Units of coordinate increment and reference value.
  for (int i = 0; i < naxis; i++) {
    if (wcs->cunit[i][0] == '\0') continue;

    snprintf(keyvalue, 96, "'%s'", wcs->cunit[i]);
    wcshdo_util(ctrl, "CUNIT", "CUN", WCSHDO_CRPXna, "CUNI", i+1, 0, 0, alt,
      colnum, colax, keyvalue, "Units of coordinate increment and value",
      nkeyrec, header, &status);
  }

  // May need to alter ctype for particular distortions so do basic checks
  // now.  Note that SIP, TPV, DSS, TNX, and ZPX are restricted to exactly
  // two axes and cannot coexist with other distortion types.
  char   tpdsrc[24];
  int    dosip = 0, dotpd = 0, dotpv = 0;
  struct disprm *dis;

  if ((dis = wcs->lin.dispre)) {
    for (int i = 0; i < naxis; i++) {
      if (strcmp(dis->dtype[i], "SIP") == 0) {
        // Simple Imaging Polynomial (SIP).  Write it in its native form
        // if possible, unless specifically requested to write it as TPD.
        dotpd = (dis->iparm[i][I_DTYPE] & DIS_DOTPD);

        if (!dotpd) {
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
            // Must have been read as a 'SIP' distortion, CPDISja = 'SIP'.
            // Cannot be written as native SIP so write it as TPD.
            dotpd = DIS_DOTPD;
          }

          if (strncmp(wcs->ctype[0], "RA---TAN", 8) != 0 ||
              strncmp(wcs->ctype[1], "DEC--TAN", 8) != 0) {
            // Must have been permuted by wcssub().
            // Native SIP doesn't have axis mapping so write it as TPD.
            dotpd = DIS_DOTPD;
          }

          if (dotpd) {
            strncpy(tpdsrc, "SIP coordinates", 24);
          } else {
            dosip = 1;
          }
        }

        break;
      }
    }
  }

  if ((dis = wcs->lin.disseq)) {
    for (int i = 0; i < naxis; i++) {
      if (strcmp(dis->dtype[i], "TPV") == 0) {
        // TPV "projection".  Write it in its native form if possible,
        // unless specifically requested to write it as TPD.
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
            // Must have been read as a 'TPV' distortion, CPDISja = 'TPV'.
            // Cannot be written as native TPV so write it as TPD.
            dotpd = DIS_DOTPD;
          }

          if (dotpd) {
            strncpy(tpdsrc, "TPV \"projection\"", 24);
          } else {
            dotpv = 1;
          }
        }

        break;

      } else if (strcmp(dis->dtype[i], "DSS") == 0) {
        // Always written as TPD.
        dotpd = DIS_DOTPD;
        strncpy(tpdsrc, dis->dtype[i], 24);

      } else if (strncmp(dis->dtype[i], "WAT", 3) == 0) {
        // Always written as TPD.
        dotpd = DIS_DOTPD;
        strncpy(tpdsrc, dis->dtype[i]+4, 24);

        if (strcmp(dis->dtype[i], "DSS") == 0) {
          strncpy(tpdsrc, wcs->wcsname, 24);
        } else {
          char *tp = tpdsrc + strlen(tpdsrc);
          size_t tsize = (tpdsrc + 24) - tp;
          strncpy(tp, " \"projection\"", tsize);
        }

        break;
      }
    }
  }

  // Coordinate type.
  for (int i = 0; i < naxis; i++) {
    if (wcs->ctype[i][0] == '\0') continue;

    snprintf(keyvalue, 96, "'%s'", wcs->ctype[i]);
    strncpy(comment, "Coordinate type code", 128);

    char *ctypei = keyvalue + 1;
    if (i == wcs->lng || i == wcs->lat) {
      // Alter ctype for particular distortions.
      if (dosip) {
        // It could have come in as CPDISja = 'SIP'.
        strncpy(ctypei+8, "-SIP'", 6);
      } else if (dotpv) {
        // Reinstate projection code edited by wcsset().
        strncpy(ctypei+4, "-TPV'", 6);
      }

      if (strncmp(ctypei+8, "-SIP", 4) == 0) {
        strncpy(comment, "TAN (gnomonic) projection + SIP distortions", 128);

      } else if (strncmp(ctypei+4, "-TPV", 4) == 0) {
        strncpy(comment, "TAN (gnomonic) projection + distortions", 128);

      } else {
        if (strncmp(ctypei, "RA--", 4) == 0) {
          strncpy(comment, "Right ascension, ", 128);

        } else if (strncmp(ctypei, "DEC-", 4) == 0) {
          strncpy(comment, "Declination, ", 128);

        } else if (strncmp(ctypei+1, "LON", 3) == 0 ||
                   strncmp(ctypei+1, "LAT", 3) == 0) {
          ctypei[0] = (char)toupper((unsigned char)ctypei[0]);

          switch (ctypei[0]) {
          case 'G':
            strncpy(comment, "Galactic l", 128);
            break;
          case 'E':
            strncpy(comment, "Ecliptic l", 128);
            break;
          case 'H':
            strncpy(comment, "Helioecliptic l", 128);
            break;
          case 'S':
            strncpy(comment, "Supergalactic l", 128);
            break;
          default:
            // User-defined coordinate system.
            strncpy(comment, "L", 128);
            break;
          }

          char *cm = comment + strlen(comment);
          size_t csize = cme - cm;
          if (i == wcs->lng) {
            strncpy(cm, "ongitude, ", csize);
          } else {
            strncpy(cm, "atitude, ", csize);
          }
        }

        char *cm = comment + strlen(comment);
        size_t csize = cme - cm;
        strncpy(cm, wcs->cel.prj.name, csize);

        cm += strlen(cm);
        csize = cme - cm;
        strncpy(cm, " projection", csize);
      }

    } else if (i == wcs->spec) {
      char ptype, xtype;
      spctyp(wcs->ctype[i], 0x0, 0x0, comment, 0x0, &ptype, &xtype, 0x0);

      char *cm = comment + strlen(comment);
      size_t csize = cme - cm;

      if (ptype == xtype) {
        strncpy(cm, " (linear)", csize);
      } else {
        switch (xtype) {
        case 'F':
          strncpy(cm, " (linear in frequency)", csize);
          break;
        case 'V':
          strncpy(cm, " (linear in velocity)", csize);
          break;
        case 'W':
          strncpy(cm, " (linear in wavelength)", csize);
          break;
        default:
          strncpy(cm, " (not recognised)", csize);
          break;
        }
      }
    }

    wcshdo_util(ctrl, "CTYPE", "CTY", WCSHDO_CRPXna, "CTYP", i+1, 0, 0, alt,
      colnum, colax, keyvalue, comment, nkeyrec, header, &status);
  }

  // Coordinate value at reference point.
  for (int i = 0; i < naxis; i++) {
    if (dofmt) wcshdo_format('G', 1, wcs->crval+i, format);
    wcsutil_double2str(keyvalue, 96, format, wcs->crval[i]);
    comment[0] = '\0';
    if (wcs->cunit[i][0]) snprintf(comment, 128, "[%s] ", wcs->cunit[i]);
    char *cm = comment + strlen(comment);
    size_t csize = cme - cm;
    strncpy(cm, "Coordinate value at reference point", csize);
    wcshdo_util(ctrl, "CRVAL", "CRV", WCSHDO_CRPXna, "CRVL", i+1, 0, 0, alt,
      colnum, colax, keyvalue, comment, nkeyrec, header, &status);
  }

  // Parameter values.
  if (dofmt) strncpy(format, "%20.12G", 16);
  for (int k = 0; k < wcs->npv; k++) {
    wcsutil_double2str(keyvalue, 96, format, (wcs->pv[k]).value);
    if ((wcs->pv[k]).i == (wcs->lng + 1)) {
      switch ((wcs->pv[k]).m) {
      case 1:
        strncpy(comment, "[deg] Native longitude of the reference point",
          128);
        break;
      case 2:
        strncpy(comment, "[deg] Native latitude  of the reference point",
          128);
        break;
      case 3:
        if (primage) {
          snprintf(keyword, 16, "LONPOLE%c", alt);
        } else if (bintab) {
          snprintf(keyword, 16, "LONP%d%c", colnum, alt);
        } else {
          snprintf(keyword, 16, "LONP%d%c", colax[(wcs->pv[k]).i - 1], alt);
        }
        snprintf(comment, 128, "[deg] alias for %s (has precedence)",
          keyword);
        break;
      case 4:
        if (primage) {
          snprintf(keyword, 16, "LATPOLE%c", alt);
        } else if (bintab) {
          snprintf(keyword, 16, "LATP%d%c", colnum, alt);
        } else {
          snprintf(keyword, 16, "LATP%d%c", colax[(wcs->pv[k]).i - 1], alt);
        }
        snprintf(comment, 128, "[deg] alias for %s (has precedence)",
          keyword);
        break;
      default:;
        char *cm = comment + strlen(comment);
        size_t csize = cme - cm;
        strncpy(cm, " (not recognised)", csize);
        break;
      }
    } else if ((wcs->pv[k]).i == (wcs->lat + 1)) {
      snprintf(comment, 128, "%s projection parameter", wcs->cel.prj.code);
    } else {
      strncpy(comment, "Coordinate transformation parameter", 128);
    }

    wcshdo_util(ctrl, "PV", "V", WCSHDO_PVn_ma, "PV", wcs->pv[k].i, -1,
      wcs->pv[k].m, alt, colnum, colax, keyvalue, comment, nkeyrec, header,
      &status);
  }

  for (int k = 0; k < wcs->nps; k++) {
    snprintf(keyvalue, 96, "'%s'", (wcs->ps[k]).value);
    wcshdo_util(ctrl, "PS", "S", WCSHDO_PVn_ma, "PS", wcs->ps[k].i, -1,
      wcs->ps[k].m, alt, colnum, colax, keyvalue,
      "Coordinate transformation parameter", nkeyrec, header, &status);
  }

  // Celestial and spectral transformation parameters.
  if (!undefined(wcs->lonpole)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->lonpole);
    wcshdo_util(ctrl, "LONPOLE", "LONP", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "[deg] Native longitude of celestial pole", nkeyrec, header,
      &status);
  }

  if (!undefined(wcs->latpole)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->latpole);
    wcshdo_util(ctrl, "LATPOLE", "LATP", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "[deg] Native latitude of celestial pole", nkeyrec, header,
      &status);
  }

  if (wcs->restfrq != 0.0) {
    wcsutil_double2str(keyvalue, 96, format, wcs->restfrq);
    wcshdo_util(ctrl, "RESTFRQ", "RFRQ", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "[Hz] Line rest frequency", nkeyrec, header, &status);
  }

  if (wcs->restwav != 0.0) {
    wcsutil_double2str(keyvalue, 96, format, wcs->restwav);
    wcshdo_util(ctrl, "RESTWAV", "RWAV", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "[Hz] Line rest wavelength", nkeyrec, header, &status);
  }

  // - - - - - - - - - - - - - - - - -  Auxiliary coordinate axis information.
  char timeunit[16];
  snprintf(timeunit, 16, "%.15s", wcs->timeunit[0] ? wcs->timeunit : "s");

  // Coordinate axis title.
  if (wcs->cname) {
    for (int i = 0; i < naxis; i++) {
      if (wcs->cname[i][0] == '\0') continue;

      snprintf(keyvalue, 96, "'%s'", wcs->cname[i]);
      wcshdo_util(ctrl, "CNAME", "CNA", WCSHDO_CNAMna, "CNAM", i+1, 0, 0,
        alt, colnum, colax, keyvalue, "Axis name for labelling purposes",
        nkeyrec, header, &status);
    }
  }

  // Random error in coordinate.
  if (wcs->crder) {
    for (int i = 0; i < naxis; i++) {
      if (undefined(wcs->crder[i])) continue;

      wcsutil_double2str(keyvalue, 96, format, wcs->crder[i]);
      comment[0] = '\0';
      if (wcs->cunit[i][0]) snprintf(comment, 128, "[%s] ", wcs->cunit[i]);
      char *cm = comment + strlen(comment);
      size_t csize = cme - cm;
      strncpy(cm, "Random error in coordinate", csize);
      wcshdo_util(ctrl, "CRDER", "CRD", WCSHDO_CNAMna, "CRDE", i+1, 0, 0,
        alt, colnum, colax, keyvalue, comment, nkeyrec, header, &status);
    }
  }

  // Systematic error in coordinate.
  if (wcs->csyer) {
    for (int i = 0; i < naxis; i++) {
      if (undefined(wcs->csyer[i])) continue;

      wcsutil_double2str(keyvalue, 96, format, wcs->csyer[i]);
      comment[0] = '\0';
      if (wcs->cunit[i][0]) snprintf(comment, 128, "[%s] ", wcs->cunit[i]);
      char *cm = comment + strlen(comment);
      size_t csize = cme - cm;
      strncpy(cm, "Systematic error in coordinate", csize);
      wcshdo_util(ctrl, "CSYER", "CSY", WCSHDO_CNAMna, "CSYE", i+1, 0, 0,
        alt, colnum, colax, keyvalue, comment, nkeyrec, header, &status);
    }
  }

  // Time at zero point of phase axis.
  if (wcs->czphs) {
    for (int i = 0; i < naxis; i++) {
      if (undefined(wcs->czphs[i])) continue;

      wcsutil_double2str(keyvalue, 96, format, wcs->czphs[i]);
      snprintf(comment, 128, "[%s] Time at zero point of phase axis",
        timeunit);
      wcshdo_util(ctrl, "CZPHS", "CZP", WCSHDO_CNAMna, "CZPH", i+1, 0, 0,
        alt, colnum, colax, keyvalue, comment, nkeyrec, header, &status);
    }
  }

  // Period of phase axis.
  if (wcs->cperi) {
    for (int i = 0; i < naxis; i++) {
      if (undefined(wcs->cperi[i])) continue;

      wcsutil_double2str(keyvalue, 96, format, wcs->cperi[i]);
      snprintf(comment, 128, "[%s] Period of phase axis", timeunit);
      wcshdo_util(ctrl, "CPERI", "CPR", WCSHDO_CNAMna, "CPER", i+1, 0, 0,
        alt, colnum, colax, keyvalue, comment, nkeyrec, header, &status);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - -  Coordinate system title.

  // Coordinate system title.
  if (wcs->wcsname[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->wcsname);
    if (bintab) {
      wcshdo_util(ctrl, "WCSNAME", "WCSN", 0, 0x0, 0, 0, 0, alt, colnum,
        colax, keyvalue, "Coordinate system title", nkeyrec, header, &status);
    } else {
      // TWCS was a mistake.
      wcshdo_util(ctrl, "WCSNAME", "TWCS", WCSHDO_WCSNna, "WCSN", 0, 0, 0,
        alt, colnum, colax, keyvalue, "Coordinate system title",
        nkeyrec, header, &status);
    }
  }

  // - - - - - - - - - - - - - - - - -  Time reference system and measurement.

  // Time scale.
  if (wcs->timesys[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->timesys);
    wcshdo_util(ctrl, "TIMESYS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "Time scale", nkeyrec, header, &status);
  }

  // Time reference position.
  if (wcs->trefpos[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->trefpos);
    wcshdo_util(ctrl, "TREFPOS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "Time reference position", nkeyrec, header, &status);
  }

  // Time reference direction.
  if (wcs->trefdir[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->trefdir);
    wcshdo_util(ctrl, "TREFDIR", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "Time reference direction", nkeyrec, header, &status);
  }

  // Ephemerides used for pathlength delay calculation.
  if (wcs->plephem[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->plephem);
    wcshdo_util(ctrl, "PLEPHEM", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "Ephemerides used for pathlength delays", nkeyrec, header,
      &status);
  }

  // Time units.
  if (wcs->timeunit[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->timeunit);
    wcshdo_util(ctrl, "TIMEUNIT", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "Time units", nkeyrec, header, &status);
  }

  // Fiducial (reference) time.
  if (wcs->mjdref[0] == 0.0 && wcs->mjdref[1] == 0.0) {
    // MJD of fiducial time (simplified if it takes its default value).
    wcsutil_double2str(keyvalue, 96, format, 0.0);
    wcshdo_util(ctrl, "MJDREF", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "[d] MJD of fiducial time", nkeyrec, header, &status);

  } else {
    // ISO-8601 fiducial time.
    if (wcs->dateref[0]) {
      snprintf(keyvalue, 96, "'%s'", wcs->dateref);
      wcshdo_util(ctrl, "DATEREF", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, "ISO-8601 fiducial time", nkeyrec, header, &status);
    }

    if (wcs->mjdref[1] == 0.0) {
      // MJD of fiducial time (no fractional part).
      if (!undefined(wcs->mjdref[0])) {
        wcsutil_double2str(keyvalue, 96, format, wcs->mjdref[0]);
        wcshdo_util(ctrl, "MJDREF", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
          keyvalue, "[d] MJD of fiducial time", nkeyrec, header, &status);
      }

    } else {
      // MJD of fiducial time, integer part.
      if (!undefined(wcs->mjdref[0])) {
        wcsutil_double2str(keyvalue, 96, format, wcs->mjdref[0]);
        wcshdo_util(ctrl, "MJDREFI", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
          keyvalue, "[d] MJD of fiducial time, integer part", nkeyrec,
          header, &status);
      }

      // MJD of fiducial time, fractional part.
      if (!undefined(wcs->mjdref[1])) {
        wcsutil_double2str(keyvalue, 96, format, wcs->mjdref[1]);
        wcshdo_util(ctrl, "MJDREFF", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
          keyvalue, "[d] MJD of fiducial time, fractional part", nkeyrec,
          header, &status);
      }
    }
  }

  // Clock correction.
  if (!undefined(wcs->timeoffs)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->timeoffs);
    snprintf(comment, 128, "[%s] Clock correction", timeunit);
    wcshdo_util(ctrl, "TIMEOFFS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, comment, nkeyrec, header, &status);
  }

  // - - - - - - - - - - - - - - - - - - - - -  Data timestamps and durations.

  // ISO-8601 time of observation.
  if (wcs->dateobs[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->dateobs);
    strncpy(comment, "ISO-8601 time of observation", 128);

    if (ctrl & 1) {
      // Allow DOBSn.
      wcshdo_util(ctrl, "DATE-OBS", "DOBS", WCSHDO_DOBSn, 0x0, 0, 0, 0, ' ',
        colnum, colax, keyvalue, comment, nkeyrec, header, &status);
    } else {
      // Force DATE-OBS.
      wcshdo_util(ctrl, "DATE-OBS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, comment, nkeyrec, header, &status);
    }
  }

  // MJD of observation.
  if (!undefined(wcs->mjdobs)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->mjdobs);
    wcshdo_util(ctrl, "MJD-OBS", "MJDOB", 0, 0x0, 0, 0, 0, ' ', colnum, colax,
      keyvalue, "[d] MJD of observation", nkeyrec, header, &status);
  }

  // Julian epoch of observation.
  if (!undefined(wcs->jepoch)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->jepoch);
    wcshdo_util(ctrl, "JEPOCH", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "[a] Julian epoch of observation", nkeyrec, header, &status);
  }

  // Besselian epoch of observation.
  if (!undefined(wcs->bepoch)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->bepoch);
    wcshdo_util(ctrl, "BEPOCH", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "[a] Besselian epoch of observation", nkeyrec, header,
      &status);
  }

  // ISO-8601 time at start of observation.
  if (wcs->datebeg[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->datebeg);
    wcshdo_util(ctrl, "DATE-BEG", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "ISO-8601 time at start of observation", nkeyrec, header,
      &status);
  }

  // MJD at start of observation.
  if (!undefined(wcs->mjdbeg)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->mjdbeg);
    wcshdo_util(ctrl, "MJD-BEG", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "[d] MJD at start of observation", nkeyrec, header, &status);
  }

  // Time elapsed at start since fiducial time.
  if (!undefined(wcs->tstart)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->tstart);
    snprintf(comment, 128, "[%s] Time elapsed since fiducial time at start",
      timeunit);
    wcshdo_util(ctrl, "TSTART", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, comment, nkeyrec, header, &status);
  }

  // ISO-8601 time at midpoint of observation.
  if (wcs->dateavg[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->dateavg);
    wcshdo_util(ctrl, "DATE-AVG", "DAVG", 0, 0x0, 0, 0, 0, ' ', colnum, colax,
      keyvalue, "ISO-8601 time at midpoint of observation", nkeyrec, header,
      &status);
  }

  // MJD at midpoint of observation.
  if (!undefined(wcs->mjdavg)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->mjdavg);
    wcshdo_util(ctrl, "MJD-AVG", "MJDA", 0, 0x0, 0, 0, 0, ' ', colnum, colax,
      keyvalue, "[d] MJD at midpoint of observation", nkeyrec, header,
      &status);
  }

  // ISO-8601 time at end of observation.
  if (wcs->dateend[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->dateend);
    wcshdo_util(ctrl, "DATE-END", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "ISO-8601 time at end of observation", nkeyrec, header,
      &status);
  }

  // MJD at end of observation.
  if (!undefined(wcs->mjdend)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->mjdend);
    wcshdo_util(ctrl, "MJD-END", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "[d] MJD at end of observation", nkeyrec, header, &status);
  }

  // Time elapsed at end since fiducial time.
  if (!undefined(wcs->tstop)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->tstop);
    snprintf(comment, 128, "[%s] Time elapsed since fiducial time at end",
      timeunit);
    wcshdo_util(ctrl, "TSTOP", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, comment, nkeyrec, header, &status);
  }

  // Exposure (integration) time.
  if (!undefined(wcs->xposure)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->xposure);
    snprintf(comment, 128, "[%s] Exposure (integration) time", timeunit);
    wcshdo_util(ctrl, "XPOSURE", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, comment, nkeyrec, header, &status);
  }

  // Elapsed time (start to stop).
  if (!undefined(wcs->telapse)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->telapse);
    snprintf(comment, 128, "[%s] Elapsed time (start to stop)", timeunit);
    wcshdo_util(ctrl, "TELAPSE", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, comment, nkeyrec, header, &status);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - -  Timing accuracy.

  // Systematic error in time measurements.
  if (!undefined(wcs->timsyer)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->timsyer);
    snprintf(comment, 128, "[%s] Systematic error in time measurements",
      timeunit);
    wcshdo_util(ctrl, "TIMSYER", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, comment, nkeyrec, header, &status);
  }

  // Relative error in time measurements.
  if (!undefined(wcs->timrder)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->timrder);
    snprintf(comment, 128, "[%s] Relative error in time measurements",
      timeunit);
    wcshdo_util(ctrl, "TIMRDER", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, comment, nkeyrec, header, &status);
  }

  // Time resolution.
  if (!undefined(wcs->timedel)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->timedel);
    snprintf(comment, 128, "[%s] Time resolution", timeunit);
    wcshdo_util(ctrl, "TIMEDEL", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, comment, nkeyrec, header, &status);
  }

  // Reference position of timestamp in binned data.
  if (!undefined(wcs->timepixr)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->timepixr);
    wcshdo_util(ctrl, "TIMEPIXR", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "Reference position of timestamp in binned data", nkeyrec,
      header, &status);
  }

  // - - - - - - - - - - - - - - - - - -  Spatial & celestial reference frame.

  // Observatory coordinates.
  if (!undefined(wcs->obsgeo[0]) &&
      !undefined(wcs->obsgeo[1]) &&
      !undefined(wcs->obsgeo[2])) {

    char obsgeo[16] = "OBSGEO-?", obsg[8] = "OBSG?", xyz[] = "XYZ";
    for (int k = 0; k < 3; k++) {
      wcsutil_double2str(keyvalue, 96, format, wcs->obsgeo[k]);
      snprintf(comment, 128, "[m] observatory %c-coordinate", xyz[k]);
      obsgeo[7] = xyz[k];
      obsg[4]   = xyz[k];
      wcshdo_util(ctrl, obsgeo, obsg, 0, 0x0, 0, 0, 0, ' ', colnum, colax,
        keyvalue, comment, nkeyrec, header, &status);
    }

  } else if (
      !undefined(wcs->obsgeo[3]) &&
      !undefined(wcs->obsgeo[4]) &&
      !undefined(wcs->obsgeo[5])) {

    wcsutil_double2str(keyvalue, 96, format, wcs->obsgeo[3]);
    wcshdo_util(ctrl, "OBSGEO-L", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "[deg] IAU(1976) observatory longitude", nkeyrec, header,
      &status);

    wcsutil_double2str(keyvalue, 96, format, wcs->obsgeo[4]);
    wcshdo_util(ctrl, "OBSGEO-B", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "[deg] IAU(1976) observatory latitude", nkeyrec, header,
      &status);

    wcsutil_double2str(keyvalue, 96, format, wcs->obsgeo[5]);
    wcshdo_util(ctrl, "OBSGEO-L", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "[m]   IAU(1976) observatory height", nkeyrec, header,
      &status);
  }

  // Spacecraft orbit ephemeris file.
  if (wcs->obsorbit[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->obsorbit);
    wcshdo_util(ctrl, "OBSORBIT", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
      keyvalue, "Spacecraft orbit ephemeris file", nkeyrec, header, &status);
  }

  // Equatorial coordinate system type.
  if (wcs->radesys[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->radesys);
    wcshdo_util(ctrl, "RADESYS", "RADE", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "Equatorial coordinate system", nkeyrec, header, &status);
  }

  // Equinox of equatorial coordinate system.
  if (!undefined(wcs->equinox)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->equinox);
    wcshdo_util(ctrl, "EQUINOX", "EQUI", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "[yr] Equinox of equatorial coordinates", nkeyrec, header,
      &status);
  }

  // Reference frame of spectral coordinates.
  if (wcs->specsys[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->specsys);
    wcshdo_util(ctrl, "SPECSYS", "SPEC", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "Reference frame of spectral coordinates", nkeyrec, header,
      &status);
  }

  // Reference frame of spectral observation.
  if (wcs->ssysobs[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->ssysobs);
    wcshdo_util(ctrl, "SSYSOBS", "SOBS", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "Reference frame of spectral observation", nkeyrec, header,
      &status);
  }

  // Observer's velocity towards source.
  if (!undefined(wcs->velosys)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->velosys);
    wcshdo_util(ctrl, "VELOSYS", "VSYS", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "[m/s] Velocity towards source", nkeyrec, header, &status);
  }

  // Redshift of the source.
  if (!undefined(wcs->zsource)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->zsource);
    wcshdo_util(ctrl, "ZSOURCE", "ZSOU", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "Redshift of the source", nkeyrec, header, &status);
  }

  // Reference frame of source redshift.
  if (wcs->ssyssrc[0]) {
    snprintf(keyvalue, 96, "'%s'", wcs->ssyssrc);
    wcshdo_util(ctrl, "SSYSSRC", "SSRC", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "Reference frame of source redshift", nkeyrec, header,
      &status);
  }

  // Velocity orientation angle.
  if (!undefined(wcs->velangl)) {
    wcsutil_double2str(keyvalue, 96, format, wcs->velangl);
    wcshdo_util(ctrl, "VELANGL", "VANG", 0, 0x0, 0, 0, 0, alt, colnum, colax,
      keyvalue, "[deg] Velocity orientation angle", nkeyrec, header, &status);
  }

  // - - - - - - - - - - - - - - - - - - - -  Additional auxiliary parameters.

  struct auxprm *aux;
  if ((aux = wcs->aux)) {
    if (!undefined(aux->rsun_ref)) {
      wcsutil_double2str(keyvalue, 96, format, aux->rsun_ref);
      wcshdo_util(ctrl, "RSUN_REF", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, "[m] Solar radius", nkeyrec, header, &status);
    }

    if (!undefined(aux->dsun_obs)) {
      wcsutil_double2str(keyvalue, 96, format, aux->dsun_obs);
      wcshdo_util(ctrl, "DSUN_OBS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, "[m] Distance from centre of Sun to observer", nkeyrec,
        header, &status);
    }

    if (!undefined(aux->crln_obs)) {
      wcsutil_double2str(keyvalue, 96, format, aux->crln_obs);
      wcshdo_util(ctrl, "CRLN_OBS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, "[deg] Carrington heliographic lng of observer", nkeyrec,
        header, &status);

      if (!undefined(aux->hglt_obs)) {
        wcsutil_double2str(keyvalue, 96, format, aux->hglt_obs);
        wcshdo_util(ctrl, "CRLT_OBS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
          keyvalue, "[deg] Heliographic latitude of observer", nkeyrec,
          header, &status);
      }
    }

    if (!undefined(aux->hgln_obs)) {
      wcsutil_double2str(keyvalue, 96, format, aux->hgln_obs);
      wcshdo_util(ctrl, "HGLN_OBS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, "[deg] Stonyhurst heliographic lng of observer", nkeyrec,
        header, &status);

      if (!undefined(aux->hglt_obs)) {
        wcsutil_double2str(keyvalue, 96, format, aux->hglt_obs);
        wcshdo_util(ctrl, "HGLT_OBS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
          keyvalue, "[deg] Heliographic latitude of observer", nkeyrec,
          header, &status);
      }
    }

    if (!undefined(aux->a_radius)) {
      wcsutil_double2str(keyvalue, 96, format, aux->a_radius);
      wcshdo_util(ctrl, "A_RADIUS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, "[m] Object ellipsoid semi-major axis", nkeyrec,
        header, &status);
    }

    if (!undefined(aux->b_radius)) {
      wcsutil_double2str(keyvalue, 96, format, aux->b_radius);
      wcshdo_util(ctrl, "B_RADIUS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, "[m] Object ellipsoid semi-intermediate axis", nkeyrec,
        header, &status);
    }

    if (!undefined(aux->c_radius)) {
      wcsutil_double2str(keyvalue, 96, format, aux->c_radius);
      wcshdo_util(ctrl, "C_RADIUS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, "[m] Object ellipsoid semi-minor axis", nkeyrec,
        header, &status);
    }

    if (!undefined(aux->blon_obs)) {
      wcsutil_double2str(keyvalue, 96, format, aux->blon_obs);
      wcshdo_util(ctrl, "BLON_OBS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, "[deg] Bodycentric longitude of observer", nkeyrec,
        header, &status);
    }

    if (!undefined(aux->blat_obs)) {
      wcsutil_double2str(keyvalue, 96, format, aux->blat_obs);
      wcshdo_util(ctrl, "BLAT_OBS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, "[deg] Bodycentric latitude of observer", nkeyrec,
        header, &status);
    }

    if (!undefined(aux->bdis_obs)) {
      wcsutil_double2str(keyvalue, 96, format, aux->bdis_obs);
      wcshdo_util(ctrl, "BDIS_OBS", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, "[m] Bodycentric distance of observer", nkeyrec,
        header, &status);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - Distortion function parameters.

  char term[16];
  double *dparm, keyval;

  if (dosip) {
    // Simple Imaging Polynomial (SIP) is handled by translating its dpkey
    // records.  Determine a suitable numerical precision for the
    // polynomial coefficients to avoid trailing zeroes common to all of
    // them.
    dis = wcs->lin.dispre;
    if (dofmt) {
      struct dpkey *ikeyp = dis->dp;
      int kp0 = 2;
      for (int idp = 0; idp < dis->ndp; idp++, ikeyp++) {
        cp = strchr(ikeyp->field, '.') + 1;
        if (strncmp(cp, "SIP.", 4) != 0) continue;
        wcsutil_double2str(keyvalue, 96, "%20.13E", dpkeyd(ikeyp));

        int kpi = 15;
        while (kp0 < kpi && keyvalue[kpi] == '0') kpi--;
        kp0 = kpi;
      }

      int precision = kp0 - 2;
      if (precision < 1)  precision = 1;
      if (13 < precision) precision = 13;
      snprintf(format, 16, "%%20.%dE", precision);
    }

    // Ensure the coefficients are written in a human-readable sequence.
    for (int j = 0; j <= 1; j++) {
      // Distortion function polynomial coefficients.
      wcshdo_util(ctrl, "", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        "", "", nkeyrec, header, &status);

      if (j == 0) {
        strncpy(keyword, "A_", 16);
      } else {
        strncpy(keyword, "B_", 16);
      }

      int degree;
      int ncoeff = dis->iparm[j][I_TPDNCO];
      for (degree = 0; degree <= 9; degree++) {
        if (ncoeff <= nTPD[degree]) break;
      }

      strncpy(keyword+2, "ORDER", 14);
      snprintf(keyvalue, 96, "%20d", degree);
      snprintf(comment, 128, "SIP polynomial degree, axis %d, pixel-to-sky",
        j+1);
      wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, comment, nkeyrec, header, &status);

      struct dpkey *ikeyp = dis->dp;
      for (int idp = 0; idp < dis->ndp; idp++, ikeyp++) {
        if (ikeyp->j != j+1) continue;
        if ((keyval = dpkeyd(ikeyp)) == 0.0) continue;

        cp = strchr(ikeyp->field, '.') + 1;
        if (strncmp(cp, "SIP.FWD.", 8) != 0) continue;
        cp += 8;
        strncpy(keyword+2, cp, 14);
        int p, q;
        sscanf(cp, "%d_%d", &p, &q);
        strncpy(term, "xxxxxxxxx", p);
        strncpy(term+p, "yyyyyyyyy", q);
        term[p+q] = '\0';

        wcsutil_double2str(keyvalue, 96, format, keyval);
        snprintf(comment, 128, "SIP distortion coefficient: %s", term);
        wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
          keyvalue, comment, nkeyrec, header, &status);
      }

      if (dis->maxdis[j] != 0.0) {
        strncpy(keyword+2, "DMAX", 14);
        wcsutil_double2str(keyvalue, 96, "%20.3f", dis->maxdis[j]);
        wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
          keyvalue, "Maximum value of distortion function", nkeyrec,
          header, &status);
      }

      // Inverse distortion function polynomial coefficients.
      if (dis->disx2p == 0x0) continue;

      wcshdo_util(ctrl, "", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        "", "", nkeyrec, header, &status);

      if (j == 0) {
        strncpy(keyword, "AP_", 16);
      } else {
        strncpy(keyword, "BP_", 16);
      }

      ncoeff = dis->iparm[j][I_NDPARM] - dis->iparm[j][I_TPDNCO];
      for (degree = 0; degree <= 9; degree++) {
        if (ncoeff <= nTPD[degree]) break;
      }

      strncpy(keyword+3, "ORDER", 13);
      snprintf(keyvalue, 96, "%20d", degree);
      snprintf(comment, 128, "SIP polynomial degree, axis %d, sky-to-pixel",
        j+1);
      wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
        keyvalue, comment, nkeyrec, header, &status);

      ikeyp = dis->dp;
      for (int idp = 0; idp < dis->ndp; idp++, ikeyp++) {
        if (ikeyp->j != j+1) continue;
        if ((keyval = dpkeyd(ikeyp)) == 0.0) continue;

        cp = strchr(ikeyp->field, '.') + 1;
        if (strncmp(cp, "SIP.REV.", 8) != 0) continue;
        cp += 8;
        strncpy(keyword+3, cp, 13);
        int p, q;
        sscanf(cp, "%d_%d", &p, &q);
        strncpy(term, "xxxxxxxxx", p);
        strncpy(term+p, "yyyyyyyyy", q);
        term[p+q] = '\0';

        wcsutil_double2str(keyvalue, 96, format, keyval);
        snprintf(comment, 128, "SIP inverse coefficient: %s", term);
        wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
          keyvalue, comment, nkeyrec, header, &status);
      }
    }
  }

  for (int idis = 0; idis < 2; idis++) {
    if (idis == 0 && (dis = wcs->lin.dispre) == 0x0) continue;
    if (idis == 1 && (dis = wcs->lin.disseq) == 0x0) continue;

    for (int j = 0; j < naxis; j++) {
      if (dis->disp2x[j] == 0x0) continue;

      int *iparm = dis->iparm[j];
      dparm = dis->dparm[j];

      // Identify the distortion type.
      if (dotpv) {
        // TPV "projection" is handled by translating its dpkey records,
        // which were originally translated from PVi_ma by wcsset(), or
        // possibly input directly as a CQDISia = 'TPV' distortion type.
        // Determine a suitable numerical precision for the polynomial
        // coefficients to avoid trailing zeroes common to all of them.
        if (dofmt) wcshdo_format('E', iparm[I_NDPARM], dparm, format);
        char fmt01[8];
        snprintf(fmt01, 8, "%.3ss", format);

        wcshdo_util(ctrl, "", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
          "", "", nkeyrec, header, &status);

        // Distortion function polynomial coefficients.
        snprintf(keyword, 16, "PV%d_", j+1);
        char *kw = keyword + strlen(keyword);

        struct dpkey *ikeyp = dis->dp;
        for (int idp = 0; idp < dis->ndp; idp++, ikeyp++) {
          if (ikeyp->j != j+1) continue;
          if ((keyval = dpkeyd(ikeyp)) == 0.0) continue;

          cp = strchr(ikeyp->field, '.') + 1;
          if (strncmp(cp, "TPV.", 4) != 0) continue;

          size_t ksize = kwe - kw;
          strncpy(kw, cp+4, ksize);

          // Identify the term of the TPV polynomial for human readers.
          int m;
          sscanf(cp+4, "%d", &m);
          wcshdo_tpdterm(m, j == wcs->lng, term);
          snprintf(comment, 128, "TPV coefficient: %s", term);

          if (keyval == 1.0) {
            snprintf(keyvalue, 96, fmt01, "1.0");
          } else {
            wcsutil_double2str(keyvalue, 96, format, keyval);
          }
          wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
            keyvalue, comment, nkeyrec, header, &status);
        }

      } else if (strcmp(dis->dtype[j], "TPD") == 0 || dotpd ||
                 strcmp(dis->dtype[j], "Polynomial")  == 0 ||
                 strcmp(dis->dtype[j], "Polynomial*") == 0) {
        // One of the Paper IV type polynomial distortions.
        wcshdo_util(ctrl, "", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
          "", "", nkeyrec, header, &status);

        if (strcmp(dis->dtype[j], "TPD") == 0) {
          // Pure TPD.
          dotpd = 1;
        } else if (strncmp(dis->dtype[j], "Polynomial", 10) == 0) {
          // Polynomial distortion.  Write it as TPD by request?
          dotpd = (iparm[I_DTYPE] & DIS_DOTPD);
          strncpy(tpdsrc, "Polynomial distortion", 24);
        }

        char pq = idis ? 'Q' : 'P';
        int Nhat = dis->Nhat[j];

        // CPDISja/CQDISia
        snprintf(keyword, 16, "C%cDIS%d", pq, j+1);
        if (idis == 0) {
          strncpy(comment, "P = prior, ", 128);
        } else {
          strncpy(comment, "Q = sequent, ", 128);
        }

        int direct = 0, doaux = 0;
        char *cm = comment + strlen(comment);
        size_t csize = cme - cm;
        if (dotpd) {
          strncpy(keyvalue, "'TPD'", 96);
          strncpy(cm, "Template Polynomial Distortion", csize);

          // For identifying terms of the TPD polynomial.
          int *axmap = dis->axmap[j];
          direct = 1;
          doaux  = iparm[I_TPDAUX];
          if (Nhat == 2) {
            // Associate x with longitude, y with latitude.
            if (axmap[0] == wcs->lng && axmap[1] == wcs->lat) {
              direct = 1;
            } else if (axmap[0] == wcs->lat && axmap[1] == wcs->lng) {
              direct = 0;
            } else {
              // Non-celestial.
              direct = (axmap[0] < axmap[1]);
            }
          }
        } else {
          strncpy(keyvalue, "'Polynomial'", 96);
          strncpy(cm, "general Polynomial distortion", csize);
        }

        wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
          keyvalue, comment, nkeyrec, header, &status);

        // NAXES.
        snprintf(keyword, 16,  "D%c%d", pq, j+1);
        snprintf(keyvalue, 96, "'NAXES:  %d'", Nhat);
        if (Nhat == 1) {
          strncpy(comment,  "One independent variable", 128);
        } else if (Nhat == 2) {
          strncpy(comment,  "Two independent variables", 128);
        } else {
          strncpy(comment,  "Number of independent variables", 128);
        }
        wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
          keyvalue, comment, nkeyrec, header, &status);

        // AXIS.jhat
        for (int jhat = 0; jhat < Nhat; jhat++) {
          int *axmap = dis->axmap[j];
          snprintf(keyvalue, 96, "'AXIS.%d: %d'", jhat+1, axmap[jhat]+1);
          if (jhat == 0) {
            strncpy(comment, "1st", 128);
          } else if (jhat == 1) {
            strncpy(comment, "2nd", 128);
          } else if (jhat == 2) {
            strncpy(comment, "3rd", 128);
          } else {
            snprintf(comment, 128, "%dth", jhat+1);
          }

          snprintf(comment, 128+strlen(comment),
            " independent variable: axis %d",
            axmap[jhat]+1);
          if (dotpd) {
            // axid is "xyxuvu".
            cp = axid;
            if (!direct) cp++;
            if (doaux) cp += 3;
            snprintf(comment, 128-strlen(comment), " (= %c)", cp[jhat]);
          }

          wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
            keyvalue, comment, nkeyrec, header, &status);
        }

        char ctemp[32];

        // OFFSET.jhat
        if (dofmt) wcshdo_format('f', Nhat, dis->offset[j], format);
        for (int jhat = 0; jhat < Nhat; jhat++) {
          if (dis->offset[j][jhat] == 0.0) continue;

          wcsutil_double2str(ctemp, 32, format, dis->offset[j][jhat]);
          snprintf(keyvalue, 96, "'OFFSET.%d: %s'", jhat+1, ctemp);
          snprintf(comment, 128, "Variable %d renormalization offset",
            jhat+1);

          wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
            keyvalue, comment, nkeyrec, header, &status);
        }

        // SCALE.jhat
        if (dofmt) wcshdo_format('f', Nhat, dis->scale[j], format);
        for (int jhat = 0; jhat < Nhat; jhat++) {
          if (dis->scale[j][jhat] == 1.0) continue;

          wcsutil_double2str(ctemp, 32, format, dis->scale[j][jhat]);
          snprintf(keyvalue, 96, "'SCALE.%d: %s'", jhat+1, ctemp);
          snprintf(comment, 128, "Variable %d renormalization scale", jhat+1);

          wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
            keyvalue, comment, nkeyrec, header, &status);
        }

        // Does the distortion function compute a correction?
        if (dis->docorr[j]) {
          wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
            "'DOCORR: 1'", "Distortion function computes a correction",
            nkeyrec, header, &status);
        } else {
          wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
            "'DOCORR: 0'", "Distortion function computes coordinates",
            nkeyrec, header, &status);
        }

        if (dotpd) {
          // Template Polynomial Distortion (TPD).  As it may have been
          // translated from SIP, TPV, DSS, TNX, ZPX, or perhaps
          // Polynomial, the dpkey records may not relate to TPD.
          // Output is therefore handled via dparm.
          if (dofmt) wcshdo_format('E', iparm[I_NDPARM], dparm, format);
          char fmt01[8];
          snprintf(fmt01, 8, "%.3ss", format);

          // AUX.jhat.COEFF.m
          if (doaux) {
            for (int idp = 0; idp < 6; idp++) {
              if (dparm[idp] == 0.0) {
                snprintf(ctemp, 32, fmt01, "0.0");
              } else if (dparm[idp] == 1.0) {
                snprintf(ctemp, 32, fmt01, "1.0");
              } else {
                wcsutil_double2str(ctemp, 32, format, dparm[idp]);
              }

              if (idp < 3) {
                snprintf(keyvalue, 96, "'AUX.1.COEFF.%d: %s'", idp%3, ctemp);
                snprintf(comment, 128, "TPD: x = c0 + c1*u + c2*v");
              } else {
                snprintf(keyvalue, 96, "'AUX.2.COEFF.%d: %s'", idp%3, ctemp);
                snprintf(comment, 128, "TPD: y = d0 + d1*u + d2*v");
              }

              wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
                keyvalue, comment, nkeyrec, header, &status);

            }

            dparm += 6;
          }

          // TPD.FWD.m
          for (int idp = 0; idp < iparm[I_TPDNCO]; idp++) {
            if (dparm[idp] == 0.0) continue;

            if (dparm[idp] == 1.0) {
              snprintf(ctemp, 32, fmt01, "1.0");
            } else {
              wcsutil_double2str(ctemp, 32, format, dparm[idp]);
            }

            int m = idp;
            snprintf(keyvalue, 96, "'TPD.FWD.%d:%s %s'", m, (m<10)?" ":"",
              ctemp);
            wcshdo_tpdterm(m, direct, term);
            snprintf(comment, 128, "TPD coefficient: %s", term);

            wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
              keyvalue, comment, nkeyrec, header, &status);
          }

          // CPERRja/CQERRia
          if (dis->maxdis[j] != 0.0) {
            snprintf(keyword, 16,  "C%cERR%d", pq, j+1);
            snprintf(keyvalue, 96, "%20.2f", dis->maxdis[j]);
            snprintf(comment, 128, "%sMaximum absolute value of distortion",
              idis?"":"[pix] ");
            wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
              keyvalue, comment, nkeyrec, header, &status);
          }

          // Inverse distortion function polynomial coefficients.
          if (dis->disx2p[j] == 0x0) continue;

          wcshdo_util(ctrl, "", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
            "", "", nkeyrec, header, &status);

          // TPD.REV.m
          snprintf(keyword, 16,  "D%c%d", pq, j+1);
          for (int idp = iparm[I_TPDNCO]; idp < iparm[I_NDPARM]; idp++) {
            if (dparm[idp] == 0.0) continue;

            wcsutil_double2str(ctemp, 32, format, dparm[idp]);
            int m = idp - iparm[I_TPDNCO];
            snprintf(keyvalue, 96, "'TPD.REV.%d:%s %s'", m, (m<10)?" ":"",
              ctemp);
            wcshdo_tpdterm(m, direct, term);
            snprintf(comment, 128, "TPD coefficient: %s", term);

            wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
              keyvalue, comment, nkeyrec, header, &status);
          }

        } else {
          // General polynomial distortion, handled via its dpkey records
          // since iparm and dparm may hold a translation to TPD.

          // Do auxiliary variables first.
          struct dpkey *ikeyp = dis->dp;
          for (int idp = 0; idp < dis->ndp; idp++, ikeyp++) {
            if (ikeyp->j != j+1) continue;

            cp = strchr(ikeyp->field, '.') + 1;
            if (strncmp(cp, "NAUX", 4) != 0) continue;

            snprintf(keyvalue, 96, "'%s: %d'", cp, dpkeyi(ikeyp));
            wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
              keyvalue, "Number of auxiliary variables", nkeyrec, header,
              &status);

            struct dpkey *jkeyp = dis->dp;
            for (int jdp = 0; jdp < dis->ndp; jdp++, jkeyp++) {
              if (jkeyp->j != j+1) continue;

              keyval = dpkeyd(jkeyp);

              cp = strchr(jkeyp->field, '.') + 1;
              if (strncmp(cp, "AUX.", 4) != 0) continue;

              int m;
              sscanf(cp+4, "%d", &m);
              snprintf(keyvalue, 96, "'%s:", cp);

              cp = strchr(cp+4, '.') + 1;

              char *kv = keyvalue + strlen(keyvalue);
              size_t ksize = kve - kv;
              if ((double)((int)keyval) == keyval) {
                snprintf(kv, ksize, "%4d'", (int)keyval);
              } else if (keyval == 0.5) {
                strncpy(kv, " 0.5'", ksize);
              } else {
                wcsutil_double2str(kv, ksize, "%21.13E", keyval);
                kv += strlen(kv);
                ksize = kve - kv;
                strncpy(kv, "'", ksize);
              }

              int p;
              sscanf(cp+6, "%d", &p);
              if (strncmp(cp, "POWER.", 4) == 0) {
                if (p) {
                  snprintf(comment, 128, "Aux %d: var %d power", m, p);
                } else {
                  snprintf(comment, 128, "Aux %d: power of sum of terms", m);
                }
              } else {
                if (p) {
                  snprintf(comment, 128, "Aux %d: var %d coefficient", m, p);
                } else {
                  snprintf(comment, 128, "Aux %d: offset term", m);
                }
              }

              wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
                keyvalue, comment, nkeyrec, header, &status);
            }

            break;
          }

          // Do polynomial terms.
          ikeyp = dis->dp;
          for (int idp = 0; idp < dis->ndp; idp++, ikeyp++) {
            if (ikeyp->j != j+1) continue;

            cp = strchr(ikeyp->field, '.') + 1;
            if (strncmp(cp, "NTERMS", 6) != 0) continue;

            snprintf(keyvalue, 96, "'%s: %d'", cp, dpkeyi(ikeyp));
            wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
              keyvalue, "Number of terms in the polynomial", nkeyrec, header,
              &status);
          }

          ikeyp = dis->dp;
          for (int idp = 0; idp < dis->ndp; idp++, ikeyp++) {
            if (ikeyp->j != j+1) continue;

            if ((keyval = dpkeyd(ikeyp)) == 0.0) continue;

            cp = strchr(ikeyp->field, '.') + 1;
            if (strncmp(cp, "TERM.", 5) != 0) continue;

            int m;
            sscanf(cp+5, "%d", &m);
            snprintf(keyvalue, 96, "'%s:%s ", cp, (m<10)?" ":"");

            cp = strchr(cp+5, '.') + 1;
            char *kv = keyvalue + strlen(keyvalue);
            size_t ksize = kve - kv;
            if (strncmp(cp, "VAR.", 4) == 0) {
              if ((double)((int)keyval) == keyval) {
                snprintf(kv, ksize, "%20d", (int)keyval);
              } else {
                wcsutil_double2str(kv, ksize, "%20.13f", keyval);
              }

              int p;
              sscanf(cp+4, "%d", &p);
              if (p <= Nhat) {
                snprintf(comment, 128, "Poly term %d: var %d power", m, p);
              } else {
                snprintf(comment, 128, "Poly term %d: aux %d power", m,
                  p-Nhat);
              }

            } else {
              wcsutil_double2str(kv, ksize, "%20.13E", keyval);
              snprintf(comment, 128, "Poly term %d: coefficient", m);
            }

            kv += strlen(kv);
            ksize = kve - kv;
            strncpy(kv, "'", ksize);

            wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
              keyvalue, comment, nkeyrec, header, &status);
          }

          // CPERRja/CQERRia
          if (dis->maxdis[j] != 0.0) {
            snprintf(keyword, 16,  "C%cERR%d", pq, j+1);
            snprintf(keyvalue, 96, "%20.2f", dis->maxdis[j]);
            snprintf(comment, 128, "%sMaximum absolute value of distortion",
              idis?"":"[pix] ");
            wcshdo_util(ctrl, keyword, 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
              keyvalue, comment, nkeyrec, header, &status);
          }
        }
      }
    }

    // DVERRa
    if (dis->totdis != 0.0) {
      snprintf(keyvalue, 96, "%20.2f", dis->totdis);
      snprintf(comment, 128, "Maximum combined distortion");
      wcshdo_util(ctrl, "DVERR", 0x0, 0, 0x0, 0, 0, 0, alt, 0, 0x0,
        keyvalue, comment, nkeyrec, header, &status);
    }
  }


  // Add identification.
  wcshdo_util(ctrl, "", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
    "", "", nkeyrec, header, &status);

  if (dotpd == DIS_DOTPD) {
    // TPD by translation.
    snprintf(comment, 128, "Translated from %s to TPD by WCSLIB %s", tpdsrc,
      wcslib_version(0x0));
  } else {
    snprintf(comment, 128, "WCS header keyrecords produced by WCSLIB %s",
      wcslib_version(0x0));
  }

  wcshdo_util(ctrl, "COMMENT", 0x0, 0, 0x0, 0, 0, 0, ' ', 0, 0x0,
    "", comment, nkeyrec, header, &status);


  if (status == WCSHDRERR_MEMORY) {
    wcserr_set(WCSHDR_ERRMSG(status));
  }
  return status;
}

//----------------------------------------------------------------------------
// Determine a suitable floating point format for a set of parameters.

static void wcshdo_format(
  int fmt,
  int nval,
  const double val[],
  char format[16])

{
  int emax = -999;
  int emin = +999;
  int precision = 0;
  for (int i = 0; i < nval; i++) {
    // Double precision has at least 15 significant digits, and up to 17:
    // http://en.wikipedia.org/wiki/Double-precision_floating-point_format
    char cval[24];
    wcsutil_double2str(cval, 24, "%21.14E", val[i]);

    int cpi = 16;
    while (2 < cpi && cval[cpi] == '0') cpi--;

    // Precision for 'E' format.
    cpi -= 2;
    if (precision < cpi) precision = cpi;

    // Range of significant digits for 'f' format.
    int expon;
    sscanf(cval+18, "%d", &expon);

    if (emax < expon) emax = expon;
    expon -= cpi;
    if (expon < emin) emin = expon;
  }


  if (fmt == 'G') {
    // Because e.g. writing 1e4 as 10000 requires an extra digit.
    emax++;

    if (emin < -15 || 15 < emax || 15 < (emax - emin)) {
      fmt = 'E';
    } else {
      fmt = 'f';
    }
  }

  if (fmt == 'f') {
    precision = -emin;
    if (precision < 1)  precision = 1;
    if (17 < precision) precision = 17;
    snprintf(format, 16, "%%20.%df", precision);

  } else {
    if (precision < 1)  precision = 1;
    if (14 < precision) precision = 14;
    if (precision < 14) {
      snprintf(format, 16, "%%20.%dE", precision);
    } else {
      snprintf(format, 16, "%%21.%dE", precision);
    }
  }
}

//----------------------------------------------------------------------------
// Construct a string that identifies the term of a TPD or TPV polynomial.

static void wcshdo_tpdterm(
  int m,
  int direct,
  char term[16])

{
  const int nTPD[] = {1, 4, 7, 12, 17, 24, 31, 40, 49, 60};

  int degree;
  for (degree = 0; degree <= 9; degree++) {
    if (m < nTPD[degree]) break;
  }

  if (degree == 0) {
    strncpy(term, "1", 16);

  } else {
    int k = degree - (m - nTPD[degree-1]);

    if (k < 0) {
      memcpy(term, "rrrrrrrrr", degree);
    } else if (direct) {
      memcpy(term, "xxxxxxxxx", k);
      memcpy(term+k, "yyyyyyyyy", degree-k);
    } else {
      memcpy(term, "yyyyyyyyy", k);
      memcpy(term+k, "xxxxxxxxx", degree-k);
    }

    term[degree] = '\0';
  }
}

//----------------------------------------------------------------------------
// Construct a keyrecord from the components given.

static void wcshdo_util(
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
  if (*status) return;

  // Reallocate memory in blocks of 2880 bytes.
  char *hptr;
  if ((*nkeyrec)%32 == 0) {
    int nbyte = ((*nkeyrec)/32 + 1) * 2880 + 1;
    if (!(hptr = realloc(*header, nbyte))) {
      *status = WCSHDRERR_MEMORY;
      return;
    }

    *header = hptr;
  }

  // Construct the keyword.
  char keyword[32];
  if (alt == ' ') alt = '\0';
  if (btcol) {
    // Binary table image array.
    if (i > 0 && j) {
      if (j > 0) {
        snprintf(keyword, 32, "%d%d%s%d%c", i, j, tbkey, btcol, alt);
      } else {
        snprintf(keyword, 32, "%d%s%d_%d%c", i, tbkey, btcol, m, alt);
      }
    } else if (i > 0) {
      snprintf(keyword, 32, "%d%s%d%c", i, tbkey, btcol, alt);
    } else if (j > 0) {
      snprintf(keyword, 32, "%d%s%d%c", j, tbkey, btcol, alt);
    } else {
      snprintf(keyword, 32, "%s%d%c", tbkey, btcol, alt);
    }

    if ((strlen(keyword) < 8) && tlkey && (relax & level)) {
      // Use the long form.
      if (i > 0 && j) {
        if (j > 0) {
          snprintf(keyword, 32, "%d%d%s%d%c", i, j, tlkey, btcol, alt);
        } else {
          snprintf(keyword, 32, "%d%s%d_%d%c", i, tlkey, btcol, m, alt);
        }
      } else if (i > 0) {
        snprintf(keyword, 32, "%d%s%d%c", i, tlkey, btcol, alt);
      } else if (j > 0) {
        snprintf(keyword, 32, "%d%s%d%c", j, tlkey, btcol, alt);
      } else {
        snprintf(keyword, 32, "%s%d%c", tlkey, btcol, alt);
      }
    }

  } else if (plcol && plcol[0]) {
    // Pixel list.
    if (i > 0 && j) {
      if (j > 0) {
        snprintf(keyword, 32, "T%s%d_%d%c", tbkey, plcol[i-1], plcol[j-1],
          alt);
      } else {
        snprintf(keyword, 32, "T%s%d_%d%c", tbkey, plcol[i-1], m, alt);
      }
    } else if (i > 0) {
      snprintf(keyword, 32, "T%s%d%c", tbkey, plcol[i-1], alt);
    } else if (j > 0) {
      snprintf(keyword, 32, "T%s%d%c", tbkey, plcol[j-1], alt);
    } else {
      snprintf(keyword, 32, "%s%d%c", tbkey, plcol[0], alt);
    }

    if ((strlen(keyword) < 8) && tlkey && (relax & level)) {
      // Use the long form.
      if (i > 0 && j) {
        if (j > 0) {
          snprintf(keyword, 32, "T%s%d_%d%c", tlkey, plcol[i-1], plcol[j-1],
            alt);
        } else {
          snprintf(keyword, 32, "T%s%d_%d%c", tlkey, plcol[i-1], m, alt);
        }
      } else if (i > 0) {
        snprintf(keyword, 32, "T%s%d%c", tlkey, plcol[i-1], alt);
      } else if (j > 0) {
        snprintf(keyword, 32, "T%s%d%c", tlkey, plcol[j-1], alt);
      } else {
        snprintf(keyword, 32, "%s%d%c", tlkey, plcol[0], alt);
      }
    }
  } else {
    if (i > 0 && j) {
      if (j > 0) {
        snprintf(keyword, 32, "%s%d_%d%c", pikey, i, j, alt);
      } else {
        snprintf(keyword, 32, "%s%d_%d%c", pikey, i, m, alt);
      }
    } else if (i > 0) {
      snprintf(keyword, 32, "%s%d%c", pikey, i, alt);
    } else if (j > 0) {
      snprintf(keyword, 32, "%s%d%c", pikey, j, alt);
    } else {
      snprintf(keyword, 32, "%s%c", pikey, alt);
    }
  }

  // Double-up single-quotes in string keyvalues.
  if (*keyvalue == '\'') {
    hptr = keyvalue + 1;
    while (*hptr) {
      if (*hptr == '\'') {
        char *kptr = hptr++;
        if (*hptr) {
          char ch0 = *kptr;
          while (*kptr) {
            char ch1 = *(++kptr);
            *kptr = ch0;
            ch0 = ch1;
          }
        } else {
          break;
        }
      }

      hptr++;
    }

    // Check length.
    if (strlen(keyvalue) > 70) {
      // Truncate.
      keyvalue[69] = '\'';
      keyvalue[70] = '\0';
    }

  } else {
    // Check length.
    if (strlen(keyvalue) > 70) {
      // Truncate.
      keyvalue[70] = '\0';
    }
  }

  int nc = 47, nv;
  if ((nv = (int)strlen(keyvalue)) > 20) {
    // Rob the keycomment to make space for the keyvalue.
    nc -= (nv - 20);
  }

  hptr = *header + (ptrdiff_t)(80*((*nkeyrec)++));
  if (*keyword == '\0') {
    snprintf(hptr, 81, "%80.80s", " ");
  } else if (strcmp(keyword, "COMMENT") == 0) {
    snprintf(hptr, 81, "%-8.8s %-71.71s", keyword, keycomment);
  } else {
    snprintf(hptr, 81, "%-8.8s= %-20s / %-*.*s", keyword, keyvalue, nc, nc,
      keycomment);
  }
}
