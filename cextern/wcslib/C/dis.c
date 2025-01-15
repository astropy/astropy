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
  $Id: dis.c,v 7.12 2022/09/09 04:57:58 mcalabre Exp $
*===========================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcserr.h"
#include "wcsprintf.h"
#include "wcsutil.h"
#include "dis.h"

const int DISSET = 137;

const int DIS_TPD        =    1;
const int DIS_POLYNOMIAL =    2;
const int DIS_DOTPD      = 1024;

// Maximum number of DPja or DQia keywords.
int NDPMAX = 256;

// Map status return value to message.
const char *dis_errmsg[] = {
  "Success",
  "Null disprm pointer passed",
  "Memory allocation failed",
  "Invalid parameter value",
  "Distort error",
  "De-distort error"};

// Convenience macro for invoking wcserr_set().
#define DIS_ERRMSG(status) WCSERR_SET(status), dis_errmsg[status]

// Internal helper functions, not for general use.
static int polyset(int j, struct disprm *dis);
static int tpdset(int j, struct disprm *dis);

static int pol2tpd(int j, struct disprm *dis);
static int tpvset(int j, struct disprm *dis);
static int sipset(int j, struct disprm *dis);
static int dssset(int j, struct disprm *dis);
static int watset(int j, struct disprm *dis);
static int cheleg(int type, int m, int n, double coeffm[], double coeffn[]);

static int dispoly(DISP2X_ARGS);
static int tpd1(DISP2X_ARGS);
static int tpd2(DISP2X_ARGS);
static int tpd3(DISP2X_ARGS);
static int tpd4(DISP2X_ARGS);
static int tpd5(DISP2X_ARGS);
static int tpd6(DISP2X_ARGS);
static int tpd7(DISP2X_ARGS);
static int tpd8(DISP2X_ARGS);
static int tpd9(DISP2X_ARGS);

// The first three iparm indices have meanings common to all distortion
// functions.  They are used by disp2x(), disx2p(), disprt(), and dishdo().
#define I_DTYPE   0	// Distortion type code.
#define I_NIPARM  1	// Full (allocated) length of iparm[].
#define I_NDPARM  2	// No. of parameters in dparm[], excl. work space.

//----------------------------------------------------------------------------

int disndp(int ndpmax) { if (ndpmax >= 0) NDPMAX = ndpmax; return NDPMAX; }

//----------------------------------------------------------------------------

int dpfill(
  struct dpkey *dp,
  const char *keyword,
  const char *field,
  int j,
  int type,
  int i,
  double f)

{
  if (keyword) {
    if (field) {
      if (j && 2 <= strlen(keyword)) {
        // Fill in the axis number from the value given.
        if (keyword[2] == '\0') {
          sprintf(dp->field, "%s%d.%s", keyword, j, field);
        } else {
          // Take care not to overwrite any alternate code.
          char axno[8];
          sprintf(dp->field, "%s.%s", keyword, field);
          sprintf(axno, "%d", j);
          dp->field[2] = axno[0];
        }

      } else {
        sprintf(dp->field, "%s.%s", keyword, field);
      }
    } else {
      strcpy(dp->field, keyword);
    }
  } else if (field) {
    strcpy(dp->field, field);
  }

  if (j) {
    dp->j = j;
  } else {
    // The field name must either be given or preset.
    char *cp;
    if ((cp = strpbrk(dp->field, "0123456789")) != 0x0) {
      sscanf(cp, "%d.", &(dp->j));
    }
  }

  if ((dp->type = type)) {
    dp->value.f = f;
  } else {
    dp->value.i = i;
  }

  return 0;
}

//----------------------------------------------------------------------------

int dpkeyi(const struct dpkey *dp)

{
  if (dp->type != 0) {
    return (int)dp->value.f;
  }

  return dp->value.i;
}

//----------------------------------------------------------------------------

double dpkeyd(const struct dpkey *dp)

{
  if (dp->type == 0) {
    return (double)dp->value.i;
  }

  return dp->value.f;
}

//----------------------------------------------------------------------------

int disini(int alloc, int naxis, struct disprm *dis)

{
  return disinit(alloc, naxis, dis, -1);
}

//----------------------------------------------------------------------------

int disinit(int alloc, int naxis, struct disprm *dis, int ndpmax)

{
  static const char *function = "disinit";

  // Check inputs.
  if (dis == 0x0) return DISERR_NULL_POINTER;

  if (ndpmax < 0) ndpmax = disndp(-1);


  // Initialize error message handling.
  if (dis->flag == -1) {
    dis->err = 0x0;
  }
  struct wcserr **err = &(dis->err);
  wcserr_clear(err);


  // Initialize pointers.
  if (dis->flag == -1 || dis->m_flag != DISSET) {
    if (dis->flag == -1) {
      dis->docorr = 0x0;
      dis->Nhat   = 0x0;

      dis->axmap  = 0x0;
      dis->offset = 0x0;
      dis->scale  = 0x0;
      dis->iparm  = 0x0;
      dis->dparm  = 0x0;

      dis->disp2x = 0x0;
      dis->disx2p = 0x0;
      dis->tmpmem = 0x0;

      dis->i_naxis = 0;
    }

    // Initialize memory management.
    dis->m_flag   = 0;
    dis->m_naxis  = 0;
    dis->m_dtype  = 0x0;
    dis->m_dp     = 0x0;
    dis->m_maxdis = 0x0;
  }

  if (naxis < 0) {
    return wcserr_set(WCSERR_SET(DISERR_MEMORY),
      "naxis must not be negative (got %d)", naxis);
  }


  // Allocate memory for arrays if required.
  if (alloc ||
      dis->dtype  == 0x0 ||
      (ndpmax && dis->dp == 0x0) ||
      dis->maxdis == 0x0) {

    // Was sufficient allocated previously?
    if (dis->m_flag == DISSET &&
       (dis->m_naxis < naxis  ||
        dis->ndpmax  < ndpmax)) {
      // No, free it.
      disfree(dis);
    }

    if (alloc || dis->dtype == 0x0) {
      if (dis->m_dtype) {
        // In case the caller fiddled with it.
        dis->dtype = dis->m_dtype;

      } else {
        if ((dis->dtype = calloc(naxis, sizeof(char [72]))) == 0x0) {
          disfree(dis);
          return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
        }

        dis->m_flag  = DISSET;
        dis->m_naxis = naxis;
        dis->m_dtype = dis->dtype;
      }
    }

    if (alloc || dis->dp == 0x0) {
      if (dis->m_dp) {
        // In case the caller fiddled with it.
        dis->dp = dis->m_dp;

      } else {
        if (ndpmax) {
          if ((dis->dp = calloc(ndpmax, sizeof(struct dpkey))) == 0x0) {
            disfree(dis);
            return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
          }
        } else {
          dis->dp = 0x0;
        }

        dis->ndpmax  = ndpmax;

        dis->m_flag  = DISSET;
        dis->m_naxis = naxis;
        dis->m_dp    = dis->dp;
      }
    }

    if (alloc || dis->maxdis == 0x0) {
      if (dis->m_maxdis) {
        // In case the caller fiddled with it.
        dis->maxdis = dis->m_maxdis;

      } else {
        if ((dis->maxdis = calloc(naxis, sizeof(double))) == 0x0) {
          disfree(dis);
          return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
        }

        dis->m_flag  = DISSET;
        dis->m_naxis = naxis;
        dis->m_maxdis = dis->maxdis;
      }
    }
  }


  // Set defaults.
  dis->flag  = 0;
  dis->naxis = naxis;

  if (naxis) {
    memset(dis->dtype, 0, naxis*sizeof(char [72]));
  }

  dis->ndp = 0;
  if (ndpmax) {
    memset(dis->dp, 0, ndpmax*sizeof(struct dpkey));
  }

  if (naxis) {
    memset(dis->maxdis, 0, naxis*sizeof(double));
  }
  dis->totdis = 0.0;

  return 0;
}


//----------------------------------------------------------------------------

int discpy(int alloc, const struct disprm *dissrc, struct disprm *disdst)

{
  static const char *function = "discpy";

  if (dissrc == 0x0) return DISERR_NULL_POINTER;
  if (disdst == 0x0) return DISERR_NULL_POINTER;
  struct wcserr **err = &(disdst->err);

  int naxis = dissrc->naxis;
  if (naxis < 1) {
    return wcserr_set(WCSERR_SET(DISERR_MEMORY),
      "naxis must be positive (got %d)", naxis);
  }

  int status;
  if ((status = disinit(alloc, naxis, disdst, dissrc->ndpmax))) {
    return status;
  }

  memcpy(disdst->dtype, dissrc->dtype, naxis*sizeof(char [72]));

  disdst->ndp = dissrc->ndp;
  memcpy(disdst->dp, dissrc->dp, dissrc->ndpmax*sizeof(struct dpkey));

  memcpy(disdst->maxdis, dissrc->maxdis, naxis*sizeof(double));
  disdst->totdis = dissrc->totdis;

  return 0;
}

//----------------------------------------------------------------------------

int disfree(struct disprm *dis)

{
  if (dis == 0x0) return DISERR_NULL_POINTER;

  if (dis->flag != -1) {
    // Optionally allocated by disinit() for given parameters.
    if (dis->m_flag == DISSET) {
      if (dis->dtype  == dis->m_dtype)  dis->dtype  = 0x0;
      if (dis->dp     == dis->m_dp)     dis->dp     = 0x0;
      if (dis->maxdis == dis->m_maxdis) dis->maxdis = 0x0;

      if (dis->m_dtype)  free(dis->m_dtype);
      if (dis->m_dp)     free(dis->m_dp);
      if (dis->m_maxdis) free(dis->m_maxdis);
    }

    // The remainder were allocated by disset().
    if (dis->docorr) free(dis->docorr);
    if (dis->Nhat)   free(dis->Nhat);

    // Recall that axmap, offset, and scale were allocated in bulk.
    if (dis->axmap  && dis->axmap[0])  free(dis->axmap[0]);
    if (dis->offset && dis->offset[0]) free(dis->offset[0]);
    if (dis->scale  && dis->scale[0])  free(dis->scale[0]);

    if (dis->axmap)  free(dis->axmap);
    if (dis->offset) free(dis->offset);
    if (dis->scale)  free(dis->scale);

    if (dis->iparm) {
      for (int j = 0; j < dis->i_naxis; j++) {
        if (dis->iparm[j]) free(dis->iparm[j]);
      }
      free(dis->iparm);
    }

    if (dis->dparm) {
      for (int j = 0; j < dis->i_naxis; j++) {
        if (dis->dparm[j]) free(dis->dparm[j]);
      }
      free(dis->dparm);
    }

    if (dis->disp2x) free(dis->disp2x);
    if (dis->disx2p) free(dis->disx2p);
    if (dis->tmpmem) free(dis->tmpmem);
  }

  dis->m_flag   = 0;
  dis->m_naxis  = 0;
  dis->m_dtype  = 0x0;
  dis->m_dp     = 0x0;
  dis->m_maxdis = 0x0;

  dis->docorr = 0x0;
  dis->Nhat   = 0x0;
  dis->axmap  = 0x0;
  dis->offset = 0x0;
  dis->scale  = 0x0;
  dis->iparm  = 0x0;
  dis->dparm  = 0x0;
  dis->disp2x = 0x0;
  dis->disx2p = 0x0;
  dis->tmpmem = 0x0;

  wcserr_clear(&(dis->err));

  dis->flag = 0;

  return 0;
}

//----------------------------------------------------------------------------

int dissize(const struct disprm *dis, int sizes[2])

{
  if (dis == 0x0) {
    sizes[0] = sizes[1] = 0;
    return DISERR_NULL_POINTER;
  }

  // Base size, in bytes.
  sizes[0] = sizeof(struct disprm);

  // Total size of allocated memory, in bytes.
  sizes[1] = 0;

  int naxis = dis->naxis;

  // disprm::dtype[].
  sizes[1] += naxis * sizeof(char [72]);

  // disprm::dp[].
  sizes[1] += dis->ndpmax * sizeof(struct dpkey);

  // disprm::maxdis[].
  sizes[1] += naxis * sizeof(double);

  // dis::err[].
  int exsizes[2];
  wcserr_size(dis->err, exsizes);
  sizes[1] += exsizes[0] + exsizes[1];

  // The remaining arrays are allocated by disset().
  if (dis->flag != DISSET) {
    return DISERR_SUCCESS;
  }

  // dis::docorr[].
  sizes[1] += naxis * sizeof(int *);

  // dis::Nhat[].
  sizes[1] += naxis * sizeof(int *);

  // dis::axmap[][].
  sizes[1] += naxis * sizeof(int *);
  sizes[1] += naxis*naxis * sizeof(int);

  // dis::offset[][].
  sizes[1] += naxis * sizeof(double *);
  sizes[1] += naxis*naxis * sizeof(double);

  // dis::scale[][].
  sizes[1] += naxis * sizeof(double *);
  sizes[1] += naxis*naxis * sizeof(double);

  // dis::iparm[][].
  sizes[1] += naxis * sizeof(int *);
  for (int j = 0; j < naxis; j++) {
    if (dis->iparm[j]) {
      sizes[1] += dis->iparm[j][I_NIPARM] * sizeof(int);
    }
  }

  // dis::dparm[][].
  sizes[1] += naxis * sizeof(double *);
  for (int j = 0; j < naxis; j++) {
    if (dis->dparm[j]) {
      sizes[1] += dis->dparm[j][I_NDPARM] * sizeof(double);
    }
  }

  // dis::disp2x[].
  sizes[1] += naxis * sizeof(int (*)(DISP2X_ARGS));

  // dis::disx2p[].
  sizes[1] += naxis * sizeof(int (*)(DISX2P_ARGS));

  // dis::tmpmem[].
  sizes[1] += 5*naxis * sizeof(double);

  return DISERR_SUCCESS;
}

//----------------------------------------------------------------------------

int disprt(const struct disprm *dis)

{
  if (dis == 0x0) return DISERR_NULL_POINTER;

  if (dis->flag != DISSET) {
    wcsprintf("The disprm struct is UNINITIALIZED.\n");
    return 0;
  }

  int naxis = dis->naxis;


  wcsprintf("       flag: %d\n", dis->flag);

  // Parameters supplied.
  wcsprintf("      naxis: %d\n", naxis);

  WCSPRINTF_PTR("      dtype: ", dis->dtype, "\n");
  for (int j = 0; j < naxis; j++) {
    wcsprintf("             \"%s\"\n", dis->dtype[j]);
  }

  wcsprintf("        ndp: %d\n", dis->ndp);
  wcsprintf("     ndpmax: %d\n", dis->ndpmax);
  WCSPRINTF_PTR("         dp: ", dis->dp, "\n");
  for (int i = 0; i < dis->ndp; i++) {
    if (dis->dp[i].type) {
      wcsprintf("             %3d%3d  %#- 11.5g  %.32s\n",
        dis->dp[i].j, dis->dp[i].type, dis->dp[i].value.f, dis->dp[i].field);
    } else {
      wcsprintf("             %3d%3d  %11d  %.32s\n",
        dis->dp[i].j, dis->dp[i].type, dis->dp[i].value.i, dis->dp[i].field);
    }
  }

  WCSPRINTF_PTR("     maxdis: ", dis->maxdis, "\n");
  wcsprintf("            ");
  for (int j = 0; j < naxis; j++) {
    wcsprintf("  %#- 11.5g", dis->maxdis[j]);
  }
  wcsprintf("\n");

  wcsprintf("     totdis:  %#- 11.5g\n", dis->totdis);

  // Derived values.
  WCSPRINTF_PTR("     docorr: ", dis->docorr, "\n");
  wcsprintf("            ");
  for (int j = 0; j < naxis; j++) {
    wcsprintf("%6d", dis->docorr[j]);
  }
  wcsprintf("\n");

  WCSPRINTF_PTR("       Nhat: ", dis->Nhat, "\n");
  wcsprintf("            ");
  for (int j = 0; j < naxis; j++) {
    wcsprintf("%6d", dis->Nhat[j]);
  }
  wcsprintf("\n");

  WCSPRINTF_PTR("      axmap: ", dis->axmap, "\n");
  for (int j = 0; j < naxis; j++) {
    wcsprintf(" axmap[%d][]:", j);
    for (int jhat = 0; jhat < naxis; jhat++) {
      wcsprintf("%6d", dis->axmap[j][jhat]);
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("     offset: ", dis->offset, "\n");
  for (int j = 0; j < naxis; j++) {
    wcsprintf("offset[%d][]:", j);
    for (int jhat = 0; jhat < naxis; jhat++) {
      wcsprintf("  %#- 11.5g", dis->offset[j][jhat]);
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("      scale: ", dis->scale, "\n");
  for (int j = 0; j < naxis; j++) {
    wcsprintf(" scale[%d][]:", j);
    for (int jhat = 0; jhat < naxis; jhat++) {
      wcsprintf("  %#- 11.5g", dis->scale[j][jhat]);
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("      iparm: ", dis->iparm, "\n");
  for (int j = 0; j < naxis; j++) {
    wcsprintf(" iparm[%d]  : ", j);
    WCSPRINTF_PTR("", dis->iparm[j], "\n");

    if (dis->iparm[j]) {
      wcsprintf(" iparm[%d][]:", j);
      for (int k = 0; k < dis->iparm[j][I_NIPARM]; k++) {
        if (k && k%5 == 0) {
          wcsprintf("\n            ");
        }
        wcsprintf("  %11d", dis->iparm[j][k]);
      }
      wcsprintf("\n");
    }
  }

  WCSPRINTF_PTR("      dparm: ", dis->dparm, "\n");
  for (int j = 0; j < naxis; j++) {
    wcsprintf(" dparm[%d]  : ", j);
    WCSPRINTF_PTR("", dis->dparm[j], "\n");

    if (dis->dparm[j]) {
      wcsprintf(" dparm[%d][]:", j);
      for (int k = 0; k < dis->iparm[j][I_NDPARM]; k++) {
        if (k && k%5 == 0) {
          wcsprintf("\n            ");
        }
        wcsprintf("  %#- 11.5g", dis->dparm[j][k]);
      }
      wcsprintf("\n");
    }
  }

  wcsprintf("    i_naxis: %d\n", dis->i_naxis);
  wcsprintf("       ndis: %d\n", dis->ndis);

  // Error handling.
  WCSPRINTF_PTR("        err: ", dis->err, "\n");
  if (dis->err) {
    wcserr_prt(dis->err, "             ");
  }

  // Work arrays.
  char hext[32];
  WCSPRINTF_PTR("     disp2x: ", dis->disp2x, "\n");
  for (int j = 0; j < naxis; j++) {
    wcsprintf("  disp2x[%d]: %s", j,
      wcsutil_fptr2str((void (*)(void))dis->disp2x[j], hext));
    if (dis->disp2x[j] == dispoly) {
      wcsprintf("  (= dispoly)\n");
    } else if (dis->disp2x[j] == tpd1) {
      wcsprintf("  (= tpd1)\n");
    } else if (dis->disp2x[j] == tpd2) {
      wcsprintf("  (= tpd2)\n");
    } else if (dis->disp2x[j] == tpd3) {
      wcsprintf("  (= tpd3)\n");
    } else if (dis->disp2x[j] == tpd4) {
      wcsprintf("  (= tpd4)\n");
    } else if (dis->disp2x[j] == tpd5) {
      wcsprintf("  (= tpd5)\n");
    } else if (dis->disp2x[j] == tpd6) {
      wcsprintf("  (= tpd6)\n");
    } else if (dis->disp2x[j] == tpd7) {
      wcsprintf("  (= tpd7)\n");
    } else if (dis->disp2x[j] == tpd8) {
      wcsprintf("  (= tpd8)\n");
    } else if (dis->disp2x[j] == tpd9) {
      wcsprintf("  (= tpd9)\n");
    } else {
      wcsprintf("\n");
    }
  }
  WCSPRINTF_PTR("     disx2p: ", dis->disx2p, "\n");
  for (int j = 0; j < naxis; j++) {
    wcsprintf("  disx2p[%d]: %s\n", j,
      wcsutil_fptr2str((void (*)(void))dis->disx2p[j], hext));
  }
  WCSPRINTF_PTR("     tmpmem: ", dis->tmpmem, "\n");

  // Memory management.
  wcsprintf("     m_flag: %d\n", dis->m_flag);
  wcsprintf("    m_naxis: %d\n", dis->m_naxis);
  WCSPRINTF_PTR("    m_dtype: ", dis->m_dtype, "");
  if (dis->m_dtype  == dis->dtype)  wcsprintf("  (= dtype)");
  wcsprintf("\n");
  WCSPRINTF_PTR("       m_dp: ", dis->m_dp, "");
  if (dis->m_dp     == dis->dp)     wcsprintf("  (= dp)");
  wcsprintf("\n");
  WCSPRINTF_PTR("   m_maxdis: ", dis->m_maxdis, "");
  if (dis->m_maxdis == dis->maxdis) wcsprintf("  (= maxdis)");
  wcsprintf("\n");

  return 0;
}

//----------------------------------------------------------------------------

int disperr(const struct disprm *dis, const char *prefix)

{
  if (dis == 0x0) return DISERR_NULL_POINTER;

  if (dis->err) {
    wcserr_prt(dis->err, prefix);
  }

  return 0;
}

//----------------------------------------------------------------------------

int dishdo(struct disprm *dis)

{
  static const char *function = "dishdo";

  if (dis == 0x0) return DISERR_NULL_POINTER;
  struct wcserr **err = &(dis->err);

  int status = 0;
  for (int j = 0; j < dis->naxis; j++) {
    if (dis->iparm[j][I_DTYPE]) {
      if (dis->iparm[j][I_DTYPE] == DIS_TPD) {
        // Implemented as TPD...
        if (strcmp(dis->dtype[j], "TPD") != 0) {
          // ... but isn't TPD.
          dis->iparm[j][I_DTYPE] |= DIS_DOTPD;
        }
      } else {
        // Must be a Polynomial that can't be implemented as TPD.
        status = wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Translation of %s to TPD is not possible", dis->dtype[j]);
      }
    }
  }

  return status;
}

//----------------------------------------------------------------------------

int disset(struct disprm *dis)

{
  static const char *function = "disset";

  if (dis == 0x0) return DISERR_NULL_POINTER;
  struct wcserr **err = &(dis->err);

  int naxis = dis->naxis;


  // Do basic checks.
  if (dis->ndp < 0) {
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "disprm::ndp is negative (%d)", dis->ndp);
  }

  int ndis = 0;
  for (int j = 0; j < naxis; j++) {
    if (strlen(dis->dtype[j])) {
      ndis++;
      break;
    }
  }

  char *dpq;
  if (dis->ndp) {
    // Is it prior or sequent?
    if (dis->dp[0].field[1] == 'P') {
      dpq = "DPja";
    } else if (dis->dp[0].field[1] == 'Q') {
      dpq = "DQia";
    } else {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "disprm::dp[0].field (%s) is invalid", dis->dp[0].field);
    }

  } else {
    if (ndis) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "No DPja or DQia keywords, NAXES at least is required for each "
        "distortion");
    }

    // A Clayton's distortion.  Avert compiler warnings about possible use of
    // uninitialized variables.
    dpq = 0x0;
  }


  // Free memory allocated separately for each axis.
  for (int j = 0; j < dis->i_naxis; j++) {
    if (dis->iparm[j]) free(dis->iparm[j]);
    if (dis->dparm[j]) free(dis->dparm[j]);
    dis->iparm[j] = 0x0;
    dis->dparm[j] = 0x0;
  }

  // Allocate or reallocate memory, if necessary, for derived parameter and
  // work arrays sized according to the number of axes.
  if (dis->i_naxis < naxis) {
    if (dis->i_naxis) {
      free(dis->docorr);
      free(dis->Nhat);

      // Noting that axmap, offset, and scale are allocated in bulk.
      free(dis->axmap[0]);
      free(dis->axmap);
      free(dis->offset[0]);
      free(dis->offset);
      free(dis->scale[0]);
      free(dis->scale);

      free(dis->iparm);
      free(dis->dparm);

      free(dis->disp2x);
      free(dis->disx2p);

      free(dis->tmpmem);
    }

    if ((dis->docorr = calloc(naxis, sizeof(int *))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    if ((dis->Nhat = calloc(naxis, sizeof(int *))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    // Allocate axmap[][] in bulk and then carve it up.
    if ((dis->axmap = calloc(naxis, sizeof(int *))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    if ((dis->axmap[0] = calloc(naxis*naxis, sizeof(int))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    for (int j = 1; j < naxis; j++) {
      dis->axmap[j] = dis->axmap[j-1] + naxis;
    }

    // Allocate offset[][] in bulk and then carve it up.
    if ((dis->offset = calloc(naxis, sizeof(double *))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    if ((dis->offset[0] = calloc(naxis*naxis, sizeof(double))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    for (int j = 1; j < naxis; j++) {
      dis->offset[j] = dis->offset[j-1] + naxis;
    }

    // Allocate scale[][] in bulk and then carve it up.
    if ((dis->scale = calloc(naxis, sizeof(double *))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    if ((dis->scale[0] = calloc(naxis*naxis, sizeof(double))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    for (int j = 1; j < naxis; j++) {
      dis->scale[j] = dis->scale[j-1] + naxis;
    }

    if ((dis->iparm = calloc(naxis, sizeof(int *))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    if ((dis->dparm = calloc(naxis, sizeof(double *))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    if ((dis->disp2x = calloc(naxis, sizeof(int (*)(DISP2X_ARGS)))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    if ((dis->disx2p = calloc(naxis, sizeof(int (*)(DISX2P_ARGS)))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    if ((dis->tmpmem = calloc(5*naxis, sizeof(double))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    dis->i_naxis = naxis;
  }

  // Start with a clean slate.
  for (int j = 0; j < naxis; j++) {
    dis->docorr[j] = 1;
  }

  memset(dis->Nhat, 0, naxis*sizeof(int));

  for (int jhat = 0; jhat < naxis*naxis; jhat++) {
    dis->axmap[0][jhat] = -1;
  }

  memset(dis->offset[0], 0, naxis*naxis*sizeof(double));

  for (int jhat = 0; jhat < naxis*naxis; jhat++) {
    dis->scale[0][jhat] = 1.0;
  }

  // polyset() etc. must look after iparm[][] and dparm[][].

  dis->i_naxis = naxis;
  dis->ndis    = 0;

  memset(dis->disp2x, 0, naxis*sizeof(int (*)(DISP2X_ARGS)));
  memset(dis->disx2p, 0, naxis*sizeof(int (*)(DISX2P_ARGS)));
  memset(dis->tmpmem, 0, naxis*sizeof(double));


  // Handle DPja or DQia keywords common to all distortions.
  struct dpkey *keyp = dis->dp;
  for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
    // Check that they're all one kind or the other.
    if (keyp->field[1] != dpq[1]) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "disprm::dp appears to contain a mix of DPja and DQia keys");
    }

    int j = keyp->j;

    if (j < 1 || naxis < j) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid axis number (%d) in %s", j, keyp->field);
    }

    char *fp;
    if ((fp = strchr(keyp->field, '.')) == 0x0) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid record field name: %s", j, keyp->field);
    }
    fp++;

    // Convert to 0-relative axis number.
    j--;

    if (strncmp(fp, "DOCORR", 7) == 0) {
      if (dpkeyi(keyp) == 0) {
        dis->docorr[j] = 0;
      }

    } else if (strncmp(fp, "NAXES", 6) == 0) {
      int Nhat = dpkeyi(keyp);
      if (Nhat < 0 || naxis < Nhat) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid value of Nhat for %s distortion in %s: %d", dis->dtype[j],
          keyp->field, Nhat);
      }

      dis->Nhat[j] = Nhat;

    } else if (strncmp(fp, "AXIS.", 5) == 0) {
      int jhat;
      sscanf(fp+5, "%d", &jhat);
      if (jhat < 1 || naxis < jhat) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid axis in axis map for %s distortion in %s: %d",
          dis->dtype[j], keyp->field, jhat);
      }

      // N.B. axis numbers in the map are 0-relative.
      dis->axmap[j][jhat-1] = dpkeyi(keyp) - 1;

    } else if (strncmp(fp, "OFFSET.", 7) == 0) {
      int jhat;
      sscanf(fp+7, "%d", &jhat);
      dis->offset[j][jhat-1] = dpkeyd(keyp);

    } else if (strncmp(fp, "SCALE.", 6) == 0) {
      int jhat;
      sscanf(fp+6, "%d", &jhat);
      dis->scale[j][jhat-1] = dpkeyd(keyp);
    }
  }

  // Set defaults and do sanity checks on axmap[][].
  for (int j = 0; j < naxis; j++) {
    if (strlen(dis->dtype[j]) == 0) {
      // No distortion on this axis, check that there are no parameters.
      keyp = dis->dp;
      for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
        if (keyp->j == j+1) {
          return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
            "No distortion type, yet %s keyvalues are present for axis %d",
            dpq, j+1);
        }
      }

      continue;
    }

    // N.B. NAXES (Nhat) has no default value.
    if (dis->Nhat[j] <= 0) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "%s.NAXES was not set (or bad) for %s distortion on axis %d",
        dpq, dis->dtype[j], j+1);
    }

    // Set defaults for axmap[][].
    int Nhat = dis->Nhat[j];
    for (int jhat = 0; jhat < Nhat; jhat++) {
      if (dis->axmap[j][jhat] == -1) {
        dis->axmap[j][jhat] = jhat;
      }
    }

    // Sanity check on the length of the axis map.
    Nhat = 0;
    for (int jhat = 0; jhat < naxis; jhat++) {
      if (dis->axmap[j][jhat] != -1) Nhat = jhat+1;
    }

    if (Nhat != dis->Nhat[j]) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Mismatch in length of axis map for %s distortion on axis %d",
        dis->dtype[j], j+1);
    }

    // Check uniqueness of entries in the axis map.
    for (int jhat = 0; jhat < Nhat; jhat++) {
      for (int k = 0; k < jhat; k++) {
        if (dis->axmap[j][jhat] == dis->axmap[j][k]) {
          return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
            "Duplicated entry in axis map for %s distortion on axis %d",
            dis->dtype[j], j+1);
        }
      }
    }
  }


  // Identify the distortion functions.
  ndis = 0;
  for (int j = 0; j < naxis; j++) {
    if (strlen(dis->dtype[j]) == 0) {
      // No distortion on this axis.
      continue;
    }

    if (dis->Nhat[j] == 0) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Empty axis map for %s distortion on axis %d", dis->dtype[j], j+1);
    }

    // Invoke the specific setup functions for each distortion.
    int status;
    if (strcmp(dis->dtype[j], "TPD") == 0) {
      // Template Polynomial Distortion.
      if ((status = tpdset(j, dis))) {
        // (Preserve the error message set by tpdset().)
        return status;
      }

    } else if (strcmp(dis->dtype[j], "TPV") == 0) {
      // TPV "projection".
      if ((status = tpvset(j, dis))) {
        // (Preserve the error message set by tpvset().)
        return status;
      }

    } else if (strcmp(dis->dtype[j], "SIP") == 0) {
      // Simple Imaging Polynomial (SIP).
      if ((status = sipset(j, dis))) {
        // (Preserve the error message set by sipset().)
        return status;
      }

    } else if (strcmp(dis->dtype[j], "DSS") == 0) {
      // Digitized Sky Survey (DSS).
      if ((status = dssset(j, dis))) {
        // (Preserve the error message set by dssset().)
        return status;
      }

    } else if (strncmp(dis->dtype[j], "WAT", 3) == 0) {
      // WAT (TNX or ZPX "projections").
      if ((status = watset(j, dis))) {
        // (Preserve the error message set by watset().)
        return status;
      }

    } else if (strcmp(dis->dtype[j], "Polynomial")  == 0 ||
               strcmp(dis->dtype[j], "Polynomial*") == 0) {
      // General polynomial distortion.
      if ((status = polyset(j, dis))) {
        // (Preserve the error message set by polyset().)
        return status;
      }

    } else {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Unrecognized/unimplemented distortion function: %s", dis->dtype[j]);
    }

    ndis++;
  }

  dis->ndis = ndis;
  dis->flag = DISSET;

  return 0;
}

//----------------------------------------------------------------------------

int disp2x(
  struct disprm *dis,
  const double rawcrd[],
  double discrd[])

{
  static const char *function = "disp2x";

  // Initialize.
  if (dis == 0x0) return DISERR_NULL_POINTER;
  struct wcserr **err = &(dis->err);

  if (dis->flag != DISSET) {
    int status;
    if ((status = disset(dis))) return status;
  }

  int naxis = dis->naxis;


  // Invoke the distortion functions for each axis.
  double *tmpcrd = dis->tmpmem;
  for (int j = 0; j < naxis; j++) {
    if (dis->disp2x[j]) {
      double *offset = dis->offset[j];
      double *scale  = dis->scale[j];

      int Nhat = dis->Nhat[j];
      for (int jhat = 0; jhat < Nhat; jhat++) {
        int axisj = dis->axmap[j][jhat];
        tmpcrd[jhat] = (rawcrd[axisj] - offset[jhat])*scale[jhat];
      }

      double dtmp;
      if ((dis->disp2x[j])(0, dis->iparm[j], dis->dparm[j], Nhat, tmpcrd,
                           &dtmp)) {
        return wcserr_set(DIS_ERRMSG(DISERR_DISTORT));
      }

      if (dis->docorr[j]) {
        // Distortion function computes a correction to be applied.
        discrd[j] = rawcrd[j] + dtmp;
      } else {
        // Distortion function computes corrected coordinates directly.
        discrd[j] = dtmp;
      }

    } else {
      discrd[j] = rawcrd[j];
    }
  }

  return 0;
}

//----------------------------------------------------------------------------
// This function is intended for debugging purposes only.
// No documentation or prototype is provided in dis.h.

int disitermax(int itermax)

{
  static int ITERMAX = 30;

  if (itermax >= 0) {
    ITERMAX = itermax;
  }

  return ITERMAX;
}

//----------------------------------------------------------------------------

int disx2p(
  struct disprm *dis,
  const double discrd[],
  double rawcrd[])

{
  static const char *function = "disx2p";

  const double TOL = 1.0e-13;

  int status;

  // Initialize.
  if (dis == 0x0) return DISERR_NULL_POINTER;
  struct wcserr **err = &(dis->err);

  int naxis = dis->naxis;

  // Carve up working memory, noting that disp2x() gets to it first.
  double *dcrd0 = dis->tmpmem + naxis;
  double *dcrd1 = dcrd0 + naxis;
  double *rcrd1 = dcrd1 + naxis;
  double *delta = rcrd1 + naxis;


  // Zeroth approximation.  The assumption here and below is that the
  // distortion is small so that, to first order in the neighbourhood of
  // the solution, discrd[j] ~= a + b*rawcrd[j], i.e. independent of
  // rawcrd[i], where i != j.  This is effectively equivalent to assuming
  // that the distortion functions are separable to first order.
  // Furthermore, a is assumed to be small, and b close to unity.
  memcpy(rawcrd, discrd, naxis*sizeof(double));

  // If available, use disprm::disx2p to improve the zeroth approximation.
  for (int j = 0; j < naxis; j++) {
    if (dis->disx2p[j]) {
      double *offset = dis->offset[j];
      double *scale  = dis->scale[j];
      double *tmpcrd = dis->tmpmem;

      int Nhat = dis->Nhat[j];
      for (int jhat = 0; jhat < Nhat; jhat++) {
        int axisj = dis->axmap[j][jhat];
        tmpcrd[jhat] = (discrd[axisj] - offset[jhat])*scale[jhat];
      }

      double rtmp;
      if ((status = (dis->disx2p[j])(1, dis->iparm[j], dis->dparm[j], Nhat,
                                     tmpcrd, &rtmp))) {
        return wcserr_set(DIS_ERRMSG(DISERR_DEDISTORT));
      }

      if (dis->docorr[j]) {
        // Inverse distortion function computes a correction to be applied.
        rawcrd[j] = discrd[j] + rtmp;
      } else {
        // Inverse distortion function computes corrected coordinates directly.
        rawcrd[j] = rtmp;
      }
    }
  }

  // Quick return debugging hook, assumes inverse functions were defined.
  int itermax;
  if ((itermax = disitermax(-1)) == 0) {
    return 0;
  }


  // Iteratively invert the (well-behaved!) distortion function.
  int convergence, iter;
  for (iter = 0; iter < itermax; iter++) {
    if ((status = disp2x(dis, rawcrd, dcrd0))) {
      return wcserr_set(DIS_ERRMSG(status));
    }

    // Check for convergence.
    convergence = 1;
    for (int j = 0; j < naxis; j++) {
      delta[j] = discrd[j] - dcrd0[j];

      double dd;
      if (fabs(discrd[j]) < 1.0) {
        dd = delta[j];
      } else {
        // TOL may be below the precision achievable from floating point
        // subtraction, so switch to a fractional tolerance.
        dd = delta[j] / discrd[j];
      }

      if (TOL < fabs(dd)) {
        // No convergence yet on this axis.
        convergence = 0;
      }
    }

    if (convergence) break;

    // Determine a suitable test point for computing the gradient.
    for (int j = 0; j < naxis; j++) {
      // Constrain the displacement.
      delta[j] /= 2.0;
      if (fabs(delta[j]) < 1.0e-6) {
        if (delta[j] < 0.0) {
          delta[j] = -1.0e-6;
        } else {
          delta[j] =  1.0e-6;
        }
      } else if (1.0 < fabs(delta[j])) {
        if (delta[j] < 0.0) {
          delta[j] = -1.0;
        } else {
          delta[j] =  1.0;
        }
      }
    }

    if (iter < itermax/2) {
      // With the assumption of small distortions (as above), the gradient
      // of discrd[j] should be dominated by the partial derivative with
      // respect to rawcrd[j], and we can neglect partials with respect
      // to rawcrd[i], where i != j.  Thus only one test point is needed,
      // not one for each axis.
      for (int j = 0; j < naxis; j++) {
        rcrd1[j] = rawcrd[j] + delta[j];
      }

      // Compute discrd[] at the test point.
      if ((status = disp2x(dis, rcrd1, dcrd1))) {
        return wcserr_set(DIS_ERRMSG(status));
      }

      // Compute the next approximation.
      for (int j = 0; j < naxis; j++) {
        rawcrd[j] += (discrd[j] - dcrd0[j]) *
                        (delta[j]/(dcrd1[j] - dcrd0[j]));
      }

    } else {
      // Convergence should not take more than seven or so iterations.  As
      // it is slow, try computing the gradient in full.
      memcpy(rcrd1, rawcrd, naxis*sizeof(double));

      for (int j = 0; j < naxis; j++) {
        rcrd1[j] += delta[j];

        // Compute discrd[] at the test point.
        if ((status = disp2x(dis, rcrd1, dcrd1))) {
          return wcserr_set(DIS_ERRMSG(status));
        }

        // Compute the next approximation.
        rawcrd[j] += (discrd[j] - dcrd0[j]) *
                       (delta[j]/(dcrd1[j] - dcrd0[j]));

        rcrd1[j] -= delta[j];
      }
    }
  }


  if (!convergence) {
    double residual = 0.0;
    for (int j = 0; j < naxis; j++) {
      double dd = discrd[j] - dcrd0[j] ;
      residual += dd*dd;
    }
    residual = sqrt(residual);

    return wcserr_set(WCSERR_SET(DISERR_DEDISTORT),
      "Convergence not achieved after %d iterations, residual %#7.2g", iter,
        residual);
  }


  return 0;
}

//----------------------------------------------------------------------------

int diswarp(
  struct disprm *dis,
  const double pixblc[],
  const double pixtrc[],
  const double pixsamp[],
  int    *nsamp,
  double maxdis[],
  double *maxtot,
  double avgdis[],
  double *avgtot,
  double rmsdis[],
  double *rmstot)

{
  static const char *function = "diswarp";

  int status = 0;

  // Initialize.
  if (dis == 0x0) return DISERR_NULL_POINTER;
  struct wcserr **err = &(dis->err);

  int naxis = dis->naxis;

  if (nsamp) *nsamp = 0;
  for (int j = 0; j < naxis; j++) {
    if (maxdis) maxdis[j] = 0.0;
    if (avgdis) avgdis[j] = 0.0;
    if (rmsdis) rmsdis[j] = 0.0;
  }
  if (maxtot) *maxtot = 0.0;
  if (avgtot) *avgtot = 0.0;
  if (rmstot) *rmstot = 0.0;

  // Quick return if no distortions.
  if (dis->ndis == 0) return 0;

  // Carve up working memory, noting that disp2x() gets to it first.
  double *pixinc = dis->tmpmem + naxis;
  double *pixend = pixinc + naxis;
  double *sumdis = pixend + naxis;
  double *ssqdis = sumdis + naxis;

  // Work out increments on each axis.
  for (int j = 0; j < naxis; j++) {
    double pixspan = pixtrc[j] - (pixblc ? pixblc[j] : 1.0);

    if (pixsamp == 0x0) {
      pixinc[j] = 1.0;
    } else if (pixsamp[j] == 0.0) {
      pixinc[j] = 1.0;
    } else if (pixsamp[j] > 0.0) {
      pixinc[j] = pixsamp[j];
    } else if (pixsamp[j] > -1.5) {
      pixinc[j] = 2.0*pixspan;
    } else {
      pixinc[j] = pixspan / ((int)(-pixsamp[j] - 0.5));
    }
  }

  // Get some more memory for coordinate vectors.
  double *pix0, *pix1;
  if ((pix0 = calloc(2*naxis, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  pix1 = pix0 + naxis;


  // Set up the array of pixel coordinates.
  for (int j = 0; j < naxis; j++) {
    pix0[j] = pixblc ? pixblc[j] : 1.0;
    pixend[j] = pixtrc[j] + 0.5*pixinc[j];
  }

  // Initialize accumulators.
  for (int j = 0; j < naxis; j++) {
    sumdis[j] = 0.0;
    ssqdis[j] = 0.0;
  }
  double sumtot = 0.0;
  double ssqtot = 0.0;


  // Loop over N dimensions.
  int carry = 0;
  while (carry == 0) {
    if ((status = disp2x(dis, pix0, pix1))) {
      // (Preserve the error message set by disp2x().)
      goto cleanup;
    }

    // Accumulate statistics.
    (*nsamp)++;

    double dssq = 0.0;
    for (int j = 0; j < naxis; j++) {
      double dpix = pix1[j] - pix0[j];
      double dpx2 = dpix*dpix;

      sumdis[j] += dpix;
      ssqdis[j] += dpx2;

      if (maxdis && (dpix = fabs(dpix)) > maxdis[j]) {
        maxdis[j] = dpix;
      }

      dssq += dpx2;
    }

    double totdis = sqrt(dssq);
    sumtot += totdis;
    ssqtot += totdis*totdis;

    if (maxtot && *maxtot < totdis) {
      *maxtot = totdis;
    }

    // Next pixel.
    for (int j = 0; j < naxis; j++) {
      pix0[j] += pixinc[j];
      if (pix0[j] < pixend[j]) {
        carry = 0;
        break;
      }

      pix0[j] = pixblc ? pixblc[j] : 1.0;
      carry = 1;
    }
  }


  // Compute the means and RMSs.
  for (int j = 0; j < naxis; j++) {
    ssqdis[j] /= *nsamp;
    sumdis[j] /= *nsamp;
    if (avgdis) avgdis[j] = sumdis[j];
    if (rmsdis) rmsdis[j] = sqrt(ssqdis[j] - sumdis[j]*sumdis[j]);
  }

  ssqtot /= *nsamp;
  sumtot /= *nsamp;
  if (avgtot) *avgtot = sumtot;
  if (rmstot) *rmstot = sqrt(ssqtot - sumtot*sumtot);


cleanup:
  free(pix0);

  return status;
}

//----------------------------------------------------------------------------

int polyset(int j, struct disprm *dis)

{
  static const char *function = "polyset";

  // Initialize.
  if (dis == 0x0) return DISERR_NULL_POINTER;
  struct wcserr **err = &(dis->err);

  int naxis = dis->naxis;

  char   id[32];
  sprintf(id, "Polynomial on axis %d", j+1);


  // Find the number of auxiliary variables and terms.
  int K = 0;
  int M = 0;
  struct dpkey *keyp = dis->dp;
  for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    char *fp;
    if ((fp = strchr(keyp->field, '.')) == 0x0) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid field name for %s: %s", id, keyp->field);
    }
    fp++;

    if (strcmp(fp, "NAUX") == 0) {
      K = dpkeyi(keyp);
    } else if (strcmp(fp, "NTERMS") == 0) {
      M = dpkeyi(keyp);
    }
  }

  if (K < 0) {
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "Invalid number of auxiliaries (%d) for %s", K, id);
  }

  if (M <= 0) {
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "Invalid number of terms (%d) for %s", M, id);
  }

  int Nhat = dis->Nhat[j];
  int nKparm = 2*(Nhat + 1);
  int nVar   = Nhat + K;
  int nTparm = 1 + nVar;
  int ndparm = K*nKparm + M*nTparm;

// These iparm indices are specific to Polynomial.
#define I_NIDX    3	// No. of indexes in iparm[].
#define I_LENDP   4	// Full (allocated) length of dparm[].
#define I_K       5	// No. of auxiliary variables.
#define I_M       6	// No. of terms in the polynomial.
#define I_NKPARM  7	// No. of parameters used to define each auxiliary.
#define I_NTPARM  8	// No. of parameters used to define each term.
#define I_NVAR    9	// No. of independent + auxiliary variables.
#define I_MNVAR  10	// No. of powers (exponents) in the polynomial.
#define I_DPOLY  11	// dparm offset for polynomial coefficients.
#define I_DAUX   12	// dparm offset for auxiliary coefficients.
#define I_DVPOW  13	// dparm offset for integral powers of variables.
#define I_MAXPOW 14	// iparm offset for max powers.
#define I_DPOFF  15	// iparm offset for dparm offsets.
#define I_FLAGS  16	// iparm offset for flags.
#define I_IPOW   17	// iparm offset for integral powers.
#define I_NPOLY  18

  // Add extra for handling integer exponents.  See "Optimization" below.
  int niparm = I_NPOLY + (2 + 2*M)*nVar;

  // Add extra memory for temporaries.
  int lendp = ndparm + K;

  // Allocate memory for the indexes and parameter array.
  if ((dis->iparm[j] = calloc(niparm, sizeof(int))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  if ((dis->dparm[j] = calloc(lendp, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  // These help a bit to stop the code from turning into hieroglyphics.
  int    *iparm = dis->iparm[j];
  double *dparm = dis->dparm[j];


  // Record the indexing parameters.  The first three are more widely used.
  iparm[I_DTYPE]  = DIS_POLYNOMIAL;
  iparm[I_NIPARM] = niparm;
  iparm[I_NDPARM] = ndparm;

  iparm[I_NIDX]   = I_NPOLY;
  iparm[I_LENDP]  = lendp;
  iparm[I_K]      = K;
  iparm[I_M]      = M;
  iparm[I_NKPARM] = nKparm;
  iparm[I_NTPARM] = nTparm;
  iparm[I_NVAR]   = nVar;
  iparm[I_MNVAR]  = M*nVar;
  iparm[I_DPOLY]  = K*nKparm;
  iparm[I_DAUX]   = ndparm;
  iparm[I_DVPOW]  = ndparm + K;
  iparm[I_MAXPOW] = iparm[I_NIDX];
  iparm[I_DPOFF]  = iparm[I_MAXPOW] + nVar;
  iparm[I_FLAGS]  = iparm[I_DPOFF]  + nVar;
  iparm[I_IPOW]   = iparm[I_FLAGS]  + M*nVar;

  // Set default values of POWER for the auxiliary variables.
  double *dptr = dparm + (1 + Nhat);
  for (int k = 0; k < K; k++) {
    for (int jhat = 0; jhat <= Nhat; jhat++) {
      dptr[jhat] = 1.0;
    }
    dptr += nKparm;
  }

  // Set default values of COEFF for the independent variables.
  dptr = dparm + iparm[I_DPOLY];
  for (int m = 0; m < M; m++) {
    *dptr = 1.0;
    dptr += nTparm;
  }

  // Extract parameter values from DPja or DQia.
  int i, k, m;
  k = m = 0;
  keyp = dis->dp;
  for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
    // N.B. keyp->j is 1-relative, but j is 0-relative.
    if (keyp->j-1 != j) continue;

    char *fp = strchr(keyp->field, '.') + 1;

    if (strncmp(fp, "AUX.", 4) == 0) {
      // N.B. k here is 1-relative.
      fp += 4;
      sscanf(fp, "%d", &k);
      if (k < 1 || K < k) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Bad auxiliary variable (%d) for %s: %s", k, id, keyp->field);
      }

      if ((fp = strchr(fp, '.')) == 0x0) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid field name for %s: %s", id, keyp->field);
      }
      fp++;

      int offset;
      if (strncmp(fp, "COEFF.", 6) == 0) {
        offset = 0;

      } else if (strncmp(fp, "POWER.", 6) == 0) {
        offset = 1 + Nhat;

      } else {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Unrecognized field name for %s: %s", id, keyp->field);
      }

      fp += 6;
      int jhat;
      sscanf(fp, "%d", &jhat);
      if (jhat < 0 || naxis < jhat) {
        // N.B. jhat == 0 is ok.
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid axis number (%d) for %s: %s", jhat, id, keyp->field);
      }

      i = (k-1)*nKparm + offset + jhat;
      dparm[i] = dpkeyd(keyp);

    } else if (strncmp(fp, "TERM.", 5) == 0) {
      // N.B. m here is 1-relative.
      fp += 5;
      sscanf(fp, "%d", &m);
      if (m < 1 || M < m) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Bad term (%d) for %s: %s", m, id, keyp->field);
      }

      if ((fp = strchr(fp, '.')) == 0x0) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid field name for %s: %s", id, keyp->field);
      }
      fp++;

      if (strcmp(fp, "COEFF") == 0) {
        i = iparm[I_DPOLY] + (m-1)*nTparm;
        dparm[i] = dpkeyd(keyp);

      } else if (strncmp(fp, "VAR.", 4) == 0) {
        // N.B. jhat here is 1-relative.
        fp += 4;
        int jhat;
        sscanf(fp, "%d", &jhat);
        if (jhat < 1 || naxis < jhat) {
          return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid axis number (%d) for %s: %s", jhat, id, keyp->field);
        }

        i = iparm[I_DPOLY] + (m-1)*nTparm + 1 + (jhat-1);
        double power = dpkeyd(keyp);
        dparm[i] = power;

      } else if (strncmp(fp, "AUX.", 4) == 0) {
        // N.B. k here is 1-relative.
        fp += 4;
        sscanf(fp, "%d", &k);
        if (k < 1 || K < k) {
          return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
            "Bad auxiliary variable (%d) for %s: %s", k, id, keyp->field);
        }

        i = iparm[I_DPOLY] + (m-1)*nTparm + 1 + Nhat + (k-1);
        double power = dpkeyd(keyp);
        dparm[i] = power;

      } else {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Unrecognized field name for %s: %s", id, keyp->field);
      }

    } else if (strcmp(fp, "DOCORR") &&
               strcmp(fp, "NAXES")  &&
              strncmp(fp, "AXIS.",   5) &&
              strncmp(fp, "OFFSET.", 7) &&
              strncmp(fp, "SCALE.",  6) &&
               strcmp(fp, "NAUX")   &&
               strcmp(fp, "NTERMS")) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Unrecognized field name for %s: %s", id, keyp->field);
    }
  }


  // Optimization: when the power is integral, it is faster to multiply
  // ------------  repeatedly than call pow().  iparm[] is constructed as
  //               follows:
  //  I_NPOLY indexing parameters, as above,
  //     nVar elements record the largest integral power for each variable,
  //     nVar elements record offsets into dparm for each variable,
  //   M*nVar flags to signal whether the power is integral,
  //   M*nVar integral powers.
  for (int ivar = 0; ivar < nVar; ivar++) {
    // Want at least the first degree power for all variables.
    i = iparm[I_MAXPOW] + ivar;
    iparm[i] = 1;
  }

  for (int ivar = 0; ivar < nVar; ivar++) {
    for (m = 0; m < M; m++) {
      i = iparm[I_DPOLY] + m*nTparm + 1 + ivar;
      double power = dparm[i];

      // Is it integral?  (Positive, negative, or zero.)
      int ipow = (int)power;
      if (power == (double)ipow) {
        // Signal that the power is integral.
        i = iparm[I_FLAGS] + m*nVar + ivar;
        if (ipow == 0) {
          iparm[i] = 3;
        } else {
          iparm[i] = 1;
        }

        // The integral power itself.
        i = iparm[I_IPOW] + m*nVar + ivar;
        iparm[i] = ipow;
      }

      // Record the largest integral power for each variable.
      i = iparm[I_MAXPOW] + ivar;
      if (iparm[i] < abs(ipow)) {
        iparm[i] = abs(ipow);
      }
    }
  }

  // How many of all powers of each variable will there be?
  int npow = 0;
  for (int ivar = 0; ivar < nVar; ivar++) {
    // Offset into dparm.
    i = iparm[I_DPOFF] + ivar;
    iparm[i] = lendp + npow;

    i = iparm[I_MAXPOW] + ivar;
    npow += iparm[i];
  }

  // Expand dparm to store the extra powers.
  if (npow) {
    lendp += npow;
    iparm[I_LENDP] = lendp;
    if ((dis->dparm[j] = realloc(dparm, lendp*sizeof(double))) == 0x0) {
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }
  }

  // No specialist de-distortions.
  dis->disp2x[j] = dispoly;
  dis->disx2p[j] = 0x0;

  // Translate Polynomial to TPD if possible, it's much faster.
  // However don't do it if the name was given as "Polynomial*".
  if (strcmp(dis->dtype[j], "Polynomial") == 0) {
    pol2tpd(j, dis);
  }

  return 0;
}

//----------------------------------------------------------------------------

int tpdset(int j, struct disprm *dis)

{
  static const char *function = "tpdset";

  if (dis == 0x0) return DISERR_NULL_POINTER;
  struct wcserr **err = &(dis->err);

  char id[32];
  sprintf(id, "TPD on axis %d", j+1);


  // TPD distortion.
  if (dis->Nhat[j] < 1 || 2 < dis->Nhat[j]) {
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "Axis map for %s must contain 1 or 2 entries, not %d", id,
      dis->Nhat[j]);
  }

  // Find the number of parameters.
  int ncoeff[2] = {0, 0};
  int doaux     = 0;
  int doradial  = 0;
  struct dpkey *keyp = dis->dp;
  for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    char *fp = strchr(keyp->field, '.') + 1;

    if (strncmp(fp, "TPD.", 4) == 0) {
      fp += 4;
      int idis;
      if (strncmp(fp, "FWD.", 4) == 0) {
        idis = 0;

      } else if (strncmp(fp, "REV.", 4) == 0) {
        // TPD may provide a polynomial approximation for the inverse.
        idis = 1;

      } else {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Unrecognized field name for %s: %s", id, keyp->field);
      }

      int k;
      sscanf(fp+4, "%d", &k);
      if (0 <= k && k <= 59) {
        if (ncoeff[idis] < k+1) ncoeff[idis] = k+1;

        // Any radial terms?
        if (k == 3 || k == 11 || k == 23 || k == 39 || k == 59) {
          doradial = 1;
        }

      } else {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid parameter number (%d) for %s: %s", k, id, keyp->field);
      }

    } else if (strncmp(fp, "AUX.", 4) == 0) {
      // Flag usage of auxiliary variables.
      doaux = 1;

    } else if (strcmp(fp, "DOCORR") &&
               strcmp(fp, "NAXES")  &&
              strncmp(fp, "AXIS.",   5) &&
              strncmp(fp, "OFFSET.", 7) &&
              strncmp(fp, "SCALE.",  6)) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Unrecognized field name for %s: %s", id, keyp->field);
    }
  }

  int (*(distpd[2]))(DISP2X_ARGS) = {0x0, 0x0};
  for (int idis = 0; idis < 2; idis++) {
    if (ncoeff[idis] <= 4) {
      if (idis) {
        // No inverse polynomial.
        break;
      }

      // First degree.
      ncoeff[idis] = 4;
      distpd[idis] = tpd1;
    } else if (ncoeff[idis] <= 7) {
      // Second degree.
      ncoeff[idis] = 7;
      distpd[idis] = tpd2;
    } else if (ncoeff[idis] <= 12) {
      // Third degree.
      ncoeff[idis] = 12;
      distpd[idis] = tpd3;
    } else if (ncoeff[idis] <= 17) {
      // Fourth degree.
      ncoeff[idis] = 17;
      distpd[idis] = tpd4;
    } else if (ncoeff[idis] <= 24) {
      // Fifth degree.
      ncoeff[idis] = 24;
      distpd[idis] = tpd5;
    } else if (ncoeff[idis] <= 31) {
      // Sixth degree.
      ncoeff[idis] = 31;
      distpd[idis] = tpd6;
    } else if (ncoeff[idis] <= 40) {
      // Seventh degree.
      ncoeff[idis] = 40;
      distpd[idis] = tpd7;
    } else if (ncoeff[idis] <= 49) {
      // Eighth degree.
      ncoeff[idis] = 49;
      distpd[idis] = tpd8;
    } else if (ncoeff[idis] <= 60) {
      // Ninth degree.
      ncoeff[idis] = 60;
      distpd[idis] = tpd9;
    } else {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid number of parameters (%d) for %s", ncoeff[idis], id);
    }
  }

  // disx2p() only uses the inverse TPD, if present, to provide a better
  // zeroth approximation.
  dis->disp2x[j] = distpd[0];
  dis->disx2p[j] = distpd[1];


// These iparm indices are specific to TPD (matching definitions in wcshdr.c).
#define I_TPDNCO  3	// No. of TPD coefficients, forward...
#define I_TPDINV  4	// ...and inverse.
#define I_TPDAUX  5	// True if auxiliary variables are used.
#define I_TPDRAD  6	// True if the radial variable is used.
#define I_NTPD    7

  // Record indexing parameters.
  int niparm = I_NTPD;
  if ((dis->iparm[j] = calloc(niparm, sizeof(int))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  int ndparm = (doaux?6:0) + ncoeff[0] + ncoeff[1];

  // The first three are more widely used.
  dis->iparm[j][I_DTYPE]  = DIS_TPD;
  dis->iparm[j][I_NIPARM] = niparm;
  dis->iparm[j][I_NDPARM] = ndparm;

  // Number of TPD coefficients.
  dis->iparm[j][I_TPDNCO] = ncoeff[0];
  dis->iparm[j][I_TPDINV] = ncoeff[1];

  // Flag for presence of auxiliary variables.
  dis->iparm[j][I_TPDAUX] = doaux;

  // Flag for presence of radial terms.
  dis->iparm[j][I_TPDRAD] = doradial;


  // Allocate memory for the polynomial coefficients and fill it.
  if ((dis->dparm[j] = calloc(ndparm, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  // Set default auxiliary coefficients.
  if (doaux) {
    dis->dparm[j][1] = 1.0;
    dis->dparm[j][5] = 1.0;
  }

  keyp = dis->dp;
  for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    char *fp = strchr(keyp->field, '.') + 1;

    if (strncmp(fp, "AUX.", 4) == 0) {
      // Auxiliary variables.
      fp += 4;
      int k;
      sscanf(fp, "%d", &k);
      if (k < 1 || 2 < k) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Bad auxiliary variable (%d) for %s: %s", k, id, keyp->field);
      }

      if ((fp = strchr(fp, '.')) == 0x0) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid field name for %s: %s", id, keyp->field);
      }
      fp++;

      if (strncmp(fp, "COEFF.", 6) != 0) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Unrecognized field name for %s: %s", id, keyp->field);
      }

      fp += 6;
      int m;
      sscanf(fp, "%d", &m);
      if (m < 0 || 2 < m) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid coefficient number (%d) for %s: %s", m, id, keyp->field);
      }

      int idis = 3*(k-1) + m;
      dis->dparm[j][idis] = dpkeyd(keyp);

    } else if (strncmp(fp, "TPD.", 4) == 0) {
      fp += 4;
      int idis = (doaux?6:0);
      if (strncmp(fp, "REV.", 4) == 0) {
        idis += ncoeff[0];
      }

      int k;
      sscanf(fp+4, "%d", &k);
      dis->dparm[j][idis+k] = dpkeyd(keyp);
    }
  }

  return 0;
}

//----------------------------------------------------------------------------

int pol2tpd(int j, struct disprm *dis)

{
  static const char *function = "pol2tpd";

  static const int map[][10] = {{ 0,  2,  6, 10, 16, 22, 30, 38, 48, 58},
                                { 1,  5,  9, 15, 21, 29, 37, 47, 57, -1},
                                { 4,  8, 14, 20, 28, 36, 46, 56, -1, -1},
                                { 7, 13, 19, 27, 35, 45, 55, -1, -1, -1},
                                {12, 18, 26, 34, 44, 54, -1, -1, -1, -1},
                                {17, 25, 33, 43, 53, -1, -1, -1, -1, -1},
                                {24, 32, 42, 52, -1, -1, -1, -1, -1, -1},
                                {31, 41, 51, -1, -1, -1, -1, -1, -1, -1},
                                {40, 50, -1, -1, -1, -1, -1, -1, -1, -1},
                                {49, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

  // Initialize.
  if (dis == 0x0) return DISERR_NULL_POINTER;
  struct wcserr **err = &(dis->err);

  int    *iparm = dis->iparm[j];
  double *dparm = dis->dparm[j];


  // Check the number of independent variables, no more than two.
  int Nhat = dis->Nhat[j];
  if (2 < Nhat) return -1;

  // Check auxiliaries: only one is allowed...
  int K = iparm[I_K];
  if (1 < K) return -1;
  if (K) {
    // ...and it must be radial.
    if (dparm[0] != 0.0) return -1;
    if (dparm[1] != 1.0) return -1;
    if (dparm[2] != 1.0) return -1;
    if (dparm[3] != 0.5) return -1;
    if (dparm[4] != 2.0) return -1;
    if (dparm[5] != 2.0) return -1;
  }

  // Check powers...
  int *iflgp = iparm + iparm[I_FLAGS];
  int *ipowp = iparm + iparm[I_IPOW];
  int degree = 0;
  for (int m = 0; m < iparm[I_M]; m++) {
    int deg = 0;
    for (int jhat = 0; jhat < Nhat; jhat++) {
      // ...they must be positive integral.
      if (*iflgp == 0)  return -1;
      if (*ipowp < 0)   return -1;
      deg += *ipowp;
      iflgp++;
      ipowp++;
    }

    // The polynomial degree can't be greater than 9.
    if (9 < deg) return -1;

    if (K) {
      // Likewise for the radial variable.
      if (*iflgp == 0)  return -1;
      if (*ipowp) {
        if (*ipowp < 0) return -1;
        if (9 < *ipowp) return -1;

        // Can't mix the radial and other terms.
        if (deg)        return -1;

        // Can't have even powers of the radial variable.
        deg = *ipowp;
        if (!(deg%2))   return -1;
      }
      iflgp++;
      ipowp++;
    }

    if (degree < deg) degree = deg;
  }


  // OK, it ticks all the boxes.  Now translate it.
  int ndparm = 0;
  if (degree == 1) {
    ndparm = 4;
    dis->disp2x[j] = tpd1;
  } else if (degree == 2) {
    ndparm = 7;
    dis->disp2x[j] = tpd2;
  } else if (degree == 3) {
    ndparm = 12;
    dis->disp2x[j] = tpd3;
  } else if (degree == 4) {
    ndparm = 17;
    dis->disp2x[j] = tpd4;
  } else if (degree == 5) {
    ndparm = 24;
    dis->disp2x[j] = tpd5;
  } else if (degree == 6) {
    ndparm = 31;
    dis->disp2x[j] = tpd6;
  } else if (degree == 7) {
    ndparm = 40;
    dis->disp2x[j] = tpd7;
  } else if (degree == 8) {
    ndparm = 49;
    dis->disp2x[j] = tpd8;
  } else if (degree == 9) {
    ndparm = 60;
    dis->disp2x[j] = tpd9;
  }

  // No specialist de-distortions.
  dis->disx2p[j] = 0x0;

  // Record indexing parameters.
  int niparm = I_NTPD;
  int *tpd_iparm;
  if ((tpd_iparm = calloc(niparm, sizeof(int))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  // The first three are more widely used.
  tpd_iparm[I_DTYPE]  = DIS_TPD;
  tpd_iparm[I_NIPARM] = niparm;
  tpd_iparm[I_NDPARM] = ndparm;

  // Number of TPD coefficients.
  tpd_iparm[I_TPDNCO] = ndparm;
  tpd_iparm[I_TPDINV] = 0;

  // No auxiliary variables yet.
  tpd_iparm[I_TPDAUX] = 0;

  // Flag for presence of radial terms.
  tpd_iparm[I_TPDRAD] = K;


  // Allocate memory for the polynomial coefficients and fill it.
  double *tpd_dparm;
  if ((tpd_dparm = calloc(ndparm, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  ipowp = iparm + iparm[I_IPOW];
  double *dpolp = dparm + iparm[I_DPOLY];
  for (int m = 0; m < iparm[I_M]; m++) {
    if (K && ipowp[Nhat]) {
      // The radial variable.
      switch (ipowp[Nhat]) {
      case 1:
        tpd_dparm[3]  = *dpolp;
        break;
      case 3:
        tpd_dparm[11] = *dpolp;
        break;
      case 5:
        tpd_dparm[23] = *dpolp;
        break;
      case 7:
        tpd_dparm[39] = *dpolp;
        break;
      case 9:
        tpd_dparm[59] = *dpolp;
        break;
      }

    } else {
      // The independent variables.
      int p[] = {0, 0};
      for (int jhat = 0; jhat < Nhat; jhat++) {
        p[jhat] = ipowp[jhat];
      }

      int n = map[p[0]][p[1]];
      tpd_dparm[n] = *dpolp;
    }


    ipowp += iparm[I_NVAR];
    dpolp += iparm[I_NVAR] + 1;
  }


  // Switch from Polynomial to TPD.
  free(iparm);
  free(dparm);
  dis->iparm[j] = tpd_iparm;
  dis->dparm[j] = tpd_dparm;

  return 0;
}

//----------------------------------------------------------------------------

int tpvset(int j, struct disprm *dis)

{
  static const char *function = "tpvset";

  // Initialize.
  if (dis == 0x0) return DISERR_NULL_POINTER;
  struct wcserr **err = &(dis->err);

  // TPV "projection".
  char id[32];
  sprintf(id, "TPV on axis %d", j+1);

  // TPV is a sequent distortion, applied to intermediate world coordinates
  // (normally used with CDi_ja).  It computes corrected coordinates directly.
  dis->docorr[j] = 0;

  if (dis->Nhat[j] != 2) {
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "Axis map for %s must contain 2 entries, not %d", id, dis->Nhat[j]);
  }

  // Find the number of parameters.
  int ndparm   = 0;
  int doradial = 0;
  struct dpkey *keyp = dis->dp;
  for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    char *fp = strchr(keyp->field, '.') + 1;

    if (strncmp(fp, "TPV.", 4) == 0) {
      int k;
      sscanf(fp+4, "%d", &k);
      if (0 <= k && k <= 39) {
        if (ndparm < k+1) ndparm = k+1;

        // Any radial terms?
        if (k == 3 || k == 11 || k == 23 || k == 39 || k == 59) {
          doradial = 1;
        }

      } else {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid parameter number (%d) for %s: %s", k, id, keyp->field);
      }

    } else if (strcmp(fp, "NAXES")  &&
              strncmp(fp, "AXIS.",   5) &&
              strncmp(fp, "OFFSET.", 7) &&
              strncmp(fp, "SCALE.",  6)) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Unrecognized field name for %s: %s", id, keyp->field);
    }
  }

  // TPD is going to do the dirty work.
  if (ndparm <= 4) {
    // First degree.
    ndparm = 4;
    dis->disp2x[j] = tpd1;
  } else if (ndparm <= 7) {
    // Second degree.
    ndparm = 7;
    dis->disp2x[j] = tpd2;
  } else if (ndparm <= 12) {
    // Third degree.
    ndparm = 12;
    dis->disp2x[j] = tpd3;
  } else if (ndparm <= 17) {
    // Fourth degree.
    ndparm = 17;
    dis->disp2x[j] = tpd4;
  } else if (ndparm <= 24) {
    // Fifth degree.
    ndparm = 24;
    dis->disp2x[j] = tpd5;
  } else if (ndparm <= 31) {
    // Sixth degree.
    ndparm = 31;
    dis->disp2x[j] = tpd6;
  } else if (ndparm <= 40) {
    // Seventh degree.
    ndparm = 40;
    dis->disp2x[j] = tpd7;
  } else {
    // Could go to ninth degree, but that wouldn't be legit.
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "Invalid number of parameters (%d) for %s", ndparm, id);
  }

  // No specialist de-distortions.
  dis->disx2p[j] = 0x0;

  // Record indexing parameters.
  int niparm = I_NTPD;
  if ((dis->iparm[j] = calloc(niparm, sizeof(int))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  // The first three are more widely used.
  dis->iparm[j][I_DTYPE]  = DIS_TPD;
  dis->iparm[j][I_NIPARM] = niparm;
  dis->iparm[j][I_NDPARM] = ndparm;

  // Number of TPD coefficients.
  dis->iparm[j][I_TPDNCO] = ndparm;
  dis->iparm[j][I_TPDINV] = 0;

  // TPV never needs auxiliary variables.
  dis->iparm[j][I_TPDAUX] = 0;

  // Flag for presence of radial terms.
  dis->iparm[j][I_TPDRAD] = doradial;


  // Allocate memory for the polynomial coefficients and fill it.
  if ((dis->dparm[j] = calloc(ndparm, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  keyp = dis->dp;
  for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    char *fp = strchr(keyp->field, '.') + 1;

    // One-to-one correspondence between TPV and TPD coefficients.
    if (strncmp(fp, "TPV.", 4) == 0) {
      int k;
      sscanf(fp+4, "%d", &k);
      dis->dparm[j][k] = dpkeyd(keyp);
    }
  }

  return 0;
}

//----------------------------------------------------------------------------

int sipset(int j, struct disprm *dis)

{
  static const char *function = "sipset";

  static const int map[][10] = {{ 0,  2,  6, 10, 16, 22, 30, 38, 48, 58},
                                { 1,  5,  9, 15, 21, 29, 37, 47, 57, -1},
                                { 4,  8, 14, 20, 28, 36, 46, 56, -1, -1},
                                { 7, 13, 19, 27, 35, 45, 55, -1, -1, -1},
                                {12, 18, 26, 34, 44, 54, -1, -1, -1, -1},
                                {17, 25, 33, 43, 53, -1, -1, -1, -1, -1},
                                {24, 32, 42, 52, -1, -1, -1, -1, -1, -1},
                                {31, 41, 51, -1, -1, -1, -1, -1, -1, -1},
                                {40, 50, -1, -1, -1, -1, -1, -1, -1, -1},
                                {49, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

  // Initialize.
  if (dis == 0x0) return DISERR_NULL_POINTER;
  struct wcserr **err = &(dis->err);

  // Simple Imaging Polynomial.
  char id[32];
  sprintf(id, "SIP on axis %d", j+1);


  // SIP is a prior distortion that computes an additive correction.
  dis->docorr[j] = 1;

  if (dis->Nhat[j] != 2) {
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "Axis map for %s must contain 2 entries, not %d", id, dis->Nhat[j]);
  }

  // Find the polynomial degree, at least 1 for the forward function.
  int degree[2] = {1, -1};
  struct dpkey *keyp = dis->dp;
  for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    char *fp = strchr(keyp->field, '.') + 1;

    if (strncmp(fp, "SIP.", 4) == 0) {
      fp += 4;
      int idis;
      if (strncmp(fp, "FWD.", 4) == 0) {
        idis = 0;

      } else if (strncmp(fp, "REV.", 4) == 0) {
        // SIP uses a polynomial approximation for the inverse.
        idis = 1;

      } else {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Unrecognized field name for %s: %s", id, keyp->field);
      }

      fp += 4;
      int p, q;
      sscanf(fp, "%d_%d", &p, &q);
      int deg = p + q;
      if (p < 0 || 9 < p || q < 0 || 9 < q || 9 < deg) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid powers (%d, %d) for %s: %s", p, q, id, keyp->field);
      }

      if (degree[idis] < deg) degree[idis] = deg;

    } else if (strcmp(fp, "NAXES")  &&
              strncmp(fp, "AXIS.",   5) &&
              strncmp(fp, "OFFSET.", 7) &&
              strncmp(fp, "SCALE.",  6)) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Unrecognized field name for %s: %s", id, keyp->field);
    }
  }

  if (degree[1] == 0 ) degree[1] = 1;

  // TPD is going to do the dirty work.
  int (*(distpd[2]))(DISP2X_ARGS) = {0x0, 0x0}, ncoeff[2];
  for (int idis = 0; idis < 2; idis++) {
    ncoeff[idis] = 0;
    if (degree[idis] == 1) {
      ncoeff[idis] = 4;
      distpd[idis] = tpd1;
    } else if (degree[idis] == 2) {
      ncoeff[idis] = 7;
      distpd[idis] = tpd2;
    } else if (degree[idis] == 3) {
      ncoeff[idis] = 12;
      distpd[idis] = tpd3;
    } else if (degree[idis] == 4) {
      ncoeff[idis] = 17;
      distpd[idis] = tpd4;
    } else if (degree[idis] == 5) {
      ncoeff[idis] = 24;
      distpd[idis] = tpd5;
    } else if (degree[idis] == 6) {
      ncoeff[idis] = 31;
      distpd[idis] = tpd6;
    } else if (degree[idis] == 7) {
      ncoeff[idis] = 40;
      distpd[idis] = tpd7;
    } else if (degree[idis] == 8) {
      ncoeff[idis] = 49;
      distpd[idis] = tpd8;
    } else if (degree[idis] == 9) {
      ncoeff[idis] = 60;
      distpd[idis] = tpd9;
    }
  }

  // SIP uses a polynomial approximation to the inverse.  It's not very
  // accurate but may provide disx2p() with a better zeroth approximation.
  dis->disp2x[j] = distpd[0];
  dis->disx2p[j] = distpd[1];


  // Record indexing parameters.
  int niparm = I_NTPD;
  if ((dis->iparm[j] = calloc(niparm, sizeof(int))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  int ndparm = ncoeff[0] + ncoeff[1];

  // The first three are more widely used.
  dis->iparm[j][I_DTYPE]  = DIS_TPD;
  dis->iparm[j][I_NIPARM] = niparm;
  dis->iparm[j][I_NDPARM] = ndparm;

  // Number of TPD coefficients.
  dis->iparm[j][I_TPDNCO] = ncoeff[0];
  dis->iparm[j][I_TPDINV] = ncoeff[1];

  // SIP never needs auxiliary variables.
  dis->iparm[j][I_TPDAUX] = 0;

  // SIP never needs the radial terms.
  dis->iparm[j][I_TPDRAD] = 0;


  // Allocate memory for the polynomial coefficients and fill it.
  if ((dis->dparm[j] = calloc(ndparm, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  keyp = dis->dp;
  for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    char *fp = strchr(keyp->field, '.') + 1;

    if (strncmp(fp, "SIP.", 4) == 0) {
      fp += 4;
      int idis;
      if (strncmp(fp, "FWD.", 4) == 0) {
        idis = 0;
      } else {
        idis = ncoeff[0];
      }

      int p, q;
      sscanf(fp+4, "%d_%d", &p, &q);

      // Map to TPD coefficient number.
      idis += map[p][q];

      dis->dparm[j][idis] = dpkeyd(keyp);
    }
  }


  return 0;
}

//----------------------------------------------------------------------------

int dssset(int j, struct disprm *dis)

{
  static const char *function = "dssset";

  // Initialize.
  if (dis == 0x0) return DISERR_NULL_POINTER;
  struct wcserr **err = &(dis->err);

  // Digitized Sky Survey.
  char id[32];
  sprintf(id, "DSS on axis %d", j+1);


  // DSS is translated into a sequent distortion, applied to intermediate
  // pixel coordinates.  It computes corrected coordinates directly.
  dis->docorr[j] = 0;

  if (dis->Nhat[j] != 2) {
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "Axis map for %s must contain 2 entries, not %d", id, dis->Nhat[j]);
  }

  // Safe to assume the polynomial degree is 5 (or less).
  int ncoeff = 24;
  dis->disp2x[j] = tpd5;

  // No specialist de-distortions.
  dis->disx2p[j] = 0x0;


  // Record indexing parameters.
  int niparm = I_NTPD;
  if ((dis->iparm[j] = calloc(niparm, sizeof(int))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  int ndparm = 6 + ncoeff;

  // The first three are more widely used.
  dis->iparm[j][I_DTYPE]  = DIS_TPD;
  dis->iparm[j][I_NIPARM] = niparm;
  dis->iparm[j][I_NDPARM] = ndparm;

  // Number of TPD coefficients.
  dis->iparm[j][I_TPDNCO] = ncoeff;
  dis->iparm[j][I_TPDINV] = 0;

  // DSS always needs auxiliary variables.
  dis->iparm[j][I_TPDAUX] = 1;

  // DSS never needs the radial terms.
  dis->iparm[j][I_TPDRAD] = 0;


  // Allocate memory for the polynomial coefficients and fill it.
  if ((dis->dparm[j] = calloc(ndparm, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  // This translation follows WCS Paper IV, Sect. 5.2 using the same
  // variable names.  Find A1, A2, A3, B1, B2, and B3.
  double A1, A2, A3, B1, B2, B3;
  A1 = A2 = A3 = 0.0;
  B1 = B2 = B3 = 0.0;
  struct dpkey *keyp = dis->dp;
  for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
    char *fp = strchr(keyp->field, '.') + 1;
    if (strncmp(fp, "DSS.AMD.", 8) == 0) {
      fp += 8;
      int m;
      sscanf(fp, "%d", &m);

      if (m == 1) {
        if (keyp->j == 1) {
          A1 = dpkeyd(keyp);
        } else {
          B1 = dpkeyd(keyp);
        }
      } else if (m == 2) {
        if (keyp->j == 1) {
          A2 = dpkeyd(keyp);
        } else {
          B2 = dpkeyd(keyp);
        }
      } else if (m == 3) {
        if (keyp->j == 1) {
          A3 = dpkeyd(keyp);
        } else {
          B3 = dpkeyd(keyp);
        }
      }
    }
  }

  double X0 = (A2*B3 - A3*B1) / (A1*B1 - A2*B2);
  double Y0 = (A3*B2 - A1*B3) / (A1*B1 - A2*B2);

  double S = sqrt(fabs(A1*B1 - A2*B2));
  if (S == 0.0) {
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "Coefficient scale for %s is zero.", id);
  }

  // Coefficients for the auxiliary variables.
  double *dparm = dis->dparm[j];
  if (j == 0) {
    dparm[0] =  X0;
    dparm[1] = -B1/S;
    dparm[2] = -A2/S;
    dparm[3] =  Y0;
    dparm[4] =  B2/S;
    dparm[5] =  A1/S;

    // Change the sign of S for scaling the A coefficients.
    S *= -1.0;

  } else {
    dparm[0] =  Y0;
    dparm[1] =  B2/S;
    dparm[2] =  A1/S;
    dparm[3] =  X0;
    dparm[4] = -B1/S;
    dparm[5] = -A2/S;
  }

  // Translate DSS coefficients to TPD.
  dparm += 6;
  int degree = 3;
  keyp = dis->dp;
  for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    char *fp = strchr(keyp->field, '.') + 1;

    if (strncmp(fp, "DSS.AMD.", 8) == 0) {
      // Skip zero coefficients.
      double coeff = dpkeyd(keyp);
      if (coeff == 0.0) continue;

      fp += 8;
      int m;
      sscanf(fp, "%d", &m);

      // Apply the coefficient scale factor.
      coeff /= S;

      if (m == 1) {
        dparm[1]  = coeff;
      } else if (m == 2) {
        dparm[2]  = coeff;
      } else if (m == 3) {
        dparm[0]  = coeff;
      } else if (m == 4) {
        dparm[4] += coeff;
      } else if (m == 5) {
        dparm[5]  = coeff;
      } else if (m == 6) {
        dparm[6] += coeff;
      } else if (m == 7) {
        dparm[4] += coeff;
        dparm[6] += coeff;
      } else if (m == 8) {
        dparm[7] += coeff;
      } else if (m == 9) {
        dparm[8]  = coeff;
      } else if (m == 10) {
        dparm[9] += coeff;
      } else if (m == 11) {
        dparm[10] = coeff;
      } else if (m == 12) {
        dparm[7] += coeff;
        dparm[9] += coeff;
      } else if (m == 13) {
        dparm[17] = coeff;
        dparm[19] = coeff * 2.0;
        dparm[21] = coeff;
        degree = 5;
      } else if (coeff != 0.0) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid parameter for %s: %s", m, id, keyp->field);
      }

    } else if (strcmp(fp, "NAXES")  &&
              strncmp(fp, "AXIS.",   5) &&
              strncmp(fp, "OFFSET.", 7) &&
              strncmp(fp, "SCALE.",  6)) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Unrecognized field name for %s: %s", id, keyp->field);
    }
  }

  // The DSS polynomial doesn't have 4th degree terms, and the 5th degree
  // coefficient is often zero.
  if (degree == 3) {
    dis->iparm[j][I_TPDNCO] = 12;
    dis->disp2x[j] = tpd3;
  }

  return 0;
}

//----------------------------------------------------------------------------

#define CHEBYSHEV 1
#define LEGENDRE  2
#define MONOMIAL  3

int watset(int j, struct disprm *dis)

{
  static const char *function = "watset";

  static const int map[][10] = {{ 0,  2,  6, 10, 16, 22, 30, 38, 48, 58},
                                { 1,  5,  9, 15, 21, 29, 37, 47, 57, -1},
                                { 4,  8, 14, 20, 28, 36, 46, 56, -1, -1},
                                { 7, 13, 19, 27, 35, 45, 55, -1, -1, -1},
                                {12, 18, 26, 34, 44, 54, -1, -1, -1, -1},
                                {17, 25, 33, 43, 53, -1, -1, -1, -1, -1},
                                {24, 32, 42, 52, -1, -1, -1, -1, -1, -1},
                                {31, 41, 51, -1, -1, -1, -1, -1, -1, -1},
                                {40, 50, -1, -1, -1, -1, -1, -1, -1, -1},
                                {49, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

  // Initialize.
  if (dis == 0x0) return DISERR_NULL_POINTER;
  struct wcserr **err = &(dis->err);

  // WAT (TNX or ZPX) Polynomial.
  char id[32];
  sprintf(id, "WAT (%s) on axis %d", dis->dtype[0]+4, j+1);


  // WAT is a sequent distortion, applied to intermediate world coordinates
  // (normally used with CDi_ja).  It computes an additive correction.
  dis->docorr[j] = 1;

  if (dis->Nhat[j] != 2) {
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "Axis map for %s must contain 2 entries, not %d", id, dis->Nhat[j]);
  }

  // Find the polynomial degree (at least 1), kind, and domain.
  int degree = 1;
  int kind = 0;
  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
  struct dpkey *keyp = dis->dp;
  for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    char *fp = strchr(keyp->field, '.') + 1;

    if (strncmp(fp, "WAT.", 4) == 0) {
      fp += 4;
      if (strncmp(fp, "CHBY.", 5) == 0 ||
          strncmp(fp, "LEGR.", 5) == 0 ||
          strncmp(fp, "MONO.", 5) == 0) {

        fp += 5;
        int m, n;
        sscanf(fp, "%d_%d", &m, &n);
        int deg = m + n;
        if (m < 0 || 9 < m || n < 0 || 9 < n || 9 < deg) {
          return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid powers (%d, %d) for %s: %s", m, n, id, keyp->field);
        }

        if (degree < deg) degree = deg;

      } else if (strcmp(fp, "POLY") == 0) {
        kind = dpkeyi(keyp);

      } else if (strcmp(fp, "XMIN") == 0) {
        xmin = dpkeyd(keyp);

      } else if (strcmp(fp, "XMAX") == 0) {
        xmax = dpkeyd(keyp);

      } else if (strcmp(fp, "YMIN") == 0) {
        ymin = dpkeyd(keyp);

      } else if (strcmp(fp, "YMAX") == 0) {
        ymax = dpkeyd(keyp);
      }

    } else if (strcmp(fp, "NAXES")  &&
              strncmp(fp, "AXIS.",   5) &&
              strncmp(fp, "OFFSET.", 7) &&
              strncmp(fp, "SCALE.",  6)) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Unrecognized field name for %s: %s", id, keyp->field);
    }
  }

  int doaux = (kind == 1 || kind == 2);

  // TPD is going to do the dirty work.
  int ncoeff = 0;
  if (degree == 1) {
    // First degree.
    ncoeff = 4;
    dis->disp2x[j] = tpd1;
  } else if (degree == 2) {
    // Second degree.
    ncoeff = 7;
    dis->disp2x[j] = tpd2;
  } else if (degree == 3) {
    // Third degree.
    ncoeff = 12;
    dis->disp2x[j] = tpd3;
  } else if (degree == 4) {
    // Fourth degree.
    ncoeff = 17;
    dis->disp2x[j] = tpd4;
  } else if (degree == 5) {
    // Fifth degree.
    ncoeff = 24;
    dis->disp2x[j] = tpd5;
  } else if (degree == 6) {
    // Sixth degree.
    ncoeff = 31;
    dis->disp2x[j] = tpd6;
  } else if (degree == 7) {
    // Seventh degree.
    ncoeff = 40;
    dis->disp2x[j] = tpd7;
  } else if (degree == 8) {
    // Eighth degree.
    ncoeff = 49;
    dis->disp2x[j] = tpd8;
  } else if (degree == 9) {
    // Ninth degree.
    ncoeff = 60;
    dis->disp2x[j] = tpd9;
  }

  // No specialist de-distortions.
  dis->disx2p[j] = 0x0;


  // Record indexing parameters.
  int niparm = I_NTPD;
  if ((dis->iparm[j] = calloc(niparm, sizeof(int))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  int *iparm = dis->iparm[j];

  int ndparm = 6 + ncoeff;

  // The first three are more widely used.
  iparm[I_DTYPE]  = DIS_TPD;
  iparm[I_NIPARM] = niparm;
  iparm[I_NDPARM] = ndparm;

  // Number of TPD coefficients.
  iparm[I_TPDNCO] = ncoeff;
  iparm[I_TPDINV] = 0;

  // The Chebyshev and Legendre polynomials use auxiliary variables.
  iparm[I_TPDAUX] = doaux;

  // WAT never needs the radial terms.
  iparm[I_TPDRAD] = 0;


  // Allocate memory for the polynomial coefficients and fill it.
  if ((dis->dparm[j] = calloc(ndparm, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  double *dparm = dis->dparm[j];


  // Coefficients for the auxiliary variables.
  if (doaux) {
    double x0 = (xmax + xmin)/2.0;
    double dx = (xmax - xmin)/2.0;
    if (dx == 0.0) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "X-span for %s is zero", id);
    }

    dparm[0] = -x0/dx;
    dparm[1] = 1.0/dx;
    dparm[2] = 0.0;

    double y0 = (ymax + ymin)/2.0;
    double dy = (ymax - ymin)/2.0;
    if (dy == 0.0) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Y-span for %s is zero", id);
    }

    dparm[3] = -y0/dy;
    dparm[4] = 0.0;
    dparm[5] = 1.0/dy;

    dparm += 6;
  }


  // Unpack the polynomial coefficients.
  keyp = dis->dp;
  for (int idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    char *fp = strchr(keyp->field, '.') + 1;

    if ((kind == CHEBYSHEV && strncmp(fp, "WAT.CHBY.", 9) == 0) ||
        (kind == LEGENDRE  && strncmp(fp, "WAT.LEGR.", 9) == 0) ||
        (kind == MONOMIAL  && strncmp(fp, "WAT.MONO.", 9) == 0)) {
      fp += 9;

      int m, n;
      sscanf(fp, "%d_%d", &m, &n);

      if (kind == MONOMIAL) {
        // Monomial coefficient, maps simply to TPD coefficient number.
        int idis = map[m][n];
        dparm[idis] = dpkeyd(keyp);

      } else {
        // Coefficient of the product of two Chebyshev or two Legendre
        // polynomials.  Find the corresponding monomial coefficients.
        double coeff = dpkeyd(keyp);

        double coeffm[10], coeffn[10];
        cheleg(kind, m, n, coeffm, coeffn);
        for (int im = 0; im <= m; im++) {
          if (coeffm[im] == 0.0) continue;

          for (int in = 0; in <= n; in++) {
            if (coeffn[in] == 0.0) continue;

            int idis = map[im][in];
            dparm[idis] += coeff*coeffm[im]*coeffn[in];
          }
        }
      }
    }
  }

  return 0;
}

//----------------------------------------------------------------------------
// Compute the coefficients of Chebyshev or Legendre polynomials of degree
// m and n.

int cheleg(int kind, int m, int n, double coeffm[], double coeffn[])

{
  int N = (m > n) ? m : n;

  // Allocate work arrays.
  double *coeff[3];
  coeff[0] = calloc(3*(N+1), sizeof(double));
  coeff[1] = coeff[0] + (N+1);
  coeff[2] = coeff[1] + (N+1);

  for (int j = 0; j <= N; j++) {
    int j0 =  j%3;

    if (j == 0) {
      coeff[0][0] = 1.0;

    } else if (j == 1) {
      coeff[1][1] = 1.0;

    } else {
      // Cyclic buffer indices.
      int j1 = (j-1)%3;
      int j2 = (j-2)%3;

      memset(coeff[j0], 0, (N+1)*sizeof(double));

      double d = (double)j;
      for (int k = 0; k < N; k++) {
        if (kind == CHEBYSHEV) {
          coeff[j0][k+1] = 2.0 * coeff[j1][k];
          coeff[j0][k]  -=       coeff[j2][k];
        } else if (kind == LEGENDRE) {
          coeff[j0][k+1] = ((2.0*d - 1.0) * coeff[j1][k]) / d;
          coeff[j0][k]  -=     ((d - 1.0) * coeff[j2][k]) / d;
        }
      }
    }

    if (j == m) memcpy(coeffm, coeff[j0], (m+1)*sizeof(double));
    if (j == n) memcpy(coeffn, coeff[j0], (n+1)*sizeof(double));
  }

  free(coeff[0]);

  return 0;
}

//----------------------------------------------------------------------------

int dispoly(
  int dummy,
  const int iparm[],
  const double dparm[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  // Avert nuisance compiler warnings about unused parameters.
  (void)dummy;

  // Check for zeroes.
  for (int jhat = 0; jhat < Nhat; jhat++) {
    if (rawcrd[jhat] == 0.0) {
      *discrd = 0.0;
      return 0;
    }
  }

  // Working memory for auxiliaries &c. was allocated at the end of p[].
  double *aux = (double *)(dparm + iparm[I_DAUX]);

  // Compute the auxiliary variables.
  for (int k = 0; k < iparm[I_K]; k++) {
    const double *cptr = dparm + k*iparm[I_NKPARM];
    const double *pptr = cptr + (1+Nhat);

    aux[k] = *(cptr++);
    double auxp0  = *(pptr++);

    for (int jhat = 0; jhat < Nhat; jhat++) {
      aux[k] += *(cptr++)*pow(rawcrd[jhat], *(pptr++));
    }

    aux[k] = pow(aux[k], auxp0);

    // Check for zeroes.
    if (aux[k] == 0.0) {
      *discrd = 0.0;
      return 0;
    }
  }


  // Compute all required integral powers of the variables.
  const int *imaxpow = iparm + iparm[I_MAXPOW];
  double *dvarpow = (double *)(dparm + iparm[I_DVPOW]);

  const int *imaxp = imaxpow;
  double *dpowp = dvarpow;
  for (int jhat = 0; jhat < Nhat; jhat++, imaxp++) {
    double var = 1.0;
    for (int ip = 0; ip < *imaxp; ip++, dpowp++) {
      var *= rawcrd[jhat];
      *dpowp = var;
    }
  }

  for (int k = 0; k < iparm[I_K]; k++, imaxp++) {
    double var = 1.0;
    for (int ip = 0; ip < *imaxp; ip++, dpowp++) {
      var *= aux[k];
      *dpowp = var;
    }
  }

  // Loop for each term of the polynomial.
  *discrd = 0.0;
  const int    *iflgp = iparm + iparm[I_FLAGS];
  const int    *ipowp = iparm + iparm[I_IPOW];
  const double *dpolp = dparm + iparm[I_DPOLY];
  for (int m = 0; m < iparm[I_M]; m++) {
    double term = *(dpolp++);

    // Loop over all variables.
    imaxp = imaxpow;
    dpowp = dvarpow - 1;
    for (int ivar = 0; ivar < iparm[I_NVAR]; ivar++) {
      if (*iflgp & 2) {
        // Nothing (zero power).

      } else if (*iflgp) {
        // Integral power.
        if (*ipowp < 0) {
          // Negative.
          term /= dpowp[*ipowp];
        } else {
          // Positive.
          term *= dpowp[*ipowp];
        }

      } else {
        // Fractional power.
        term *= pow(dpowp[0], *dpolp);
      }

      iflgp++;
      ipowp++;
      dpolp++;

      dpowp += *imaxp;
      imaxp++;
    }

    *discrd += term;
  }

  return 0;
}

//----------------------------------------------------------------------------

int tpd1(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  if (i[I_TPDNCO+inverse] != 4 || 2 < Nhat) {
    return 1;
  }

  double r, s;
  double u = rawcrd[0];
  double v = rawcrd[1];

  // Auxiliary variables?
  if (i[I_TPDAUX]) {
    r = p[0] + p[1]*u + p[2]*v;
    v = p[3] + p[4]*u + p[5]*v;
    u = r;
    p += 6;
  }

  if (inverse) p += i[I_TPDNCO];

  // First degree.
  *discrd = p[0] + u*p[1];

  if (Nhat == 1) return 0;

  *discrd += v*p[2];

  // Radial terms?
  if (i[I_TPDRAD]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*p[3];
  }

  return 0;
}

//----------------------------------------------------------------------------

int tpd2(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  if (i[I_TPDNCO+inverse] != 7 || 2 < Nhat) {
    return 1;
  }

  double r, s;
  double u = rawcrd[0];
  double v = rawcrd[1];

  // Auxiliary variables?
  if (i[I_TPDAUX]) {
    r = p[0] + p[1]*u + p[2]*v;
    v = p[3] + p[4]*u + p[5]*v;
    u = r;
    p += 6;
  }

  if (inverse) p += i[I_TPDNCO];

  // Second degree.
  *discrd = p[0] + u*(p[1] + u*(p[4]));

  if (Nhat == 1) return 0;

  *discrd +=
      v*(p[2]  + v*(p[6]))
    + u*(p[5])*v;

  // Radial terms?
  if (i[I_TPDRAD]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*p[3];
  }

  return 0;
}

//----------------------------------------------------------------------------

int tpd3(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  if (i[I_TPDNCO+inverse] != 12 || 2 < Nhat) {
    return 1;
  }

  double r, s;
  double u = rawcrd[0];
  double v = rawcrd[1];

  // Auxiliary variables?
  if (i[I_TPDAUX]) {
    r = p[0] + p[1]*u + p[2]*v;
    v = p[3] + p[4]*u + p[5]*v;
    u = r;
    p += 6;
  }

  if (inverse) p += i[I_TPDNCO];

  // Third degree.
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7])));

  if (Nhat == 1) return 0;

  *discrd +=
      v*(p[2]  + v*(p[6]  + v*(p[10])))
    + u*(p[5]  + v*(p[9])
    + u*(p[8]))*v;

  // Radial terms?
  if (i[I_TPDRAD]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11]));
  }

  return 0;
}

//----------------------------------------------------------------------------

int tpd4(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  if (i[I_TPDNCO+inverse] != 17 || 2 < Nhat) {
    return 1;
  }

  double r, s;
  double u = rawcrd[0];
  double v = rawcrd[1];

  // Auxiliary variables?
  if (i[I_TPDAUX]) {
    r = p[0] + p[1]*u + p[2]*v;
    v = p[3] + p[4]*u + p[5]*v;
    u = r;
    p += 6;
  }

  if (inverse) p += i[I_TPDNCO];

  // Fourth degree.
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7] + u*(p[12]))));

  if (Nhat == 1) return 0;

  *discrd +=
      v*(p[2]  + v*(p[6]  + v*(p[10] + v*(p[16]))))
    + u*(p[5]  + v*(p[9]  + v*(p[15]))
    + u*(p[8]  + v*(p[14])
    + u*(p[13])))*v;

  // Radial terms?
  if (i[I_TPDRAD]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11]));
  }

  return 0;
}

//----------------------------------------------------------------------------

int tpd5(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  if (i[I_TPDNCO+inverse] != 24 || 2 < Nhat) {
    return 1;
  }

  double r, s;
  double u = rawcrd[0];
  double v = rawcrd[1];

  // Auxiliary variables?
  if (i[I_TPDAUX]) {
    r = p[0] + p[1]*u + p[2]*v;
    v = p[3] + p[4]*u + p[5]*v;
    u = r;
    p += 6;
  }

  if (inverse) p += i[I_TPDNCO];

  // Fifth degree.
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7] + u*(p[12] + u*(p[17])))));

  if (Nhat == 1) return 0;

  *discrd +=
      v*(p[2]  + v*(p[6]  + v*(p[10] + v*(p[16] + v*(p[22])))))
    + u*(p[5]  + v*(p[9]  + v*(p[15] + v*(p[21])))
    + u*(p[8]  + v*(p[14] + v*(p[20]))
    + u*(p[13] + v*(p[19])
    + u*(p[18]))))*v;

  // Radial terms?
  if (i[I_TPDRAD]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11] + s*(p[23])));
  }

  return 0;
}

//----------------------------------------------------------------------------

int tpd6(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  if (i[I_TPDNCO+inverse] != 31 || 2 < Nhat) {
    return 1;
  }

  double r, s;
  double u = rawcrd[0];
  double v = rawcrd[1];

  // Auxiliary variables?
  if (i[I_TPDAUX]) {
    r = p[0] + p[1]*u + p[2]*v;
    v = p[3] + p[4]*u + p[5]*v;
    u = r;
    p += 6;
  }

  if (inverse) p += i[I_TPDNCO];

  // Sixth degree.
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7] + u*(p[12] + u*(p[17] + u*(p[24]))))));

  if (Nhat == 1) return 0;

  *discrd +=
      v*(p[2]  + v*(p[6]  + v*(p[10] + v*(p[16] + v*(p[22] + v*(p[30]))))))
    + u*(p[5]  + v*(p[9]  + v*(p[15] + v*(p[21] + v*(p[29]))))
    + u*(p[8]  + v*(p[14] + v*(p[20] + v*(p[28])))
    + u*(p[13] + v*(p[19] + v*(p[27]))
    + u*(p[18] + v*(p[26])
    + u*(p[25])))))*v;

  // Radial terms?
  if (i[I_TPDRAD]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11] + s*(p[23])));
  }

  return 0;
}

//----------------------------------------------------------------------------

int tpd7(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  if (i[I_TPDNCO+inverse] != 40 || 2 < Nhat) {
    return 1;
  }

  double r, s;
  double u = rawcrd[0];
  double v = rawcrd[1];

  // Auxiliary variables?
  if (i[I_TPDAUX]) {
    r = p[0] + p[1]*u + p[2]*v;
    v = p[3] + p[4]*u + p[5]*v;
    u = r;
    p += 6;
  }

  if (inverse) p += i[I_TPDNCO];

  // Seventh degree.
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7] + u*(p[12] + u*(p[17] + u*(p[24] + u*(p[31])))))));

  if (Nhat == 1) return 0;

  *discrd +=
      v*(p[2]  + v*(p[6]  + v*(p[10] + v*(p[16] + v*(p[22] + v*(p[30] + v*(p[38])))))))
    + u*(p[5]  + v*(p[9]  + v*(p[15] + v*(p[21] + v*(p[29] + v*(p[37])))))
    + u*(p[8]  + v*(p[14] + v*(p[20] + v*(p[28] + v*(p[36]))))
    + u*(p[13] + v*(p[19] + v*(p[27] + v*(p[35])))
    + u*(p[18] + v*(p[26] + v*(p[34]))
    + u*(p[25] + v*(p[33])
    + u*(p[32]))))))*v;

  // Radial terms?
  if (i[I_TPDRAD]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11] + s*(p[23] + s*(p[39]))));
  }

  return 0;
}

//----------------------------------------------------------------------------

int tpd8(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  if (i[I_TPDNCO+inverse] != 49 || 2 < Nhat) {
    return 1;
  }

  double r, s;
  double u = rawcrd[0];
  double v = rawcrd[1];

  // Auxiliary variables?
  if (i[I_TPDAUX]) {
    r = p[0] + p[1]*u + p[2]*v;
    v = p[3] + p[4]*u + p[5]*v;
    u = r;
    p += 6;
  }

  if (inverse) p += i[I_TPDNCO];

  // Eighth degree.
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7] + u*(p[12] + u*(p[17] + u*(p[24] + u*(p[31] + u*(p[40]))))))));

  if (Nhat == 1) return 0;

  *discrd +=
      v*(p[2]  + v*(p[6]  + v*(p[10] + v*(p[16] + v*(p[22] + v*(p[30] + v*(p[38] + v*(p[48]))))))))
    + u*(p[5]  + v*(p[9]  + v*(p[15] + v*(p[21] + v*(p[29] + v*(p[37] + v*(p[47]))))))
    + u*(p[8]  + v*(p[14] + v*(p[20] + v*(p[28] + v*(p[36] + v*(p[46])))))
    + u*(p[13] + v*(p[19] + v*(p[27] + v*(p[35] + v*(p[45]))))
    + u*(p[18] + v*(p[26] + v*(p[34] + v*(p[44])))
    + u*(p[25] + v*(p[33] + v*(p[43]))
    + u*(p[32] + v*(p[42])
    + u*(p[41])))))))*v;

  // Radial terms?
  if (i[I_TPDRAD]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11] + s*(p[23] + s*(p[39]))));
  }

  return 0;
}

//----------------------------------------------------------------------------

int tpd9(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  if (i[I_TPDNCO+inverse] != 60 || 2 < Nhat) {
    return 1;
  }

  double r, s;
  double u = rawcrd[0];
  double v = rawcrd[1];

  // Auxiliary variables?
  if (i[I_TPDAUX]) {
    r = p[0] + p[1]*u + p[2]*v;
    v = p[3] + p[4]*u + p[5]*v;
    u = r;
    p += 6;
  }

  if (inverse) p += i[I_TPDNCO];

  // Ninth degree.
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7] + u*(p[12] + u*(p[17] + u*(p[24] + u*(p[31] + u*(p[40] + u*(p[49])))))))));

  if (Nhat == 1) return 0;

  *discrd +=
      v*(p[2]  + v*(p[6]  + v*(p[10] + v*(p[16] + v*(p[22] + v*(p[30] + v*(p[38] + v*(p[48] + v*(p[58])))))))))
    + u*(p[5]  + v*(p[9]  + v*(p[15] + v*(p[21] + v*(p[29] + v*(p[37] + v*(p[47] + v*(p[57])))))))
    + u*(p[8]  + v*(p[14] + v*(p[20] + v*(p[28] + v*(p[36] + v*(p[46] + v*(p[56]))))))
    + u*(p[13] + v*(p[19] + v*(p[27] + v*(p[35] + v*(p[45] + v*(p[55])))))
    + u*(p[18] + v*(p[26] + v*(p[34] + v*(p[44] + v*(p[54]))))
    + u*(p[25] + v*(p[33] + v*(p[43] + v*(p[53])))
    + u*(p[32] + v*(p[42] + v*(p[52]))
    + u*(p[41] + v*(p[51])
    + u*(p[50]))))))))*v;

  // Radial terms?
  if (i[I_TPDRAD]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11] + s*(p[23] + s*(p[39] + s*(p[59])))));
  }

  return 0;
}
