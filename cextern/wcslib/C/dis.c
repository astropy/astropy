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
  $Id: dis.c,v 5.5.1.1 2015/05/05 13:21:00 mcalabre Exp mcalabre $
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

const int POLYNOMIAL = 0;
const int TPD        = 1;

/* Maximum number of DPja or DQia keywords. */
int NDPMAX = 256;

/* Map status return value to message. */
const char *dis_errmsg[] = {
  "Success",
  "Null disprm pointer passed",
  "Memory allocation failed",
  "Failed to initialize distortion functions",
  "Distort error",
  "De-distort error"};

/* Convenience macro for invoking wcserr_set(). */
#define DIS_ERRMSG(status) WCSERR_SET(status), dis_errmsg[status]

/*--------------------------------------------------------------------------*/

int disndp(int ndpmax) { if (ndpmax >= 0) NDPMAX = ndpmax; return NDPMAX; }

/*--------------------------------------------------------------------------*/

int dpfill(
  struct dpkey *dp,
  const char *keyword,
  const char *field,
  int j,
  int type,
  int ival,
  double fval)

{
  char *cp;

  if (keyword) {
    if (field) {
      sprintf(dp->field, "%s.%s", keyword, field);
    } else {
      strcpy(dp->field, keyword);
    }
  } else if (field) {
    strcpy(dp->field, field);
  }

  if (j) {
    dp->j = j;
  } else {
    /* The field name must either be given or preset. */
    if ((cp = strpbrk(dp->field, "0123456789")) != 0x0) {
      sscanf(cp, "%d.", &(dp->j));
    }
  }

  if ((dp->type = type)) {
    dp->value.f = fval;
  } else {
    dp->value.i = ival;
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int disini(int alloc, int naxis, struct disprm *dis)

{
  static const char *function = "disini";

  struct wcserr **err;

  if (dis == 0x0) return DISERR_NULL_POINTER;

  /* Initialize error message handling. */
  err = &(dis->err);
  if (dis->flag != -1) {
    if (dis->err) free(dis->err);
  }
  dis->err = 0x0;


  /* Initialize pointers. */
  if (dis->flag == -1 || dis->m_flag != DISSET) {
    if (dis->flag == -1) {
      dis->axmap  = 0x0;
      dis->Nhat   = 0x0;
      dis->offset = 0x0;
      dis->scale  = 0x0;
      dis->iparm  = 0x0;
      dis->dparm  = 0x0;

      dis->disp2x = 0x0;
      dis->disx2p = 0x0;
      dis->tmpmem = 0x0;

      dis->i_naxis = 0;
    }

    /* Initialize memory management. */
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


  /* Allocate memory for arrays if required. */
  if (alloc ||
      dis->dtype  == 0x0 ||
      (NDPMAX && dis->dp == 0x0) ||
      dis->maxdis == 0x0) {

    /* Was sufficient allocated previously? */
    if (dis->m_flag == DISSET &&
       (dis->m_naxis < naxis  ||
        dis->ndpmax  < NDPMAX)) {
      /* No, free it. */
      disfree(dis);
    }

    if (alloc || dis->dtype == 0x0) {
      if (dis->m_dtype) {
        /* In case the caller fiddled with it. */
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
        /* In case the caller fiddled with it. */
        dis->dp = dis->m_dp;

      } else {
        if (NDPMAX) {
          if ((dis->dp = calloc(NDPMAX, sizeof(struct dpkey))) == 0x0) {
            disfree(dis);
            return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
          }
        } else {
          dis->dp = 0x0;
        }

        dis->ndpmax  = NDPMAX;

        dis->m_flag  = DISSET;
        dis->m_naxis = naxis;
        dis->m_dp    = dis->dp;
      }
    }

    if (alloc || dis->maxdis == 0x0) {
      if (dis->m_maxdis) {
        /* In case the caller fiddled with it. */
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


  /* Set defaults. */
  dis->flag  = 0;
  dis->naxis = naxis;

  memset(dis->dtype,  0, naxis*sizeof(char [72]));
  dis->ndp = 0;
  memset(dis->dp,     0, NDPMAX*sizeof(struct dpkey));
  memset(dis->maxdis, 0, naxis*sizeof(double));
  dis->totdis = 0.0;

  return 0;
}


/*--------------------------------------------------------------------------*/

int discpy(int alloc, const struct disprm *dissrc, struct disprm *disdst)

{
  static const char *function = "discpy";

  int naxis, ndp, status;
  struct wcserr **err;

  if (dissrc == 0x0) return DISERR_NULL_POINTER;
  if (disdst == 0x0) return DISERR_NULL_POINTER;
  err = &(disdst->err);

  naxis = dissrc->naxis;
  if (naxis < 1) {
    return wcserr_set(WCSERR_SET(DISERR_MEMORY),
      "naxis must be positive (got %d)", naxis);
  }

  ndp = NDPMAX;
  NDPMAX = dissrc->ndpmax;

  if ((status = disini(alloc, naxis, disdst))) {
    return status;
  }

  NDPMAX = ndp;

  memcpy(disdst->dtype, dissrc->dtype, naxis*sizeof(char [72]));

  disdst->ndp = dissrc->ndp;
  disdst->ndpmax = dissrc->ndpmax;
  memcpy(disdst->dp, dissrc->dp, dissrc->ndpmax*sizeof(struct dpkey));

  memcpy(disdst->maxdis, dissrc->maxdis, naxis*sizeof(double));
  disdst->totdis = dissrc->totdis;

  return 0;
}

/*--------------------------------------------------------------------------*/

int disfree(struct disprm *dis)

{
  int j;

  if (dis == 0x0) return DISERR_NULL_POINTER;

  if (dis->flag != -1) {
    /* Optionally allocated by disini() for given parameters. */
    if (dis->m_flag == DISSET) {
      if (dis->dtype  == dis->m_dtype)  dis->dtype  = 0x0;
      if (dis->dp     == dis->m_dp)     dis->dp     = 0x0;
      if (dis->maxdis == dis->m_maxdis) dis->maxdis = 0x0;

      if (dis->m_dtype)  free(dis->m_dtype);
      if (dis->m_dp)     free(dis->m_dp);
      if (dis->m_maxdis) free(dis->m_maxdis);
    }

    /* Recall that these were allocated in bulk by disset(). */
    if (dis->axmap  && dis->axmap[0])  free(dis->axmap[0]);
    if (dis->offset && dis->offset[0]) free(dis->offset[0]);
    if (dis->scale  && dis->scale[0])  free(dis->scale[0]);

    if (dis->axmap)  free(dis->axmap);
    if (dis->Nhat)   free(dis->Nhat);
    if (dis->offset) free(dis->offset);
    if (dis->scale)  free(dis->scale);
    for (j = 0; j < dis->i_naxis; j++) {
      if (dis->iparm[j]) free(dis->iparm[j]);
      if (dis->dparm[j]) free(dis->dparm[j]);
    }
    if (dis->iparm)  free(dis->iparm);
    if (dis->dparm)  free(dis->dparm);

    if (dis->disp2x) free(dis->disp2x);
    if (dis->disx2p) free(dis->disx2p);
    if (dis->tmpmem) free(dis->tmpmem);

    if (dis->err) free(dis->err);
  }

  dis->m_flag   = 0;
  dis->m_naxis  = 0;
  dis->m_dtype  = 0x0;
  dis->m_dp     = 0x0;
  dis->m_maxdis = 0x0;

  dis->axmap  = 0x0;
  dis->Nhat   = 0x0;
  dis->offset = 0x0;
  dis->scale  = 0x0;
  dis->iparm  = 0x0;
  dis->dparm  = 0x0;
  dis->disp2x = 0x0;
  dis->disx2p = 0x0;
  dis->tmpmem = 0x0;

  dis->err  = 0x0;

  dis->flag = 0;

  return 0;
}

/*--------------------------------------------------------------------------*/

int disprt(const struct disprm *dis)

{
  char hext[32];
  int i, j, k, naxis;

  if (dis == 0x0) return DISERR_NULL_POINTER;

  if (dis->flag != DISSET) {
    wcsprintf("The disprm struct is UNINITIALIZED.\n");
    return 0;
  }

  naxis = dis->naxis;


  wcsprintf("       flag: %d\n", dis->flag);

  /* Parameters supplied. */
  wcsprintf("      naxis: %d\n", naxis);

  WCSPRINTF_PTR("      dtype: ", dis->dtype, "\n");
  for (j = 0; j < naxis; j++) {
    wcsprintf("             \"%s\"\n", dis->dtype[j]);
  }

  wcsprintf("        ndp: %d\n", dis->ndp);
  wcsprintf("     ndpmax: %d\n", dis->ndpmax);
  WCSPRINTF_PTR("         dp: ", dis->dp, "\n");
  for (i = 0; i < dis->ndp; i++) {
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
  for (j = 0; j < naxis; j++) {
    wcsprintf("  %#- 11.5g", dis->maxdis[j]);
  }
  wcsprintf("\n");

  wcsprintf("     totdis:  %#- 11.5g\n", dis->totdis);

  /* Derived values. */
  WCSPRINTF_PTR("      axmap: ", dis->axmap, "\n");
  for (j = 0; j < naxis; j++) {
    wcsprintf(" axmap[%d][]:", j);
    for (i = 0; i < naxis; i++) {
      wcsprintf("%6d", dis->axmap[j][i]);
    }
    wcsprintf("\n            ");
    for (i = naxis; i < 2*naxis; i++) {
      wcsprintf("%6d", dis->axmap[j][i]);
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("       Nhat: ", dis->Nhat, "\n");
  wcsprintf("            ");
  for (j = 0; j < naxis; j++) {
    wcsprintf("%6d", dis->Nhat[j]);
  }
  wcsprintf("\n");

  WCSPRINTF_PTR("     offset: ", dis->offset, "\n");
  for (j = 0; j < naxis; j++) {
    wcsprintf("offset[%d][]:", j);
    for (i = 0; i < naxis; i++) {
      wcsprintf("  %#- 11.5g", dis->offset[j][i]);
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("      scale: ", dis->scale, "\n");
  for (j = 0; j < naxis; j++) {
    wcsprintf(" scale[%d][]:", j);
    for (i = 0; i < naxis; i++) {
      wcsprintf("  %#- 11.5g", dis->scale[j][i]);
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("      iparm: ", dis->iparm, "\n");
  for (j = 0; j < naxis; j++) {
    wcsprintf(" iparm[%d][]:", j);
    for (k = 0; k < dis->iparm[j][1]; k++) {
      if (k && k%5 == 0) {
        wcsprintf("\n            ");
      }
      wcsprintf("  %11d", dis->iparm[j][k]);
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("      dparm: ", dis->dparm, "\n");
  for (j = 0; j < naxis; j++) {
    wcsprintf(" dparm[%d][]:", j);
    for (k = 0; k < dis->iparm[j][2]; k++) {
      if (k && k%5 == 0) {
        wcsprintf("\n            ");
      }
      wcsprintf("  %#- 11.5g", dis->dparm[j][k]);
    }
    wcsprintf("\n");
  }

  wcsprintf("    i_naxis: %d\n", dis->i_naxis);
  wcsprintf("       ndis: %d\n", dis->ndis);

  /* Error handling. */
  WCSPRINTF_PTR("        err: ", dis->err, "\n");
  if (dis->err) {
    wcserr_prt(dis->err, "             ");
  }

  /* Work arrays. */
  WCSPRINTF_PTR("     disp2x: ", dis->disp2x, "\n");
  for (j = 0; j < naxis; j++) {
    wcsprintf("  disp2x[%d]: %s", j,
      wcsutil_fptr2str((int (*)(void))dis->disp2x[j], hext));
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
  for (j = 0; j < naxis; j++) {
    wcsprintf("  disx2p[%d]: %s\n", j,
      wcsutil_fptr2str((int (*)(void))dis->disx2p[j], hext));
  }
  WCSPRINTF_PTR("     tmpmem: ", dis->tmpmem, "\n");

  /* Memory management. */
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

/*--------------------------------------------------------------------------*/

int disset(struct disprm *dis)

{
  static const char *function = "disset";

  char   *dpq, *fp;
  int    idp, j, jhat, k, naxis, ndis, Nhat, status;
  struct dpkey *keyp;
  struct wcserr **err;

  if (dis == 0x0) return DISERR_NULL_POINTER;
  err = &(dis->err);

  naxis = dis->naxis;


  /* Do basic checks. */
  if (dis->ndp < 0) {
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "disprm::ndp is negative (%d)", dis->ndp);
  }

  ndis = 0;
  for (j = 0; j < naxis; j++) {
    if (strlen(dis->dtype[j])) {
      ndis++;
      break;
    }
  }

  if (dis->ndp) {
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
        "No DPja or DQia keywords, NAXES at least is required for each"
        "distortion");
    }
  }


  /* Allocate memory for derived parameters and work arrays. */
  if (dis->i_naxis < naxis) {
    if (dis->i_naxis) {
      /* Recall that axmap, offset, and scale are allocated in bulk. */
      free(dis->axmap[0]);
      free(dis->axmap);
      free(dis->Nhat);
      free(dis->offset[0]);
      free(dis->offset);
      free(dis->scale[0]);
      free(dis->scale);

      for (j = 0; j < dis->i_naxis; j++) {
        /* Memory allocated separately for each axis. */
        if (dis->iparm[j]) free(dis->iparm[j]);
        if (dis->dparm[j]) free(dis->dparm[j]);
      }
      free(dis->iparm);
      free(dis->dparm);

      free(dis->disp2x);
      free(dis->disx2p);

      free(dis->tmpmem);
    }

    if ((dis->axmap = calloc(naxis, sizeof(int *))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    /* Allocate axmap[][] in bulk and then carve it up. */
    if ((dis->axmap[0] = calloc(2*naxis*naxis, sizeof(int))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    for (j = 1; j < naxis; j++) {
      dis->axmap[j] = dis->axmap[j-1] + 2*naxis;
    }

    if ((dis->Nhat = calloc(naxis, sizeof(int *))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    if ((dis->offset = calloc(naxis, sizeof(double *))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    /* Allocate offset[][] in bulk and then carve it up. */
    if ((dis->offset[0] = calloc(naxis*naxis, sizeof(double))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    for (j = 1; j < naxis; j++) {
      dis->offset[j] = dis->offset[j-1] + naxis;
    }

    if ((dis->scale = calloc(naxis, sizeof(double *))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    /* Allocate scale[][] in bulk and then carve it up. */
    if ((dis->scale[0] = calloc(naxis*naxis, sizeof(double))) == 0x0) {
      disfree(dis);
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }

    for (j = 1; j < naxis; j++) {
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

  /* Start with a clean slate. */
  for (j = 0; j < 2*naxis*naxis; j++) {
    dis->axmap[0][j] = -1;
  }

  memset(dis->Nhat,      0, naxis*sizeof(int));
  memset(dis->offset[0], 0, naxis*naxis*sizeof(double));

  for (j = 0; j < naxis*naxis; j++) {
    dis->scale[0][j] = 1.0;
  }

  /* polyset() etc. must look after iparm[][] and dparm[][]. */

  dis->i_naxis = naxis;
  dis->ndis    = 0;

  memset(dis->disp2x, 0, naxis*sizeof(int (*)(DISP2X_ARGS)));
  memset(dis->disx2p, 0, naxis*sizeof(int (*)(DISX2P_ARGS)));
  memset(dis->tmpmem, 0, naxis*sizeof(double));


  /* Handle DPja or DQia keywords common to all distortions. */
  keyp = dis->dp;
  for (idp = 0; idp < dis->ndp; idp++, keyp++) {
    /* Check that they're all one kind or the other. */
    if (dis->dp[0].field[1] != dpq[1]) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "disprm::dp appears to contain a mix of DPja and DQia keys");
    }

    j = keyp->j;

    if (j < 1 || naxis < j) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid axis number (%d) in %s", j, keyp->field);
    }

    if ((fp = strpbrk(keyp->field, ".")) == 0x0) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid record field name: %s", j, keyp->field);
    }
    fp++;

    j--;
    if (strncmp(fp, "NAXES", 6) == 0) {
      Nhat = wcsutil_dpkey_int(keyp);
      if (Nhat < 0 || naxis < Nhat) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid value of Nhat for %s distortion in %s: %d", dis->dtype[j],
          keyp->field, Nhat);
      }

      dis->Nhat[j] = Nhat;

    } else if (strncmp(fp, "AXIS.", 5) == 0) {
      sscanf(fp+5, "%d", &jhat);
      if (jhat < 0 || naxis < jhat) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid axis in axis map for %s distortion in %s: %d",
          dis->dtype[j], keyp->field, jhat);
      }

      /* N.B. axis numbers in the map are 0-relative. */
      dis->axmap[j][jhat-1] = wcsutil_dpkey_int(keyp) - 1;

    } else if (strncmp(fp, "OFFSET.", 7) == 0) {
      sscanf(fp+7, "%d", &jhat);
      dis->offset[j][jhat-1] = wcsutil_dpkey_double(keyp);

    } else if (strncmp(fp, "SCALE.", 6) == 0) {
      sscanf(fp+6, "%d", &jhat);
      dis->scale[j][jhat-1] = wcsutil_dpkey_double(keyp);
    }
  }

  /* Set defaults and do sanity checks on axmap[][].  */
  for (j = 0; j < naxis; j++) {
    if (strlen(dis->dtype[j]) == 0) {
      /* No distortion on this axis, check that there are no parameters. */
      keyp = dis->dp;
      for (idp = 0; idp < dis->ndp; idp++, keyp++) {
        if (keyp->j == j+1) {
          return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
            "No distortion type, yet %s keyvalues are present for axis %d",
            dpq, j+1);
        }
      }

      continue;
    }

    /* N.B. NAXES (Nhat) has no default value. */
    if (dis->Nhat[j] <= 0) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "%s.NAXES was not set (or bad) for %s distortion on axis %d",
        dpq, dis->dtype[j], j+1);
    }

    /* Set defaults for axmap[][]. */
    Nhat = dis->Nhat[j];
    for (jhat = 0; jhat < Nhat; jhat++) {
      if (dis->axmap[j][jhat] == -1) {
        dis->axmap[j][jhat] = jhat;
      }
    }

    /* Sanity check on the length of the axis map. */
    Nhat = 0;
    for (jhat = 0; jhat < naxis; jhat++) {
      if (dis->axmap[j][jhat] != -1) Nhat = jhat+1;
    }

    if (Nhat != dis->Nhat[j]) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Mismatch in length of axis map for %s distortion on axis %d",
        dis->dtype[j], j+1);
    }

    /* Check uniqueness of entries in the axis map. */
    for (jhat = 0; jhat < Nhat; jhat++) {
      for (k = 0; k < jhat; k++) {
        if (dis->axmap[j][jhat] == dis->axmap[j][k]) {
          return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
            "Duplicated entry in axis map for %s distortion on axis %d",
            dis->dtype[j], j+1);
        }
      }
    }

    /* Construct the inverse axis map. */
    for (jhat = 0; jhat < Nhat; jhat++) {
      k = naxis + dis->axmap[j][jhat];
      dis->axmap[j][k] = jhat;
    }
  }


  /* Identify the distortion functions. */
  ndis = 0;
  for (j = 0; j < naxis; j++) {
    if (strlen(dis->dtype[j]) == 0) {
      /* No distortion on this axis. */
      continue;
    }

    if (dis->Nhat[j] == 0) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Empty axis map for %s distortion on axis %d", dis->dtype[j], j+1);
    }

    /* Invoke the specific setup functions for each distortion. */
    if (strcmp(dis->dtype[j], "TPD") == 0) {
      /* Template Polynomial Distortion". */
      if ((status = tpdset(j, dis))) {
        /* (Preserve the error message set by tpdset().) */
        return status;
      }

    } else if (strcmp(dis->dtype[j], "TPV") == 0) {
      /* TPV "projection". */
      if ((status = tpvset(j, dis))) {
        /* (Preserve the error message set by tpvset().) */
        return status;
      }

    } else if (strcmp(dis->dtype[j], "SIP") == 0) {
      /* Simple Imaging Polynomial (SIP). */
      if ((status = sipset(j, dis))) {
        /* (Preserve the error message set by sipset().) */
        return status;
      }

    } else if (strcmp(dis->dtype[j], "Polynomial")  == 0 ||
               strcmp(dis->dtype[j], "Polynomial*") == 0) {
      /* General polynomial distortion. */
      if ((status = polyset(j, dis))) {
        /* (Preserve the error message set by polyset().) */
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

/*--------------------------------------------------------------------------*/

int disp2x(
  struct disprm *dis,
  const double rawcrd[],
  double discrd[])

{
  static const char *function = "disp2x";

  int axisj, j, jhat, naxis, Nhat, status;
  double dtmp, *offset, *scale, *tmpcrd;
  struct wcserr **err;


  /* Initialize. */
  if (dis == 0x0) return DISERR_NULL_POINTER;
  err = &(dis->err);

  if (dis->flag != DISSET) {
    if ((status = disset(dis))) return status;
  }

  naxis = dis->naxis;


  /* Invoke the distortion functions for each axis. */
  tmpcrd = dis->tmpmem;
  for (j = 0; j < naxis; j++) {
    offset = dis->offset[j];
    scale  = dis->scale[j];

    if (dis->disp2x[j]) {
      Nhat = dis->Nhat[j];
      for (jhat = 0; jhat < Nhat; jhat++) {
        axisj = dis->axmap[j][jhat];
        tmpcrd[jhat] = (rawcrd[axisj] - offset[jhat])*scale[jhat];
      }

      if ((status = (dis->disp2x[j])(0, dis->iparm[j], dis->dparm[j], Nhat,
                                     tmpcrd, &dtmp))) {
        return wcserr_set(DIS_ERRMSG(DISERR_DISTORT));
      }

      discrd[j] = dtmp;

    } else {
      discrd[j] = rawcrd[j];
    }
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int disx2p(
  struct disprm *dis,
  const double discrd[],
  double rawcrd[])

{
  static const char *function = "disx2p";

  const int ITERMAX = 30;
  const double TOL = 1.0e-13;

  int convergence, iter, j, jhat, naxis, Nhat, status;
  double dd, *dcrd0, *dcrd1, *delta, residual, *rcrd1, rtmp;
  struct wcserr **err;


  /* Initialize. */
  if (dis == 0x0) return DISERR_NULL_POINTER;
  err = &(dis->err);

  naxis = dis->naxis;

  /* Carve up working memory, noting that disp2x() gets to it first. */
  dcrd0 = dis->tmpmem + naxis;
  dcrd1 = dcrd0 + naxis;
  rcrd1 = dcrd1 + naxis;
  delta = rcrd1 + naxis;


  /* Zeroth approximation.  The assumption here and below is that the     */
  /* distortion is small so that, to first order in the neighbourhood of  */
  /* the solution, discrd[j] ~= a + b*rawcrd[j], i.e. independent of      */
  /* rawcrd[i], where i != j.  This is effectively equivalent to assuming */
  /* that the distortion functions are separable to first order.          */
  /* Furthermore, a is assumed to be small, and b close to unity.         */
  memcpy(rawcrd, discrd, naxis*sizeof(double));

  /* If available, use disprm::disx2p to improve the zeroth approximation. */
  for (j = 0; j < naxis; j++) {
    if (dis->disx2p[j]) {
      Nhat = dis->Nhat[j];
      if ((status = (dis->disx2p[j])(1, dis->iparm[j], dis->dparm[j], Nhat,
                                     discrd, &rtmp))) {
        return wcserr_set(DIS_ERRMSG(DISERR_DEDISTORT));
      }

      /* Using the inverse axis map. */
      jhat = dis->axmap[j][naxis+j];
      rawcrd[j] = dis->offset[j][jhat] + rtmp/dis->scale[j][jhat];
    }
  }

  return 0;

  /* Iteratively invert the (well-behaved!) distortion function. */
  for (iter = 0; iter < ITERMAX; iter++) {
    if ((status = disp2x(dis, rawcrd, dcrd0))) {
      return wcserr_set(DIS_ERRMSG(status));
    }

    /* Check for convergence. */
    convergence = 1;
    for (j = 0; j < naxis; j++) {
      delta[j] = discrd[j] - dcrd0[j];

      if (fabs(discrd[j]) < 1.0) {
        dd = delta[j];
      } else {
        /* TOL may be below the precision achievable from floating point */
        /* subtraction, so switch to a fractional tolerance.             */
        dd = delta[j] / discrd[j];
      }

      if (TOL < fabs(dd)) {
        /* No convergence yet on this axis. */
        convergence = 0;
      }
    }

    if (convergence) break;

    /* Determine a suitable test point for computing the gradient. */
    for (j = 0; j < naxis; j++) {
      /* Constrain the displacement. */
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

    if (iter < ITERMAX/2) {
      /* With the assumption of small distortions (as above), the gradient */
      /* of discrd[j] should be dominated by the partial derivative with   */
      /* respect  to rawcrd[j], and we can neglect partials with respect   */
      /* to rawcrd[i], where i != j.  Thus only one test point is needed,  */
      /* not one for each axis.                                            */
      for (j = 0; j < naxis; j++) {
        rcrd1[j] = rawcrd[j] + delta[j];
      }

      /* Compute discrd[] at the test point. */
      if ((status = disp2x(dis, rcrd1, dcrd1))) {
        return wcserr_set(DIS_ERRMSG(status));
      }

      /* Compute the next approximation. */
      for (j = 0; j < naxis; j++) {
        rawcrd[j] += (discrd[j] - dcrd0[j]) *
                        (delta[j]/(dcrd1[j] - dcrd0[j]));
      }

    } else {
      /* Convergence should not take more than seven or so iterations.  As */
      /* it is slow, try computing the gradient in full.                   */
      memcpy(rcrd1, rawcrd, naxis*sizeof(double));

      for (j = 0; j < naxis; j++) {
        rcrd1[j] += delta[j];

        /* Compute discrd[] at the test point. */
        if ((status = disp2x(dis, rcrd1, dcrd1))) {
          return wcserr_set(DIS_ERRMSG(status));
        }

        /* Compute the next approximation. */
        rawcrd[j] += (discrd[j] - dcrd0[j]) *
                       (delta[j]/(dcrd1[j] - dcrd0[j]));

        rcrd1[j] -= delta[j];
      }
    }
  }


  if (!convergence) {
    residual = 0.0;
    for (j = 0; j < naxis; j++) {
      dd = discrd[j] - dcrd0[j] ;
      residual += dd*dd;
    }
    residual = sqrt(residual);

    return wcserr_set(WCSERR_SET(DISERR_DEDISTORT),
      "Convergence not achieved after %d iterations, residual %#7.2g", iter,
        residual);
  }


  return 0;
}

/*--------------------------------------------------------------------------*/

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

  int carry, j, naxis, status = 0;
  double dpix, dpx2, dssq, *pix0, *pix1, *pixend, *pixinc, pixspan, *ssqdis,
         ssqtot, *sumdis, sumtot, totdis;
  struct wcserr **err;


  /* Initialize. */
  if (dis == 0x0) return DISERR_NULL_POINTER;
  err = &(dis->err);

  naxis = dis->naxis;

  if (nsamp) *nsamp = 0;
  for (j = 0; j < naxis; j++) {
    if (maxdis) maxdis[j] = 0.0;
    if (avgdis) avgdis[j] = 0.0;
    if (rmsdis) rmsdis[j] = 0.0;
  }
  if (maxtot) *maxtot = 0.0;
  if (avgtot) *avgtot = 0.0;
  if (rmstot) *rmstot = 0.0;

  /* Quick return if no distortions. */
  if (dis->ndis == 0) return 0;

  /* Carve up working memory, noting that disp2x() gets to it first. */
  pixinc = dis->tmpmem + naxis;
  pixend = pixinc + naxis;
  sumdis = pixend + naxis;
  ssqdis = sumdis + naxis;

  /* Work out increments on each axis. */
  for (j = 0; j < naxis; j++) {
    pixspan = pixtrc[j] - (pixblc ? pixblc[j] : 1.0);

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

  /* Get some more memory for coordinate vectors. */
  if ((pix0 = calloc(2*naxis, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  pix1 = pix0 + naxis;


  /* Set up the array of pixel coordinates. */
  for (j = 0; j < naxis; j++) {
    pix0[j] = pixblc ? pixblc[j] : 1.0;
    pixend[j] = pixtrc[j] + 0.5*pixinc[j];
  }

  /* Initialize accumulators. */
  for (j = 0; j < naxis; j++) {
    sumdis[j] = 0.0;
    ssqdis[j] = 0.0;
  }
  sumtot = 0.0;
  ssqtot = 0.0;


  /* Loop over N dimensions. */
  carry = 0;
  while (carry == 0) {
    if ((status = disp2x(dis, pix0, pix1))) {
      /* (Preserve the error message set by disp2x().) */
      goto cleanup;
    }

    /* Accumulate statistics. */
    (*nsamp)++;

    dssq = 0.0;
    for (j = 0; j < naxis; j++) {
      dpix = pix1[j] - pix0[j];
      dpx2 = dpix*dpix;

      sumdis[j] += dpix;
      ssqdis[j] += dpx2;

      if (maxdis && (dpix = fabs(dpix)) > maxdis[j]) {
        maxdis[j] = dpix;
      }

      dssq += dpx2;
    }

    totdis = sqrt(dssq);
    sumtot += totdis;
    ssqtot += totdis*totdis;

    if (maxtot && *maxtot < totdis) {
      *maxtot = totdis;
    }

    /* Next pixel. */
    for (j = 0; j < naxis; j++) {
      pix0[j] += pixinc[j];
      if (pix0[j] < pixend[j]) {
        carry = 0;
        break;
      }

      pix0[j] = pixblc ? pixblc[j] : 1.0;
      carry = 1;
    }
  }


  /* Compute the means and RMSs. */
  for (j = 0; j < naxis; j++) {
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

/*--------------------------------------------------------------------------*/

int polyset(int j, struct disprm *dis)

{
  static const char *function = "polyset";

  char   *fp, id[32];
  int    i, idp, *iparm, ipow, ivar, jhat, k, K, lendp, m, M, naxis, ndparm,
         Nhat, niparm, nKparm, npow, nTparm, nVar, offset;
  double *dparm, *dptr, power;
  struct dpkey *keyp;
  struct wcserr **err;


  /* Initialize. */
  if (dis == 0x0) return DISERR_NULL_POINTER;
  err = &(dis->err);

  naxis = dis->naxis;
  sprintf(id, "Polynomial on axis %d", j+1);


  /* Find the number of auxiliary variables and terms. */
  K = 0;
  M = 0;
  keyp = dis->dp;
  for (idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    if ((fp = strpbrk(keyp->field, ".")) == 0x0) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid field name for %s: %s", id, keyp->field);
    }
    fp++;

    if (strcmp(fp, "NAUX") == 0) {
      K = wcsutil_dpkey_int(keyp);
    } else if (strcmp(fp, "NTERMS") == 0) {
      M = wcsutil_dpkey_int(keyp);
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

  Nhat = dis->Nhat[j];
  nKparm = 2*(Nhat + 1);
  nVar   = Nhat + K;
  nTparm = 1 + nVar;
  ndparm = K*nKparm + M*nTparm;

#define I_DTYPE   0	/* Distortion type code.                            */
#define I_NIPARM  1	/* Full (allocated) length of iparm[].              */
#define I_NDPARM  2	/* No. of parameters in dparm[], excl. work space.  */
#define I_NIDX    3	/* No. of indexes in iparm[].                       */
#define I_LENDP   4	/* Full (allocated) length of dparm[].              */
#define I_K       5	/* No. of independent variables.                    */
#define I_M       6	/* No. of terms in the polynomial.                  */
#define I_NKPARM  7	/* No. of parameters used to define each auxiliary. */
#define I_NTPARM  8	/* No. of parameters used to define each term.      */
#define I_NVAR    9	/* No. of independent + auxiliary variables.        */
#define I_MNVAR  10	/* No. of powers (exponents) in the polynomial.     */
#define I_DPOLY  11	/* dparm offset for polynomial coefficients.        */
#define I_DAUX   12	/* dparm offset for auxiliary coefficients.         */
#define I_DVPOW  13	/* dparm offset for integral powers of variables.   */
#define I_MAXPOW 14	/* iparm offset for max powers.                     */
#define I_DPOFF  15	/* iparm offset for dparm offsets.                  */
#define I_FLAGS  16	/* iparm offset for flags.                          */
#define I_IPOW   17	/* iparm offset for integral powers.                */
#define NIDX     18

  /* Add extra for integer exponents.  See "Optimization" below. */
  niparm = NIDX + (2 + 2*M)*nVar;

  /* Add extra memory for temporaries. */
  lendp = ndparm + K;

  /* Allocate memory for the indexes and parameter array. */
  if ((dis->iparm[j] = calloc(niparm, sizeof(int))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  if ((dis->dparm[j] = calloc(lendp, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  /* These help a bit to stop the code from turning into hieroglyphics. */
  iparm = dis->iparm[j];
  dparm = dis->dparm[j];


  /* Record the indexing parameters.  The first two are used by disprt(). */
  iparm[I_DTYPE]  = POLYNOMIAL;
  iparm[I_NIPARM] = niparm;
  iparm[I_NDPARM] = ndparm;
  iparm[I_NIDX]   = NIDX;
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

  /* Set default values of POWER for the auxiliary variables. */
  dptr = dparm + (1 + Nhat);
  for (k = 0; k < K; k++) {
    for (jhat = 0; jhat <= Nhat; jhat++) {
      dptr[jhat] = 1.0;
    }
    dptr += nKparm;
  }

  /* Set default values of COEFF for the independent variables. */
  dptr = dparm + iparm[I_DPOLY];
  for (m = 0; m < M; m++) {
    *dptr = 1.0;
    dptr += nTparm;
  }

  /* Extract parameter values from DPja or DQia. */
  k = m = 0;
  keyp = dis->dp;
  for (idp = 0; idp < dis->ndp; idp++, keyp++) {
    /* N.B. keyp->j is 1-relative, but j is 0-relative. */
    if (keyp->j-1 != j) continue;

    fp = strpbrk(keyp->field, ".") + 1;

    if (strncmp(fp, "AUX.", 4) == 0) {
      /* N.B. k here is 1-relative. */
      fp += 4;
      sscanf(fp, "%d", &k);
      if (k < 1 || K < k) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Bad auxiliary variable (%d) for %s: %s", k, id, keyp->field);
      }

      if ((fp = strpbrk(fp, ".")) == 0x0) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid field name for %s: %s", id, keyp->field);
      }
      fp++;

      if (strncmp(fp, "COEFF.", 6) == 0) {
        offset = 0;

      } else if (strncmp(fp, "POWER.", 6) == 0) {
        offset = 1 + Nhat;

      } else {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Unrecognized field name for %s: %s", id, keyp->field);
      }

      fp += 6;
      sscanf(fp, "%d", &jhat);
      if (jhat < 0 || naxis < jhat) {
        /* N.B. jhat == 0 is ok. */
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid axis number (%d) for %s: %s", jhat, id, keyp->field);
      }

      i = (k-1)*nKparm + offset + jhat;
      dparm[i] = wcsutil_dpkey_double(keyp);

    } else if (strncmp(fp, "TERM.", 5) == 0) {
      /* N.B. m here is 1-relative. */
      fp += 5;
      sscanf(fp, "%d", &m);
      if (m < 1 || M < m) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Bad term (%d) for %s: %s", m, id, keyp->field);
      }

      if ((fp = strpbrk(fp, ".")) == 0x0) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid field name for %s: %s", id, keyp->field);
      }
      fp++;

      if (strcmp(fp, "COEFF") == 0) {
        i = iparm[I_DPOLY] + (m-1)*nTparm;
        dparm[i] = wcsutil_dpkey_double(keyp);

      } else if (strncmp(fp, "VAR.", 4) == 0) {
        /* N.B. jhat here is 1-relative. */
        fp += 4;
        sscanf(fp, "%d", &jhat);
        if (jhat < 1 || naxis < jhat) {
          return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Invalid axis number (%d) for %s: %s", jhat, id, keyp->field);
        }

        i = iparm[I_DPOLY] + (m-1)*nTparm + 1 + (jhat-1);
        power = wcsutil_dpkey_double(keyp);
        dparm[i] = power;

      } else if (strncmp(fp, "AUX.", 4) == 0) {
        /* N.B. k here is 1-relative. */
        fp += 4;
        sscanf(fp, "%d", &k);
        if (k < 1 || K < k) {
          return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
            "Bad auxiliary variable (%d) for %s: %s", k, id, keyp->field);
        }

        i = iparm[I_DPOLY] + (m-1)*nTparm + 1 + Nhat + (k-1);
        power = wcsutil_dpkey_double(keyp);
        dparm[i] = power;

      } else {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Unrecognized field name for %s: %s", id, keyp->field);
      }

    } else if (strcmp(fp, "NAXES")  &&
              strncmp(fp, "AXIS.",   5) &&
              strncmp(fp, "OFFSET.", 7) &&
              strncmp(fp, "SCALE.",  6) &&
               strcmp(fp, "NAUX")   &&
               strcmp(fp, "NTERMS")) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Unrecognized field name for %s: %s", id, keyp->field);
    }
  }


  /* Optimization: when the power is integral, it is faster to multiply     */
  /* ------------  repeatedly than call pow().  iparm[] is constructed as   */
  /*               follows:                                                 */
  /*     NIDX indexing parameters, as above,                                */
  /*     nVar elements record the largest integral power for each variable, */
  /*     nVar elements record offsets into dparm for each variable,         */
  /*   M*nVar flags to signal whether the power is integral,                */
  /*   M*nVar integral powers.                                              */
  for (ivar = 0; ivar < nVar; ivar++) {
    /* Want at least the first degree power for all variables. */
    i = iparm[I_MAXPOW] + ivar;
    iparm[i] = 1;
  }

  for (ivar = 0; ivar < nVar; ivar++) {
    for (m = 0; m < M; m++) {
      i = iparm[I_DPOLY] + m*nTparm + 1 + ivar;
      power = dparm[i];

      /* Is it integral?  (Positive, negative, or zero.) */
      ipow = (int)power;
      if (power == (double)ipow) {
        /* Signal that the power is integral. */
        i = iparm[I_FLAGS] + m*nVar + ivar;
        if (ipow == 0) {
          iparm[i] = 3;
        } else {
          iparm[i] = 1;
        }

        /* The integral power itself. */
        i = iparm[I_IPOW] + m*nVar + ivar;
        iparm[i] = ipow;
      }

      /* Record the largest integral power for each variable. */
      i = iparm[I_MAXPOW] + ivar;
      if (iparm[i] < abs(ipow)) {
        iparm[i] = abs(ipow);
      }
    }
  }

  /* How many of all powers of each variable will there be? */
  npow = 0;
  for (ivar = 0; ivar < nVar; ivar++) {
    /* Offset into dparm. */
    i = iparm[I_DPOFF] + ivar;
    iparm[i] = lendp + npow;

    i = iparm[I_MAXPOW] + ivar;
    npow += iparm[i];
  }

  /* Expand dparm to store the extra powers. */
  if (npow) {
    lendp += npow;
    iparm[I_LENDP] = lendp;
    if ((dis->dparm[j] = realloc(dparm, lendp*sizeof(double))) == 0x0) {
      return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
    }
  }

  /* No specialist de-distortions. */
  dis->disp2x[j] = dispoly;
  dis->disx2p[j] = 0x0;

  /* Translate Polynomial to TPD if possible, it's much faster.  */
  /* However don't do it if the name was given as "Polynomial*". */
  if (strcmp(dis->dtype[j], "Polynomial") == 0) {
    pol2tpd(j, dis);
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

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

  int deg, degree, *iflgp, *iparm, *ipowp, jhat, K, m, n, ndparm, Nhat,
      niparm, p[2], *tpd_iparm;
  double *dparm, *dpolp, *tpd_dparm;
  struct wcserr **err;

  /* Initialize. */
  if (dis == 0x0) return DISERR_NULL_POINTER;
  err = &(dis->err);

  iparm = dis->iparm[j];
  dparm = dis->dparm[j];


  /* Check the number of independent variables, no more than two. */
  Nhat = dis->Nhat[j];
  if (2 < Nhat) return -1;

  /* Check auxiliaries: only one is allowed... */
  K = iparm[I_K];
  if (1 < K) return -1;
  if (K) {
    /* ...and it must be radial. */
    if (dparm[0] != 0.0) return -1;
    if (dparm[1] != 1.0) return -1;
    if (dparm[2] != 1.0) return -1;
    if (dparm[3] != 0.5) return -1;
    if (dparm[4] != 2.0) return -1;
    if (dparm[5] != 2.0) return -1;
  }

  /* Check powers... */
  iflgp = iparm + iparm[I_FLAGS];
  ipowp = iparm + iparm[I_IPOW];
  degree = 0;
  for (m = 0; m < iparm[I_M]; m++) {
    deg = 0;
    for (jhat = 0; jhat < Nhat; jhat++) {
      /* ...they must be positive integral. */
      if (*iflgp == 0)  return -1;
      if (*ipowp < 0)   return -1;
      deg += *ipowp;
      iflgp++;
      ipowp++;
    }

    /* The polynomial degree can't be greater than 9. */
    if (9 < deg) return -1;

    if (K) {
      /* Likewise for the radial variable. */
      if (*iflgp == 0)  return -1;
      if (*ipowp) {
        if (*ipowp < 0) return -1;
        if (9 < *ipowp) return -1;

        /* Can't mix the radial and other terms. */
        if (deg)        return -1;

        /* Can't have even powers of the radial variable. */
        deg = *ipowp;
        if (!(deg%2))   return -1;
      }
      iflgp++;
      ipowp++;
    }

    if (degree < deg) degree = deg;
  }


  /* OK, it ticks all the boxes.  Now translate it. */
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

  /* No specialist de-distortions. */
  dis->disx2p[j] = 0x0;

  /* Record indexing parameters. */
  niparm = 5;
  if ((tpd_iparm = calloc(niparm, sizeof(int))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  /* The first three are required by disprt(). */
  tpd_iparm[0] = TPD;
  tpd_iparm[1] = niparm;
  tpd_iparm[2] = ndparm;
  tpd_iparm[3] = ndparm;

  /* Flag for presence of radial terms. */
  tpd_iparm[4] = K;


  /* Allocate memory for the polynomial coefficients and fill it. */
  if ((tpd_dparm = calloc(ndparm, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  ipowp = iparm + iparm[I_IPOW];
  dpolp = dparm + iparm[I_DPOLY];
  for (m = 0; m < iparm[I_M]; m++) {
    if (K && ipowp[Nhat]) {
      /* The radial variable. */
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
      /* The independent variables. */
      p[0] = p[1] = 0;
      for (jhat = 0; jhat < Nhat; jhat++) {
        p[jhat] = ipowp[jhat];
      }

      n = map[p[0]][p[1]];
      tpd_dparm[n] = *dpolp;
    }


    ipowp += iparm[I_NVAR];
    dpolp += iparm[I_NVAR] + 1;
  }


  /* Switch from Polynomial to TPD. */
  free(iparm);
  free(dparm);
  dis->iparm[j] = tpd_iparm;
  dis->dparm[j] = tpd_dparm;

  return 0;
}

/*--------------------------------------------------------------------------*/

int dispoly(
  int inverse,
  const int iparm[],
  const double dparm[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  const int *iflgp, *imaxp, *imaxpow, *ipowp;
  int    ip, ivar, jhat, k, m;
  const double *cptr, *dpolp, *pptr;
  double *aux, auxp0, *dvarpow, *dpowp, term, var;

  /* Check for zeroes. */
  for (jhat = 0; jhat < Nhat; jhat++) {
    if (rawcrd[jhat] == 0.0) {
      *discrd = 0.0;
      return 0;
    }
  }

  /* Working memory for auxiliaries &c. was allocated at the end of p[]. */
  aux = (double *)(dparm + iparm[I_DAUX]);

  /* Compute the auxiliary variables. */
  for (k = 0; k < iparm[I_K]; k++) {
    cptr = dparm + k*iparm[I_NKPARM];
    pptr = cptr + (1+Nhat);

    aux[k] = *(cptr++);
    auxp0  = *(pptr++);

    for (jhat = 0; jhat < Nhat; jhat++) {
      aux[k] += *(cptr++)*pow(rawcrd[jhat], *(pptr++));
    }

    aux[k] = pow(aux[k], auxp0);

    /* Check for zeroes. */
    if (aux[k] == 0.0) {
      *discrd = 0.0;
      return 0;
    }
  }


  /* Compute all required integral powers of the variables. */
  imaxpow = iparm + iparm[I_MAXPOW];
  dvarpow = (double *)(dparm + iparm[I_DVPOW]);

  imaxp = imaxpow;
  dpowp = dvarpow;
  for (jhat = 0; jhat < Nhat; jhat++, imaxp++) {
    var = 1.0;
    for (ip = 0; ip < *imaxp; ip++, dpowp++) {
      var *= rawcrd[jhat];
      *dpowp = var;
    }
  }

  for (k = 0; k < iparm[I_K]; k++, imaxp++) {
    var = 1.0;
    for (ip = 0; ip < *imaxp; ip++, dpowp++) {
      var *= aux[k];
      *dpowp = var;
    }
  }

  /* Loop for each term of the polynomial. */
  *discrd = 0.0;
  iflgp = iparm + iparm[I_FLAGS];
  ipowp = iparm + iparm[I_IPOW];
  dpolp = dparm + iparm[I_DPOLY];
  for (m = 0; m < iparm[I_M]; m++) {
    term = *(dpolp++);

    /* Loop over all variables. */
    imaxp = imaxpow;
    dpowp = dvarpow - 1;
    for (ivar = 0; ivar < iparm[I_NVAR]; ivar++) {
      if (*iflgp & 2) {
        /* Nothing (zero power). */

      } else if (*iflgp) {
        /* Integral power. */
        if (*ipowp < 0) {
          /* Negative. */
          term /= dpowp[*ipowp];
        } else {
          /* Positive. */
          term *= dpowp[*ipowp];
        }

      } else {
        /* Fractional power. */
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

/*--------------------------------------------------------------------------*/

int tpvset(int j, struct disprm *dis)

{
  static const char *function = "tpvset";

  char   *fp, id[16];
  int    doradial, idp, k, ndparm, niparm;
  struct dpkey *keyp;
  struct wcserr **err;

  if (dis == 0x0) return DISERR_NULL_POINTER;
  err = &(dis->err);

  sprintf(id, "TPV on axis %d", j+1);


  /* TPV "projection". */
  if (dis->Nhat[j] != 2) {
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "Axis map for %s must contain 2 entries, not %d", id, dis->Nhat[j]);
  }

  /* Find the number of parameters. */
  ndparm   = 0;
  doradial = 0;
  keyp = dis->dp;
  for (idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    fp = strpbrk(keyp->field, ".") + 1;

    if (strncmp(fp, "TPV.", 4) == 0) {
      sscanf(fp+4, "%d", &k);
      if (0 <= k && k <= 39) {
        if (ndparm < k+1) ndparm = k+1;

	/* Any radial terms? */
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

  /* TPD is going to do the dirty work. */
  if (ndparm <= 4) {
    /* First degree. */
    ndparm = 4;
    dis->disp2x[j] = tpd1;
  } else if (ndparm <= 7) {
    /* Second degree. */
    ndparm = 7;
    dis->disp2x[j] = tpd2;
  } else if (ndparm <= 12) {
    /* Third degree. */
    ndparm = 12;
    dis->disp2x[j] = tpd3;
  } else if (ndparm <= 17) {
    /* Fourth degree. */
    ndparm = 17;
    dis->disp2x[j] = tpd4;
  } else if (ndparm <= 24) {
    /* Fifth degree. */
    ndparm = 24;
    dis->disp2x[j] = tpd5;
  } else if (ndparm <= 31) {
    /* Sixth degree. */
    ndparm = 31;
    dis->disp2x[j] = tpd6;
  } else if (ndparm <= 40) {
    /* Seventh degree. */
    ndparm = 40;
    dis->disp2x[j] = tpd7;
  } else {
    /* Could go to ninth degree, but that wouldn't be legit. */
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "Invalid number of parameters (%d) for %s", ndparm, id);
  }

  /* No specialist de-distortions. */
  dis->disx2p[j] = 0x0;

  /* Record indexing parameters. */
  niparm = 5;
  if ((dis->iparm[j] = calloc(niparm, sizeof(int))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  /* The first three are required by disprt(). */
  dis->iparm[j][0] = TPD;
  dis->iparm[j][1] = niparm;
  dis->iparm[j][2] = ndparm;
  dis->iparm[j][3] = ndparm;

  /* Flag for presence of radial terms. */
  dis->iparm[j][4] = doradial;


  /* Allocate memory for the polynomial coefficients and fill it. */
  if ((dis->dparm[j] = calloc(ndparm, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  keyp = dis->dp;
  for (idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    fp = strpbrk(keyp->field, ".") + 1;

    /* One-to-one correspondence between TPV and TPD coefficients. */
    if (strncmp(fp, "TPV.", 4) == 0) {
      sscanf(fp+4, "%d", &k);
      dis->dparm[j][k] = wcsutil_dpkey_double(keyp);
    }
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

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

  char   *fp, id[16];
  int    deg, degree[2], idp, m, ncoeff[2], ndparm, niparm, p, q;
  struct dpkey *keyp;
  struct wcserr **err;
  int (*(distpd[2]))(DISP2X_ARGS);

  if (dis == 0x0) return DISERR_NULL_POINTER;
  err = &(dis->err);

  sprintf(id, "SIP on axis %d", j+1);


  /* Simple Imaging Polynomial. */
  if (dis->Nhat[j] != 2) {
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "Axis map for %s must contain 2 entries, not %d", id, dis->Nhat[j]);
  }

  /* Find the polynomial degree, at least 1 for the forward function. */
  degree[0] =  1;
  degree[1] = -1;
  keyp = dis->dp;
  for (idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    fp = strpbrk(keyp->field, ".") + 1;

    if (strncmp(fp, "SIP.", 4) == 0) {
      fp += 4;
      if (strncmp(fp, "FWD.", 4) == 0) {
        m = 0;

      } else if (strncmp(fp, "REV.", 4) == 0) {
        /* SIP uses a polynomial approximation for the inverse. */
        m = 1;

      } else {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Unrecognized field name for %s: %s", id, keyp->field);
      }

      fp += 4;
      sscanf(fp, "%d_%d", &p, &q);
      deg = p + q;
      if (p < 0 || 9 < p || q < 0 || 9 < q || 9 < deg) {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid powers (%d, %d) for %s: %s", p, q, id, keyp->field);
      }

      if (degree[m] < deg) degree[m] = deg;

    } else if (strcmp(fp, "NAXES")  &&
              strncmp(fp, "AXIS.",   5) &&
              strncmp(fp, "OFFSET.", 7) &&
              strncmp(fp, "SCALE.",  6)) {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Unrecognized field name for %s: %s", id, keyp->field);
    }
  }

  if (degree[1] == 0 ) degree[1] = 1;

  /* TPD is going to do the dirty work. */
  distpd[0] = 0x0;
  distpd[1] = 0x0;
  for (m = 0; m < 2; m++) {
    ncoeff[m] = 0;
    if (degree[m] == 1) {
      ncoeff[m] = 4;
      distpd[m] = tpd1;
    } else if (degree[m] == 2) {
      ncoeff[m] = 7;
      distpd[m] = tpd2;
    } else if (degree[m] == 3) {
      ncoeff[m] = 12;
      distpd[m] = tpd3;
    } else if (degree[m] == 4) {
      ncoeff[m] = 17;
      distpd[m] = tpd4;
    } else if (degree[m] == 5) {
      ncoeff[m] = 24;
      distpd[m] = tpd5;
    } else if (degree[m] == 6) {
      ncoeff[m] = 31;
      distpd[m] = tpd6;
    } else if (degree[m] == 7) {
      ncoeff[m] = 40;
      distpd[m] = tpd7;
    } else if (degree[m] == 8) {
      ncoeff[m] = 49;
      distpd[m] = tpd8;
    } else if (degree[m] == 9) {
      ncoeff[m] = 60;
      distpd[m] = tpd9;
    }
  }

  /* SIP uses a polynomial approximation to the inverse.  It's not very    */
  /* accurate but may provide disx2p() with a better zeroth approximation. */
  dis->disp2x[j] = distpd[0];
  dis->disx2p[j] = distpd[1];


  /* Record indexing parameters. */
  niparm = 5;
  if ((dis->iparm[j] = calloc(niparm, sizeof(int))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  /* The first three are required by disprt(). */
  ndparm = ncoeff[0] + ncoeff[1];
  dis->iparm[j][0] = TPD;
  dis->iparm[j][1] = niparm;
  dis->iparm[j][2] = ndparm;
  dis->iparm[j][3] = ncoeff[0];

  /* SIP never needs the radial terms. */
  dis->iparm[j][4] = 0;


  /* Allocate memory for the polynomial coefficients and fill it. */
  if ((dis->dparm[j] = calloc(ndparm, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  keyp = dis->dp;
  for (idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    fp = strpbrk(keyp->field, ".") + 1;

    if (strncmp(fp, "SIP.", 4) == 0) {
      fp += 4;
      if (strncmp(fp, "FWD.", 4) == 0) {
        m = 0;
      } else {
        m = ncoeff[0];
      }

      sscanf(fp+4, "%d_%d", &p, &q);

      /* Map to TPD coefficient number. */
      m += map[p][q];

      dis->dparm[j][m] = wcsutil_dpkey_double(keyp);
    }
  }

  /* Account for the fact that the SIP distortion provides an additive    */
  /* correction to the offset of the pixel coordinate from CRPIX, whereas */
  /* we expect the distortion function to provide the actual value of the */
  /* distorted pixel coordinate.  Note also that SIP has no axis mapping, */
  /* i.e. jhat = j always.                                                */
  m = j + 1;
  dis->dparm[j][0]  = dis->offset[j][j];
  dis->dparm[j][m] += 1.0;
  if (degree[1] > 0) {
    dis->dparm[j][ncoeff[0]]    = dis->offset[j][j];
    dis->dparm[j][ncoeff[0]+m] += 1.0;
  }


  return 0;
}

/*--------------------------------------------------------------------------*/

int tpdset(int j, struct disprm *dis)

{
  static const char *function = "tpdset";

  char   *fp, id[16];
  int    doradial, idp, k, m, ncoeff[2], ndparm, niparm;
  struct dpkey *keyp;
  struct wcserr **err;
  int (*(distpd[2]))(DISP2X_ARGS);

  if (dis == 0x0) return DISERR_NULL_POINTER;
  err = &(dis->err);

  sprintf(id, "TPD on axis %d", j+1);


  /* TPD distortion. */
  if (dis->Nhat[j] < 1 || 2 < dis->Nhat[j]) {
    return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
      "Axis map for %s must contain 1 or 2 entries, not %d", id,
      dis->Nhat[j]);
  }

  /* Find the number of parameters. */
  ncoeff[0] = 0;
  ncoeff[1] = 0;
  doradial  = 0;
  keyp = dis->dp;
  for (idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    fp = strpbrk(keyp->field, ".") + 1;

    if (strncmp(fp, "TPD.", 4) == 0) {
      fp += 4;
      if (strncmp(fp, "FWD.", 4) == 0) {
        m = 0;

      } else if (strncmp(fp, "REV.", 4) == 0) {
        /* TPD may provide a polynomial approximation for the inverse. */
        m = 1;

      } else {
        return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
          "Unrecognized field name for %s: %s", id, keyp->field);
      }

      sscanf(fp+4, "%d", &k);
      if (0 <= k && k <= 59) {
        if (ncoeff[m] < k+1) ncoeff[m] = k+1;

        /* Any radial terms? */
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

  distpd[0] = 0x0;
  distpd[1] = 0x0;
  for (m = 0; m < 2; m++) {
    if (ncoeff[m] <= 4) {
      /* First degree. */
      ncoeff[m] = 4;
      distpd[m] = tpd1;
    } else if (ncoeff[m] <= 7) {
      /* Second degree. */
      ncoeff[m] = 7;
      distpd[m] = tpd2;
    } else if (ncoeff[m] <= 12) {
      /* Third degree. */
      ncoeff[m] = 12;
      distpd[m] = tpd3;
    } else if (ncoeff[m] <= 17) {
      /* Fourth degree. */
      ncoeff[m] = 17;
      distpd[m] = tpd4;
    } else if (ncoeff[m] <= 24) {
      /* Fifth degree. */
      ncoeff[m] = 24;
      distpd[m] = tpd5;
    } else if (ncoeff[m] <= 31) {
      /* Sixth degree. */
      ncoeff[m] = 31;
      distpd[m] = tpd6;
    } else if (ncoeff[m] <= 40) {
      /* Seventh degree. */
      ncoeff[m] = 40;
      distpd[m] = tpd7;
    } else if (ncoeff[m] <= 49) {
      /* Eighth degree. */
      ncoeff[m] = 49;
      distpd[m] = tpd8;
    } else if (ncoeff[m] <= 60) {
      /* Ninth degree. */
      ncoeff[m] = 60;
      distpd[m] = tpd9;
    } else {
      return wcserr_set(WCSERR_SET(DISERR_BAD_PARAM),
        "Invalid number of parameters (%d) for %s", ncoeff[m], id);
    }
  }

  /* disx2p() only uses the inverse TPD, if present, to provide a better */
  /* zeroth approximation. */
  dis->disp2x[j] = distpd[0];
  dis->disx2p[j] = distpd[1];


  /* Record indexing parameters. */
  niparm = 5;
  if ((dis->iparm[j] = calloc(niparm, sizeof(int))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  /* The first three are required by disprt(). */
  ndparm = ncoeff[0] + ncoeff[1];
  dis->iparm[j][0] = TPD;
  dis->iparm[j][1] = niparm;
  dis->iparm[j][2] = ndparm;
  dis->iparm[j][3] = ncoeff[0];

  /* Flag for presence of radial terms. */
  dis->iparm[j][4] = doradial;


  /* Allocate memory for the polynomial coefficients and fill it. */
  if ((dis->dparm[j] = calloc(ndparm, sizeof(double))) == 0x0) {
    return wcserr_set(DIS_ERRMSG(DISERR_MEMORY));
  }

  keyp = dis->dp;
  for (idp = 0; idp < dis->ndp; idp++, keyp++) {
    if (keyp->j-1 != j) continue;

    fp = strpbrk(keyp->field, ".") + 1;

    if (strncmp(fp, "TPD.", 4) == 0) {
      fp += 4;
      if (strncmp(fp, "FWD.", 4) == 0) {
        m = 0;
      } else {
        m = ncoeff[0];
      }

      sscanf(fp+4, "%d", &k);
      dis->dparm[j][m+k] = wcsutil_dpkey_double(keyp);
    }
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int tpd1(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  double r, s, u, v;

  if (i[3] != 4 || 2 < Nhat) {
    return 1;
  }

  if (inverse) p += i[3];

  u = rawcrd[0];

  /* First degree. */
  *discrd = p[0] + u*p[1];

  if (Nhat == 1) return 0;

  v = rawcrd[1];

  *discrd += v*p[2];

  if (i[4]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*p[3];
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int tpd2(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  double r, s, u, v;

  if (i[3] != 7 || 2 < Nhat) {
    return 1;
  }

  if (inverse) p += i[3];

  u = rawcrd[0];

  /* Second degree. */
  *discrd = p[0] + u*(p[1] + u*(p[4]));

  if (Nhat == 1) return 0;

  v = rawcrd[1];

  *discrd +=
      v*(p[2]  + v*(p[6]))
    + u*(p[5])*v;

  if (i[4]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*p[3];
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int tpd3(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  double r, s, u, v;

  if (i[3] != 12 || 2 < Nhat) {
    return 1;
  }

  if (inverse) p += i[3];

  u = rawcrd[0];

  /* Third degree. */
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7])));

  if (Nhat == 1) return 0;

  v = rawcrd[1];

  *discrd +=
      v*(p[2]  + v*(p[6]  + v*(p[10])))
    + u*(p[5]  + v*(p[9])
    + u*(p[8]))*v;

  if (i[4]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11]));
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int tpd4(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  double r, s, u, v;

  if (i[3] != 17 || 2 < Nhat) {
    return 1;
  }

  if (inverse) p += i[3];

  u = rawcrd[0];

  /* Fourth degree. */
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7] + u*(p[12]))));

  if (Nhat == 1) return 0;

  v = rawcrd[1];

  *discrd +=
      v*(p[2]  + v*(p[6]  + v*(p[10] + v*(p[16]))))
    + u*(p[5]  + v*(p[9]  + v*(p[15]))
    + u*(p[8]  + v*(p[14])
    + u*(p[13])))*v;

  if (i[4]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11]));
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int tpd5(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  double r, s, u, v;

  if (i[3] != 24 || 2 < Nhat) {
    return 1;
  }

  if (inverse) p += i[3];

  u = rawcrd[0];

  /* Fifth degree. */
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7] + u*(p[12] + u*(p[17])))));

  if (Nhat == 1) return 0;

  v = rawcrd[1];

  *discrd +=
      v*(p[2]  + v*(p[6]  + v*(p[10] + v*(p[16] + v*(p[22])))))
    + u*(p[5]  + v*(p[9]  + v*(p[15] + v*(p[21])))
    + u*(p[8]  + v*(p[14] + v*(p[20]))
    + u*(p[13] + v*(p[19])
    + u*(p[18]))))*v;

  if (i[4]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11] + s*(p[23])));
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int tpd6(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  double r, s, u, v;

  if (i[3] != 31 || 2 < Nhat) {
    return 1;
  }

  if (inverse) p += i[3];

  u = rawcrd[0];

  /* Sixth degree. */
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7] + u*(p[12] + u*(p[17] + u*(p[24]))))));

  if (Nhat == 1) return 0;

  v = rawcrd[1];

  *discrd +=
      v*(p[2]  + v*(p[6]  + v*(p[10] + v*(p[16] + v*(p[22] + v*(p[30]))))))
    + u*(p[5]  + v*(p[9]  + v*(p[15] + v*(p[21] + v*(p[29]))))
    + u*(p[8]  + v*(p[14] + v*(p[20] + v*(p[28])))
    + u*(p[13] + v*(p[19] + v*(p[27]))
    + u*(p[18] + v*(p[26])
    + u*(p[25])))))*v;

  if (i[4]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11] + s*(p[23])));
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int tpd7(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  double r, s, u, v;

  if (i[3] != 40 || 2 < Nhat) {
    return 1;
  }

  if (inverse) p += i[3];

  u = rawcrd[0];

  /* Seventh degree. */
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7] + u*(p[12] + u*(p[17] + u*(p[24] + u*(p[31])))))));

  if (Nhat == 1) return 0;

  v = rawcrd[1];

  *discrd +=
      v*(p[2]  + v*(p[6]  + v*(p[10] + v*(p[16] + v*(p[22] + v*(p[30] + v*(p[38])))))))
    + u*(p[5]  + v*(p[9]  + v*(p[15] + v*(p[21] + v*(p[29] + v*(p[37])))))
    + u*(p[8]  + v*(p[14] + v*(p[20] + v*(p[28] + v*(p[36]))))
    + u*(p[13] + v*(p[19] + v*(p[27] + v*(p[35])))
    + u*(p[18] + v*(p[26] + v*(p[34]))
    + u*(p[25] + v*(p[33])
    + u*(p[32]))))))*v;

  if (i[4]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11] + s*(p[23] + s*(p[39]))));
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int tpd8(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  double r, s, u, v;

  if (i[3] != 49 || 2 < Nhat) {
    return 1;
  }

  if (inverse) p += i[3];

  u = rawcrd[0];

  /* Eighth degree. */
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7] + u*(p[12] + u*(p[17] + u*(p[24] + u*(p[31] + u*(p[40]))))))));

  if (Nhat == 1) return 0;

  v = rawcrd[1];

  *discrd +=
      v*(p[2]  + v*(p[6]  + v*(p[10] + v*(p[16] + v*(p[22] + v*(p[30] + v*(p[38] + v*(p[48]))))))))
    + u*(p[5]  + v*(p[9]  + v*(p[15] + v*(p[21] + v*(p[29] + v*(p[37] + v*(p[47]))))))
    + u*(p[8]  + v*(p[14] + v*(p[20] + v*(p[28] + v*(p[36] + v*(p[46])))))
    + u*(p[13] + v*(p[19] + v*(p[27] + v*(p[35] + v*(p[45]))))
    + u*(p[18] + v*(p[26] + v*(p[34] + v*(p[44])))
    + u*(p[25] + v*(p[33] + v*(p[43]))
    + u*(p[32] + v*(p[42])
    + u*(p[41])))))))*v;

  if (i[4]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11] + s*(p[23] + s*(p[39]))));
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int tpd9(
  int inverse,
  const int i[],
  const double p[],
  int Nhat,
  const double rawcrd[],
  double *discrd)

{
  double r, s, u, v;

  if (i[3] != 60 || 2 < Nhat) {
    return 1;
  }

  if (inverse) p += i[3];

  u = rawcrd[0];

  /* Ninth degree. */
  *discrd = p[0] + u*(p[1] + u*(p[4] + u*(p[7] + u*(p[12] + u*(p[17] + u*(p[24] + u*(p[31] + u*(p[40] + u*(p[49])))))))));

  if (Nhat == 1) return 0;

  v = rawcrd[1];

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

  if (i[4]) {
    s = u*u + v*v;
    r = sqrt(s);

    *discrd += r*(p[3] + s*(p[11] + s*(p[23] + s*(p[39] + s*(p[59])))));
  }

  return 0;
}
