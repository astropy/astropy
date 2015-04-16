/*============================================================================

  WCSLIB 5.2 - an implementation of the FITS WCS standard.
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
  $Id: tab.c,v 5.2 2015/04/15 12:35:07 mcalabre Exp $
*===========================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcserr.h"
#include "wcsmath.h"
#include "wcsprintf.h"
#include "wcsutil.h"
#include "tab.h"

const int TABSET = 137;

/* Map status return value to message. */
const char *tab_errmsg[] = {
  "Success",
  "Null tabprm pointer passed",
  "Memory allocation failed",
  "Invalid tabular parameters",
  "One or more of the x coordinates were invalid",
  "One or more of the world coordinates were invalid"};

/* Convenience macro for invoking wcserr_set(). */
#define TAB_ERRMSG(status) WCSERR_SET(status), tab_errmsg[status]

/*--------------------------------------------------------------------------*/

int tabini(int alloc, int M, const int K[], struct tabprm *tab)

{
  static const char *function = "tabini";

  int k, m, N;
  double *dp;
  struct wcserr **err;

  if (tab == 0x0) return TABERR_NULL_POINTER;

  /* Initialize error message handling. */
  err = &(tab->err);
  if (tab->flag != -1) {
    if (tab->err) free(tab->err);
  }
  tab->err = 0x0;


  if (M <= 0) {
    return wcserr_set(WCSERR_SET(TABERR_BAD_PARAMS),
      "M must be positive, got %d", M);
  }

  /* Determine the total number of elements in the coordinate array. */
  if (K) {
    N = M;

    for (m = 0; m < M; m++) {
      if (K[m] < 0) {
        return wcserr_set(WCSERR_SET(TABERR_BAD_PARAMS),
          "Invalid tabular parameters: Each element of K must be "
          "non-negative, got %d", K[m]);
      }

      N *= K[m];
    }

  } else {
    /* Axis lengths as yet unknown. */
    N = 0;
  }


  /* Initialize memory management. */
  if (tab->flag == -1 || tab->m_flag != TABSET) {
    if (tab->flag == -1) {
      tab->sense   = 0x0;
      tab->p0      = 0x0;
      tab->delta   = 0x0;
      tab->extrema = 0x0;
      tab->set_M   = 0;
    }

    tab->m_flag  = 0;
    tab->m_M     = 0;
    tab->m_N     = 0;
    tab->m_K     = 0x0;
    tab->m_map   = 0x0;
    tab->m_crval = 0x0;
    tab->m_index = 0x0;
    tab->m_indxs = 0x0;
    tab->m_coord = 0x0;

  } else {
    /* Clear any outstanding signals set by wcstab(). */
    for (m = 0; m < tab->m_M; m++) {
      if (tab->m_indxs[m] == (double *)0x1) tab->m_indxs[m] = 0x0;
    }

    if (tab->m_coord == (double *)0x1) tab->m_coord = 0x0;
  }


  /* Allocate memory for arrays if required. */
  if (alloc ||
     tab->K == 0x0 ||
     tab->map == 0x0 ||
     tab->crval == 0x0 ||
     tab->index == 0x0 ||
     tab->coord == 0x0) {

    /* Was sufficient allocated previously? */
    if (tab->m_flag == TABSET && (tab->m_M < M || tab->m_N < N)) {
      /* No, free it. */
      tabfree(tab);
    }

    if (alloc || tab->K == 0x0) {
      if (tab->m_K) {
        /* In case the caller fiddled with it. */
        tab->K = tab->m_K;

      } else {
        if (!(tab->K = calloc(M, sizeof(int)))) {
          return wcserr_set(TAB_ERRMSG(TABERR_MEMORY));
        }

        tab->m_flag = TABSET;
        tab->m_M = M;
        tab->m_K = tab->K;
      }
    }

    if (alloc || tab->map == 0x0) {
      if (tab->m_map) {
        /* In case the caller fiddled with it. */
        tab->map = tab->m_map;

      } else {
        if (!(tab->map = calloc(M, sizeof(int)))) {
          return wcserr_set(TAB_ERRMSG(TABERR_MEMORY));
        }

        tab->m_flag = TABSET;
        tab->m_M = M;
        tab->m_map = tab->map;
      }
    }

    if (alloc || tab->crval == 0x0) {
      if (tab->m_crval) {
        /* In case the caller fiddled with it. */
        tab->crval = tab->m_crval;

      } else {
        if (!(tab->crval = calloc(M, sizeof(double)))) {
          return wcserr_set(TAB_ERRMSG(TABERR_MEMORY));
        }

        tab->m_flag = TABSET;
        tab->m_M = M;
        tab->m_crval = tab->crval;
      }
    }

    if (alloc || tab->index == 0x0) {
      if (tab->m_index) {
        /* In case the caller fiddled with it. */
        tab->index = tab->m_index;

      } else {
        if (!(tab->index = calloc(M, sizeof(double *)))) {
          return wcserr_set(TAB_ERRMSG(TABERR_MEMORY));
        }

        tab->m_flag = TABSET;
        tab->m_M = M;
        tab->m_N = N;
        tab->m_index = tab->index;

        if (!(tab->m_indxs = calloc(M, sizeof(double *)))) {
          return wcserr_set(TAB_ERRMSG(TABERR_MEMORY));
        }

        /* Recall that calloc() initializes these pointers to zero. */
        if (K) {
          for (m = 0; m < M; m++) {
            if (K[m]) {
              if (!(tab->index[m] = calloc(K[m], sizeof(double)))) {
                return wcserr_set(TAB_ERRMSG(TABERR_MEMORY));
              }

              tab->m_indxs[m] = tab->index[m];
            }
          }
        }
      }
    }

    if (alloc || tab->coord == 0x0) {
      if (tab->m_coord) {
        /* In case the caller fiddled with it. */
        tab->coord = tab->m_coord;

      } else if (N) {
        if (!(tab->coord = calloc(N, sizeof(double)))) {
          return wcserr_set(TAB_ERRMSG(TABERR_MEMORY));
        }

        tab->m_flag = TABSET;
        tab->m_M = M;
        tab->m_N = N;
        tab->m_coord = tab->coord;
      }
    }
  }

  tab->flag = 0;
  tab->M = M;

  /* Set defaults. */
  for (m = 0; m < M; m++) {
    tab->map[m] = -1;
    tab->crval[m] = 0.0;

    if (K) {
      tab->K[m] = K[m];
      if ((dp = tab->index[m])) {
        for (k = 0; k < K[m]; k++) {
          *(dp++) = k;
        }
      }
    } else {
      tab->K[m] = 0;
    }
  }

  /* Initialize the coordinate array. */
  for (dp = tab->coord; dp < tab->coord + N; dp++) {
    *dp = UNDEFINED;
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int tabmem(struct tabprm *tab)

{
  static const char *function = "tabmem";

  int m, M, N;
  struct wcserr **err;

  if (tab == 0x0) return TABERR_NULL_POINTER;
  err = &(tab->err);

  if (tab->M == 0 || tab->K == 0x0) {
    /* Should have been set by this time. */
    return wcserr_set(WCSERR_SET(TABERR_MEMORY),
      "Null pointers in tabprm struct");
  }


  N = M = tab->M;
  for (m = 0; m < M; m++) {
    if (tab->K[m] < 0) {
      return wcserr_set(WCSERR_SET(TABERR_BAD_PARAMS),
        "Invalid tabular parameters: Each element of K must be "
        "non-negative, got %d", M);
    }

    N *= tab->K[m];
  }


  if (tab->m_M == 0) {
    tab->m_M = M;
  } else if (tab->m_M < M) {
    /* Only possible if the user changed M. */
    return wcserr_set(WCSERR_SET(TABERR_MEMORY),
      "tabprm struct inconsistent");
  }

  if (tab->m_N == 0) {
    tab->m_N = N;
  } else if (tab->m_N < N) {
    /* Only possible if the user changed K[]. */
    return wcserr_set(WCSERR_SET(TABERR_MEMORY),
      "tabprm struct inconsistent");
  }

  if (tab->m_K == 0x0) {
    if ((tab->m_K = tab->K)) {
      tab->m_flag = TABSET;
    }
  }

  if (tab->m_map == 0x0) {
    if ((tab->m_map = tab->map)) {
      tab->m_flag = TABSET;
    }
  }

  if (tab->m_crval == 0x0) {
    if ((tab->m_crval = tab->crval)) {
      tab->m_flag = TABSET;
    }
  }

  if (tab->m_index == 0x0) {
    if ((tab->m_index = tab->index)) {
      tab->m_flag = TABSET;
    }
  }

  for (m = 0; m < tab->m_M; m++) {
    if (tab->m_indxs[m] == 0x0 || tab->m_indxs[m] == (double *)0x1) {
      if ((tab->m_indxs[m] = tab->index[m])) {
        tab->m_flag = TABSET;
      }
    }
  }

  if (tab->m_coord == 0x0 || tab->m_coord == (double *)0x1) {
    if ((tab->m_coord = tab->coord)) {
      tab->m_flag = TABSET;
    }
  }

  tab->flag = 0;

  return 0;
}

/*--------------------------------------------------------------------------*/

int tabcpy(int alloc, const struct tabprm *tabsrc, struct tabprm *tabdst)

{
  static const char *function = "tabcpy";

  int k, m, M, n, N, status;
  double *dstp, *srcp;
  struct wcserr **err;

  if (tabsrc == 0x0) return TABERR_NULL_POINTER;
  if (tabdst == 0x0) return TABERR_NULL_POINTER;
  err = &(tabdst->err);

  M = tabsrc->M;
  if (M <= 0) {
    return wcserr_set(WCSERR_SET(TABERR_BAD_PARAMS),
      "M must be positive, got %d", M);
  }

  if ((status = tabini(alloc, M, tabsrc->K, tabdst))) {
    return status;
  }

  N = M;
  for (m = 0; m < M; m++) {
    tabdst->map[m]   = tabsrc->map[m];
    tabdst->crval[m] = tabsrc->crval[m];
    N *= tabsrc->K[m];
  }

  for (m = 0; m < M; m++) {
    if ((srcp = tabsrc->index[m])) {
      dstp = tabdst->index[m];
      for (k = 0; k < tabsrc->K[m]; k++) {
        *(dstp++) = *(srcp++);
      }
    }
  }

  srcp = tabsrc->coord;
  dstp = tabdst->coord;
  for (n = 0; n < N; n++) {
    *(dstp++) = *(srcp++);
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int tabcmp(int cmp,
           double tol,
           const struct tabprm *tab1,
           const struct tabprm *tab2,
           int *equal)

{
  int m, M, N;

  if (tab1  == 0x0) return TABERR_NULL_POINTER;
  if (tab2  == 0x0) return TABERR_NULL_POINTER;
  if (equal == 0x0) return TABERR_NULL_POINTER;

  *equal = 0;

  if (tab1->M != tab2->M) {
    return 0;
  }

  M = tab1->M;

  if (!wcsutil_intEq(M, tab1->K, tab2->K) ||
      !wcsutil_intEq(M, tab1->map, tab2->map) ||
      !wcsutil_Eq(M, tol, tab1->crval, tab2->crval)) {
    return 0;
  }

  N = M;
  for (m = 0; m < M; m++) {
    if (!wcsutil_Eq(tab1->K[m], tol, tab1->index[m], tab2->index[m])) {
      return 0;
    }

    N *= tab1->K[m];
  }

  if (!wcsutil_Eq(N, tol, tab1->coord, tab2->coord)) {
    return 0;
  }

  *equal = 1;

  return 0;
}


/*--------------------------------------------------------------------------*/

int tabfree(struct tabprm *tab)

{
  int m;

  if (tab == 0x0) return TABERR_NULL_POINTER;

  if (tab->flag != -1) {
    /* Clear any outstanding signals set by wcstab(). */
    for (m = 0; m < tab->m_M; m++) {
      if (tab->m_indxs[m] == (double *)0x1) tab->m_indxs[m] = 0x0;
    }

    if (tab->m_coord == (double *)0x1) tab->m_coord = 0x0;

    /* Free memory allocated by tabini(). */
    if (tab->m_flag == TABSET) {
      if (tab->K     == tab->m_K)     tab->K = 0x0;
      if (tab->map   == tab->m_map)   tab->map = 0x0;
      if (tab->crval == tab->m_crval) tab->crval = 0x0;
      if (tab->index == tab->m_index) tab->index = 0x0;
      if (tab->coord == tab->m_coord) tab->coord = 0x0;

      if (tab->m_K)     free(tab->m_K);
      if (tab->m_map)   free(tab->m_map);
      if (tab->m_crval) free(tab->m_crval);

      if (tab->m_index) {
        for (m = 0; m < tab->m_M; m++) {
          if (tab->m_indxs[m]) free(tab->m_indxs[m]);
        }
        free(tab->m_index);
        free(tab->m_indxs);
      }

      if (tab->m_coord) free(tab->m_coord);
    }

    /* Free memory allocated by tabset(). */
    if (tab->sense)   free(tab->sense);
    if (tab->p0)      free(tab->p0);
    if (tab->delta)   free(tab->delta);
    if (tab->extrema) free(tab->extrema);
  }

  tab->m_flag  = 0;
  tab->m_M     = 0;
  tab->m_N     = 0;
  tab->m_K     = 0x0;
  tab->m_map   = 0x0;
  tab->m_crval = 0x0;
  tab->m_index = 0x0;
  tab->m_indxs = 0x0;
  tab->m_coord = 0x0;

  tab->sense   = 0x0;
  tab->p0      = 0x0;
  tab->delta   = 0x0;
  tab->extrema = 0x0;
  tab->set_M   = 0;

  if (tab->err) {
    free(tab->err);
    tab->err = 0x0;
  }

  tab->flag = 0;

  return 0;
}

/*--------------------------------------------------------------------------*/

int tabprt(const struct tabprm *tab)

{
  char   *cp, text[128];
  int    j, k, m, n, nd;
  double *dp;

  if (tab == 0x0) return TABERR_NULL_POINTER;

  if (tab->flag != TABSET) {
    wcsprintf("The tabprm struct is UNINITIALIZED.\n");
    return 0;
  }

  wcsprintf("       flag: %d\n", tab->flag);
  wcsprintf("          M: %d\n", tab->M);

  /* Array dimensions. */
  WCSPRINTF_PTR("          K: ", tab->K, "\n");
  wcsprintf("            ");
  for (m = 0; m < tab->M; m++) {
    wcsprintf("%6d", tab->K[m]);
  }
  wcsprintf("\n");

  /* Map vector. */
  WCSPRINTF_PTR("        map: ", tab->map, "\n");
  wcsprintf("            ");
  for (m = 0; m < tab->M; m++) {
    wcsprintf("%6d", tab->map[m]);
  }
  wcsprintf("\n");

  /* Reference index value. */
  WCSPRINTF_PTR("      crval: ", tab->crval, "\n");
  wcsprintf("            ");
  for (m = 0; m < tab->M; m++) {
    wcsprintf("  %#- 11.5g", tab->crval[m]);
  }
  wcsprintf("\n");

  /* Index vectors. */
  WCSPRINTF_PTR("      index: ", tab->index, "\n");
  for (m = 0; m < tab->M; m++) {
    wcsprintf("   index[%d]: ", m);
    WCSPRINTF_PTR("", tab->index[m], "");
    if (tab->index[m]) {
      for (k = 0; k < tab->K[m]; k++) {
        if (k%5 == 0) {
          wcsprintf("\n            ");
        }
        wcsprintf("  %#- 11.5g", tab->index[m][k]);
      }
      wcsprintf("\n");
    }
  }

  /* Coordinate array. */
  WCSPRINTF_PTR("      coord: ", tab->coord, "\n");
  dp = tab->coord;
  for (n = 0; n < tab->nc; n++) {
    /* Array index. */
    j = n;
    cp = text;
    for (m = 0; m < tab->M; m++) {
      nd = (tab->K[m] < 10) ? 1 : 2;
      sprintf(cp, ",%*d", nd, j % tab->K[m] + 1);
      j /= tab->K[m];
      cp += strlen(cp);
    }

    wcsprintf("             (*%s)", text);
    for (m = 0; m < tab->M; m++) {
      wcsprintf("  %#- 11.5g", *(dp++));
    }
    wcsprintf("\n");
  }

  wcsprintf("         nc: %d\n", tab->nc);

  WCSPRINTF_PTR("      sense: ", tab->sense, "\n");
  if (tab->sense) {
    wcsprintf("            ");
    for (m = 0; m < tab->M; m++) {
      wcsprintf("%6d", tab->sense[m]);
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("         p0: ", tab->p0, "\n");
  if (tab->p0) {
    wcsprintf("            ");
    for (m = 0; m < tab->M; m++) {
      wcsprintf("%6d", tab->p0[m]);
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("      delta: ", tab->delta, "\n");
  if (tab->delta) {
    wcsprintf("            ");
    for (m = 0; m < tab->M; m++) {
      wcsprintf("  %#- 11.5g", tab->delta[m]);
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("    extrema: ", tab->extrema, "\n");
  dp = tab->extrema;
  for (n = 0; n < tab->nc/tab->K[0]; n++) {
    /* Array index. */
    j = n;
    cp = text;
    *cp = '\0';
    for (m = 1; m < tab->M; m++) {
      nd = (tab->K[m] < 10) ? 1 : 2;
      sprintf(cp, ",%*d", nd, j % tab->K[m] + 1);
      j /= tab->K[m];
      cp += strlen(cp);
    }

    wcsprintf("             (*,*%s)", text);
    for (m = 0; m < 2*tab->M; m++) {
      if (m == tab->M) wcsprintf("->  ");
      wcsprintf("  %#- 11.5g", *(dp++));
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("        err: ", tab->err, "\n");
  if (tab->err) {
    wcserr_prt(tab->err, "             ");
  }

  /* Memory management. */
  wcsprintf("     m_flag: %d\n", tab->m_flag);
  wcsprintf("        m_M: %d\n", tab->m_M);
  wcsprintf("        m_N: %d\n", tab->m_N);

  WCSPRINTF_PTR("        m_K: ", tab->m_K, "");
  if (tab->m_K == tab->K) wcsprintf("  (= K)");
  wcsprintf("\n");

  WCSPRINTF_PTR("      m_map: ", tab->m_map, "");
  if (tab->m_map == tab->map) wcsprintf("  (= map)");
  wcsprintf("\n");

  WCSPRINTF_PTR("    m_crval: ", tab->m_crval, "");
  if (tab->m_crval == tab->crval) wcsprintf("  (= crval)");
  wcsprintf("\n");

  WCSPRINTF_PTR("    m_index: ", tab->m_index, "");
  if (tab->m_index == tab->index) wcsprintf("  (= index)");
  wcsprintf("\n");
  for (m = 0; m < tab->M; m++) {
    wcsprintf(" m_indxs[%d]: ", m);
    WCSPRINTF_PTR("", tab->m_indxs[m], "");
    if (tab->m_indxs[m] == tab->index[m]) wcsprintf("  (= index[%d])", m);
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("    m_coord: ", tab->m_coord, "");
  if (tab->m_coord == tab->coord) wcsprintf("  (= coord)");
  wcsprintf("\n");

  return 0;
}

/*--------------------------------------------------------------------------*/

int tabset(struct tabprm *tab)

{
  static const char *function = "tabset";

  int i, ic, k, *Km, m, M, ne;
  double *dcrd, *dmax, *dmin, dPsi, dval, *Psi;
  struct wcserr **err;

  if (tab == 0x0) return TABERR_NULL_POINTER;
  err = &(tab->err);

  /* Check the number of tabular coordinate axes. */
  if ((M = tab->M) < 1) {
    return wcserr_set(WCSERR_SET(TABERR_BAD_PARAMS),
      "Invalid tabular parameters: M must be positive, got %d", M);
  }

  /* Check the axis lengths. */
  if (!tab->K) {
    return wcserr_set(WCSERR_SET(TABERR_MEMORY),
      "Null pointers in tabprm struct");
  }

  tab->nc = 1;
  for (m = 0; m < M; m++) {
    if (tab->K[m] < 1) {
      return wcserr_set(WCSERR_SET(TABERR_BAD_PARAMS),
        "Invalid tabular parameters: Each element of K must be positive, "
        "got %d", tab->K[m]);
    }

    /* Number of coordinate vectors in the coordinate array. */
    tab->nc *= tab->K[m];
  }

  /* Check that the map vector is sensible. */
  if (!tab->map) {
    return wcserr_set(WCSERR_SET(TABERR_MEMORY),
      "Null pointers in tabprm struct");
  }

  for (m = 0; m < M; m++) {
    i = tab->map[m];
    if (i < 0) {
      return wcserr_set(WCSERR_SET(TABERR_BAD_PARAMS),
        "Invalid tabular parameters: Each element of map must be "
        "non-negative, got %d", i);
    }
  }

  /* Check memory allocation for the remaining vectors. */
  if (!tab->crval || !tab->index || !tab->coord) {
    return wcserr_set(WCSERR_SET(TABERR_MEMORY),
      "Null pointers in tabprm struct");
  }

  /* Take memory if signalled to by wcstab(). */
  for (m = 0; m < tab->m_M; m++) {
    if (tab->m_indxs[m] == (double *)0x1 &&
      (tab->m_indxs[m] = tab->index[m])) {
      tab->m_flag = TABSET;
    }
  }

  if (tab->m_coord == (double *)0x1 &&
    (tab->m_coord = tab->coord)) {
    tab->m_flag = TABSET;
  }


  /* Allocate memory for work vectors. */
  if (tab->flag != TABSET || tab->set_M < M) {
    /* Free memory that may have been allocated previously. */
    if (tab->sense)   free(tab->sense);
    if (tab->p0)      free(tab->p0);
    if (tab->delta)   free(tab->delta);
    if (tab->extrema) free(tab->extrema);

    /* Allocate memory for internal arrays. */
    if (!(tab->sense = calloc(M, sizeof(int)))) {
      return wcserr_set(TAB_ERRMSG(TABERR_MEMORY));
    }

    if (!(tab->p0 = calloc(M, sizeof(int)))) {
      free(tab->sense);
      return wcserr_set(TAB_ERRMSG(TABERR_MEMORY));
    }

    if (!(tab->delta = calloc(M, sizeof(double)))) {
      free(tab->sense);
      free(tab->p0);
      return wcserr_set(TAB_ERRMSG(TABERR_MEMORY));
    }

    ne = M * tab->nc * 2 / tab->K[0];
    if (!(tab->extrema = calloc(ne, sizeof(double)))) {
      free(tab->sense);
      free(tab->p0);
      free(tab->delta);
      return wcserr_set(TAB_ERRMSG(TABERR_MEMORY));
    }

    tab->set_M = M;
  }

  /* Check that the index vectors are monotonic. */
  Km = tab->K;
  for (m = 0; m < M; m++, Km++) {
    tab->sense[m] = 0;

    if (*Km > 1) {
      if ((Psi = tab->index[m]) == 0x0) {
        /* Default indexing. */
        tab->sense[m] = 1;

      } else {
        for (k = 0; k < *Km-1; k++) {
          switch (tab->sense[m]) {
          case 0:
            if (Psi[k] < Psi[k+1]) {
              /* Monotonic increasing. */
              tab->sense[m] = 1;
            } else if (Psi[k] > Psi[k+1]) {
              /* Monotonic decreasing. */
              tab->sense[m] = -1;
            }
            break;

          case 1:
            if (Psi[k] > Psi[k+1]) {
              /* Should be monotonic increasing. */
              free(tab->sense);
              free(tab->p0);
              free(tab->delta);
              free(tab->extrema);
              return wcserr_set(WCSERR_SET(TABERR_BAD_PARAMS),
                "Invalid tabular parameters: Index vectors are not "
                "monotonically increasing");
            }
            break;

          case -1:
            if (Psi[k] < Psi[k+1]) {
              /* Should be monotonic decreasing. */
              free(tab->sense);
              free(tab->p0);
              free(tab->delta);
              free(tab->extrema);
              return wcserr_set(WCSERR_SET(TABERR_BAD_PARAMS),
                "Invalid tabular parameters: Index vectors are not "
                "monotonically decreasing");
            }
            break;
          }
        }
      }

      if (tab->sense[m] == 0) {
        free(tab->sense);
        free(tab->p0);
        free(tab->delta);
        free(tab->extrema);
        return wcserr_set(WCSERR_SET(TABERR_BAD_PARAMS),
          "Invalid tabular parameters: Index vectors are not monotonic");
      }
    }
  }

  /* Find the extremal values of the coordinate elements in each row. */
  dcrd = tab->coord;
  dmin = tab->extrema;
  dmax = tab->extrema + M;
  for (ic = 0; ic < tab->nc; ic += tab->K[0]) {
    for (m = 0; m < M; m++, dcrd++) {
      if (tab->K[0] > 1) {
        /* Extrapolate a little before the start of the row. */
        Psi = tab->index[0];
        if (Psi == 0x0) {
          dPsi = 1.0;
        } else {
          dPsi = Psi[1] - Psi[0];
        }

        dval = *dcrd;
        if (dPsi != 0.0) {
          dval -= 0.5 * (*(dcrd+M) - *dcrd)/dPsi;
        }

        *(dmax+m) = *(dmin+m) = dval;
      } else {
        *(dmax+m) = *(dmin+m) = *dcrd;
      }
    }

    dcrd -= M;
    for (i = 0; i < tab->K[0]; i++) {
      for (m = 0; m < M; m++, dcrd++) {
        if (*(dmax+m) < *dcrd) *(dmax+m) = *dcrd;
        if (*(dmin+m) > *dcrd) *(dmin+m) = *dcrd;

        if (tab->K[0] > 1 && i == tab->K[0]-1) {
          /* Extrapolate a little beyond the end of the row. */
          Psi = tab->index[0];
          if (Psi == 0x0) {
            dPsi = 1.0;
          } else {
            dPsi = Psi[i] - Psi[i-1];
          }

          dval = *dcrd;
          if (dPsi != 0.0) {
            dval += 0.5 * (*dcrd - *(dcrd-M))/dPsi;
          }

          if (*(dmax+m) < dval) *(dmax+m) = dval;
          if (*(dmin+m) > dval) *(dmin+m) = dval;
        }
      }
    }

    dmin += 2*M;
    dmax += 2*M;
  }

  tab->flag = TABSET;

  return 0;
}

/*--------------------------------------------------------------------------*/

int tabx2s(
  struct tabprm *tab,
  int ncoord,
  int nelem,
  const double x[],
  double world[],
  int stat[])

{
  static const char *function = "tabx2s";

  int i, iv, k, *Km, m, M, n, nv, offset, p1, status;
  double *coord, *Psi, psi_m, upsilon, wgt;
  register int *statp;
  register const double *xp;
  register double *wp;
  struct wcserr **err;

  if (tab == 0x0) return TABERR_NULL_POINTER;
  err = &(tab->err);

  /* Initialize if required. */
  if (tab->flag != TABSET) {
    if ((status = tabset(tab))) return status;
  }

  /* This is used a lot. */
  M = tab->M;

  status = 0;
  xp = x;
  wp = world;
  statp = stat;
  for (n = 0; n < ncoord; n++) {
    /* Determine the indexes. */
    Km = tab->K;
    for (m = 0; m < M; m++, Km++) {
      /* N.B. psi_m and Upsilon_m are 1-relative FITS indexes. */
      i = tab->map[m];
      psi_m = *(xp+i) + tab->crval[m];

      Psi = tab->index[m];
      if (Psi == 0x0) {
        /* Default indexing is simple. */
        upsilon = psi_m;

      } else {
        /* To ease confusion, decrement Psi so that we can use 1-relative
           C array indexing to match the 1-relative FITS indexing. */
        Psi--;

        if (*Km == 1) {
          /* Index vector is degenerate. */
          if (Psi[1]-0.5 <= psi_m && psi_m <= Psi[1]+0.5) {
            upsilon = psi_m;
          } else {
            *statp = 1;
            status = wcserr_set(TAB_ERRMSG(TABERR_BAD_X));
            goto next;
          }

        } else {
          /* Interpolate in the indexing vector. */
          if (tab->sense[m] == 1) {
            /* Monotonic increasing index values. */
            if (psi_m < Psi[1]) {
              if (Psi[1] - 0.5*(Psi[2]-Psi[1]) <= psi_m) {
                /* Allow minor extrapolation. */
                k = 1;

              } else {
                /* Index is out of range. */
                *statp = 1;
                status = wcserr_set(TAB_ERRMSG(TABERR_BAD_X));
                goto next;
              }

            } else if (Psi[*Km] < psi_m) {
              if (psi_m <= Psi[*Km] + 0.5*(Psi[*Km]-Psi[*Km-1])) {
                /* Allow minor extrapolation. */
                k = *Km - 1;

              } else {
                /* Index is out of range. */
                *statp = 1;
                status = wcserr_set(TAB_ERRMSG(TABERR_BAD_X));
                goto next;
              }

            } else {
              for (k = 1; k < *Km; k++) {
                if (psi_m < Psi[k]) {
                  continue;
                }
                if (Psi[k] == psi_m && psi_m < Psi[k+1]) {
                  break;
                }
                if (Psi[k] < psi_m && psi_m <= Psi[k+1]) {
                  break;
                }
              }
            }

          } else {
            /* Monotonic decreasing index values. */
            if (psi_m > Psi[1]) {
              if (Psi[1] + 0.5*(Psi[1]-Psi[2]) >= psi_m) {
                /* Allow minor extrapolation. */
                k = 1;

              } else {
                /* Index is out of range. */
                *statp = 1;
                status = wcserr_set(TAB_ERRMSG(TABERR_BAD_X));
                goto next;
              }

            } else if (psi_m < Psi[*Km]) {
              if (Psi[*Km] - 0.5*(Psi[*Km-1]-Psi[*Km]) <= psi_m) {
                /* Allow minor extrapolation. */
                k = *Km - 1;

              } else {
                /* Index is out of range. */
                *statp = 1;
                status = wcserr_set(TAB_ERRMSG(TABERR_BAD_X));
                goto next;
              }

            } else {
              for (k = 1; k < *Km; k++) {
                if (psi_m > Psi[k]) {
                  continue;
                }
                if (Psi[k] == psi_m && psi_m > Psi[k+1]) {
                  break;
                }
                if (Psi[k] > psi_m && psi_m >= Psi[k+1]) {
                  break;
                }
              }
            }
          }

          upsilon = k + (psi_m - Psi[k]) / (Psi[k+1] - Psi[k]);
        }
      }

      if (upsilon < 0.5 || upsilon > *Km + 0.5) {
        /* Index out of range. */
        *statp = 1;
        status = wcserr_set(TAB_ERRMSG(TABERR_BAD_X));
        goto next;
      }

      /* Fiducial array indices and fractional offset.
         p1 is 1-relative while tab::p0 is 0-relative. */
      p1 = (int)floor(upsilon);
      tab->p0[m] = p1 - 1;
      tab->delta[m] = upsilon - p1;

      if (p1 == 0) {
        tab->p0[m] += 1;
        tab->delta[m] -= 1.0;
      } else if (p1 == *Km && *Km > 1) {
        tab->p0[m] -= 1;
        tab->delta[m] += 1.0;
      }
    }


    /* Now interpolate in the coordinate array; the M-dimensional linear  */
    /* interpolation algorithm is described in Sect. 3.4 of WCS Paper IV. */
    for (m = 0; m < M; m++) {
     i = tab->map[m];
     *(wp+i) = 0.0;
    }

    /* Loop over the 2^M vertices surrounding P. */
    nv = 1 << M;
    for (iv = 0; iv < nv; iv++) {
      /* Locate vertex in the coordinate array and compute its weight. */
      offset = 0;
      wgt = 1.0;
      for (m = M-1; m >= 0; m--) {
        offset *= tab->K[m];
        offset += tab->p0[m];
        if (iv & (1 << m)) {
          if (tab->K[m] > 1) offset++;
          wgt *= tab->delta[m];
        } else {
          wgt *= 1.0 - tab->delta[m];
        }
      }

      if (wgt == 0.0) continue;

      /* Add the contribution from this vertex to each element. */
      coord = tab->coord + offset*M;
      for (m = 0; m < M; m++) {
        i = tab->map[m];
        *(wp+i) += *(coord++) * wgt;
      }

      if (wgt == 1.0) break;
    }

    *statp = 0;

next:
    xp += nelem;
    wp += nelem;
    statp++;
  }

  return status;
}

/*--------------------------------------------------------------------------*/

int tabs2x(
  struct tabprm* tab,
  int ncoord,
  int nelem,
  const double world[],
  double x[],
  int stat[])

{
  static const char *function = "tabs2x";

  int tabedge(struct tabprm *);
  int tabrow(struct tabprm *, const double *);
  int tabvox(struct tabprm *, const double *, int, double **, unsigned int *);

  int edge, i, ic, iv, k, *Km, M, m, n, nv, offset, status;
  double *dcrd, delta, *Psi, psi_m, **tabcoord, upsilon;
  register int *statp;
  register const double *wp;
  register double *xp;
  struct wcserr **err;

  if (tab == 0x0) return TABERR_NULL_POINTER;
  err = &(tab->err);

  /* Initialize if required. */
  if (tab->flag != TABSET) {
    if ((status = tabset(tab))) return status;
  }

  /* This is used a lot. */
  M = tab->M;

  tabcoord = 0x0;
  nv = 0;
  if (M > 1) {
    nv = 1 << M;
    tabcoord = calloc(nv, sizeof(double *));
  }


  status = 0;
  wp = world;
  xp = x;
  statp = stat;
  for (n = 0; n < ncoord; n++) {
    /* Locate this coordinate in the coordinate array. */
    edge = 0;
    for (m = 0; m < M; m++) {
      tab->p0[m] = 0;
    }

    for (ic = 0; ic < tab->nc; ic++) {
      if (tab->p0[0] == 0) {
        /* New row, could it contain a solution? */
        if (edge || tabrow(tab, wp)) {
          /* No, skip it. */
          ic += tab->K[0];
          tab->p0[1]++;
          edge = tabedge(tab);

          /* Because ic will be incremented when the loop is reentered. */
          ic--;
          continue;
        }
      }

      if (M == 1) {
        /* Deal with the one-dimensional case separately for efficiency. */
        if (*wp == tab->coord[0]) {
          tab->p0[0] = 0;
          tab->delta[0] = 0.0;
          break;

        } else if (ic < tab->nc - 1) {
          if (((tab->coord[ic] <= *wp && *wp <= tab->coord[ic+1]) ||
               (tab->coord[ic] >= *wp && *wp >= tab->coord[ic+1])) &&
               (tab->index[0] == 0x0 ||
                tab->index[0][ic] != tab->index[0][ic+1])) {
            tab->p0[0] = ic;
            tab->delta[0] = (*wp - tab->coord[ic]) /
                            (tab->coord[ic+1] - tab->coord[ic]);
            break;
          }
        }

      } else {
        /* Multi-dimensional tables are harder. */
        if (!edge) {
          /* Addresses of the coordinates for each corner of the "voxel". */
          for (iv = 0; iv < nv; iv++) {
            offset = 0;
            for (m = M-1; m >= 0; m--) {
              offset *= tab->K[m];
              offset += tab->p0[m];
              if ((iv & (1 << m)) && (tab->K[m] > 1)) offset++;
            }
            tabcoord[iv] = tab->coord + offset*M;
          }

          if (tabvox(tab, wp, 0, tabcoord, 0x0) == 0) {
            /* Found a solution. */
            break;
          }
        }

        /* Next voxel. */
        tab->p0[0]++;
        edge = tabedge(tab);
      }
    }


    if (ic == tab->nc) {
      /* Coordinate not found; allow minor extrapolation. */
      if (M == 1) {
        /* Should there be a solution? */
        if (tab->extrema[0] <= *wp && *wp <= tab->extrema[1]) {
          dcrd = tab->coord;
          for (i = 0; i < 2; i++) {
            if (i) dcrd += tab->K[0] - 2;

            delta = (*wp - *dcrd) / (*(dcrd+1) - *dcrd);

            if (i == 0) {
              if (-0.5 <= delta && delta <= 0.0) {
                tab->p0[0] = 0;
                tab->delta[0] = delta;
                ic = 0;
                break;
              }
            } else {
              if (1.0 <= delta && delta <= 1.5) {
                tab->p0[0] = tab->K[0] - 1;
                tab->delta[0] = delta - 1.0;
                ic = 0;
              }
            }
          }
        }

      } else {
        /* Multi-dimensional tables. */
        /* >>> TBD <<< */
      }
    }


    if (ic == tab->nc) {
      /* Coordinate not found. */
      *statp = 1;
      status = wcserr_set(TAB_ERRMSG(TABERR_BAD_WORLD));

    } else {
      /* Determine the intermediate world coordinates. */
      Km = tab->K;
      for (m = 0; m < M; m++, Km++) {
        /* N.B. Upsilon_m and psi_m are 1-relative FITS indexes. */
        upsilon = (tab->p0[m] + 1) + tab->delta[m];

        if (upsilon < 0.5 || upsilon > *Km + 0.5) {
          /* Index out of range. */
          *statp = 1;
          status = wcserr_set(TAB_ERRMSG(TABERR_BAD_WORLD));

        } else {
          /* Do inverse lookup of the index vector. */
          Psi = tab->index[m];
          if (Psi == 0x0) {
            /* Default indexing. */
            psi_m = upsilon;

          } else {
            /* Decrement Psi and use 1-relative C array indexing to match the
               1-relative FITS indexing. */
            Psi--;

            if (*Km == 1) {
              /* Degenerate index vector. */
              psi_m = Psi[1];
            } else {
              k = (int)(upsilon);
              psi_m = Psi[k];
              if (k < *Km) {
                psi_m += (upsilon - k) * (Psi[k+1] - Psi[k]);
              }
            }
          }

          i = tab->map[m];
          xp[i] = psi_m - tab->crval[m];
        }
      }
      *statp = 0;
    }

    wp += nelem;
    xp += nelem;
    statp++;
  }

  if (tabcoord) free(tabcoord);

  return status;
}

/*----------------------------------------------------------------------------
* Convenience routine to deal with of edge effects in tabprm::p0.
*---------------------------------------------------------------------------*/

int tabedge(struct tabprm* tab)

{
  int edge, *Km, m;

  edge = 0;
  Km = tab->K;
  for (m = 0; m < tab->M; m++, Km++) {
    if (tab->p0[m] == *Km) {
      /* p0 has been incremented beyond the end of the row, point it to the
         next one. */
      tab->p0[m] = 0;
      tab->p0[m+1]++;
    } else if (tab->p0[m] == *Km - 1 && *Km > 1) {
      /* p0 is sitting at the end of a non-degenerate row. */
      edge = 1;
    }
  }

  return edge;
}

/*----------------------------------------------------------------------------
* Quick test to see whether the world coordinate indicated by wp could lie
* somewhere along (or near) the row of the image indexed by tabprm::p0.
* Return 0 if so, 1 otherwise.
*
* tabprm::p0 selects a particular row of the image, p0[0] being ignored (i.e.
* treated as zero).  Adjacent rows that delimit a row of "voxels" are formed
* by incrementing elements other than p0[0] in all binary combinations.  N.B.
* these are not the same as the voxels (pixels) that are indexed by, and
* centred on, integral pixel coordinates in FITS.
*
* To see why it is necessary to examine the adjacent rows, consider the 2-D
* case where the first world coordinate element is constant along each row.
* If the first element of wp has value 0.5, and its value in the row indexed
* by p0 has value 0, and in the next row it has value 1, then it is clear that
* the solution lies in neither row but somewhere between them.  Thus both rows
* will be involved in finding the solution.
*
* tabprm::extrema is the address of the first element of a 1-D array that
* records the minimum and maximum value of each element of the coordinate
* vector in each row of the coordinate array, treated as though it were
* defined as
*
*   double extrema[K_M]...[K_2][2][M]
*
* The minimum is recorded in the first element of the compressed K_1
* dimension, then the maximum.
*---------------------------------------------------------------------------*/

int tabrow(struct tabprm* tab, const double *wp)

{
  int iv, M, m, nv, offset;
  unsigned int eq, gt, lt;
  const double tol = 1e-10;
  double *cp, w;

  M = tab->M;

  /* The number of corners in a "voxel".  We need examine only half this
     number of rows.  The extra factor of two will be used to select between
     the minimal and maximal values in each row. */
  nv = 1 << M;

  eq = 0;
  lt = 0;
  gt = 0;
  for (iv = 0; iv < nv; iv++) {
    /* Find the index into tabprm::extrema for this row. */
    offset = 0;
    for (m = M-1; m > 0; m--) {
      offset *= tab->K[m];
      offset += tab->p0[m];

      /* Select the row. */
      if (iv & (1 << m)) {
        if (tab->K[m] > 1) offset++;
      }
    }

    /* The K_1 dimension has length 2 (see prologue). */
    offset *= 2;

    /* Select the minimum on even numbered iterations, else the maximum. */
    if (iv & 1) offset++;

    /* The last dimension has length M (see prologue). */
    offset *= M;

    /* Address of the extremal elements (min or max) for this row. */
    cp = tab->extrema + offset;

    /* For each coordinate element, we only need to find one row where its
       minimum value is less than that of wp, and one row where the maximum
       value is greater.  That doesn't mean that there is a solution, only
       that there might be. */
    for (m = 0; m < M; m++, cp++) {
      /* Apply the axis mapping. */
      w = wp[tab->map[m]];

      /* Finally the test itself; set bits in the bitmask. */
      if (fabs(*cp - w) < tol) {
        eq |= (1 << m);
      } else if (*cp < w) {
        lt |= (1 << m);
      } else if (*cp > w) {
        gt |= (1 << m);
      }
    }

    /* Have all bits been switched on? */
    if ((lt | eq) == nv-1 && (gt | eq) == nv-1) {
      /* A solution could lie within this row of voxels. */
      return 0;
    }
  }

  /* No solution in this row. */
  return 1;
}

/*----------------------------------------------------------------------------
* Does the world coordinate indicated by wp lie within the voxel indexed by
* tabprm::p0?  If so, do a binary chop of the interior of the voxel to find
* it and return 0, with tabprm::delta set to the solution.  Else return 1.
*
* As in tabrow(), a "voxel" is formed by incrementing the elements of
* tabprm::p0 in all binary combinations.  Note that these are not the same as
* the voxels (pixels) that are indexed by, and centred on, integral pixel
* coordinates in FITS.
*
* tabvox() calls itself recursively.  When called from outside, level, being
* the level of recursion, should be given as zero.  tabcoord is an array
* holding the addresses of the coordinates for each corner of the voxel.
* vox is the address of a work array (vox2) used during recursive calls to
* dissect the voxel.  It is ignored when tabvox() is called from outside
* (level == 0).
*
* It is assumed that the image dimensions are no greater than 16.
----------------------------------------------------------------------------*/

int tabvox(
  struct tabprm* tab,
  const double *wp,
  int level,
  double **tabcoord,
  unsigned int *vox)

{
  int i, iv, jv, M, m, nv;
  unsigned int eq, et, gt, lt, vox2[16];
  const double tol = 1e-10;
  double coord[16], *cp, dv, w, wgt;

  M = tab->M;

  /* The number of corners in a voxel. */
  nv = 1 << M;

  dv = 1.0;
  for (i = 0; i < level; i++) {
    dv /= 2.0;
  }

  /* Could the coordinate lie within this voxel (level == 0) or sub-voxel
     (level > 0)?  We use the fact that with linear interpolation the
     coordinate elements are extremal in a corner and test each one. */
  lt = 0;
  gt = 0;
  eq = 0;
  for (iv = 0; iv < nv; iv++) {
    /* Select a corner of the sub-voxel. */
    for (m = 0; m < M; m++) {
      coord[m] = 0.0;
      tab->delta[m] = level ? dv*vox[m] : 0.0;

      if (iv & (1 << m)) {
        tab->delta[m] += dv;
      }
    }

    /* Compute the coordinates of this corner of the sub-voxel by linear
       interpolation using the weighting algorithm described in Sect. 3.4 of
       WCS Paper IV. */
    for (jv = 0; jv < nv; jv++) {
      /* Find the weight for this corner of the parent voxel. */
      wgt = 1.0;
      for (m = 0; m < M; m++) {
        if (jv & (1 << m)) {
          wgt *= tab->delta[m];
        } else {
          wgt *= 1.0 - tab->delta[m];
        }
      }

      if (wgt == 0.0) continue;

      /* Add its contribution to each coordinate element. */
      cp = tabcoord[jv];
      for (m = 0; m < M; m++) {
        coord[m] += *(cp++) * wgt;
      }

      if (wgt == 1.0) break;
    }

    /* Coordinate elements are minimal or maximal in a corner. */
    et = 0;
    for (m = 0; m < M; m++) {
      /* Apply the axis mapping. */
      w = wp[tab->map[m]];

      /* Finally the test itself; set bits in the bitmask. */
      if (fabs(coord[m] - w) < tol) {
        et |= (1 << m);
      } else if (coord[m] < w) {
        lt |= (1 << m);
      } else if (coord[m] > w) {
        gt |= (1 << m);
      }
    }

    if (et == nv-1) {
      /* We've stumbled across a solution in this corner of the sub-voxel. */
      return 0;
    }

    eq |= et;
  }

  /* Could the coordinate lie within this sub-voxel? */
  if ((lt | eq) == nv-1 && (gt | eq) == nv-1) {
    /* Yes it could, but does it? */

    /* Is it time to stop the recursion? */
    if (level == 31) {
      /* We have a solution, squeeze out the last bit of juice. */
      dv /= 2.0;
      for (m = 0; m < M; m++) {
        tab->delta[m] = dv * (2.0*vox[m] + 1.0);
      }

      return 0;
    }

    /* Subdivide the sub-voxel and try again for each subdivision. */
    for (iv = 0; iv < nv; iv++) {
      /* Select the subdivision. */
      for (m = 0; m < M; m++) {
        vox2[m] = level ? 2*vox[m] : 0;
        if (iv & (1 << m)) {
          vox2[m]++;
        }
      }

      /* Recurse. */
      if (tabvox(tab, wp, level+1, tabcoord, vox2) == 0) {
        return 0;
      }
    }
  }

  /* No solution in this sub-voxel. */
  return 1;
}
