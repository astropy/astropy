/*============================================================================

  WCSLIB 4.8 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2011, Mark Calabretta

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
  along with WCSLIB.  If not, see <http://www.gnu.org/licenses/>.

  Correspondence concerning WCSLIB may be directed to:
    Internet email: mcalabre@atnf.csiro.au
    Postal address: Dr. Mark Calabretta
                    Australia Telescope National Facility, CSIRO
                    PO Box 76
                    Epping NSW 1710
                    AUSTRALIA

  Author: Mark Calabretta, Australia Telescope National Facility
  http://www.atnf.csiro.au/~mcalabre/index.html
  $Id: wcsutil.c,v 4.8.1.1 2011/08/15 08:07:06 cal103 Exp cal103 $
*===========================================================================*/

#include <stdio.h>
#include <string.h>

#include "wcsutil.h"

/*--------------------------------------------------------------------------*/

void wcsutil_blank_fill(int n, char c[])

{
  int k;

  for (k = strlen(c); k < n; k++) {
    c[k] = ' ';
  }

  return;
}

/*--------------------------------------------------------------------------*/

void wcsutil_null_fill(int n, char c[])

{
  int j, k;

  if (n <= 0) return;

  /* Null-fill the string. */
  *(c+n-1) = '\0';
  for (j = 0; j < n; j++) {
    if (c[j] == '\0') {
      for (k = j+1; k < n; k++) {
        c[k] = '\0';
      }
      break;
    }
  }

  for (k = j-1; k > 0; k--) {
    if (c[k] != ' ') break;
    c[k] = '\0';
  }

   return;
}

/*--------------------------------------------------------------------------*/

int wcsutil_allEq(int nvec, int nelem, const double *first)

{
  double v0;
  const double *vp;

  if (nvec <= 0 || nelem <= 0) return 0;

  v0 = *first;
  for (vp = first+nelem; vp < first + nvec*nelem; vp += nelem) {
    if (*vp != v0) return 0;
  }

  return 1;
}

/*--------------------------------------------------------------------------*/

void wcsutil_setAll(int nvec, int nelem, double *first)

{
  double v0, *vp;

  if (nvec <= 0 || nelem <= 0) return;

  v0 = *first;
  for (vp = first+nelem; vp < first + nvec*nelem; vp += nelem) {
    *vp = v0;
  }
}

/*--------------------------------------------------------------------------*/

void wcsutil_setAli(int nvec, int nelem, int *first)

{
  int v0, *vp;

  if (nvec <= 0 || nelem <= 0) return;

  v0 = *first;
  for (vp = first+nelem; vp < first + nvec*nelem; vp += nelem) {
    *vp = v0;
  }
}

/*--------------------------------------------------------------------------*/

void wcsutil_setBit(int nelem, const int *sel, int bits, int *array)

{
  int *arrp;

  if (bits == 0 || nelem <= 0) return;

  if (sel == 0x0) {
    /* All elements selected. */
    for (arrp = array; arrp < array + nelem; arrp++) {
      *arrp |= bits;
    }

  } else {
    /* Some elements selected. */
    for (arrp = array; arrp < array + nelem; arrp++) {
      if (*(sel++)) *arrp |= bits;
    }
  }
}

/*--------------------------------------------------------------------------*/

char *wcsutil_fptr2str(int (*func)(), char hext[])

{
  unsigned char *p = (unsigned char *)(&func);
  char *t = hext;
  int i, *(ip[2]), j[2], le = 1, gotone = 0;

  /* Test for little-endian addresses. */
  ip[0] = j;
  ip[1] = j + 1;
  if ((unsigned char *)ip[0] < (unsigned char *)ip[1]) {
    /* Little-endian, reverse it. */
    p += sizeof(func) - 1;
    le = -1;
  }

  sprintf(t, "0x0");
  t += 2;

  for (i = 0; i < sizeof(func); i++) {
    /* Skip leading zeroes. */
    if (*p) gotone = 1;

    if (gotone) {
      sprintf(t, "%02x", *p);
      t += 2;
    }

    p += le;
  }

  return hext;
}
