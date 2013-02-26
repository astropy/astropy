/*============================================================================

  WCSLIB 4.17 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2013, Mark Calabretta

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
  $Id: wcsutil.c,v 4.17 2013/01/29 05:29:20 cal103 Exp $
*===========================================================================*/

#include <ctype.h>
#include <locale.h>
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

char *wcsutil_fptr2str(int (*func)(void), char hext[19])

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

/*--------------------------------------------------------------------------*/

static const char *wcsutil_dot_to_locale(const char *inbuf, char *outbuf)

{
  struct lconv *locale_data = localeconv();
  const char *decimal_point = locale_data->decimal_point;

  if (decimal_point[0] != '.' || decimal_point[1] != 0) {
    char *out = outbuf;
    size_t decimal_point_len = strlen(decimal_point);

    for ( ; *inbuf; ++inbuf) {
      if (*inbuf == '.') {
        strncpy(out, decimal_point, decimal_point_len);
        out += decimal_point_len;
      } else {
        *out++ = *inbuf;
      }
    }

    *out = '\0';

    return outbuf;
  } else {
    return inbuf;
  }
}


int wcsutil_str2double(const char *buf, const char *format, double *value)

{
  char ctmp[72];
  return sscanf(wcsutil_dot_to_locale(buf, ctmp), "%lf", value) < 1;
}

/*--------------------------------------------------------------------------*/

static void wcsutil_locale_to_dot(char *buf)

{
  struct lconv *locale_data = localeconv();
  const char *decimal_point = locale_data->decimal_point;

  if (decimal_point[0] != '.' || decimal_point[1] != 0) {
    size_t decimal_point_len = strlen(decimal_point);
    char *inbuf = buf;
    char *outbuf = buf;

    for ( ; *inbuf; ++inbuf) {
      if (strncmp(inbuf, decimal_point, decimal_point_len) == 0) {
        *outbuf++ = '.';
        inbuf += decimal_point_len - 1;
      } else {
        *outbuf++ = *inbuf;
      }
    }

    *outbuf = '\0';
  }
}


void wcsutil_double2str(char *buf, const char *format, double value)

{
  sprintf(buf, format, value);
  wcsutil_locale_to_dot(buf);
}
