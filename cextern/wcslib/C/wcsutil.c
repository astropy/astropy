/*============================================================================

  WCSLIB 7.3 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2020, Mark Calabretta

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
  $Id: wcsutil.c,v 7.3 2020/06/03 03:37:02 mcalabre Exp $
*===========================================================================*/

#include <ctype.h>
#include <locale.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcsutil.h"
#include "wcsmath.h"

/*--------------------------------------------------------------------------*/

void wcsdealloc(void *ptr)

{
  free(ptr);

  return;
}

/*--------------------------------------------------------------------------*/

void wcsutil_strcvt(int n, char c, const char src[], char dst[])

{
  int j;

  if (n <= 0) return;

  if (c != '\0') c = ' ';

  if (src == 0x0) {
    if (dst) {
       memset(dst, c, n);
    }

    return;
  }

  /* Copy to the first NULL character. */
  for (j = 0; j < n; j++) {
    if ((dst[j] = src[j]) == '\0') {
      break;
    }
  }

  if (j < n) {
    /* The given string is null-terminated. */
    memset(dst+j, c, n-j);

  } else {
    /* The given string is not null-terminated. */
    if (c == '\0') {
      j = n - 1;
      dst[j] = '\0';

      j--;

      /* Work backwards, looking for the first non-blank. */
      for (; j >= 0; j--) {
        if (dst[j] != ' ') {
          break;
        }
      }

      j++;
      memset(dst+j, '\0', n-j);
    }
  }

  return;
}

/*--------------------------------------------------------------------------*/

void wcsutil_blank_fill(int n, char c[])

{
  int j;

  if (n <= 0) return;

  if (c == 0x0) {
    return;
  }

  /* Replace the terminating null and all successive characters. */
  for (j = 0; j < n; j++) {
    if (c[j] == '\0') {
      memset(c+j, ' ', n-j);
      break;
    }
  }

  return;
}

/*--------------------------------------------------------------------------*/

void wcsutil_null_fill(int n, char c[])

{
  int j;

  if (n <= 0) return;

  if (c == 0x0) {
    return;
  }

  /* Find the first NULL character. */
  for (j = 0; j < n; j++) {
    if (c[j] == '\0') {
      break;
    }
  }

  /* Ensure null-termination. */
  if (j == n) {
    j = n - 1;
    c[j] = '\0';
  }

  /* Work backwards, looking for the first non-blank. */
  j--;
  for (; j > 0; j--) {
    if (c[j] != ' ') {
      break;
    }
  }

  if (++j < n) {
    memset(c+j, '\0', n-j);
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

int wcsutil_Eq(int nelem, double tol, const double *arr1, const double *arr2)

{
  int i;

  if (nelem == 0) return 1;
  if (nelem  < 0) return 0;

  if (arr1 == 0x0 && arr2 == 0x0) return 1;
  if (arr1 == 0x0 || arr2 == 0x0) return 0;

  if (tol == 0.0) {
    /* Handled separately for speed of execution. */
    for (i = 0; i < nelem; i++, arr1++, arr2++) {
      if (*arr1 != *arr2) return 0;
    }

  } else {
    for (i = 0; i < nelem; i++, arr1++, arr2++) {
      /* Undefined values must match exactly. */
      if (*arr1 == UNDEFINED && *arr2 != UNDEFINED) return 0;
      if (*arr1 != UNDEFINED && *arr2 == UNDEFINED) return 0;

      /* Otherwise, compare within the specified tolerance. */
      if (fabs(*arr1 - *arr2) > 0.5*tol) return 0;
    }
  }

  return 1;
}

/*--------------------------------------------------------------------------*/

int wcsutil_intEq(int nelem, const int *arr1, const int *arr2)

{
  int i;

  if (nelem == 0) return 1;
  if (nelem  < 0) return 0;

  if (arr1 == 0x0 && arr2 == 0x0) return 1;
  if (arr1 == 0x0 || arr2 == 0x0) return 0;

  for (i = 0; i < nelem; i++, arr1++, arr2++) {
    if (*arr1 != *arr2) return 0;
  }

  return 1;
}

/*--------------------------------------------------------------------------*/

int wcsutil_strEq(int nelem, char (*arr1)[72], char (*arr2)[72])

{
  int i;

  if (nelem == 0) return 1;
  if (nelem  < 0) return 0;

  if (arr1 == 0x0 && arr2 == 0x0) return 1;
  if (arr1 == 0x0 || arr2 == 0x0) return 0;

  for (i = 0; i < nelem; i++, arr1++, arr2++) {
    if (strncmp(*arr1, *arr2, 72)) return 0;
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

char *wcsutil_fptr2str(void (*fptr)(void), char hext[19])

{
  unsigned char *p = (unsigned char *)(&fptr);
  char *t = hext;
  unsigned int i;
  int *(ip[2]), j[2], le = 1, gotone = 0;

  /* Test for little-endian addresses. */
  ip[0] = j;
  ip[1] = j + 1;
  if ((unsigned char *)ip[0] < (unsigned char *)ip[1]) {
    /* Little-endian, reverse it. */
    p += sizeof(fptr) - 1;
    le = -1;
  }

  sprintf(t, "0x0");
  t += 2;

  for (i = 0; i < sizeof(fptr); i++) {
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

static void wcsutil_locale_to_dot(char *buf)

{
  struct lconv *locale_data = localeconv();
  const char *decimal_point = locale_data->decimal_point;

  if (decimal_point[0] != '.' || decimal_point[1] != 0) {
    size_t decimal_point_len = strlen(decimal_point);
    char *inbuf = buf;
    char *outbuf = buf;

    for ( ; *inbuf; inbuf++) {
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
  char *bp, *cp;

  sprintf(buf, format, value);
  wcsutil_locale_to_dot(buf);

  /* Look for a decimal point or exponent. */
  bp = buf;
  while (*bp) {
    if (*bp != ' ') {
      if (*bp == '.') return;
      if (*bp == 'e') return;
      if (*bp == 'E') return;
    }
    bp++;
  }

  /* Not found, add a fractional part. */
  bp = buf;
  if (*bp == ' ') {
    cp = buf + 1;
    if (*cp == ' ') cp++;

    while (*cp) {
      *bp = *cp;
      bp++;
      cp++;
    }

    *bp = '.';
    bp++;
    if (bp < cp) *bp = '0';
  }
}

/*--------------------------------------------------------------------------*/

static const char *wcsutil_dot_to_locale(const char *inbuf, char *outbuf)

{
  struct lconv *locale_data = localeconv();
  const char *decimal_point = locale_data->decimal_point;

  if (decimal_point[0] != '.' || decimal_point[1] != 0) {
    char *out = outbuf;
    size_t decimal_point_len = strlen(decimal_point);

    for ( ; *inbuf; inbuf++) {
      if (*inbuf == '.') {
        memcpy(out, decimal_point, decimal_point_len);
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


int wcsutil_str2double(const char *buf, double *value)

{
  char ctmp[72];
  return sscanf(wcsutil_dot_to_locale(buf, ctmp), "%lf", value) < 1;
}


int wcsutil_str2double2(const char *buf, double *value)

{
  char   *cptr, ctmp[72], *dptr, *eptr, ltmp[72];
  int    exp = 0;

  value[0] = 0.0;
  value[1] = 0.0;

  /* Get the integer part. */
  if (sscanf(wcsutil_dot_to_locale(buf, ltmp), "%lf", value) < 1) {
    return 1;
  }
  value[0] = floor(value[0]);

  strcpy(ctmp, buf);

  /* Look for a decimal point. */
  dptr = strchr(ctmp, '.');

  /* Look for an exponent. */
  if ((eptr = strchr(ctmp, 'E')) == NULL) {
    if ((eptr = strchr(ctmp, 'D')) == NULL) {
      if ((eptr = strchr(ctmp, 'e')) == NULL) {
        eptr = strchr(ctmp, 'd');
      }
    }
  }

  if (eptr) {
    /* Get the exponent. */
    if (sscanf(eptr+1, "%d", &exp) < 1) {
      return 1;
    }

    if (!dptr) {
      dptr = eptr;
      eptr++;
    }

    if (dptr+exp <= ctmp) {
      /* There is only a fractional part. */
      return sscanf(wcsutil_dot_to_locale(buf, ctmp), "%lf", value+1) < 1;
    } else if (eptr <= dptr+exp+1) {
      /* There is no fractional part. */
      return 0;
    }
  }

  /* Get the fractional part. */
  if (dptr) {
    cptr = ctmp;
    while (cptr <= dptr+exp) {
      if ('0' < *cptr && *cptr <= '9') *cptr = '0';
      cptr++;
    }

    if (sscanf(wcsutil_dot_to_locale(ctmp, ltmp), "%lf", value+1) < 1) {
      return 1;
    }
  }

  return 0;
}
