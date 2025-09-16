/*============================================================================
  WCSLIB 8.4 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2024, Mark Calabretta

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
  $Id: wcsutil.c,v 8.4 2024/10/28 13:56:16 mcalabre Exp $
*===========================================================================*/

#include <ctype.h>
#include <locale.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcsutil.h"
#include "wcsmath.h"

//----------------------------------------------------------------------------

void wcsdealloc(void *ptr)

{
  free(ptr);

  return;
}

//----------------------------------------------------------------------------

void wcsutil_strcvt(int n, char c, int nt, const char src[], char dst[])

{
  if (n <= 0) return;

  if (c != '\0') c = ' ';

  if (src == 0x0) {
    if (dst) {
      memset(dst, c, n);
    }

  } else {
    // Copy to the first NULL character.
    int j;
    for (j = 0; j < n; j++) {
      if ((dst[j] = src[j]) == '\0') {
        break;
      }
    }

    if (j < n) {
      // The given string is null-terminated.
      memset(dst+j, c, n-j);

    } else {
      // The given string is not null-terminated.
      if (c == '\0') {
        // Work backwards, looking for the first non-blank.
        for (j = n - 1; j >= 0; j--) {
          if (dst[j] != ' ') {
            break;
          }
        }

        j++;
	if (j == n && !nt) {
	  dst[n-1] = '\0';
	} else {
          memset(dst+j, '\0', n-j);
	}
      }
    }
  }

  if (nt) dst[n] = '\0';

  return;
}

//----------------------------------------------------------------------------

void wcsutil_blank_fill(int n, char c[])

{
  if (n <= 0) return;

  if (c == 0x0) {
    return;
  }

  // Replace the terminating null and all successive characters.
  for (int j = 0; j < n; j++) {
    if (c[j] == '\0') {
      memset(c+j, ' ', n-j);
      break;
    }
  }

  return;
}

//----------------------------------------------------------------------------

void wcsutil_null_fill(int n, char c[])

{
  if (n <= 0) return;

  if (c == 0x0) {
    return;
  }

  // Find the first NULL character.
  int j;
  for (j = 0; j < n; j++) {
    if (c[j] == '\0') {
      break;
    }
  }

  // Ensure null-termination.
  if (j == n) {
    j = n - 1;
    c[j] = '\0';
  }

  // Work backwards, looking for the first non-blank.
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

//----------------------------------------------------------------------------

int wcsutil_all_ival(int nelem, int ival, const int iarr[])

{
  for (int i = 0; i < nelem; i++) {
    if (iarr[i] != ival) return 0;
  }

  return 1;
}

//----------------------------------------------------------------------------

int wcsutil_all_dval(int nelem, double dval, const double darr[])

{
  for (int i = 0; i < nelem; i++) {
    if (darr[i] != dval) return 0;
  }

  return 1;
}

//----------------------------------------------------------------------------

int wcsutil_all_sval(int nelem, const char *sval, const char (*sarr)[72])

{
  for (int i = 0; i < nelem; i++) {
    if (strncmp(sarr[i], sval, 72)) return 0;
  }

  return 1;
}

//----------------------------------------------------------------------------

int wcsutil_allEq(int nvec, int nelem, const double *first)

{
  if (nvec <= 0 || nelem <= 0) return 0;

  double v0 = *first;
  for (const double *vp = first+nelem; vp < first + nvec*nelem; vp += nelem) {
    if (*vp != v0) return 0;
  }

  return 1;
}

//----------------------------------------------------------------------------

int wcsutil_dblEq(
  int nelem,
  double tol,
  const double *darr1,
  const double *darr2)

{
  if (nelem == 0) return 1;
  if (nelem  < 0) return 0;

  if (darr1 == 0x0 && darr2 == 0x0) return 1;

  if (tol == 0.0) {
    // Handled separately for speed of execution.
    for (int i = 0; i < nelem; i++) {
      double dval1 = (darr1 ? darr1[i] : UNDEFINED);
      double dval2 = (darr2 ? darr2[i] : UNDEFINED);

      // Undefined values must match exactly.
      if (dval1 == UNDEFINED && dval2 != UNDEFINED) return 0;
      if (dval1 != UNDEFINED && dval2 == UNDEFINED) return 0;

      if (dval1 != dval2) return 0;
    }

  } else {
    for (int i = 0; i < nelem; i++) {
      double dval1 = (darr1 ? darr1[i] : UNDEFINED);
      double dval2 = (darr2 ? darr2[i] : UNDEFINED);

      // Undefined values must match exactly.
      if (dval1 == UNDEFINED && dval2 != UNDEFINED) return 0;
      if (dval1 != UNDEFINED && dval2 == UNDEFINED) return 0;

      // Otherwise, compare within the specified tolerance.
      if (fabs(dval1 - dval2) > 0.5*tol) return 0;
    }
  }

  return 1;
}

//----------------------------------------------------------------------------

int wcsutil_intEq(int nelem, const int *iarr1, const int *iarr2)

{
  if (nelem == 0) return 1;
  if (nelem  < 0) return 0;

  if (iarr1 == 0x0 && iarr2 == 0x0) return 1;

  for (int i = 0; i < nelem; i++) {
    int ival1 = (iarr1 ?  iarr1[i] : 0);
    int ival2 = (iarr2 ?  iarr2[i] : 0);

    if (ival1 != ival2) return 0;
  }

  return 1;
}

//----------------------------------------------------------------------------

int wcsutil_strEq(int nelem, char (*sarr1)[72], char (*sarr2)[72])

{
  if (nelem == 0) return 1;
  if (nelem  < 0) return 0;

  if (sarr1 == 0x0 && sarr2 == 0x0) return 1;

  for (int i = 0; i < nelem; i++) {
    char *sval1 = (sarr1 ?  sarr1[i] : "");
    char *sval2 = (sarr2 ?  sarr2[i] : "");

    if (strncmp(sval1, sval2, 72)) return 0;
  }

  return 1;
}

//----------------------------------------------------------------------------

void wcsutil_setAll(int nvec, int nelem, double *first)

{
  if (nvec <= 0 || nelem <= 0) return;

  double v0 = *first;
  for (double *vp = first+nelem; vp < first + nvec*nelem; vp += nelem) {
    *vp = v0;
  }
}

//----------------------------------------------------------------------------

void wcsutil_setAli(int nvec, int nelem, int *first)

{
  if (nvec <= 0 || nelem <= 0) return;

  int v0 = *first;
  for (int *vp = first+nelem; vp < first + nvec*nelem; vp += nelem) {
    *vp = v0;
  }
}

//----------------------------------------------------------------------------

void wcsutil_setBit(int nelem, const int *sel, int bits, int *array)

{
  if (bits == 0 || nelem <= 0) return;

  if (sel == 0x0) {
    // All elements selected.
    for (int *arrp = array; arrp < array + nelem; arrp++) {
      *arrp |= bits;
    }

  } else {
    // Some elements selected.
    for (int *arrp = array; arrp < array + nelem; arrp++) {
      if (*(sel++)) *arrp |= bits;
    }
  }
}

//----------------------------------------------------------------------------

char *wcsutil_fptr2str(void (*fptr)(void), char hext[19])

{
  // Test for little-endian addresses.
  int *(ip[2]), j[2], le = 1;
  ip[0] = j;
  ip[1] = j + 1;
  unsigned char *p = (unsigned char *)(&fptr);
  if ((unsigned char *)ip[0] < (unsigned char *)ip[1]) {
    // Little-endian, reverse it.
    p += sizeof(fptr) - 1;
    le = -1;
  }

  char *t = hext;
  sprintf(t, "0x0");
  t += 2;

  int gotone = 0;
  for (size_t i = 0; i < sizeof(fptr); i++) {
    // Skip leading zeroes.
    if (*p) gotone = 1;

    if (gotone) {
      sprintf(t, "%02x", *p);
      t += 2;
    }

    p += le;
  }

  return hext;
}

//----------------------------------------------------------------------------

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
  sprintf(buf, format, value);
  wcsutil_locale_to_dot(buf);

  // Look for a decimal point or exponent.
  char *bp = buf;
  while (*bp) {
    if (*bp != ' ') {
      if (*bp == '.') return;
      if (*bp == 'e') return;
      if (*bp == 'E') return;
    }
    bp++;
  }

  // Not found, add a fractional part.
  bp = buf;
  if (*bp == ' ') {
    char *cp = buf + 1;
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

//----------------------------------------------------------------------------

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
  value[0] = 0.0;
  value[1] = 0.0;

  // Get the integer part.
  char ltmp[72];
  if (sscanf(wcsutil_dot_to_locale(buf, ltmp), "%lf", value) < 1) {
    return 1;
  }
  value[0] = floor(value[0]);

  char ctmp[72];
  strcpy(ctmp, buf);

  // Look for a decimal point.
  char *dptr = strchr(ctmp, '.');

  // Look for an exponent.
  char *eptr;
  if ((eptr = strchr(ctmp, 'E')) == NULL) {
    if ((eptr = strchr(ctmp, 'D')) == NULL) {
      if ((eptr = strchr(ctmp, 'e')) == NULL) {
        eptr = strchr(ctmp, 'd');
      }
    }
  }

  int exp = 0;
  if (eptr) {
    // Get the exponent.
    if (sscanf(eptr+1, "%d", &exp) < 1) {
      return 1;
    }

    if (!dptr) {
      dptr = eptr;
      eptr++;
    }

    if (dptr+exp <= ctmp) {
      // There is only a fractional part.
      return sscanf(wcsutil_dot_to_locale(buf, ctmp), "%lf", value+1) < 1;
    } else if (eptr <= dptr+exp+1) {
      // There is no fractional part.
      return 0;
    }
  }

  // Get the fractional part.
  if (dptr) {
    char *cptr = ctmp;
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
