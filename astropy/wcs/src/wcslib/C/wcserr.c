/*============================================================================

  WCSLIB 4.10 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2012, Mark Calabretta

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
  Module author: Michael Droettboom
  http://www.atnf.csiro.au/~mcalabre/index.html
  $Id: wcserr.c,v 4.10 2012/02/05 23:41:44 cal103 Exp $
*===========================================================================*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcserr.h"
#include "wcsprintf.h"

static int wcserr_enabled = 0;

/*--------------------------------------------------------------------------*/

int wcserr_enable(int enable)

{
  return wcserr_enabled = (enable ? 1 : 0);
}

/*--------------------------------------------------------------------------*/

int wcserr_prt(
  const struct wcserr *err,
  const char *prefix)

{
  if (!wcserr_enabled) {
    wcsprintf("Error messaging is not enabled, use wcserr_enable().\n");
    return 2;
  }

  if (err == 0x0) {
    return 0;
  }

  if (err->status) {
    if (prefix == 0x0) prefix = "";

    if (err->status > 0) {
      wcsprintf("%sERROR %d in %s() at line %d of file %s:\n%s%s.\n",
        prefix, err->status, err->function, err->line_no, err->file, prefix,
        err->msg);
    } else {
      /* An informative message only. */
      wcsprintf("%sINFORMATIVE message from %s() at line %d of file "
        "%s:\n%s%s.\n", prefix, err->function, err->line_no, err->file,
        prefix, err->msg);
    }
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int wcserr_clear(
  struct wcserr **err)

{
  if (*err) free(*err);
  *err = 0x0;

  return 0;
}

/*--------------------------------------------------------------------------*/

int wcserr_set(
  struct wcserr **err,
  int status,
  const char *function,
  const char *file,
  int line_no,
  const char *format,
  ...)

{
  va_list argp;

  if (!wcserr_enabled) return status;

  if (err == 0x0) {
    return status;
  }

  if (status) {
    if (*err == 0x0) {
      *err = calloc(1, sizeof(struct wcserr));
    }

    (*err)->status   = status;
    (*err)->function = function;
    (*err)->file     = file;
    (*err)->line_no  = line_no;

    va_start(argp, format);
    vsnprintf((*err)->msg, WCSERR_MSG_LENGTH, format, argp);
    va_end(argp);

  } else {
    if (*err) free(*err);
    *err = 0x0;
  }

  return status;
}

/*--------------------------------------------------------------------------*/

int wcserr_copy(
  const struct wcserr *src,
  struct wcserr *dst)

{
  if (src == 0x0) {
    if (dst) {
      memset(dst, 0, sizeof(struct wcserr));
    }
    return 0;
  }

  if (dst) {
    memcpy(dst, src, sizeof(struct wcserr));
  }

  return src->status;
}
