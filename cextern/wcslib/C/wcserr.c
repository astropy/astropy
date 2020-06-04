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
  Module author: Michael Droettboom
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: wcserr.c,v 7.3 2020/06/03 03:37:02 mcalabre Exp $
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

int wcserr_prt(const struct wcserr *err, const char *prefix)

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

int wcserr_clear(struct wcserr **errp)

{
  if (*errp) {
    if ((*errp)->msg) {
      free((*errp)->msg);
    }
    free(*errp);
    *errp = 0x0;
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int wcserr_set(
  struct wcserr **errp,
  int status,
  const char *function,
  const char *file,
  int line_no,
  const char *format,
  ...)

{
  int  msglen;
  struct wcserr *err;
  va_list argp;

  if (!wcserr_enabled) return status;

  if (errp == 0x0) {
    return status;
  }
  err = *errp;

  if (status) {
    if (err == 0x0) {
      *errp = err = calloc(1, sizeof(struct wcserr));
    }

    if (err == 0x0) {
      return status;
    }

    err->status   = status;
    err->function = function;
    err->file     = file;
    err->line_no  = line_no;
    err->msg      = 0x0;

    /* Determine the required message buffer size. */
    va_start(argp, format);
    msglen = vsnprintf(0x0, 0, format, argp) + 1;
    va_end(argp);

    if (msglen <= 0 || (err->msg = malloc(msglen)) == 0x0) {
      wcserr_clear(errp);
      return status;
    }

    /* Write the message. */
    va_start(argp, format);
    msglen = vsnprintf(err->msg, msglen, format, argp);
    va_end(argp);

    if (msglen < 0) {
      wcserr_clear(errp);
    }
  }

  return status;
}

/*--------------------------------------------------------------------------*/

int wcserr_copy(const struct wcserr *src, struct wcserr *dst)

{
  size_t msglen;

  if (src == 0x0) {
    if (dst) {
      memset(dst, 0, sizeof(struct wcserr));
    }
    return 0;
  }

  if (dst) {
    memcpy(dst, src, sizeof(struct wcserr));

    if (src->msg) {
      msglen = strlen(src->msg) + 1;
      if ((dst->msg = malloc(msglen))) {
        strcpy(dst->msg, src->msg);
      }
    }
  }

  return src->status;
}
