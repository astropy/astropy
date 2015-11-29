/*============================================================================

  WCSLIB 5.10 - an implementation of the FITS WCS standard.
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
  $Id: wcsprintf.c,v 5.10 2015/10/09 08:19:15 mcalabre Exp $
*===========================================================================*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "wcsprintf.h"

static FILE  *wcsprintf_file = 0x0;
static char  *wcsprintf_buff = 0x0;
static char  *wcsprintf_bufp = 0x0;
static size_t wcsprintf_size = 0;

/*--------------------------------------------------------------------------*/

int wcsprintf_set(FILE *wcsout)

{
  if (wcsout != 0x0) {
    /* Output to file. */
    wcsprintf_file = wcsout;

    if (wcsprintf_buff != 0x0) {
      /* Release the buffer. */
      free(wcsprintf_buff);
      wcsprintf_buff = 0x0;
    }

  } else {
    /* Output to buffer. */
    if (wcsprintf_buff == 0x0) {
      /* Allocate a buffer. */
      wcsprintf_buff = malloc(1024);
      if (wcsprintf_buff == NULL) {
        return 1;
      }
      wcsprintf_size = 1024;
    }

    /* Reset pointer to the start of the buffer. */
    wcsprintf_bufp = wcsprintf_buff;
    *wcsprintf_bufp = '\0';
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

const char *wcsprintf_buf(void)

{
  return wcsprintf_buff;
}

/*--------------------------------------------------------------------------*/

int wcsprintf(const char *format, ...)

{
  char *realloc_buff;
  int  nbytes;
  size_t  used;
  va_list arg_list;

  if (wcsprintf_buff == 0x0 && wcsprintf_file == 0x0) {
    /* Send output to stdout if wcsprintf_set() hasn't been called. */
    wcsprintf_file = stdout;
  }

  va_start(arg_list, format);

  if (wcsprintf_file) {
    /* Output to file. */
    nbytes = vfprintf(wcsprintf_file, format, arg_list);

  } else {
    /* Output to buffer. */
    used = wcsprintf_bufp - wcsprintf_buff;
    if (wcsprintf_size - used < 128) {
      /* Expand the buffer. */
      wcsprintf_size += 1024;
      realloc_buff = realloc(wcsprintf_buff, wcsprintf_size);
      if (realloc_buff == NULL) {
        free(wcsprintf_buff);
        wcsprintf_buff = 0x0;
        return 1;
      }
      wcsprintf_buff = realloc_buff;
      wcsprintf_bufp = wcsprintf_buff + used;
    }

    nbytes = vsprintf(wcsprintf_bufp, format, arg_list);
    wcsprintf_bufp += nbytes;
  }

  va_end(arg_list);

  return nbytes;
}

/*--------------------------------------------------------------------------*/

int wcsfprintf(FILE *stream, const char *format, ...)

{
  char *realloc_buff;
  int  nbytes;
  size_t  used;
  va_list arg_list;

  if (wcsprintf_buff == 0x0 && wcsprintf_file == 0x0) {
    /* Send output to stream if wcsprintf_set() hasn't been called. */
    wcsprintf_file = stream;
  }

  va_start(arg_list, format);

  if (wcsprintf_file) {
    /* Output to file. */
    nbytes = vfprintf(wcsprintf_file, format, arg_list);

  } else {
    /* Output to buffer. */
    used = wcsprintf_bufp - wcsprintf_buff;
    if (wcsprintf_size - used < 128) {
      /* Expand the buffer. */
      wcsprintf_size += 1024;
      realloc_buff = realloc(wcsprintf_buff, wcsprintf_size);
      if (realloc_buff == NULL) {
        free(wcsprintf_buff);
        wcsprintf_buff = 0x0;
        return 1;
      }
      wcsprintf_buff = realloc_buff;
      wcsprintf_bufp = wcsprintf_buff + used;
    }

    nbytes = vsprintf(wcsprintf_bufp, format, arg_list);
    wcsprintf_bufp += nbytes;
  }

  va_end(arg_list);

  return nbytes;
}
