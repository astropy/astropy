/*============================================================================

  WCSLIB 5.7 - an implementation of the FITS WCS standard.
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
  $Id: getwcstab.c,v 5.7 2015/06/29 02:44:16 mcalabre Exp $
*===========================================================================*/

#include <stdlib.h>
#include <string.h>

#include <fitsio.h>

#include "getwcstab.h"

/*--------------------------------------------------------------------------*/

int fits_read_wcstab(
  fitsfile   *fptr,
  int  nwtb,
  wtbarr *wtb,
  int  *status)

{
  int  anynul, colnum, hdunum, iwtb, m, naxis, nostat;
  long *naxes = 0, nelem;
  wtbarr *wtbp;


  if (*status) return *status;

  if (fptr == 0) {
    return (*status = NULL_INPUT_PTR);
  }

  if (nwtb == 0) return 0;

  /* Zero the array pointers. */
  wtbp = wtb;
  for (iwtb = 0; iwtb < nwtb; iwtb++, wtbp++) {
    *wtbp->arrayp = 0x0;
  }

  /* Save HDU number so that we can move back to it later. */
  fits_get_hdu_num(fptr, &hdunum);

  wtbp = wtb;
  for (iwtb = 0; iwtb < nwtb; iwtb++, wtbp++) {
    /* Move to the required binary table extension. */
    if (fits_movnam_hdu(fptr, BINARY_TBL, (char *)(wtbp->extnam),
        wtbp->extver, status)) {
      goto cleanup;
    }

    /* Locate the table column. */
    if (fits_get_colnum(fptr, CASEINSEN, (char *)(wtbp->ttype), &colnum,
        status)) {
      goto cleanup;
    }

    /* Get the array dimensions and check for consistency. */
    if (wtbp->ndim < 1) {
      *status = NEG_AXIS;
      goto cleanup;
    }

    if (!(naxes = calloc(wtbp->ndim, sizeof(long)))) {
      *status = MEMORY_ALLOCATION;
      goto cleanup;
    }

    if (fits_read_tdim(fptr, colnum, wtbp->ndim, &naxis, naxes, status)) {
      goto cleanup;
    }

    if (naxis != wtbp->ndim) {
      if (wtbp->kind == 'c' && wtbp->ndim == 2) {
        /* Allow TDIMn to be omitted for degenerate coordinate arrays. */
        naxis = 2;
        naxes[1] = naxes[0];
        naxes[0] = 1;
      } else {
        *status = BAD_TDIM;
        goto cleanup;
      }
    }

    if (wtbp->kind == 'c') {
      /* Coordinate array; calculate the array size. */
      nelem = naxes[0];
      for (m = 0; m < naxis-1; m++) {
        *(wtbp->dimlen + m) = naxes[m+1];
        nelem *= naxes[m+1];
      }
    } else {
      /* Index vector; check length. */
      if ((nelem = naxes[0]) != *(wtbp->dimlen)) {
        /* N.B. coordinate array precedes the index vectors. */
        *status = BAD_TDIM;
        goto cleanup;
      }
    }

    free(naxes);
    naxes = 0;

    /* Allocate memory for the array. */
    if (!(*wtbp->arrayp = calloc((size_t)nelem, sizeof(double)))) {
      *status = MEMORY_ALLOCATION;
      goto cleanup;
    }

    /* Read the array from the table. */
    if (fits_read_col_dbl(fptr, colnum, wtbp->row, 1L, nelem, 0.0,
        *wtbp->arrayp, &anynul, status)) {
      goto cleanup;
    }
  }

cleanup:
  /* Move back to the starting HDU. */
  nostat = 0;
  fits_movabs_hdu(fptr, hdunum, 0, &nostat);

  /* Release allocated memory. */
  if (naxes) free(naxes);
  if (*status) {
    wtbp = wtb;
    for (iwtb = 0; iwtb < nwtb; iwtb++, wtbp++) {
      if (*wtbp->arrayp) free(*wtbp->arrayp);
    }
  }

  return *status;
}
