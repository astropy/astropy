/*============================================================================
  WCSLIB 8.3 - an implementation of the FITS WCS standard.
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
  $Id: lin.c,v 8.3 2024/05/13 16:33:00 mcalabre Exp $
*===========================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wcserr.h"
#include "wcsprintf.h"
#include "lin.h"
#include "dis.h"

// Map status return value to message.
const char *lin_errmsg[] = {
  "Success",
  "Null linprm pointer passed",
  "Memory allocation failed",
  "PCi_ja matrix is singular",
  "Failed to initialize distortion functions",
  "Distort error",
  "De-distort error"};

// Map error returns for lower-level routines.
const int lin_diserr[] = {
  LINERR_SUCCESS,		//  0: DISERR_SUCCESS
  LINERR_NULL_POINTER,		//  1: DISERR_NULL_POINTER
  LINERR_MEMORY,		//  2: DISERR_MEMORY
  LINERR_DISTORT_INIT,		//  3: DISERR_BAD_PARAM
  LINERR_DISTORT,		//  4: DISERR_DISTORT
  LINERR_DEDISTORT		//  5: DISERR_DEDISTORT
};

static const int LINSET = 137;

// Convenience macro for invoking wcserr_set().
#define LIN_ERRMSG(status) WCSERR_SET(status), lin_errmsg[status]

//----------------------------------------------------------------------------

int linini(int alloc, int naxis, struct linprm *lin)

{
  return lininit(alloc, naxis, lin, -1);
}

//----------------------------------------------------------------------------

int lininit(int alloc, int naxis, struct linprm *lin, int ndpmax)

{
  static const char *function = "lininit";

  if (lin == 0x0) return LINERR_NULL_POINTER;

  // Initialize error message handling.
  if (lin->flag == -1) {
    lin->err = 0x0;
  }
  struct wcserr **err = &(lin->err);
  wcserr_clear(err);


  // Initialize memory management.
  if (lin->flag == -1 || lin->m_flag != LINSET) {
    if (lin->flag == -1) {
      lin->dispre = 0x0;
      lin->disseq = 0x0;
      lin->tmpcrd = 0x0;
    }

    lin->m_flag   = 0;
    lin->m_naxis  = 0;
    lin->m_crpix  = 0x0;
    lin->m_pc     = 0x0;
    lin->m_cdelt  = 0x0;
    lin->m_dispre = 0x0;
    lin->m_disseq = 0x0;
  }

  if (naxis < 0) {
    return wcserr_set(WCSERR_SET(LINERR_MEMORY),
      "naxis must not be negative (got %d)", naxis);
  }


  // Allocate memory for arrays if required.
  if (alloc ||
      lin->crpix  == 0x0 ||
      lin->pc     == 0x0 ||
      lin->cdelt  == 0x0) {

    // Was sufficient allocated previously?
    if (lin->m_flag == LINSET && lin->m_naxis < naxis) {
      // No, free it.
      linfree(lin);
    }

    if (alloc || lin->crpix == 0x0) {
      if (lin->m_crpix) {
        // In case the caller fiddled with it.
        lin->crpix = lin->m_crpix;

      } else {
        if ((lin->crpix = calloc(naxis, sizeof(double))) == 0x0) {
          return wcserr_set(LIN_ERRMSG(LINERR_MEMORY));
        }

        lin->m_flag  = LINSET;
        lin->m_naxis = naxis;
        lin->m_crpix = lin->crpix;
      }
    }

    if (alloc || lin->pc == 0x0) {
      if (lin->m_pc) {
        // In case the caller fiddled with it.
        lin->pc = lin->m_pc;

      } else {
        if ((lin->pc = calloc(naxis*naxis, sizeof(double))) == 0x0) {
          linfree(lin);
          return wcserr_set(LIN_ERRMSG(LINERR_MEMORY));
        }

        lin->m_flag  = LINSET;
        lin->m_naxis = naxis;
        lin->m_pc    = lin->pc;
      }
    }

    if (alloc || lin->cdelt == 0x0) {
      if (lin->m_cdelt) {
        // In case the caller fiddled with it.
        lin->cdelt = lin->m_cdelt;

      } else {
        if ((lin->cdelt = calloc(naxis, sizeof(double))) == 0x0) {
          linfree(lin);
          return wcserr_set(LIN_ERRMSG(LINERR_MEMORY));
        }

        lin->m_flag  = LINSET;
        lin->m_naxis = naxis;
        lin->m_cdelt = lin->cdelt;
      }
    }
  }


  // Reinitialize disprm structs if we are managing them.
  if (lin->m_dispre) {
    disinit(1, naxis, lin->dispre, ndpmax);
  }

  if (lin->m_disseq) {
    disinit(1, naxis, lin->disseq, ndpmax);
  }


  // Free memory allocated by linset().
  if (abs(lin->flag) == LINSET) {
    if (lin->piximg) free(lin->piximg);
    if (lin->imgpix) free(lin->imgpix);
    if (lin->tmpcrd) free(lin->tmpcrd);
  }

  lin->piximg  = 0x0;
  lin->imgpix  = 0x0;
  lin->i_naxis = 0;
  lin->unity   = 0;
  lin->affine  = 0;
  lin->simple  = 0;
  lin->tmpcrd  = 0x0;


  lin->naxis = naxis;

  // CRPIXja defaults to 0.0.
  for (int j = 0; j < naxis; j++) {
    lin->crpix[j] = 0.0;
  }

  // PCi_ja defaults to the unit matrix.
  double *pc = lin->pc;
  for (int i = 0; i < naxis; i++) {
    for (int j = 0; j < naxis; j++) {
      if (j == i) {
        *pc = 1.0;
      } else {
        *pc = 0.0;
      }
      pc++;
    }
  }

  // CDELTia defaults to 1.0.
  for (int i = 0; i < naxis; i++) {
    lin->cdelt[i] = 1.0;
  }

  lin->flag = 0;

  return 0;
}

//----------------------------------------------------------------------------

int lindis(int sequence, struct linprm *lin, struct disprm *dis)

{
  return lindist(sequence, lin, dis, -1);
}

//----------------------------------------------------------------------------

int lindist(int sequence, struct linprm *lin, struct disprm *dis, int ndpmax)

{
  static const char *function = "lindist";

  if (lin == 0x0) return LINERR_NULL_POINTER;
  struct wcserr **err = &(lin->err);

  if (sequence == 1) {
    if (lin->m_dispre) {
      disfree(lin->m_dispre);
      free(lin->m_dispre);
    }

    lin->dispre   = dis;
    lin->m_flag   = LINSET;
    lin->m_dispre = dis;

  } else if (sequence == 2) {
    if (lin->m_disseq) {
      disfree(lin->m_disseq);
      free(lin->m_disseq);
    }

    lin->disseq   = dis;
    lin->m_flag   = LINSET;
    lin->m_disseq = dis;

  } else {
    return wcserr_set(WCSERR_SET(LINERR_DISTORT_INIT),
      "Invalid sequence (%d)", sequence);
  }

  if (dis) {
    int status = disinit(1, lin->naxis, dis, ndpmax);
    if (status) {
      return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
    }
  }

  return 0;
}

//----------------------------------------------------------------------------

int lincpy(int alloc, const struct linprm *linsrc, struct linprm *lindst)

{
  static const char *function = "lincpy";

  if (linsrc == 0x0) return LINERR_NULL_POINTER;
  if (lindst == 0x0) return LINERR_NULL_POINTER;
  struct wcserr **err = &(lindst->err);

  int naxis = linsrc->naxis;
  if (naxis < 1) {
    return wcserr_set(WCSERR_SET(LINERR_MEMORY),
      "naxis must be positive (got %d)", naxis);
  }

  int status = lininit(alloc, naxis, lindst, 0);
  if (status) {
    return status;
  }

  const double *srcp = linsrc->crpix;
  double *dstp = lindst->crpix;
  for (int j = 0; j < naxis; j++) {
    *(dstp++) = *(srcp++);
  }

  srcp = linsrc->pc;
  dstp = lindst->pc;
  for (int i = 0; i < naxis; i++) {
    for (int j = 0; j < naxis; j++) {
      *(dstp++) = *(srcp++);
    }
  }

  srcp = linsrc->cdelt;
  dstp = lindst->cdelt;
  for (int i = 0; i < naxis; i++) {
    *(dstp++) = *(srcp++);
  }

  if (linsrc->dispre) {
    if (!lindst->dispre) {
      if ((lindst->dispre = calloc(1, sizeof(struct disprm))) == 0x0) {
        return wcserr_set(LIN_ERRMSG(LINERR_MEMORY));
      }

      lindst->m_dispre = lindst->dispre;
    }

    if ((status = discpy(alloc, linsrc->dispre, lindst->dispre))) {
      status = wcserr_set(LIN_ERRMSG(lin_diserr[status]));
      goto cleanup;
    }
  }

  if (linsrc->disseq) {
    if (!lindst->disseq) {
      if ((lindst->disseq = calloc(1, sizeof(struct disprm))) == 0x0) {
        return wcserr_set(LIN_ERRMSG(LINERR_MEMORY));
      }

      lindst->m_disseq = lindst->disseq;
    }

    if ((status = discpy(alloc, linsrc->disseq, lindst->disseq))) {
      status = wcserr_set(LIN_ERRMSG(lin_diserr[status]));
      goto cleanup;
    }
  }

cleanup:
  if (status) {
    if (lindst->m_dispre) {
      disfree(lindst->m_dispre);
      free(lindst->m_dispre);
      lindst->m_dispre = 0x0;
      lindst->dispre = 0x0;
    }

    if (lindst->m_disseq) {
      disfree(lindst->m_disseq);
      free(lindst->m_disseq);
      lindst->m_disseq = 0x0;
      lindst->disseq = 0x0;
    }
  }

  return status;
}

//----------------------------------------------------------------------------

int linfree(struct linprm *lin)

{
  if (lin == 0x0) return LINERR_NULL_POINTER;

  if (lin->flag != -1) {
    // Optionally allocated by lininit() for given parameters.
    if (lin->m_flag == LINSET) {
      if (lin->crpix  == lin->m_crpix)  lin->crpix  = 0x0;
      if (lin->pc     == lin->m_pc)     lin->pc     = 0x0;
      if (lin->cdelt  == lin->m_cdelt)  lin->cdelt  = 0x0;
      if (lin->dispre == lin->m_dispre) lin->dispre = 0x0;
      if (lin->disseq == lin->m_disseq) lin->disseq = 0x0;

      if (lin->m_crpix)  free(lin->m_crpix);
      if (lin->m_pc)     free(lin->m_pc);
      if (lin->m_cdelt)  free(lin->m_cdelt);

      if (lin->m_dispre) {
        disfree(lin->m_dispre);
        free(lin->m_dispre);
      }

      if (lin->m_disseq) {
        disfree(lin->m_disseq);
        free(lin->m_disseq);
      }
    }

    // Allocated unconditionally by linset().
    if (lin->piximg) free(lin->piximg);
    if (lin->imgpix) free(lin->imgpix);
    if (lin->tmpcrd) free(lin->tmpcrd);
  }


  lin->m_flag   = 0;
  lin->m_naxis  = 0;
  lin->m_crpix  = 0x0;
  lin->m_pc     = 0x0;
  lin->m_cdelt  = 0x0;
  lin->m_dispre = 0x0;
  lin->m_disseq = 0x0;

  lin->piximg   = 0x0;
  lin->imgpix   = 0x0;
  lin->i_naxis  = 0;

  lin->tmpcrd   = 0x0;

  wcserr_clear(&(lin->err));

  lin->flag = 0;

  return 0;
}

//----------------------------------------------------------------------------

int linsize(const struct linprm *lin, int sizes[2])

{
  if (lin == 0x0) {
    sizes[0] = sizes[1] = 0;
    return 0;
  }

  // Base size, in bytes.
  sizes[0] = sizeof(struct linprm);

  // Total size of allocated memory, in bytes.
  sizes[1] = 0;

  int naxis = lin->naxis;

  // linprm::crpix[].
  sizes[1] += naxis * sizeof(double);

  // linprm::pc[].
  sizes[1] += naxis*naxis * sizeof(double);

  // linprm::cdelt[].
  sizes[1] += naxis * sizeof(double);

  // linprm::dispre[].
  int exsizes[2];
  dissize(lin->dispre, exsizes);
  sizes[1] += exsizes[0] + exsizes[1];

  // linprm::disseq[].
  dissize(lin->disseq, exsizes);
  sizes[1] += exsizes[0] + exsizes[1];

  // linprm::err[].
  wcserr_size(lin->err, exsizes);
  sizes[1] += exsizes[0] + exsizes[1];

  // The remaining arrays are allocated unconditionally by linset().
  if (abs(lin->flag) != LINSET) {
    return 0;
  }

  // linprm::piximg[].
  sizes[1] += naxis*naxis * sizeof(double);

  // linprm::imgpix[].
  sizes[1] += naxis*naxis * sizeof(double);

  // linprm::tmpcrd[].
  sizes[1] += naxis * sizeof(double);

  return 0;
}

//----------------------------------------------------------------------------

int linenq(const struct linprm *lin, int enquiry)

{
  // Initialize.
  if (lin == 0x0) return LINERR_NULL_POINTER;

  int answer = 0;

  if (enquiry & LINENQ_MEM) {
    if (lin->m_flag != LINSET) return 0;
    answer = 1;
  }

  if (enquiry & LINENQ_SET) {
    if (abs(lin->flag) != LINSET) return 0;
    answer = 1;
  }

  if (enquiry & LINENQ_BYP) {
    if (lin->flag != 1 && lin->flag != -LINSET) return 0;
    answer = 1;
  }

  return answer;
}

//----------------------------------------------------------------------------

int linprt(const struct linprm *lin)

{
  if (lin == 0x0) return LINERR_NULL_POINTER;

  if (abs(lin->flag) != LINSET) {
    wcsprintf("The linprm struct is UNINITIALIZED.\n");
    return 0;
  }

  // Parameters supplied.
  wcsprintf("       flag: %d\n", lin->flag);
  wcsprintf("      naxis: %d\n", lin->naxis);

  WCSPRINTF_PTR("      crpix: ", lin->crpix, "\n");
  wcsprintf("            ");
  for (int j = 0; j < lin->naxis; j++) {
    wcsprintf("  %#- 11.5g", lin->crpix[j]);
  }
  wcsprintf("\n");

  int k = 0;
  WCSPRINTF_PTR("         pc: ", lin->pc, "\n");
  for (int i = 0; i < lin->naxis; i++) {
    wcsprintf("    pc[%d][]:", i);
    for (int j = 0; j < lin->naxis; j++) {
      wcsprintf("  %#- 11.5g", lin->pc[k++]);
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("      cdelt: ", lin->cdelt, "\n");
  wcsprintf("            ");
  for (int i = 0; i < lin->naxis; i++) {
    wcsprintf("  %#- 11.5g", lin->cdelt[i]);
  }
  wcsprintf("\n");

  WCSPRINTF_PTR("     dispre: ", lin->dispre, "");
  if (lin->dispre != 0x0) wcsprintf("  (see below)");
  wcsprintf("\n");
  WCSPRINTF_PTR("     disseq: ", lin->disseq, "");
  if (lin->disseq != 0x0) wcsprintf("  (see below)");
  wcsprintf("\n");

  // Derived values.
  if (lin->piximg == 0x0) {
    wcsprintf("     piximg: (nil)\n");
  } else {
    int k = 0;
    for (int i = 0; i < lin->naxis; i++) {
      wcsprintf("piximg[%d][]:", i);
      for (int j = 0; j < lin->naxis; j++) {
        wcsprintf("  %#- 11.5g", lin->piximg[k++]);
      }
      wcsprintf("\n");
    }
  }

  if (lin->imgpix == 0x0) {
    wcsprintf("     imgpix: (nil)\n");
  } else {
    int k = 0;
    for (int i = 0; i < lin->naxis; i++) {
      wcsprintf("imgpix[%d][]:", i);
      for (int j = 0; j < lin->naxis; j++) {
        wcsprintf("  %#- 11.5g", lin->imgpix[k++]);
      }
      wcsprintf("\n");
    }
  }

  wcsprintf("    i_naxis: %d\n", lin->i_naxis);
  wcsprintf("      unity: %d\n", lin->unity);
  wcsprintf("     affine: %d\n", lin->affine);
  wcsprintf("     simple: %d\n", lin->simple);

  // Error handling.
  WCSPRINTF_PTR("        err: ", lin->err, "\n");
  if (lin->err) {
    wcserr_prt(lin->err, "             ");
  }

  // Work arrays.
  WCSPRINTF_PTR("     tmpcrd: ", lin->tmpcrd, "\n");

  // Memory management.
  wcsprintf("     m_flag: %d\n", lin->m_flag);
  wcsprintf("    m_naxis: %d\n", lin->m_naxis);
  WCSPRINTF_PTR("    m_crpix: ", lin->m_crpix, "");
  if (lin->m_crpix == lin->crpix) wcsprintf("  (= crpix)");
  wcsprintf("\n");
  WCSPRINTF_PTR("       m_pc: ", lin->m_pc, "");
  if (lin->m_pc == lin->pc) wcsprintf("  (= pc)");
  wcsprintf("\n");
  WCSPRINTF_PTR("    m_cdelt: ", lin->m_cdelt, "");
  if (lin->m_cdelt == lin->cdelt) wcsprintf("  (= cdelt)");
  wcsprintf("\n");
  WCSPRINTF_PTR("   m_dispre: ", lin->m_dispre, "");
  if (lin->dispre && lin->m_dispre == lin->dispre) wcsprintf("  (= dispre)");
  wcsprintf("\n");
  WCSPRINTF_PTR("   m_disseq: ", lin->m_disseq, "");
  if (lin->disseq && lin->m_disseq == lin->disseq) wcsprintf("  (= disseq)");
  wcsprintf("\n");

  // Distortion parameters (from above).
  if (lin->dispre) {
    wcsprintf("\n");
    wcsprintf("dispre.*\n");
    disprt(lin->dispre);
  }

  if (lin->disseq) {
    wcsprintf("\n");
    wcsprintf("disseq.*\n");
    disprt(lin->disseq);
  }

  return 0;
}

//----------------------------------------------------------------------------

int linperr(const struct linprm *lin, const char *prefix)

{
  if (lin == 0x0) return LINERR_NULL_POINTER;

  if (lin->err && wcserr_prt(lin->err, prefix) == 0) {
    if (lin->dispre) wcserr_prt(lin->dispre->err, prefix);
    if (lin->disseq) wcserr_prt(lin->disseq->err, prefix);
  }

  return 0;
}

//----------------------------------------------------------------------------

int linset(struct linprm *lin)

{
  static const char *function = "linset";

  if (lin == 0x0) return LINERR_NULL_POINTER;
  if (lin->flag == -LINSET) return 0;
  struct wcserr **err = &(lin->err);

  int naxis = lin->naxis;

  // Check for a unit matrix.
  lin->unity = 1;
  double *pc = lin->pc;
  for (int i = 0; i < naxis; i++) {
    for (int j = 0; j < naxis; j++) {
      if (j == i) {
        if (*(pc++) != 1.0) {
          lin->unity = 0;
          break;
        }
      } else {
        if (*(pc++) != 0.0) {
          lin->unity = 0;
          break;
        }
      }
    }
  }


  if (lin->unity) {
    if (abs(lin->flag) == LINSET) {
      // Free memory that may have been allocated previously.
      if (lin->piximg) free(lin->piximg);
      if (lin->imgpix) free(lin->imgpix);
    }

    lin->piximg  = 0x0;
    lin->imgpix  = 0x0;
    lin->i_naxis = 0;

    // Check cdelt.
    for (int i = 0; i < naxis; i++) {
      if (lin->cdelt[i] == 0.0) {
        return wcserr_set(LIN_ERRMSG(LINERR_SINGULAR_MTX));
      }
    }

  } else {
    if (abs(lin->flag) != LINSET || lin->i_naxis < naxis) {
      if (abs(lin->flag) == LINSET) {
        // Free memory that may have been allocated previously.
        if (lin->piximg) free(lin->piximg);
        if (lin->imgpix) free(lin->imgpix);
      }

      // Allocate memory for internal arrays.
      if ((lin->piximg = calloc(naxis*naxis, sizeof(double))) == 0x0) {
        return wcserr_set(LIN_ERRMSG(LINERR_MEMORY));
      }

      if ((lin->imgpix = calloc(naxis*naxis, sizeof(double))) == 0x0) {
        free(lin->piximg);
        return wcserr_set(LIN_ERRMSG(LINERR_MEMORY));
      }

      lin->i_naxis = naxis;
    }

    // Compute the pixel-to-image transformation matrix.
    pc = lin->pc;
    double *piximg = lin->piximg;
    for (int i = 0; i < naxis; i++) {
      if (lin->disseq == 0x0) {
        // No sequent distortions.  Incorporate cdelt into piximg.
        for (int j = 0; j < naxis; j++) {
          *(piximg++) = lin->cdelt[i] * (*(pc++));
        }
      } else {
        for (int j = 0; j < naxis; j++) {
          *(piximg++) = *(pc++);
        }
      }
    }

    // Compute the image-to-pixel transformation matrix.
    int status = matinv(naxis, lin->piximg, lin->imgpix);
    if (status) {
      return wcserr_set(LIN_ERRMSG(status));
    }
  }


  // Set up the distortion functions.
  lin->affine = 1;
  if (lin->dispre) {
    (lin->dispre)->flag = 0;
    int status = disset(lin->dispre);
    if (status) {
      return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
    }

    lin->affine = 0;
  }

  if (lin->disseq) {
    (lin->disseq)->flag = 0;
    int status = disset(lin->disseq);
    if (status) {
      return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
    }

    lin->affine = 0;
  }

  lin->simple = lin->unity && lin->affine;


  // Create work arrays.
  if (lin->tmpcrd) free(lin->tmpcrd);
  if ((lin->tmpcrd = calloc(naxis, sizeof(double))) == 0x0) {
    linfree(lin);
    return wcserr_set(LIN_ERRMSG(LINERR_MEMORY));
  }


  lin->flag = (lin->flag == 1) ? -LINSET : LINSET;

  return 0;
}

//----------------------------------------------------------------------------

int linp2x(
  struct linprm *lin,
  int ncoord,
  int nelem,
  const double pixcrd[],
  double imgcrd[])

{
  static const char *function = "linp2x";

  // Initialize.
  if (lin == 0x0) return LINERR_NULL_POINTER;
  struct wcserr **err = &(lin->err);

  if (abs(lin->flag) != LINSET) {
    int status = linset(lin);
    if (status) {
      return status;
    }
  }

  int naxis = lin->naxis;


  // Convert pixel coordinates to intermediate world coordinates.
  const double *pix = pixcrd;
  double *img = imgcrd;

  if (lin->simple) {
    // Handle the simplest and most common case with maximum efficiency.
    int nelemn = nelem - naxis;
    for (int k = 0; k < ncoord; k++) {
      for (int i = 0; i < naxis; i++) {
        *(img++) = lin->cdelt[i] * (*(pix++) - lin->crpix[i]);
      }

      pix += nelemn;
      img += nelemn;
    }

  } else if (lin->affine) {
    // No distortions.
    int ndbl   = naxis * sizeof(double);
    int nelemn = nelem - naxis;
    for (int k = 0; k < ncoord; k++) {
      memset(img, 0, ndbl);

      for (int j = 0; j < naxis; j++) {
        // cdelt will have been incorporated into piximg.
        double *piximg = lin->piximg + j;

        // Column-wise multiplication allows this to be cached.
        double temp = *(pix++) - lin->crpix[j];
        for (int i = 0; i < naxis; i++, piximg += naxis) {
          img[i] += *piximg * temp;
        }
      }

      pix += nelemn;
      img += nelem;
    }

  } else {
    // Distortions are present.
    int ndbl = naxis * sizeof(double);
    double *tmp  = lin->tmpcrd;

    for (int k = 0; k < ncoord; k++) {
      if (lin->dispre) {
        int status = disp2x(lin->dispre, pix, tmp);
        if (status) {
          return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
        }
      } else {
        memcpy(tmp, pix, ndbl);
      }

      if (lin->unity) {
        for (int i = 0; i < naxis; i++) {
          img[i] = tmp[i] - lin->crpix[i];
        }

      } else {
        for (int j = 0; j < naxis; j++) {
          tmp[j] -= lin->crpix[j];
        }

        double *piximg = lin->piximg;
        for (int i = 0; i < naxis; i++) {
          img[i] = 0.0;
          for (int j = 0; j < naxis; j++) {
            img[i] += *(piximg++) * tmp[j];
          }
        }
      }

      if (lin->disseq) {
        int status = disp2x(lin->disseq, img, tmp);
        if (status) {
          return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
        }

        // With sequent distortions, cdelt is not incorporated into piximg...
        for (int i = 0; i < naxis; i++) {
          img[i] = lin->cdelt[i] * tmp[i];
        }

      } else if (lin->unity) {
        // ...nor if the matrix is unity.
        for (int i = 0; i < naxis; i++) {
          img[i] *= lin->cdelt[i];
        }
      }

      pix += nelem;
      img += nelem;
    }
  }

  return 0;
}

//----------------------------------------------------------------------------

int linx2p(
  struct linprm *lin,
  int ncoord,
  int nelem,
  const double imgcrd[],
  double pixcrd[])

{
  static const char *function = "linx2p";

  // Initialize.
  if (lin == 0x0) return LINERR_NULL_POINTER;
  struct wcserr **err = &(lin->err);

  if (abs(lin->flag) != LINSET) {
    int status = linset(lin);
    if (status) {
      return status;
    }
  }

  int naxis = lin->naxis;


  // Convert intermediate world coordinates to pixel coordinates.
  const double *img = imgcrd;
  double *pix = pixcrd;

  if (lin->simple) {
    // Handle the simplest and most common case with maximum efficiency.
    int nelemn = nelem - naxis;
    for (int k = 0; k < ncoord; k++) {
      for (int j = 0; j < naxis; j++) {
        *(pix++) = (*(img++) / lin->cdelt[j]) + lin->crpix[j];
      }

      img += nelemn;
      pix += nelemn;
    }

  } else if (lin->affine) {
    // No distortions.
    int nelemn = nelem - naxis;
    for (int k = 0; k < ncoord; k++) {
      // cdelt will have been incorporated into imgpix.
      double *imgpix = lin->imgpix;

      for (int j = 0; j < naxis; j++) {
        *pix = 0.0;
        for (int i = 0; i < naxis; i++) {
          *pix += *imgpix * img[i];
          imgpix++;
        }

        *(pix++) += lin->crpix[j];
      }

      img += nelem;
      pix += nelemn;
    }

  } else {
    // Distortions are present.
    int ndbl = naxis * sizeof(double);
    double *tmp  = lin->tmpcrd;

    for (int k = 0; k < ncoord; k++) {
      if (lin->disseq) {
        // With sequent distortions, cdelt is not incorporated into imgpix...
        for (int i = 0; i < naxis; i++) {
          tmp[i] = img[i] / lin->cdelt[i];
        }

        int status = disx2p(lin->disseq, tmp, pix);
        if (status) {
          return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
        }

        memcpy(tmp, pix, ndbl);

      } else if (lin->unity) {
        // ...nor if the matrix is unity.
        for (int i = 0; i < naxis; i++) {
          tmp[i] = img[i] / lin->cdelt[i];
        }

      } else {
        // cdelt will have been incorporated into imgpix.
        memcpy(tmp, img, ndbl);
      }

      if (lin->unity) {
        for (int j = 0; j < naxis; j++) {
          pix[j] = tmp[j] + lin->crpix[j];
        }

      } else {
        double *imgpix = lin->imgpix;
        for (int j = 0; j < naxis; j++) {
          pix[j] = lin->crpix[j];
          for (int i = 0; i < naxis; i++) {
            pix[j] += *(imgpix++) * tmp[i];
          }
        }
      }

      if (lin->dispre) {
        memcpy(tmp, pix, ndbl);

        int status = disx2p(lin->dispre, tmp, pix);
        if (status) {
          return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
        }
      }

      img += nelem;
      pix += nelem;
    }
  }

  return 0;
}

//----------------------------------------------------------------------------

int linwarp(
  struct linprm *lin,
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
  static const char *function = "linwarp";

  // Initialize.
  if (lin == 0x0) return LINERR_NULL_POINTER;
  struct wcserr **err = &(lin->err);

  int naxis = lin->naxis;

  if (nsamp) *nsamp = 0;
  for (int j = 0; j < naxis; j++) {
    if (maxdis) maxdis[j] = 0.0;
    if (avgdis) avgdis[j] = 0.0;
    if (rmsdis) rmsdis[j] = 0.0;
  }
  if (maxtot) *maxtot = 0.0;
  if (avgtot) *avgtot = 0.0;
  if (rmstot) *rmstot = 0.0;

  // Quick return if no distortions.
  if (lin->affine) return 0;

  // It's easier if there are no sequent distortions!
  if (lin->disseq == 0x0) {
    int status = diswarp(lin->dispre, pixblc, pixtrc, pixsamp, nsamp,
                         maxdis, maxtot, avgdis, avgtot, rmsdis, rmstot);
    return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
  }

  // Make a reference copy of lin without distortions.
  struct linprm affine;
  affine.flag = -1;

  int status = lincpy(1, lin, &affine) ||
               lindist(1, &affine, 0x0, 0) ||
               lindist(2, &affine, 0x0, 0) ||
               linset(&affine);
  if (status) {
    return wcserr_set(LIN_ERRMSG(status));
  }

  // Work out increments on each axis.
  int ncoord = 0;
  for (int j = 0; j < naxis; j++) {
    double *pixinc = lin->tmpcrd;
    double pixspan = pixtrc[j] - (pixblc ? pixblc[j] : 1.0);

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

    if (j == 0) {
      // Number of samples on axis 1.
      ncoord = 1 + (int)((pixspan/pixinc[0]) + 0.5);
    }
  }

  // Allocate memory in bulk for processing the image row by row.
  double *pix0 = calloc((3*ncoord+4)*naxis, sizeof(double));
  if (pix0 == 0x0) {
    return wcserr_set(LIN_ERRMSG(LINERR_MEMORY));
  }

  // Carve up the allocated memory.
  double *img    = pix0 + naxis*ncoord;
  double *pix1   = img  + naxis*ncoord;
  double *pixinc = pix1 + naxis*ncoord;
  double *pixend = pixinc + naxis;
  double *sumdis = pixend + naxis;
  double *ssqdis = sumdis + naxis;


  // Copy tmpcrd since linp2x() will overwrite it.
  memcpy(pixinc, lin->tmpcrd, naxis*sizeof(double));

  // Set up the array of pixel coordinates.
  for (int j = 0; j < naxis; j++) {
    pix0[j] = pixblc ? pixblc[j] : 1.0;
    pixend[j] = pixtrc[j] + 0.5*pixinc[j];
  }

  double *pix0p = pix0 + naxis;
  for (int i = 1; i < ncoord; i++) {
    *(pix0p++) = pix0[0] + i*pixinc[0];

    for (int j = 1; j < naxis; j++) {
      *(pix0p++) = pix0[j];
    }
  }

  // Initialize accumulators.
  for (int j = 0; j < naxis; j++) {
    sumdis[j] = 0.0;
    ssqdis[j] = 0.0;
  }
  double sumtot = 0.0;
  double ssqtot = 0.0;


  // Loop over N dimensions.
  int carry = 0;
  while (carry == 0) {
    int status;
    if ((status = linp2x(lin, ncoord, naxis, pix0, img))) {
      // (Preserve the error message set by linp2x().)
      goto cleanup;
    }

    if ((status = linx2p(&affine, ncoord, naxis, img, pix1))) {
      // (Preserve the error message set by linx2p().)
      goto cleanup;
    }

    // Accumulate statistics.
    double *pix0p = pix0;
    double *pix1p = pix1;
    for (int i = 0; i < ncoord; i++) {
      (*nsamp)++;

      double dssq = 0.0;
      for (int j = 0; j < naxis; j++) {
        double dpix = *(pix1p++) - *(pix0p++);
        double dpx2 = dpix*dpix;

        sumdis[j] += dpix;
        ssqdis[j] += dpx2;

        if (maxdis && (dpix = fabs(dpix)) > maxdis[j]) maxdis[j] = dpix;

        dssq += dpx2;
      }

      double totdis = sqrt(dssq);
      sumtot += totdis;
      ssqtot += totdis*totdis;

      if (maxtot && *maxtot < totdis) *maxtot = totdis;
    }

    // Next array of pixel coordinates.
    for (int j = 1; j < naxis; j++) {
      pix0[j] += pixinc[j];
      if ((carry = (pix0[j] > pixend[j]))) {
        pix0[j] = pixblc ? pixblc[j] : 1.0;
      }

      pix0p = pix0 + naxis + j;
      for (int i = 1; i < ncoord; i++) {
        *pix0p = pix0[j];
        pix0p += naxis;
      }

      if (carry == 0) break;
    }
  }


  // Compute the means and RMSs.
  for (int j = 0; j < naxis; j++) {
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
  linfree(&affine);
  free(pix0);

  return status;
}

//----------------------------------------------------------------------------

int matinv(int n, const double mat[], double inv[])

{
  // Allocate memory for internal arrays.
  int *mxl = calloc(n, sizeof(int));
  if (mxl == 0x0) {
    return LINERR_MEMORY;
  }

  int *lxm = calloc(n, sizeof(int));
  if (lxm == 0x0) {
    free(mxl);
    return LINERR_MEMORY;
  }

  double *rowmax = calloc(n, sizeof(double));
  if (rowmax == 0x0) {
    free(mxl);
    free(lxm);
    return LINERR_MEMORY;
  }

  double *lu = calloc(n*n, sizeof(double));
  if (lu == 0x0) {
    free(mxl);
    free(lxm);
    free(rowmax);
    return LINERR_MEMORY;
  }


  // Initialize arrays.
  int ij = 0;
  for (int i = 0; i < n; i++) {
    // Vector that records row interchanges.
    mxl[i] = i;

    rowmax[i] = 0.0;

    for (int j = 0; j < n; j++, ij++) {
      double dtemp = fabs(mat[ij]);
      if (dtemp > rowmax[i]) rowmax[i] = dtemp;

      lu[ij] = mat[ij];
    }

    // A row of zeroes indicates a singular matrix.
    if (rowmax[i] == 0.0) {
      free(mxl);
      free(lxm);
      free(rowmax);
      free(lu);
      return LINERR_SINGULAR_MTX;
    }
  }


  // Form the LU triangular factorization using scaled partial pivoting.
  for (int k = 0; k < n; k++) {
    // Decide whether to pivot.
    int pivot = k;
    double colmax = fabs(lu[k*n+k]) / rowmax[k];

    for (int i = k+1; i < n; i++) {
      int ik = i*n + k;
      double dtemp = fabs(lu[ik]) / rowmax[i];
      if (dtemp > colmax) {
        colmax = dtemp;
        pivot = i;
      }
    }

    if (pivot > k) {
      // We must pivot, interchange the rows of the design matrix.
      int kj = k*n;
      int pj = pivot*n;
      for (int j = 0; j < n; j++, pj++, kj++) {
        double dtemp = lu[pj];
        lu[pj] = lu[kj];
        lu[kj] = dtemp;
      }

      // Amend the vector of row maxima.
      double dtemp = rowmax[pivot];
      rowmax[pivot] = rowmax[k];
      rowmax[k] = dtemp;

      // Record the interchange for later use.
      int itemp = mxl[pivot];
      mxl[pivot] = mxl[k];
      mxl[k] = itemp;
    }

    // Gaussian elimination.
    for (int i = k+1; i < n; i++) {
      int ik = i*n + k;

      // Nothing to do if lu[ik] is zero.
      if (lu[ik] != 0.0) {
        // Save the scaling factor.
        lu[ik] /= lu[k*n+k];

        // Subtract rows.
        for (int j = k+1; j < n; j++) {
          lu[i*n+j] -= lu[ik]*lu[k*n+j];
        }
      }
    }
  }


  // mxl[i] records which row of mat corresponds to row i of lu.
  // lxm[i] records which row of lu  corresponds to row i of mat.
  for (int i = 0; i < n; i++) {
    lxm[mxl[i]] = i;
  }


  // Determine the inverse matrix.
  ij = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++, ij++) {
      inv[ij] = 0.0;
    }
  }

  for (int k = 0; k < n; k++) {
    inv[lxm[k]*n+k] = 1.0;

    // Forward substitution.
    for (int i = lxm[k]+1; i < n; i++) {
      for (int j = lxm[k]; j < i; j++) {
        inv[i*n+k] -= lu[i*n+j]*inv[j*n+k];
      }
    }

    // Backward substitution.
    for (int i = n-1; i >= 0; i--) {
      for (int j = i+1; j < n; j++) {
        inv[i*n+k] -= lu[i*n+j]*inv[j*n+k];
      }
      inv[i*n+k] /= lu[i*n+i];
    }
  }

  free(mxl);
  free(lxm);
  free(rowmax);
  free(lu);

  return 0;
}
