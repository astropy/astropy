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
  $Id: lin.c,v 7.3 2020/06/03 03:37:02 mcalabre Exp $
*===========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "wcserr.h"
#include "wcsprintf.h"
#include "lin.h"
#include "dis.h"

const int LINSET = 137;

/* Map status return value to message. */
const char *lin_errmsg[] = {
  "Success",
  "Null linprm pointer passed",
  "Memory allocation failed",
  "PCi_ja matrix is singular",
  "Failed to initialize distortion functions",
  "Distort error",
  "De-distort error"};

/* Map error returns for lower-level routines. */
const int lin_diserr[] = {
  LINERR_SUCCESS,		/*  0: DISERR_SUCCESS         */
  LINERR_NULL_POINTER,		/*  1: DISERR_NULL_POINTER    */
  LINERR_MEMORY,		/*  2: DISERR_MEMORY          */
  LINERR_DISTORT_INIT,		/*  3: DISERR_BAD_PARAM       */
  LINERR_DISTORT,		/*  4: DISERR_DISTORT         */
  LINERR_DEDISTORT		/*  5: DISERR_DEDISTORT       */
};

/* Convenience macro for invoking wcserr_set(). */
#define LIN_ERRMSG(status) WCSERR_SET(status), lin_errmsg[status]

/*--------------------------------------------------------------------------*/

int linini(int alloc, int naxis, struct linprm *lin)

{
  return lininit(alloc, naxis, lin, -1);
}

/*--------------------------------------------------------------------------*/

int lininit(int alloc, int naxis, struct linprm *lin, int ndpmax)

{
  static const char *function = "lininit";

  int i, j;
  double *pc;
  struct wcserr **err;

  if (lin == 0x0) return LINERR_NULL_POINTER;

  /* Initialize error message handling. */
  if (lin->flag == -1) {
    lin->err = 0x0;
  }
  err = &(lin->err);
  wcserr_clear(err);


  /* Initialize memory management. */
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


  /* Allocate memory for arrays if required. */
  if (alloc ||
      lin->crpix  == 0x0 ||
      lin->pc     == 0x0 ||
      lin->cdelt  == 0x0) {

    /* Was sufficient allocated previously? */
    if (lin->m_flag == LINSET && lin->m_naxis < naxis) {
      /* No, free it. */
      linfree(lin);
    }

    if (alloc || lin->crpix == 0x0) {
      if (lin->m_crpix) {
        /* In case the caller fiddled with it. */
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
        /* In case the caller fiddled with it. */
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
        /* In case the caller fiddled with it. */
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


  /* Reinitialize disprm structs if we are managing them. */
  if (lin->m_dispre) {
    disinit(1, naxis, lin->dispre, ndpmax);
  }

  if (lin->m_disseq) {
    disinit(1, naxis, lin->disseq, ndpmax);
  }


  /* Free memory allocated by linset(). */
  if (lin->flag == LINSET) {
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


  lin->flag  = 0;
  lin->naxis = naxis;


  /* CRPIXja defaults to 0.0. */
  for (j = 0; j < naxis; j++) {
    lin->crpix[j] = 0.0;
  }

  /* PCi_ja defaults to the unit matrix. */
  pc = lin->pc;
  for (i = 0; i < naxis; i++) {
    for (j = 0; j < naxis; j++) {
      if (j == i) {
        *pc = 1.0;
      } else {
        *pc = 0.0;
      }
      pc++;
    }
  }

  /* CDELTia defaults to 1.0. */
  for (i = 0; i < naxis; i++) {
    lin->cdelt[i] = 1.0;
  }


  return 0;
}

/*--------------------------------------------------------------------------*/

int lindis(int sequence, struct linprm *lin, struct disprm *dis)

{
  return lindist(sequence, lin, dis, -1);
}

/*--------------------------------------------------------------------------*/

int lindist(int sequence, struct linprm *lin, struct disprm *dis, int ndpmax)

{
  static const char *function = "lindist";

  int status;
  struct wcserr **err;

  if (lin == 0x0) return LINERR_NULL_POINTER;
  err = &(lin->err);

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
    if ((status = disinit(1, lin->naxis, dis, ndpmax))) {
      return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
    }
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int lincpy(int alloc, const struct linprm *linsrc, struct linprm *lindst)

{
  static const char *function = "lincpy";

  int i, j, naxis, status;
  const double *srcp;
  double *dstp;
  struct wcserr **err;

  if (linsrc == 0x0) return LINERR_NULL_POINTER;
  if (lindst == 0x0) return LINERR_NULL_POINTER;
  err = &(lindst->err);

  naxis = linsrc->naxis;
  if (naxis < 1) {
    return wcserr_set(WCSERR_SET(LINERR_MEMORY),
      "naxis must be positive (got %d)", naxis);
  }

  if ((status = lininit(alloc, naxis, lindst, 0))) {
    return status;
  }

  srcp = linsrc->crpix;
  dstp = lindst->crpix;
  for (j = 0; j < naxis; j++) {
    *(dstp++) = *(srcp++);
  }

  srcp = linsrc->pc;
  dstp = lindst->pc;
  for (i = 0; i < naxis; i++) {
    for (j = 0; j < naxis; j++) {
      *(dstp++) = *(srcp++);
    }
  }

  srcp = linsrc->cdelt;
  dstp = lindst->cdelt;
  for (i = 0; i < naxis; i++) {
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

/*--------------------------------------------------------------------------*/

int linfree(struct linprm *lin)

{
  if (lin == 0x0) return LINERR_NULL_POINTER;

  if (lin->flag != -1) {
    /* Optionally allocated by lininit() for given parameters. */
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

    /* Allocated unconditionally by linset(). */
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

/*--------------------------------------------------------------------------*/

int linprt(const struct linprm *lin)

{
  int i, j, k;

  if (lin == 0x0) return LINERR_NULL_POINTER;

  if (lin->flag != LINSET) {
    wcsprintf("The linprm struct is UNINITIALIZED.\n");
    return 0;
  }
  wcsprintf("       flag: %d\n", lin->flag);

  /* Parameters supplied. */
  wcsprintf("      naxis: %d\n", lin->naxis);

  WCSPRINTF_PTR("      crpix: ", lin->crpix, "\n");
  wcsprintf("            ");
  for (j = 0; j < lin->naxis; j++) {
    wcsprintf("  %#- 11.5g", lin->crpix[j]);
  }
  wcsprintf("\n");

  k = 0;
  WCSPRINTF_PTR("         pc: ", lin->pc, "\n");
  for (i = 0; i < lin->naxis; i++) {
    wcsprintf("    pc[%d][]:", i);
    for (j = 0; j < lin->naxis; j++) {
      wcsprintf("  %#- 11.5g", lin->pc[k++]);
    }
    wcsprintf("\n");
  }

  WCSPRINTF_PTR("      cdelt: ", lin->cdelt, "\n");
  wcsprintf("            ");
  for (i = 0; i < lin->naxis; i++) {
    wcsprintf("  %#- 11.5g", lin->cdelt[i]);
  }
  wcsprintf("\n");

  WCSPRINTF_PTR("     dispre: ", lin->dispre, "");
  if (lin->dispre != 0x0) wcsprintf("  (see below)");
  wcsprintf("\n");
  WCSPRINTF_PTR("     disseq: ", lin->disseq, "");
  if (lin->disseq != 0x0) wcsprintf("  (see below)");
  wcsprintf("\n");

  /* Derived values. */
  if (lin->piximg == 0x0) {
    wcsprintf("     piximg: (nil)\n");
  } else {
    k = 0;
    for (i = 0; i < lin->naxis; i++) {
      wcsprintf("piximg[%d][]:", i);
      for (j = 0; j < lin->naxis; j++) {
        wcsprintf("  %#- 11.5g", lin->piximg[k++]);
      }
      wcsprintf("\n");
    }
  }

  if (lin->imgpix == 0x0) {
    wcsprintf("     imgpix: (nil)\n");
  } else {
    k = 0;
    for (i = 0; i < lin->naxis; i++) {
      wcsprintf("imgpix[%d][]:", i);
      for (j = 0; j < lin->naxis; j++) {
        wcsprintf("  %#- 11.5g", lin->imgpix[k++]);
      }
      wcsprintf("\n");
    }
  }

  wcsprintf("    i_naxis: %d\n", lin->i_naxis);
  wcsprintf("      unity: %d\n", lin->unity);
  wcsprintf("     affine: %d\n", lin->affine);
  wcsprintf("     simple: %d\n", lin->simple);

  /* Error handling. */
  WCSPRINTF_PTR("        err: ", lin->err, "\n");
  if (lin->err) {
    wcserr_prt(lin->err, "             ");
  }

  /* Work arrays. */
  WCSPRINTF_PTR("     tmpcrd: ", lin->tmpcrd, "\n");

  /* Memory management. */
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

  /* Distortion parameters (from above). */
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

/*--------------------------------------------------------------------------*/

int linperr(const struct linprm *lin, const char *prefix)

{
  if (lin == 0x0) return LINERR_NULL_POINTER;

  if (lin->err && wcserr_prt(lin->err, prefix) == 0) {
    if (lin->dispre) wcserr_prt(lin->dispre->err, prefix);
    if (lin->disseq) wcserr_prt(lin->disseq->err, prefix);
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int linset(struct linprm *lin)

{
  static const char *function = "linset";

  int i, j, naxis, status;
  double *pc, *piximg;
  struct wcserr **err;

  if (lin == 0x0) return LINERR_NULL_POINTER;
  err = &(lin->err);

  naxis = lin->naxis;

  /* Check for a unit matrix. */
  lin->unity = 1;
  pc = lin->pc;
  for (i = 0; i < naxis; i++) {
    for (j = 0; j < naxis; j++) {
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
    if (lin->flag == LINSET) {
      /* Free memory that may have been allocated previously. */
      if (lin->piximg) free(lin->piximg);
      if (lin->imgpix) free(lin->imgpix);
    }

    lin->piximg  = 0x0;
    lin->imgpix  = 0x0;
    lin->i_naxis = 0;

    /* Check cdelt. */
    for (i = 0; i < naxis; i++) {
      if (lin->cdelt[i] == 0.0) {
        return wcserr_set(LIN_ERRMSG(LINERR_SINGULAR_MTX));
      }
    }

  } else {
    if (lin->flag != LINSET || lin->i_naxis < naxis) {
      if (lin->flag == LINSET) {
        /* Free memory that may have been allocated previously. */
        if (lin->piximg) free(lin->piximg);
        if (lin->imgpix) free(lin->imgpix);
      }

      /* Allocate memory for internal arrays. */
      if ((lin->piximg = calloc(naxis*naxis, sizeof(double))) == 0x0) {
        return wcserr_set(LIN_ERRMSG(LINERR_MEMORY));
      }

      if ((lin->imgpix = calloc(naxis*naxis, sizeof(double))) == 0x0) {
        free(lin->piximg);
        return wcserr_set(LIN_ERRMSG(LINERR_MEMORY));
      }

      lin->i_naxis = naxis;
    }

    /* Compute the pixel-to-image transformation matrix. */
    pc     = lin->pc;
    piximg = lin->piximg;
    for (i = 0; i < naxis; i++) {
      for (j = 0; j < naxis; j++) {
        if (lin->disseq == 0x0) {
          /* No sequent distortions, incorporate cdelt into piximg. */
          *(piximg++) = lin->cdelt[i] * (*(pc++));
        } else {
          *(piximg++) = *(pc++);
        }
      }
    }

    /* Compute the image-to-pixel transformation matrix. */
    if ((status = matinv(naxis, lin->piximg, lin->imgpix))) {
      return wcserr_set(LIN_ERRMSG(status));
    }
  }


  /* Set up the distortion functions. */
  lin->affine = 1;
  if (lin->dispre) {
    if ((status = disset(lin->dispre))) {
      return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
    }

    lin->affine = 0;
  }

  if (lin->disseq) {
    if ((status = disset(lin->disseq))) {
      return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
    }

    lin->affine = 0;
  }

  lin->simple = lin->unity && lin->affine;


  /* Create work arrays. */
  if (lin->tmpcrd) free(lin->tmpcrd);
  if ((lin->tmpcrd = calloc(naxis, sizeof(double))) == 0x0) {
    linfree(lin);
    return wcserr_set(LIN_ERRMSG(LINERR_MEMORY));
  }


  lin->flag = LINSET;

  return 0;
}

/*--------------------------------------------------------------------------*/

int linp2x(
  struct linprm *lin,
  int ncoord,
  int nelem,
  const double pixcrd[],
  double imgcrd[])

{
  static const char *function = "linp2x";

  int i, j, k, n, ndbl, nelemn, status;
  double temp;
  register const double *pix;
  register double *img, *piximg, *tmp;
  struct wcserr **err;


  /* Initialize. */
  if (lin == 0x0) return LINERR_NULL_POINTER;
  err = &(lin->err);

  if (lin->flag != LINSET) {
    if ((status = linset(lin))) return status;
  }

  n = lin->naxis;


  /* Convert pixel coordinates to intermediate world coordinates. */
  pix = pixcrd;
  img = imgcrd;

  if (lin->simple) {
    /* Handle the simplest and most common case with maximum efficiency. */
    nelemn = nelem - n;
    for (k = 0; k < ncoord; k++) {
      for (i = 0; i < n; i++) {
        *(img++) = lin->cdelt[i] * (*(pix++) - lin->crpix[i]);
      }

      pix += nelemn;
      img += nelemn;
    }

  } else if (lin->affine) {
    /* No distortions. */
    ndbl   = n * sizeof(double);
    nelemn = nelem - n;
    for (k = 0; k < ncoord; k++) {
      memset(img, 0, ndbl);

      for (j = 0; j < n; j++) {
        /* cdelt will have been incorporated into piximg. */
        piximg = lin->piximg + j;

        /* Column-wise multiplication allows this to be cached. */
        temp = *(pix++) - lin->crpix[j];
        for (i = 0; i < n; i++, piximg += n) {
          img[i] += *piximg * temp;
        }
      }

      pix += nelemn;
      img += nelem;
    }

  } else {
    /* Distortions are present. */
    ndbl = n * sizeof(double);
    tmp  = lin->tmpcrd;

    for (k = 0; k < ncoord; k++) {
      if (lin->dispre) {
        if ((status = disp2x(lin->dispre, pix, tmp))) {
          return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
        }
      } else {
        memcpy(tmp, pix, ndbl);
      }

      if (lin->unity) {
        for (i = 0; i < n; i++) {
          img[i] = tmp[i] - lin->crpix[i];
        }

      } else {
        for (j = 0; j < n; j++) {
          tmp[j] -= lin->crpix[j];
        }

        piximg = lin->piximg;
        for (i = 0; i < n; i++) {
          img[i] = 0.0;
          for (j = 0; j < n; j++) {
            img[i] += *(piximg++) * tmp[j];
          }
        }
      }

      if (lin->disseq) {
        if ((status = disp2x(lin->disseq, img, tmp))) {
          return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
        }

        /* With sequent distortions, cdelt is not incorporated into piximg. */
        for (i = 0; i < n; i++) {
          img[i] = lin->cdelt[i] * tmp[i];
        }

      } else if (lin->unity) {
        /* ...nor if the matrix is unity. */
        for (i = 0; i < n; i++) {
          img[i] *= lin->cdelt[i];
        }
      }

      pix += nelem;
      img += nelem;
    }
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int linx2p(
  struct linprm *lin,
  int ncoord,
  int nelem,
  const double imgcrd[],
  double pixcrd[])

{
  static const char *function = "linx2p";

  int i, j, k, n, ndbl, nelemn, status;
  register const double *img;
  register double *imgpix, *pix, *tmp;
  struct wcserr **err;


  /* Initialize. */
  if (lin == 0x0) return LINERR_NULL_POINTER;
  err = &(lin->err);

  if (lin->flag != LINSET) {
    if ((status = linset(lin))) return status;
  }

  n = lin->naxis;


  /* Convert intermediate world coordinates to pixel coordinates. */
  img = imgcrd;
  pix = pixcrd;

  if (lin->simple) {
    /* Handle the simplest and most common case with maximum efficiency. */
    nelemn = nelem - n;
    for (k = 0; k < ncoord; k++) {
      for (j = 0; j < n; j++) {
        *(pix++) = (*(img++) / lin->cdelt[j]) + lin->crpix[j];
      }

      img += nelemn;
      pix += nelemn;
    }

  } else if (lin->affine) {
    /* No distortions. */
    nelemn = nelem - n;
    for (k = 0; k < ncoord; k++) {
      /* cdelt will have been incorporated into imgpix. */
      imgpix = lin->imgpix;

      for (j = 0; j < n; j++) {
        *pix = 0.0;
        for (i = 0; i < n; i++) {
          *pix += *imgpix * img[i];
          imgpix++;
        }

        *(pix++) += lin->crpix[j];
      }

      img += nelem;
      pix += nelemn;
    }

  } else {
    /* Distortions are present. */
    ndbl = n * sizeof(double);
    tmp  = lin->tmpcrd;

    for (k = 0; k < ncoord; k++) {
      if (lin->disseq) {
        /* With sequent distortions, cdelt is not incorporated into imgpix. */
        for (i = 0; i < n; i++) {
          tmp[i] = img[i] / lin->cdelt[i];
        }

        if ((status = disx2p(lin->disseq, tmp, pix))) {
          return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
        }

        memcpy(tmp, pix, ndbl);

      } else if (lin->unity) {
        /* ...nor if the matrix is unity. */
        for (i = 0; i < n; i++) {
          tmp[i] = img[i] / lin->cdelt[i];
        }

      } else {
        /* cdelt will have been incorporated into imgpix. */
        memcpy(tmp, img, ndbl);
      }

      if (lin->unity) {
        for (j = 0; j < n; j++) {
          pix[j] = tmp[j] + lin->crpix[j];
        }

      } else {
        imgpix = lin->imgpix;
        for (j = 0; j < n; j++) {
          pix[j] = lin->crpix[j];
          for (i = 0; i < n; i++) {
            pix[j] += *(imgpix++) * tmp[i];
          }
        }
      }

      if (lin->dispre) {
        memcpy(tmp, pix, ndbl);

        if ((status = disx2p(lin->dispre, tmp, pix))) {
          return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
        }
      }

      img += nelem;
      pix += nelem;
    }
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

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

  int carry, i, j, naxis, ncoord, status = 0;
  double dpix, dpx2, dssq, *img, *pix0, *pix0p, *pix1, *pix1p, *pixend,
         *pixinc, pixspan, *ssqdis, ssqtot, *sumdis, sumtot, totdis;
  struct linprm affine;
  struct wcserr **err;


  /* Initialize. */
  if (lin == 0x0) return LINERR_NULL_POINTER;
  err = &(lin->err);

  naxis = lin->naxis;

  if (nsamp) *nsamp = 0;
  for (j = 0; j < naxis; j++) {
    if (maxdis) maxdis[j] = 0.0;
    if (avgdis) avgdis[j] = 0.0;
    if (rmsdis) rmsdis[j] = 0.0;
  }
  if (maxtot) *maxtot = 0.0;
  if (avgtot) *avgtot = 0.0;
  if (rmstot) *rmstot = 0.0;

  /* Quick return if no distortions. */
  if (lin->affine) return 0;

  /* It's easier if there are no sequent distortions! */
  if (lin->disseq == 0x0) {
    status = diswarp(lin->dispre, pixblc, pixtrc, pixsamp, nsamp,
                     maxdis, maxtot, avgdis, avgtot, rmsdis, rmstot);
    return wcserr_set(LIN_ERRMSG(lin_diserr[status]));
  }

  /* Make a reference copy of lin without distortions. */
  affine.flag = -1;
  if ((status = (lincpy(1, lin, &affine) ||
                 lindist(1, &affine, 0x0, 0) ||
                 lindist(2, &affine, 0x0, 0) ||
                 linset(&affine)))) {
    return wcserr_set(LIN_ERRMSG(status));
  }

  /* Work out increments on each axis. */
  pixinc = lin->tmpcrd;
  ncoord = 0;
  for (j = 0; j < naxis; j++) {
    pixspan = pixtrc[j] - (pixblc ? pixblc[j] : 1.0);

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
      /* Number of samples on axis 1. */
      ncoord = 1 + (int)((pixspan/pixinc[0]) + 0.5);
    }
  }

  /* Get memory for processing the image row by row. */
  if ((pix0 = calloc((3*ncoord+4)*naxis, sizeof(double))) == 0x0) {
    return wcserr_set(LIN_ERRMSG(LINERR_MEMORY));
  }

  img    = pix0 + naxis*ncoord;
  pix1   = img  + naxis*ncoord;
  pixinc = pix1 + naxis*ncoord;
  pixend = pixinc + naxis;
  sumdis = pixend + naxis;
  ssqdis = sumdis + naxis;


  /* Copy tmpcrd since linp2x() will overwrite it. */
  memcpy(pixinc, lin->tmpcrd, naxis*sizeof(double));

  /* Set up the array of pixel coordinates. */
  for (j = 0; j < naxis; j++) {
    pix0[j] = pixblc ? pixblc[j] : 1.0;
    pixend[j] = pixtrc[j] + 0.5*pixinc[j];
  }

  pix0p = pix0 + naxis;
  for (i = 1; i < ncoord; i++) {
    *(pix0p++) = pix0[0] + i*pixinc[0];

    for (j = 1; j < naxis; j++) {
      *(pix0p++) = pix0[j];
    }
  }

  /* Initialize accumulators. */
  for (j = 0; j < naxis; j++) {
    sumdis[j] = 0.0;
    ssqdis[j] = 0.0;
  }
  sumtot = 0.0;
  ssqtot = 0.0;


  /* Loop over N dimensions. */
  carry = 0;
  while (carry == 0) {
    if ((status = linp2x(lin, ncoord, naxis, pix0, img))) {
      /* (Preserve the error message set by linp2x().) */
      goto cleanup;
    }

    if ((status = linx2p(&affine, ncoord, naxis, img, pix1))) {
      /* (Preserve the error message set by linx2p().) */
      goto cleanup;
    }

    /* Accumulate statistics. */
    pix0p = pix0;
    pix1p = pix1;
    for (i = 0; i < ncoord; i++) {
      (*nsamp)++;

      dssq = 0.0;
      for (j = 0; j < naxis; j++) {
        dpix = *(pix1p++) - *(pix0p++);
        dpx2 = dpix*dpix;

        sumdis[j] += dpix;
        ssqdis[j] += dpx2;

        if (maxdis && (dpix = fabs(dpix)) > maxdis[j]) maxdis[j] = dpix;

        dssq += dpx2;
      }

      totdis = sqrt(dssq);
      sumtot += totdis;
      ssqtot += totdis*totdis;

      if (maxtot && *maxtot < totdis) *maxtot = totdis;
    }

    /* Next array of pixel coordinates. */
    for (j = 1; j < naxis; j++) {
      pix0[j] += pixinc[j];
      if ((carry = (pix0[j] > pixend[j]))) {
        pix0[j] = pixblc ? pixblc[j] : 1.0;
      }

      pix0p = pix0 + naxis + j;
      for (i = 1; i < ncoord; i++) {
        *pix0p = pix0[j];
        pix0p += naxis;
      }

      if (carry == 0) break;
    }
  }


  /* Compute the means and RMSs. */
  for (j = 0; j < naxis; j++) {
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

/*--------------------------------------------------------------------------*/

int matinv(int n, const double mat[], double inv[])

{
  register int i, ij, ik, j, k, kj, pj;
  int    itemp, *mxl, *lxm, pivot;
  double colmax, *lu, *rowmax, dtemp;


  /* Allocate memory for internal arrays. */
  if ((mxl = calloc(n, sizeof(int))) == 0x0) {
    return LINERR_MEMORY;
  }
  if ((lxm = calloc(n, sizeof(int))) == 0x0) {
    free(mxl);
    return LINERR_MEMORY;
  }

  if ((rowmax = calloc(n, sizeof(double))) == 0x0) {
    free(mxl);
    free(lxm);
    return LINERR_MEMORY;
  }

  if ((lu = calloc(n*n, sizeof(double))) == 0x0) {
    free(mxl);
    free(lxm);
    free(rowmax);
    return LINERR_MEMORY;
  }


  /* Initialize arrays. */
  for (i = 0, ij = 0; i < n; i++) {
    /* Vector that records row interchanges. */
    mxl[i] = i;

    rowmax[i] = 0.0;

    for (j = 0; j < n; j++, ij++) {
      dtemp = fabs(mat[ij]);
      if (dtemp > rowmax[i]) rowmax[i] = dtemp;

      lu[ij] = mat[ij];
    }

    /* A row of zeroes indicates a singular matrix. */
    if (rowmax[i] == 0.0) {
      free(mxl);
      free(lxm);
      free(rowmax);
      free(lu);
      return LINERR_SINGULAR_MTX;
    }
  }


  /* Form the LU triangular factorization using scaled partial pivoting. */
  for (k = 0; k < n; k++) {
    /* Decide whether to pivot. */
    colmax = fabs(lu[k*n+k]) / rowmax[k];
    pivot = k;

    for (i = k+1; i < n; i++) {
      ik = i*n + k;
      dtemp = fabs(lu[ik]) / rowmax[i];
      if (dtemp > colmax) {
        colmax = dtemp;
        pivot = i;
      }
    }

    if (pivot > k) {
      /* We must pivot, interchange the rows of the design matrix. */
      for (j = 0, pj = pivot*n, kj = k*n; j < n; j++, pj++, kj++) {
        dtemp = lu[pj];
        lu[pj] = lu[kj];
        lu[kj] = dtemp;
      }

      /* Amend the vector of row maxima. */
      dtemp = rowmax[pivot];
      rowmax[pivot] = rowmax[k];
      rowmax[k] = dtemp;

      /* Record the interchange for later use. */
      itemp = mxl[pivot];
      mxl[pivot] = mxl[k];
      mxl[k] = itemp;
    }

    /* Gaussian elimination. */
    for (i = k+1; i < n; i++) {
      ik = i*n + k;

      /* Nothing to do if lu[ik] is zero. */
      if (lu[ik] != 0.0) {
        /* Save the scaling factor. */
        lu[ik] /= lu[k*n+k];

        /* Subtract rows. */
        for (j = k+1; j < n; j++) {
          lu[i*n+j] -= lu[ik]*lu[k*n+j];
        }
      }
    }
  }


  /* mxl[i] records which row of mat corresponds to row i of lu.  */
  /* lxm[i] records which row of lu  corresponds to row i of mat. */
  for (i = 0; i < n; i++) {
    lxm[mxl[i]] = i;
  }


  /* Determine the inverse matrix. */
  for (i = 0, ij = 0; i < n; i++) {
    for (j = 0; j < n; j++, ij++) {
      inv[ij] = 0.0;
    }
  }

  for (k = 0; k < n; k++) {
    inv[lxm[k]*n+k] = 1.0;

    /* Forward substitution. */
    for (i = lxm[k]+1; i < n; i++) {
      for (j = lxm[k]; j < i; j++) {
        inv[i*n+k] -= lu[i*n+j]*inv[j*n+k];
      }
    }

    /* Backward substitution. */
    for (i = n-1; i >= 0; i--) {
      for (j = i+1; j < n; j++) {
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
