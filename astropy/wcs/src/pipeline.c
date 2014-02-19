/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#include "astropy_wcs/pipeline.h"
#include "astropy_wcs/util.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#define PIP_ERRMSG(status) WCSERR_SET(status)

void
pipeline_clear(
    pipeline_t* pipeline) {

  pipeline->det2im[0] = NULL;
  pipeline->det2im[1] = NULL;
  pipeline->sip = NULL;
  pipeline->cpdis[0] = NULL;
  pipeline->cpdis[1] = NULL;
  pipeline->wcs = NULL;
  pipeline->err = NULL;
}

void
pipeline_init(
    pipeline_t* pipeline,
    /*@shared@*/ distortion_lookup_t** det2im /* [2] */,
    /*@shared@*/ sip_t* sip,
    /*@shared@*/ distortion_lookup_t** cpdis /* [2] */,
    /*@shared@*/ struct wcsprm* wcs) {

  pipeline->det2im[0] = det2im[0];
  pipeline->det2im[1] = det2im[1];
  pipeline->sip = sip;
  pipeline->cpdis[0] = cpdis[0];
  pipeline->cpdis[1] = cpdis[1];
  pipeline->wcs = wcs;
  pipeline->err = NULL;
}

void
pipeline_free(
    pipeline_t* pipeline) {

  free(pipeline->err);
  pipeline->err = NULL;
}

int
pipeline_all_pixel2world(
    pipeline_t* pipeline,
    const unsigned int ncoord,
    const unsigned int nelem,
    const double* const pixcrd /* [ncoord][nelem] */,
    double* world /* [ncoord][nelem] */) {

  static const char* function = "pipeline_all_pixel2world";

  const double*   wcs_input  = NULL;
  double*         wcs_output = NULL;
  int             has_det2im;
  int             has_sip;
  int             has_p4;
  int             has_wcs;
  int             status     = 1;
  struct wcserr **err;

  /* Temporary buffer for performing WCS calculations */
  unsigned char*     buffer = NULL;
  unsigned char*     mem = NULL;
  /*@null@*/ double* tmp;
  /*@null@*/ double* imgcrd;
  /*@null@*/ double* phi;
  /*@null@*/ double* theta;
  /*@null@*/ int*    stat;

  if (pipeline == NULL || pixcrd == NULL || world == NULL) {
    return WCSERR_NULL_POINTER;
  }

  err = &(pipeline->err);

  has_det2im = pipeline->det2im[0] != NULL || pipeline->det2im[1] != NULL;
  has_sip    = pipeline->sip != NULL;
  has_p4     = pipeline->cpdis[0] != NULL || pipeline->cpdis[1] != NULL;
  has_wcs    = pipeline->wcs != NULL;

  if (has_det2im || has_sip || has_p4) {
    if (nelem != 2) {
      status = wcserr_set(
        PIP_ERRMSG(WCSERR_BAD_COORD_TRANS),
        "Data must be 2-dimensional when Paper IV lookup table or SIP transform is present.");
      goto exit;
    }
  }

  if (has_wcs) {
    buffer = mem = malloc(
        ncoord * nelem * sizeof(double) + /* imgcrd */
        ncoord * sizeof(double) +         /* phi */
        ncoord * sizeof(double) +         /* theta */
        ncoord * nelem * sizeof(double) + /* tmp */
        ncoord * nelem * sizeof(int)      /* stat */
        );

    if (buffer == NULL) {
      status = wcserr_set(
        PIP_ERRMSG(WCSERR_MEMORY), "Memory allocation failed");
      goto exit;
    }

    imgcrd = (double *)mem;
    mem += ncoord * nelem * sizeof(double);

    phi = (double *)mem;
    mem += ncoord * sizeof(double);

    theta = (double *)mem;
    mem += ncoord * sizeof(double);

    tmp = (double *)mem;
    mem += ncoord * nelem * sizeof(double);

    stat = (int *)mem;
    /* mem += ncoord * nelem * sizeof(int); */

    if (has_det2im || has_sip || has_p4) {
      status = pipeline_pix2foc(pipeline, ncoord, nelem, pixcrd, tmp);
      if (status != 0) {
        goto exit;
      }

      wcs_input = tmp;
      wcs_output = world;
    } else {
      wcs_input = pixcrd;
      wcs_output = world;
    }

    if ((status = wcsp2s(pipeline->wcs, (int)ncoord, (int)nelem, wcs_input, imgcrd,
                         phi, theta, wcs_output, stat))) {
      if (pipeline->err == NULL) {
        pipeline->err = calloc(1, sizeof(struct wcserr));
      }
      wcserr_copy(pipeline->wcs->err, pipeline->err);
    }

    if (status == 8) {
      set_invalid_to_nan((int)ncoord, (int)nelem, wcs_output, stat);
    }
  } else {
    if (has_det2im || has_sip || has_p4) {
      status = pipeline_pix2foc(pipeline, ncoord, nelem, pixcrd, world);
    }
  }

 exit:
  free(buffer);

  return status;
}

int pipeline_pix2foc(
    pipeline_t* pipeline,
    const unsigned int ncoord,
    const unsigned int nelem,
    const double* const pixcrd /* [ncoord][nelem] */,
    double* foc /* [ncoord][nelem] */) {

  static const char* function = "pipeline_pix2foc";

  int              has_det2im;
  int              has_sip;
  int              has_p4;
  const double *   input  = NULL;
  double *         tmp    = NULL;
  int              status = 1;
  struct wcserr  **err;

  assert(nelem == 2);
  assert(pixcrd != foc);

  if (pipeline == NULL || pixcrd == NULL || foc == NULL) {
    return WCSERR_NULL_POINTER;
  }

  err = &(pipeline->err);

  has_det2im = pipeline->det2im[0] != NULL || pipeline->det2im[1] != NULL;
  has_sip    = pipeline->sip != NULL;
  has_p4     = pipeline->cpdis[0] != NULL || pipeline->cpdis[1] != NULL;

  if (has_det2im) {
    if (has_sip || has_p4) {
      tmp = malloc(ncoord * nelem * sizeof(double));
      if (tmp == NULL) {
        status = wcserr_set(
          PIP_ERRMSG(WCSERR_MEMORY), "Memory allocation failed");
        goto exit;
      }

      memcpy(tmp, pixcrd, sizeof(double) * ncoord * nelem);

      status = p4_pix2deltas(2, (void*)pipeline->det2im, ncoord, pixcrd, tmp);
      if (status) {
        wcserr_set(PIP_ERRMSG(WCSERR_NULL_POINTER), "NULL pointer passed");
        goto exit;
      }

      input = tmp;
      memcpy(foc, input, sizeof(double) * ncoord * nelem);
    } else {
      memcpy(foc, pixcrd, sizeof(double) * ncoord * nelem);

      status = p4_pix2deltas(2, (void*)pipeline->det2im, ncoord, pixcrd, foc);
      if (status) {
        wcserr_set(PIP_ERRMSG(WCSERR_NULL_POINTER), "NULL pointer passed");
        goto exit;
      }
    }
  } else {
    /* Copy pixcrd to foc as a starting point.  The "deltas" functions
       below will undistort from there */
    memcpy(foc, pixcrd, sizeof(double) * ncoord * nelem);
    input = pixcrd;
  }

  if (has_sip) {
    status = sip_pix2deltas(pipeline->sip, 2, ncoord, input, foc);
    if (status) {
      if (pipeline->err == NULL) {
        pipeline->err = calloc(1, sizeof(struct wcserr));
      }
      wcserr_copy(pipeline->sip->err, pipeline->err);
      goto exit;
    }
  }

  if (has_p4) {
    status = p4_pix2deltas(2, (void*)pipeline->cpdis, ncoord, input, foc);
    if (status) {
      wcserr_set(PIP_ERRMSG(WCSERR_NULL_POINTER), "NULL pointer passed");
      goto exit;
    }
  }

  status = 0;

 exit:
  free(tmp);

  return status;
}
