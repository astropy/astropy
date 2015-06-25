/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __PIPELINE_H__
#define __PIPELINE_H__

#include "sip.h"
#include "distortion.h"
#include "wcs.h"

typedef struct {
  distortion_lookup_t*                   det2im[2];
#if !defined(HAVE_WCSLIB_VERSION)
  /*@shared@*/ /*@null@*/ sip_t*         sip;
#endif
  distortion_lookup_t*                   cpdis[2];
  /*@shared@*/ /*@null@*/ struct wcsprm* wcs;
  struct wcserr*                         err;
} pipeline_t;

/**
Initialize all the values in a pipeline_t to NULL.
*/
void
pipeline_clear(
    pipeline_t* pipeline);

/**
Set all the values of a pipeline_t.
*/
void
pipeline_init(
    pipeline_t* pipeline,
    /*@shared@*/ distortion_lookup_t** det2im /* [2] */,
    /*@shared@*/ sip_t* sip,
    /*@shared@*/ distortion_lookup_t** cpdis /* [2] */,
    /*@shared@*/ struct wcsprm* wcs);

/**
Free all the temporary buffers of a pipeline_t.  It does not free
the underlying sip_t, distortion_lookup_t or wcsprm objects.
*/
void
pipeline_free(
    pipeline_t* pipeline);

/**
Perform the entire pipeline from pixel coordinates to world
coordinates, in the following order:

    - Detector to image plane correction (optionally)

    - SIP distortion correction (optionally)

    - FITS WCS distortion paper correction (optionally)

    - wcslib WCS transformation

@param ncoord:

@param nelem:

@param pixcrd [in]: Array of pixel coordinates.

@param world [out]: Array of world coordinates (output).

@return: A wcslib error code.
*/
int
pipeline_all_pixel2world(
    pipeline_t* pipeline,
    const unsigned int ncoord,
    const unsigned int nelem,
    const double* const pixcrd /* [ncoord][nelem] */,
    double* world /* [ncoord][nelem] */);

/**
Perform just the distortion correction part of the pipeline from pixel
coordinates to focal plane coordinates.

    - Detector to image plane correction (optionally)

    - SIP distortion correction (optionally)

    - FITS WCS distortion paper correction (optionally)

@param ncoord:

@param nelem:

@param pixcrd [in]: Array of pixel coordinates.

@param foc [out]: Array of focal plane coordinates.

@return: A wcslib error code.
*/
int
pipeline_pix2foc(
    pipeline_t* pipeline,
    const unsigned int ncoord,
    const unsigned int nelem,
    const double* const pixcrd /* [ncoord][nelem] */,
    double* foc /* [ncoord][nelem] */);

#endif
