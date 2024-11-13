/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#include "astropy_wcs/distortion.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* TODO: n-dimensional support */

int
distortion_lookup_t_init(
    distortion_lookup_t* lookup) {

  unsigned int i;

  for (i = 0; i < NAXES; ++i) {
    lookup->naxis[i] = 0;
    lookup->crpix[i] = 0.0;
    lookup->crval[i] = 0.0;
    lookup->cdelt[i] = 1.0;
  }

  lookup->data = NULL;

  return 0;
}

void
distortion_lookup_t_free(
    /*@unused@*/ distortion_lookup_t* lookup) {

  /*@empty@*/
}

/**
 * Get a value at a specific integral location in the lookup table.
 * (This is nothing more special than an array lookup with range
 * checking.)
 */
static INLINE float
get_dist_clamp(
    const float* const data,
    const unsigned int* const naxis,
    const int x,
    const int y) {

  return data[
    ((naxis[0] * CLAMP(y, 0, (long)naxis[1] - 1)) +
     CLAMP(x, 0, (long)naxis[0] - 1))];
}

static INLINE float
get_dist(
    const float* const data,
    const unsigned int* const naxis,
    const int x,
    const int y) {

  return data[(naxis[0] * y) + x];
}

/**
 * Converts a pixel coordinate to a fractional coordinate in the
 * lookup table on a single axis
 */
static INLINE double
image_coord_to_distortion_coord(
    const distortion_lookup_t * const lookup,
    const unsigned int axis,
    const double img) {

  double result;

  assert(lookup != NULL);
  assert(axis < NAXES);

  /* The "- 1" is here because the input coordinates are 1-based,
     but this is a C-array underneath */
  result = (
      ((img - lookup->crval[axis]) / lookup->cdelt[axis]) +
      lookup->crpix[axis]) - 1.0;

  return CLAMP(result, 0.0, (double)(lookup->naxis[axis] - 1));
}

/**
 * Converts a pixel coordinate to a fractional coordinate in the
 * lookup table.
 */
static INLINE void
image_coords_to_distortion_coords(
    const distortion_lookup_t * const lookup,
    const double * const img /* [NAXES] */,
    /* Output parameters */
    /*@out@*/ double *dist /* [NAXES] */) {

  unsigned int i;

  assert(lookup != NULL);
  assert(img != NULL);
  assert(dist != NULL);

  for (i = 0; i < NAXES; ++i) {
    dist[i] = image_coord_to_distortion_coord(lookup, i, img[i]);
  }
}

INLINE double
get_distortion_offset(
    const distortion_lookup_t * const lookup,
    const double * const img /*[NAXES]*/) {

  double              dist[NAXES];
  double              dist_floor[NAXES];
  int                 dist_ifloor[NAXES];
  double              dist_weight[NAXES];
  double              dist_iweight[NAXES];
  double              result;
  const unsigned int* naxis = lookup->naxis;
  const float*        data  = lookup->data;
  unsigned int        i;

  assert(lookup != NULL);
  assert(img != NULL);

  image_coords_to_distortion_coords(lookup, img, dist);

  for (i = 0; i < NAXES; ++i) {
    dist_floor[i] = floor(dist[i]);
    dist_ifloor[i] = (int)dist_floor[i];
    dist_weight[i] = dist[i] - dist_floor[i];
    dist_iweight[i] = 1.0 - dist_weight[i];
  }

  /* If we may need to clamp the lookups, use this slower approach */
  if (dist_ifloor[0] < 0 ||
      dist_ifloor[1] < 0 ||
      dist_ifloor[0] >= (long)lookup->naxis[0] - 1 ||
      dist_ifloor[1] >= (long)lookup->naxis[1] - 1) {
    result =
      (double)get_dist_clamp(data, naxis, dist_ifloor[0],     dist_ifloor[1])     * dist_iweight[0] * dist_iweight[1] +
      (double)get_dist_clamp(data, naxis, dist_ifloor[0],     dist_ifloor[1] + 1) * dist_iweight[0] * dist_weight[1] +
      (double)get_dist_clamp(data, naxis, dist_ifloor[0] + 1, dist_ifloor[1])     * dist_weight[0] * dist_iweight[1] +
      (double)get_dist_clamp(data, naxis, dist_ifloor[0] + 1, dist_ifloor[1] + 1) * dist_weight[0] * dist_weight[1];
  /* Else, we don't need to clamp 4 times for each pixel */
  } else {
    result =
      (double)get_dist(data, naxis, dist_ifloor[0],     dist_ifloor[1])     * dist_iweight[0] * dist_iweight[1] +
      (double)get_dist(data, naxis, dist_ifloor[0],     dist_ifloor[1] + 1) * dist_iweight[0] * dist_weight[1] +
      (double)get_dist(data, naxis, dist_ifloor[0] + 1, dist_ifloor[1])     * dist_weight[0] * dist_iweight[1] +
      (double)get_dist(data, naxis, dist_ifloor[0] + 1, dist_ifloor[1] + 1) * dist_weight[0] * dist_weight[1];
  }

  return result;
}

int
p4_pix2deltas(
    const unsigned int naxes,
    const distortion_lookup_t **lookup, /* [NAXES] */
    const unsigned int nelem,
    const double* pix, /* [NAXES][nelem] */
    double *foc /* [NAXES][nelem] */) {

  int i;
  double* foc0;
  const double* pix0;
  const double* pixend;

#ifndef NDEBUG
  unsigned int k;
#endif

  assert(naxes == NAXES);
  assert(lookup != NULL);
  assert(pix != NULL);
  assert(foc != NULL);

#ifndef NDEBUG
  for (k = 0; k < naxes; ++k) {
    if (lookup[k] != NULL) {
      assert(lookup[k]->data != NULL);
    }
  }
#endif

  if (pix == NULL || foc == NULL) {
    return 1;
  }

  pixend = pix + nelem * NAXES;
  /* This can't be parallelized, because pix may be equal to foc */
  /* For the same reason, i needs to be in the inner loop */
  for (pix0 = pix, foc0 = foc; pix0 < pixend; pix0 += NAXES, foc0 += NAXES) {
    for (i = 0; i < NAXES; ++i) {
      if (lookup[i]) {
        foc0[i] += get_distortion_offset(lookup[i], pix0);
      }
    }
  }

  return 0;
}

int
p4_pix2foc(
    const unsigned int naxes,
    const distortion_lookup_t **lookup, /* [NAXES] */
    const unsigned int nelem,
    const double* pix, /* [NAXES][nelem] */
    double *foc /* [NAXES][nelem] */) {

  assert(pix);
  assert(foc);

  if (pix != foc) {
    memcpy(foc, pix, sizeof(double) * naxes * nelem);
  }

  return p4_pix2deltas(naxes, lookup, nelem, pix, foc);
}
