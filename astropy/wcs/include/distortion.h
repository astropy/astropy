/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __DISTORTION_H__
#define __DISTORTION_H__

#include "util.h"

/* TODO: This is all two-dimensional.  Should be made
   multi-dimensional in the future. */
#define NAXES 2

#define MAXAXES 6

/**
A structure to contain the information for a single distortion lookup table
 */
typedef struct {
  unsigned int                   naxis[NAXES]; /* size of distortion image */
  double                         crpix[NAXES];
  double                         crval[NAXES];
  double                         cdelt[NAXES];
  /* The data is not "owned" by this structure.  It is the user's
     responsibility to free it. */
  /*@shared@*/ /*@null@*/ float *data;
} distortion_lookup_t;

/**
Initialize a lookup table to reasonable default values.
 */
int
distortion_lookup_t_init(distortion_lookup_t* lookup);

/**
Cleanup after a lookup table.  Currently does nothing, but may do
something in the future, so please call it when you are done with
the lookup table.  It does not free the data pointed to be the
lookup table -- it is the user's responsibility to free that array.
 */
void
distortion_lookup_t_free(distortion_lookup_t* lookup);

/**
Lookup the distortion offset for a particular pixel coordinate in
the lookup table.

@param lookup A lookup table object

@param A coordinate pair

@return The offset as determined by binlinear interpolation in the
lookup table
*/
double
get_distortion_offset(
    const distortion_lookup_t * const lookup,
    const double * const img /* [NAXES] */);

/**
Perform just the distortion table part of Paper IV.

@param naxes

@param lookups A pair of lookup table objects

@param nelem

@param pix [in]: An array of pixel coordinates

@param foc [out]: An array of focal plane coordinates

@return A wcslib error code
*/
int
p4_pix2foc(
    const unsigned int naxes,
    const distortion_lookup_t** lookups, /* [NAXES] */
    const unsigned int nelem,
    const double* pix, /* [NAXES][nelem] */
    double *foc /* [NAXES][nelem] */);

/**
Perform just the distortion table part of Paper IV, by adding
distortion to the values already in place in foc.

@param naxes

@param lookups A pair of lookup table objects

@param nelem

@param pix [in]: An array of pixel coordinates

@param foc [in/out]: An array of focal plane coordinates

@return A wcslib error code
*/
int
p4_pix2deltas(
    const unsigned int naxes,
    const distortion_lookup_t** lookups, /* [NAXES] */
    const unsigned int nelem,
    const double* pix, /* [NAXES][nelem] */
    double *foc /* [NAXES][nelem] */);

#endif /* __DISTORTION_H__ */
