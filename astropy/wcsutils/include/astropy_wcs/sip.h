/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __SIP_H__
#define __SIP_H__

#include "util.h"

typedef struct {
  unsigned int                    a_order;
  /*@null@*/ /*@shared@*/ double* a;
  unsigned int                    b_order;
  /*@null@*/ /*@shared@*/ double* b;
  unsigned int                    ap_order;
  /*@null@*/ /*@shared@*/ double* ap;
  unsigned int                    bp_order;
  /*@null@*/ /*@shared@*/ double* bp;
  double                          crpix[2];
  /*@null@*/ double*              scratch;
  struct wcserr*                  err;
} sip_t;

/**
Sets all the values of the sip_t structure to NULLs or zeros.
*/
void
sip_clear(sip_t* sip);

/**
Set the values of the sip_t structure.

The values expected are all exactly as defined in the FITS SIP header
keywords.

The arrays/matrices are all *copied* into the SIP struct.  To free the
memory that sip_t allocates for itself, call sip_free.

@param a_order: The order of the A_i_j matrix

@param a: The A_i_j array, which must be of size [a_order+1][a_order+1]

@param b_order: The order of the B_i_j matrix

@param b: The B_i_j array, which must be of size [b_order+1][b_order+1]

@param ap_order: The order of the AP_i_j matrix

@param ap: The AP_i_j array, which must be of size [ap_order+1][ap_order+1]

@param bp_order: The order of the BP_i_j matrix

@param bp: The BP_i_j array, which must be of size [bp_order+1][bp_order+1]

@param crpix: The position of the reference pixel
*/
int
sip_init(
    sip_t* sip,
    const unsigned int a_order, const double* a,
    const unsigned int b_order, const double* b,
    const unsigned int ap_order, const double* ap,
    const unsigned int bp_order, const double* bp,
    const double* crpix /* [2] */);

/**
Frees the memory allocated for the sip_t struct.
*/
void
sip_free(sip_t* sip);

/**
Converts pixel coordinates to focal plane coordinates using the SIP
polynomial distortion convention, and the values stored in the sip_t
struct.

@param naxes

@param nelem

@param pix [in]: An array of pixel coordinates

@param foc [out]: An array of focal plane coordinates

@return A wcslib error code
*/
int
sip_pix2foc(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* pix /* [NAXES][nelem] */,
    double* foc /* [NAXES][nelem] */);

/**
Computes the offset deltas necessary to convert pixel coordinates to
focal plane coordinates using the SIP polynomial distortion
convention, and the values stored in the sip_t struct.  The deltas are
added to the existing values in pix.

@param naxes

@param nelem

@param pix [in]: An array of pixel coordinates

@param foc [in/out]: An array of deltas, that when added to pix
results in focal plane coordinates.

@return A wcslib error code
*/
int
sip_pix2deltas(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* pix /* [NAXES][nelem] */,
    double* foc /* [NAXES][nelem] */);

/**
Adds the offset deltas necessary to convert focal plane
coordinates to pixel coordinates using the SIP polynomial distortion
convention, and the values stored in the sip_t struct.  The deltas
are added to the existing values in pix.

@param naxes

@param nelem

@param foc [in]: An array of focal plane coordinates

@param pix [in/out]: An array of pixel coordinates

@return A wcslib error code
*/
int
sip_foc2pix(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* foc /* [NAXES][nelem] */,
    double* pix /* [NAXES][nelem] */);

/**
Computes the offset deltas necessary to convert focal plane
coordinates to pixel coordinates using the SIP polynomial distortion
convention, and the values stored in the sip_t struct.  The deltas are
added to the existing values in foc.

@param naxes

@param nelem

@param foc [in]: An array of focal plane coordinates

@param foc [in/out]: An array of deltas, that when added to pix
results in focal plane coordinates.

@return A wcslib error code
*/
int
sip_foc2deltas(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* foc /* [NAXES][nelem] */,
    double* deltas /* [NAXES][nelem] */);

#endif
