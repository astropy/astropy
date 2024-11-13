/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#include "astropy_wcs/sip.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <wcserr.h>

#define SIP_ERRMSG(status) WCSERR_SET(status)

void
sip_clear(
    sip_t* sip) {

  assert(sip != NULL);

  sip->a_order = 0;
  sip->a = NULL;
  sip->b_order = 0;
  sip->b = NULL;
  sip->ap_order = 0;
  sip->ap = NULL;
  sip->bp_order = 0;
  sip->bp = NULL;
  sip->crpix[0] = 0.0;
  sip->crpix[1] = 0.0;
  sip->scratch = NULL;
  sip->err = NULL;
}

int
sip_init(
    sip_t* sip,
    const unsigned int a_order, const double* a,
    const unsigned int b_order, const double* b,
    const unsigned int ap_order, const double* ap,
    const unsigned int bp_order, const double* bp,
    const double* crpix /* [2] */) {

  unsigned int       a_size       = 0;
  unsigned int       b_size       = 0;
  unsigned int       ap_size      = 0;
  unsigned int       bp_size      = 0;
  unsigned int       scratch_size = 0;
  int                status       = 0;
  struct wcserr**    err          = NULL;
  static const char *function     = "sip_init";

  assert(sip != NULL);
  sip_clear(sip);
  err = &(sip->err);

  /* We we have one of A/B or AP/BP, we must have both. */
  if ((a == NULL) ^ (b == NULL)) {
    return wcserr_set(
      SIP_ERRMSG(WCSERR_BAD_COORD_TRANS),
      "Both A and B SIP transform must be defined");
  }

  if ((ap == NULL) ^ (bp == NULL)) {
    return wcserr_set(
      SIP_ERRMSG(WCSERR_BAD_COORD_TRANS),
      "Both AP and BP SIP transform must be defined");
  }


  if (a != NULL) {
    sip->a_order = a_order;
    a_size = (a_order + 1) * (a_order + 1) * sizeof(double);
    sip->a = malloc(a_size);
    if (sip->a == NULL) {
      sip_free(sip);
      status = wcserr_set(
        SIP_ERRMSG(WCSERR_MEMORY), "Memory allocation failed");
      goto exit;
    }
    memcpy(sip->a, a, a_size);
    if (a_order > scratch_size) {
      scratch_size = a_order;
    }

    sip->b_order = b_order;
    b_size = (b_order + 1) * (b_order + 1) * sizeof(double);
    sip->b = malloc(b_size);
    if (sip->b == NULL) {
      sip_free(sip);
      status = wcserr_set(
        SIP_ERRMSG(WCSERR_MEMORY), "Memory allocation failed");
      goto exit;
    }
    memcpy(sip->b, b, b_size);
    if (b_order > scratch_size) {
      scratch_size = b_order;
    }
  }

  if (ap != NULL) {
    sip->ap_order = ap_order;
    ap_size = (ap_order + 1) * (ap_order + 1) * sizeof(double);
    sip->ap = malloc(ap_size);
    if (sip->ap == NULL) {
      sip_free(sip);
      status = wcserr_set(
        SIP_ERRMSG(WCSERR_MEMORY), "Memory allocation failed");
      goto exit;
    }
    memcpy(sip->ap, ap, ap_size);
    if (ap_order > scratch_size) {
      scratch_size = ap_order;
    }

    sip->bp_order = bp_order;
    bp_size = (bp_order + 1) * (bp_order + 1) * sizeof(double);
    sip->bp = malloc(bp_size);
    if (sip->bp == NULL) {
      sip_free(sip);
      status = wcserr_set(
        SIP_ERRMSG(WCSERR_MEMORY), "Memory allocation failed");
      goto exit;
    }
    memcpy(sip->bp, bp, bp_size);
    if (bp_order > scratch_size) {
      scratch_size = bp_order;
    }
  }

  scratch_size = (scratch_size + 1) * sizeof(double);
  sip->scratch = malloc(scratch_size);
  if (sip->scratch == NULL) {
    sip_free(sip);
    status = wcserr_set(
      SIP_ERRMSG(WCSERR_MEMORY), "Memory allocation failed");
    goto exit;
  }

  sip->crpix[0] = crpix[0];
  sip->crpix[1] = crpix[1];

 exit:

  return status;
}

void
sip_free(sip_t* sip) {
  free(sip->a);
  sip->a = NULL;
  free(sip->b);
  sip->b = NULL;
  free(sip->ap);
  sip->ap = NULL;
  free(sip->bp);
  sip->bp = NULL;
  free(sip->scratch);
  sip->scratch = NULL;
  free(sip->err);
  sip->err = NULL;
}

static INLINE double
lu(
    const unsigned int order,
    const double* const matrix,
    const int x,
    const int y) {

  int index;
  assert(x >= 0 && x <= (int)order);
  assert(y >= 0 && y <= (int)order);

  index = x * ((int)order + 1) + y;
  assert(index >= 0 && index < ((int)order + 1) * ((int)order + 1));

  return matrix[index];
}

static int
sip_compute(
    /*@unused@*/ const unsigned int naxes,
    const unsigned int nelem,
    const unsigned int m,
    /*@null@*/ const double* a,
    const unsigned int n,
    /*@null@*/ const double* b,
    const double* crpix /* [2] */,
    /*@null@*/ double* tmp,
    /*@null@*/ const double* input /* [NAXES][nelem] */,
    /*@null@*/ double* output /* [NAXES][nelem] */) {

  unsigned int  i;
  int           j, k;
  double        x, y;
  double        sum;
  const double* input_ptr;
  double*       output_ptr;

  assert(a != NULL);
  assert(b != NULL);
  assert(crpix != NULL);
  assert(tmp != NULL);
  assert(input != NULL);
  assert(output != NULL);

  /* Avoid segfaults */
  if (input == NULL || output == NULL || tmp == NULL || crpix == NULL) {
    return 1;
  }

  /* If we have one, we must have both... */
  if ((a == NULL) ^ (b == NULL)) {
    return 6;
  }

  /* If no distortion, just return values */
  if (a == NULL /* && b == NULL ... implied */) {
    return 0;
  }

  input_ptr = input;
  output_ptr = output;
  for (i = 0; i < nelem; ++i) {
    x = *input_ptr++ - crpix[0];
    y = *input_ptr++ - crpix[1];

    for (j = 0; j <= (int)m; ++j) {
      tmp[j] = lu(m, a, (int)m-j, j);
      for (k = j-1; k >= 0; --k) {
        tmp[j] = (y * tmp[j]) + lu(m, a, (int)m-j, k);
      }
    }

    sum = tmp[0];
    for (j = (int)m; j > 0; --j) {
      sum = x * sum + tmp[(int)m - j + 1];
    }
    *output_ptr++ += sum;

    for (j = 0; j <= (int)n; ++j) {
      tmp[j] = lu(n, b, (int)n-j, j);
      for (k = j-1; k >= 0; --k) {
          tmp[j] = (y * tmp[j]) + lu(n, b, (int)n-j, k);
      }
    }

    sum = tmp[0];
    for (j = (int)n; j > 0; --j) {
      sum = x * sum + tmp[n - j + 1];
    }
    *output_ptr++ += sum;
  }

  return 0;
}

int
sip_pix2deltas(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* pix /* [NAXES][nelem] */,
    double* deltas /* [NAXES][nelem] */) {

  if (sip == NULL) {
    return 1;
  }

  return sip_compute(naxes, nelem,
                     sip->a_order, sip->a,
                     sip->b_order, sip->b,
                     sip->crpix,
                     (double *)sip->scratch,
                     pix, deltas);
}

int
sip_foc2deltas(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* foc /* [NAXES][nelem] */,
    double* deltas /* [NAXES][nelem] */) {

  if (sip == NULL) {
    return 1;
  }

  return sip_compute(naxes, nelem,
                     sip->ap_order, sip->ap,
                     sip->bp_order, sip->bp,
                     sip->crpix,
                     (double *)sip->scratch,
                     foc, deltas);
}

int
sip_pix2foc(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* pix /* [NAXES][nelem] */,
    double* foc /* [NAXES][nelem] */) {
  assert(pix);
  assert(foc);

  if (pix != foc) {
      memcpy(foc, pix, sizeof(double) * naxes * nelem);
  }

  return sip_pix2deltas(sip, naxes, nelem, pix, foc);
}

int
sip_foc2pix(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* foc /* [NAXES][nelem] */,
    double* pix /* [NAXES][nelem] */) {
  assert(pix);
  assert(foc);

  if (pix != foc) {
      memcpy(pix, foc, sizeof(double) * naxes * nelem);
  }

  return sip_foc2deltas(sip, naxes, nelem, foc, pix);
}
