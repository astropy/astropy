/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

#include "astropy_wcs/util.h"
#include <math.h>
#include <float.h>

void set_invalid_to_nan(
    const int ncoord,
    const int nelem,
    double* const data,
    const int* const stat)
{
  int i = 0;
  double* d = data;
  const int* s = stat;
  const int* s_end = stat + ncoord;
  double n;

  #ifndef NAN
    #define INF (DBL_MAX+DBL_MAX)
    #define NAN (INF-INF)
  #endif

  // Note that stat is a bit mask, so we need to mask only some of
  // the coordinates depending on the bit mask values.

  n = NAN;

  for ( ; s != s_end; ++s) {
    if (*s) {
      int bit = 1;
      for (i = 0; i < nelem; ++i) {
        if (*s & bit) {
          *d = n;
        }
        d++;
        // We don't need to worry about overflow here because the WCS
        // class cannot be used for naxis > 15 so nelem will always
        // be <=15.
        bit <<= 1;
      }
    } else {
      d += nelem;
    }
  }
}
