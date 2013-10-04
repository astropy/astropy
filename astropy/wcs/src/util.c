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

  n = NAN;

  for ( ; s != s_end; ++s) {
    if (*s) {
      for (i = 0; i < nelem; ++i) {
        *d++ = n;
      }
    } else {
      d += nelem;
    }
  }
}
