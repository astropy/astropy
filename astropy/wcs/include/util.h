/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __UTIL_H__
#define __UTIL_H__

#ifdef __SUNPRO_C
#define inline
#endif

#ifdef _MSC_VER
#define inline __inline
#endif

#include <wcs.h>
#include <wcsmath.h>

#include "isnan.h"

#undef	CLAMP
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

void set_invalid_to_nan(
    const int ncoord,
    const int nelem,
    double* const data,
    const int* const stat);

#endif /* __UTIL_H__ */
