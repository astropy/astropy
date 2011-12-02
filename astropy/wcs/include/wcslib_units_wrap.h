/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __WCSLIB_UNITS_WRAP_H__
#define __WCSLIB_UNITS_WRAP_H__

#include "pyutil.h"
#include "wcsunits.h"

extern PyTypeObject PyUnitsType;

typedef struct {
  PyObject_HEAD
  char have[80];
  char want[80];
  double scale;
  double offset;
  double power;
} PyUnits;

PyUnits*
PyUnits_cnew(
    const char* const have,
    const char* const want,
    const double scale,
    const double offset,
    const double power);

int _setup_units_type(PyObject* m);

#endif
