/*
 Author: Michael Droettboom
*/

#ifndef __WCSLIB_WRAP_H__
#define __WCSLIB_WRAP_H__

#include "pyutil.h"

extern PyObject* PyWcsprmType;

typedef struct {

  PyObject_HEAD
  struct wcsprm x;

  // WCSLIB converts units to SI (for example nanometers for a specrtal axis to
  // meters, or arcseconds for a celestial frame to degrees). There is no easy
  // way to turn off this behavior in WCSLIB itself, but since we expose WCSLIB
  // through a wrapper layer, we can optionally do conversions on-the-fly to
  // make it look to the user as if the units are the original ones. This is
  // controlled via the preserve_units option which can be 0 (default WCSLIB
  // behavior) or 1 (preserve original units). The original_cunit array is used
  // to store the original CUNIT values, while the unit_scaling array contains
  // the multiplicative scaling required to convert the units from the original
  // units (e.g. arcsec) to the internal WCSLIB units (e.g. deg)
  int    preserve_units;
  char   (*original_cunit)[72];
  double *unit_scaling;

} PyWcsprm;

int _setup_wcsprm_type(PyObject* m);

PyObject*
PyWcsprm_find_all_wcs(
    PyObject* self,
    PyObject* args,
    PyObject* kwds);

int _update_wtbarr_from_hdulist(PyObject *hdulist, struct wtbarr *wtb);

void _set_wtbarr_callback(PyObject* callback);

int PyWcsprm_cset(PyWcsprm* self, const int convert);

#endif
