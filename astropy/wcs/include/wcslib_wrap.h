/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __WCSLIB_WRAP_H__
#define __WCSLIB_WRAP_H__

#include "pyutil.h"

extern PyTypeObject PyWcsprmType;

typedef struct {
  PyObject_HEAD
  struct wcsprm x;
} PyWcsprm;

int _setup_wcsprm_type(PyObject* m);

PyObject*
PyWcsprm_find_all_wcs(
    PyObject* self,
    PyObject* args,
    PyObject* kwds);

#endif
