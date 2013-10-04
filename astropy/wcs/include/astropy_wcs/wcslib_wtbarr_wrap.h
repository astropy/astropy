/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __WCSLIB_WTBARR_WRAP_H__
#define __WCSLIB_WTBARR_WRAP_H__

#include "pyutil.h"
#include "wcs.h"

extern PyTypeObject PyWtbarrType;

typedef struct {
  PyObject_HEAD
  struct wtbarr* x;
  PyObject* owner;
} PyWtbarr;

PyWtbarr*
PyWtbarr_cnew(PyObject* wcsprm, struct wtbarr* x);

int _setup_wtbarr_type(PyObject* m);

#endif
