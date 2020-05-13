#ifndef __WCSLIB_AUXPRM_WRAP_H__
#define __WCSLIB_AUXPRM_WRAP_H__

#include "pyutil.h"
#include "wcs.h"

extern PyTypeObject PyAuxprmType;

typedef struct {
  PyObject_HEAD
  struct auxprm* x;
  PyObject* owner;
} PyAuxprm;

PyAuxprm*
PyAuxprm_cnew(PyObject* wcsprm, struct auxprm* x);

int _setup_auxprm_type(PyObject* m);

#endif
