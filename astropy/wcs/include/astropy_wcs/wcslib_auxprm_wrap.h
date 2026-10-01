#ifndef __WCSLIB_AUXPRM_WRAP_H__
#define __WCSLIB_AUXPRM_WRAP_H__

#include "pyutil.h"
#include "wcs.h"

extern PyObject* AuxprmType;

typedef struct {
#ifndef _Py_OPAQUE_PYOBJECT
  PyObject_HEAD
#endif
  struct auxprm* x;
  PyObject* owner;
} Auxprm;

Auxprm*
Auxprm_cnew(PyObject* wcsprm, struct auxprm* x);

int _setup_auxprm_type(PyObject* m);

#endif
