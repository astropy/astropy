#ifndef __WCSLIB_PRJPRM_WRAP_H__
#define __WCSLIB_PRJPRM_WRAP_H__

#include "pyutil.h"
#include "wcs.h"

extern PyObject* PrjprmType;

typedef struct {
    PyObject_HEAD
    struct prjprm* x;
    int* prefcount;
    PyObject* owner;
} Prjprm;

Prjprm* Prjprm_cnew(PyObject* celprm, struct prjprm* x, int* prefcount);

int _setup_prjprm_type(PyObject* m);

#endif
