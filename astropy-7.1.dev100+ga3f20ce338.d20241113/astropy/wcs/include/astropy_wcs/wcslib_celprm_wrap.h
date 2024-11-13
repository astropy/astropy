#ifndef __WCSLIB_CELPRM_WRAP_H__
#define __WCSLIB_CELPRM_WRAP_H__

#include "pyutil.h"
#include "wcs.h"

extern PyTypeObject PyCelprmType;

typedef struct {
    PyObject_HEAD
    struct celprm* x;
    int* prefcount;
    PyObject* owner;
} PyCelprm;

PyCelprm* PyCelprm_cnew(PyObject* wcsprm_obj, struct celprm* x, int* prefcount);

int _setup_celprm_type(PyObject* m);

#endif
