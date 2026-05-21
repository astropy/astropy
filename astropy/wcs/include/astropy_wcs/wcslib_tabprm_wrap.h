/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __WCSLIB_TABPRM_WRAP_H__
#define __WCSLIB_TABPRM_WRAP_H__

#include "pyutil.h"
#include "wcs.h"

extern PyObject* TabprmType;

typedef struct {
  PyObject_HEAD
  struct tabprm* x;
  PyObject* owner;
} Tabprm;

Tabprm*
Tabprm_cnew(PyObject* wcsprm, struct tabprm* x);

int _setup_tabprm_type(PyObject* m);

#endif
