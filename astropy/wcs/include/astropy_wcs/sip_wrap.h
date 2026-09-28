/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __SIP_WRAP_H__
#define __SIP_WRAP_H__

#include "pyutil.h"
#include "sip.h"

extern PyObject* SipType;

typedef struct {
#ifndef _Py_OPAQUE_PYOBJECT
  PyObject_HEAD
#endif
  sip_t x;
} Sip;

int
_setup_sip_type(
    PyObject* m);

#endif
