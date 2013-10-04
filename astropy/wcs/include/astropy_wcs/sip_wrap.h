/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __SIP_WRAP_H__
#define __SIP_WRAP_H__

#include "pyutil.h"
#include "sip.h"

extern PyTypeObject PySipType;

typedef struct {
  PyObject_HEAD
  sip_t x;
} PySip;

int
_setup_sip_type(
    PyObject* m);

#endif
