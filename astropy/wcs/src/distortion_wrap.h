/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __DISTORTION_WRAP_H__
#define __DISTORTION_WRAP_H__

#include "pyutil.h"
#include "distortion.h"

extern PyTypeObject PyDistLookupType;

typedef struct {
  PyObject_HEAD
  distortion_lookup_t                    x;
  /*@null@*/ /*@shared@*/ PyArrayObject* py_data;
} PyDistLookup;

int
_setup_distortion_type(
    PyObject* m);

#endif
