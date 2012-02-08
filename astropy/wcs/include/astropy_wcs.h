/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __ASTROPY_WCS_H__
#define __ASTROPY_WCS_H__

/* util.h must be imported first */
#include "pyutil.h"
#include "pipeline.h"

typedef struct {
  PyObject_HEAD
  pipeline_t x;
  /*@shared@*/ PyObject*            py_det2im[2];
  /*@null@*/ /*@shared@*/ PyObject* py_sip;
  /*@shared@*/ PyObject*            py_distortion_lookup[2];
  /*@null@*/ /*@shared@*/ PyObject* py_wcsprm;
} Wcs;

#endif /* __ASTROPY_WCS_H__ */
