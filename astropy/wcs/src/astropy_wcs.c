/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#include "astropy_wcs/astropy_wcs.h"
#include "astropy_wcs/wcslib_wrap.h"
#include "astropy_wcs/wcslib_tabprm_wrap.h"
#include "astropy_wcs/wcslib_auxprm_wrap.h"
#include "astropy_wcs/wcslib_prjprm_wrap.h"
#include "astropy_wcs/wcslib_celprm_wrap.h"
#include "astropy_wcs/wcslib_units_wrap.h"
#include "astropy_wcs/wcslib_wtbarr_wrap.h"
#include "astropy_wcs/distortion_wrap.h"
#include "astropy_wcs/sip_wrap.h"
#include "astropy_wcs/docstrings.h"
#include "astropy_wcs/astropy_wcs_api.h"
#include "astropy_wcs/unit_list_proxy.h"

#include <structmember.h> /* from Python */

#include <stdlib.h>
#include <time.h>

#include <tab.h>
#include <wtbarr.h>

/***************************************************************************
 * Wcs type
 ***************************************************************************/

static PyTypeObject WcsType;

static int _setup_wcs_type(PyObject* m);


PyObject* PyWcsprm_set_wtbarr_fitsio_callback(PyObject *dummy, PyObject *args) {
    PyObject *callback;

    if (PyArg_ParseTuple(args, "O:set_wtbarr_fitsio_callback", &callback)) {
        if (!PyCallable_Check(callback)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be callable");
            return NULL;
        }
        _set_wtbarr_callback(callback);

        Py_RETURN_NONE;
    }
    return NULL;
}


/***************************************************************************
 * PyWcs methods
 */

static int
Wcs_traverse(
    Wcs* self,
    visitproc visit,
    void* arg) {

  Py_VISIT(self->py_det2im[0]);
  Py_VISIT(self->py_det2im[1]);
  Py_VISIT(self->py_sip);
  Py_VISIT(self->py_distortion_lookup[0]);
  Py_VISIT(self->py_distortion_lookup[1]);
  Py_VISIT(self->py_wcsprm);

  return 0;
}

static int
Wcs_clear(
    Wcs* self) {

  Py_CLEAR(self->py_det2im[0]);
  Py_CLEAR(self->py_det2im[1]);
  Py_CLEAR(self->py_sip);
  Py_CLEAR(self->py_distortion_lookup[0]);
  Py_CLEAR(self->py_distortion_lookup[1]);
  Py_CLEAR(self->py_wcsprm);

  return 0;
}

static void
Wcs_dealloc(
    Wcs* self) {

  PyObject_GC_UnTrack(self);
  Wcs_clear(self);
  pipeline_free(&self->x);
  Py_TYPE(self)->tp_free((PyObject*)self);
}

/*@null@*/ static PyObject *
Wcs_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  Wcs* self;
  self = (Wcs*)type->tp_alloc(type, 0);
  if (self != NULL) {
    pipeline_clear(&self->x);
    self->py_det2im[0]            = NULL;
    self->py_det2im[1]            = NULL;
    self->py_sip                  = NULL;
    self->py_distortion_lookup[0] = NULL;
    self->py_distortion_lookup[1] = NULL;
    self->py_wcsprm               = NULL;
  }
  return (PyObject*)self;
}

static int
Wcs_init(
    Wcs* self,
    PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  size_t       i;
  PyObject*    py_sip;
  PyObject*    py_wcsprm;
  PyObject*    py_distortion_lookup[2];
  PyObject*    py_det2im[2];

  if (!PyArg_ParseTuple
      (args, "O(OO)O(OO):Wcs.__init__",
       &py_sip,
       &py_distortion_lookup[0],
       &py_distortion_lookup[1],
       &py_wcsprm,
       &py_det2im[0],
       &py_det2im[1])) {
    return -1;
  }

  /* Check and set Distortion lookup tables */
  for (i = 0; i < 2; ++i) {
    if (py_det2im[i] != NULL && py_det2im[i] != Py_None) {
      if (!PyObject_TypeCheck(py_det2im[i], &PyDistLookupType)) {
        PyErr_SetString(PyExc_TypeError,
                        "Arg 4 must be a pair of DistortionLookupTable or None objects");
        return -1;
      }

      Py_CLEAR(self->py_det2im[i]);
      self->py_det2im[i] = py_det2im[i];
      Py_INCREF(py_det2im[i]);
      self->x.det2im[i] = &(((PyDistLookup*)py_det2im[i])->x);
    }
  }

  /* Check and set SIP */
  if (py_sip != NULL && py_sip != Py_None) {
    if (!PyObject_TypeCheck(py_sip, &PySipType)) {
      PyErr_SetString(PyExc_TypeError,
                      "Arg 1 must be Sip object");
      return -1;
    }

    Py_CLEAR(self->py_sip);
    self->py_sip = py_sip;
    Py_INCREF(py_sip);
    self->x.sip = &(((PySip*)py_sip)->x);
  }

  /* Check and set Distortion lookup tables */
  for (i = 0; i < 2; ++i) {
    if (py_distortion_lookup[i] != NULL && py_distortion_lookup[i] != Py_None) {
      if (!PyObject_TypeCheck(py_distortion_lookup[i], &PyDistLookupType)) {
        PyErr_SetString(PyExc_TypeError,
                        "Arg 2 must be a pair of DistortionLookupTable or None objects");
        return -1;
      }

      Py_CLEAR(self->py_distortion_lookup[i]);
      self->py_distortion_lookup[i] = py_distortion_lookup[i];
      Py_INCREF(py_distortion_lookup[i]);
      self->x.cpdis[i] = &(((PyDistLookup*)py_distortion_lookup[i])->x);
    }
  }

  /* Set and lookup Wcsprm object */
  if (py_wcsprm != NULL && py_wcsprm != Py_None) {
    if (!PyObject_TypeCheck(py_wcsprm, &PyWcsprmType)) {
      PyErr_SetString(PyExc_TypeError,
                      "Arg 3 must be Wcsprm object");
      return -1;
    }

    Py_CLEAR(self->py_wcsprm);
    self->py_wcsprm = py_wcsprm;
    Py_INCREF(py_wcsprm);
    self->x.wcs = &(((PyWcsprm*)py_wcsprm)->x);
  }

  return 0;
}

/*@null@*/ static PyObject*
Wcs_all_pix2world(
    Wcs* self,
    PyObject* args,
    PyObject* kwds) {

  int            naxis      = 2;
  PyObject*      pixcrd_obj = NULL;
  int            origin     = 1;
  PyArrayObject* pixcrd     = NULL;
  PyArrayObject* world      = NULL;
  int            status     = -1;
  const char*    keywords[] = {
    "pixcrd", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "Oi:all_pix2world", (char **)keywords,
          &pixcrd_obj, &origin)) {
    return NULL;
  }

  naxis = self->x.wcs->naxis;

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny(pixcrd_obj, NPY_DOUBLE, 2, 2);
  if (pixcrd == NULL) {
    return NULL;
  }

  if (PyArray_DIM(pixcrd, 1) < naxis) {
    PyErr_Format(
      PyExc_RuntimeError,
      "Input array must be 2-dimensional, where the second dimension >= %d",
      naxis);
    goto exit;
  }

  world = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(pixcrd), NPY_DOUBLE);
  if (world == NULL) {
    goto exit;
  }

  /* Make the call */
  Py_BEGIN_ALLOW_THREADS
  preoffset_array(pixcrd, origin);
  wcsprm_python2c(self->x.wcs);
  status = pipeline_all_pixel2world(&self->x,
                                    (unsigned int)PyArray_DIM(pixcrd, 0),
                                    (unsigned int)PyArray_DIM(pixcrd, 1),
                                    (double*)PyArray_DATA(pixcrd),
                                    (double*)PyArray_DATA(world));
  wcsprm_c2python(self->x.wcs);
  unoffset_array(pixcrd, origin);
  Py_END_ALLOW_THREADS
  /* unoffset_array(world, origin); */

 exit:
  Py_XDECREF(pixcrd);

  if (status == 0 || status == 8) {
    return (PyObject*)world;
  } else {
    Py_XDECREF(world);
    if (status == -1) {
      PyErr_SetString(
        PyExc_ValueError,
        "Wrong number of dimensions in input array.  Expected 2.");
      return NULL;
    } else {
      if (status == -1) {
        /* exception already set */
        return NULL;
      } else {
        wcserr_to_python_exc(self->x.err);
        return NULL;
      }
    }
  }
}

/*@null@*/ static PyObject*
Wcs_p4_pix2foc(
    Wcs* self,
    PyObject* args,
    PyObject* kwds) {

  PyObject*      pixcrd_obj = NULL;
  int            origin     = 1;
  PyArrayObject* pixcrd     = NULL;
  PyArrayObject* foccrd     = NULL;
  int            status     = -1;
  const char*    keywords[] = {
    "pixcrd", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oi:p4_pix2foc", (char **)keywords,
                                   &pixcrd_obj, &origin)) {
    return NULL;
  }

  if (self->x.cpdis[0] == NULL && self->x.cpdis[1] == NULL) {
    Py_INCREF(pixcrd_obj);
    return pixcrd_obj;
  }

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny(pixcrd_obj, NPY_DOUBLE, 2, 2);
  if (pixcrd == NULL) {
    return NULL;
  }

  if (PyArray_DIM(pixcrd, 1) != NAXES) {
    PyErr_SetString(PyExc_ValueError, "Pixel array must be an Nx2 array");
    goto exit;
  }

  foccrd = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(pixcrd), NPY_DOUBLE);
  if (foccrd == NULL) {
    status = 2;
    goto exit;
  }

  Py_BEGIN_ALLOW_THREADS
  preoffset_array(pixcrd, origin);
  status = p4_pix2foc(2, (void *)self->x.cpdis,
                      (unsigned int)PyArray_DIM(pixcrd, 0),
                      (double*)PyArray_DATA(pixcrd),
                      (double*)PyArray_DATA(foccrd));
  unoffset_array(pixcrd, origin);
  unoffset_array(foccrd, origin);
  Py_END_ALLOW_THREADS

 exit:

  Py_XDECREF(pixcrd);

  if (status == 0) {
    return (PyObject*)foccrd;
  } else {
    Py_XDECREF(foccrd);
    if (status == -1) {
      /* Exception already set */
      return NULL;
    } else {
      PyErr_SetString(PyExc_MemoryError, "NULL pointer passed");
      return NULL;
    }
  }
}

/*@null@*/ static PyObject*
Wcs_det2im(
    Wcs* self,
    PyObject* args,
    PyObject* kwds) {

  PyObject*      detcrd_obj = NULL;
  int            origin     = 1;
  PyArrayObject* detcrd     = NULL;
  PyArrayObject* imcrd     = NULL;
  int            status     = -1;
  const char*    keywords[] = {
    "detcrd", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oi:det2im", (char **)keywords,
                                   &detcrd_obj, &origin)) {
    return NULL;
  }

  if (self->x.det2im[0] == NULL && self->x.det2im[1] == NULL) {
    Py_INCREF(detcrd_obj);
    return detcrd_obj;
  }

  detcrd = (PyArrayObject*)PyArray_ContiguousFromAny(detcrd_obj, NPY_DOUBLE, 2, 2);
  if (detcrd == NULL) {
    return NULL;
  }

  if (PyArray_DIM(detcrd, 1) != NAXES) {
    PyErr_SetString(PyExc_ValueError, "Pixel array must be an Nx2 array");
    goto exit;
  }

  imcrd = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(detcrd), NPY_DOUBLE);
  if (imcrd == NULL) {
    status = 2;
    goto exit;
  }

  Py_BEGIN_ALLOW_THREADS
  preoffset_array(detcrd, origin);
  status = p4_pix2foc(2, (void *)self->x.det2im,
                      (unsigned int)PyArray_DIM(detcrd, 0),
                      (double*)PyArray_DATA(detcrd),
                      (double*)PyArray_DATA(imcrd));
  unoffset_array(detcrd, origin);
  unoffset_array(imcrd, origin);
  Py_END_ALLOW_THREADS

 exit:

  Py_XDECREF(detcrd);

  if (status == 0) {
    return (PyObject*)imcrd;
  } else {
    Py_XDECREF(imcrd);
    if (status == -1) {
      /* Exception already set */
      return NULL;
    } else {
      PyErr_SetString(PyExc_MemoryError, "NULL pointer passed");
      return NULL;
    }
  }
}

/*@null@*/ static PyObject*
Wcs_pix2foc(
    Wcs* self,
    PyObject* args,
    PyObject* kwds) {

  PyObject*      pixcrd_obj = NULL;
  int            origin     = 1;
  PyArrayObject* pixcrd     = NULL;
  PyArrayObject* foccrd     = NULL;
  int            status     = -1;
  const char*    keywords[] = {
    "pixcrd", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oi:pix2foc", (char **)keywords,
                                   &pixcrd_obj, &origin)) {
    return NULL;
  }

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny(pixcrd_obj, NPY_DOUBLE, 2, 2);
  if (pixcrd == NULL) {
    return NULL;
  }

  if (PyArray_DIM(pixcrd, 1) != NAXES) {
    PyErr_SetString(PyExc_ValueError, "Pixel array must be an Nx2 array");
    goto _exit;
  }

  foccrd = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(pixcrd), NPY_DOUBLE);
  if (foccrd == NULL) {
    goto _exit;
  }

  Py_BEGIN_ALLOW_THREADS
  preoffset_array(pixcrd, origin);
  status = pipeline_pix2foc(&self->x,
                            (unsigned int)PyArray_DIM(pixcrd, 0),
                            (unsigned int)PyArray_DIM(pixcrd, 1),
                            (double*)PyArray_DATA(pixcrd),
                            (double*)PyArray_DATA(foccrd));
  unoffset_array(pixcrd, origin);
  unoffset_array(foccrd, origin);
  Py_END_ALLOW_THREADS

 _exit:

  Py_XDECREF(pixcrd);

  if (status == 0) {
    return (PyObject*)foccrd;
  } else {
    Py_XDECREF(foccrd);
    if (status == -1) {
      /* Exception already set */
      return NULL;
    } else {
      wcserr_to_python_exc(self->x.err);
      return NULL;
    }
  }
}

/*@null@*/ static PyObject*
Wcs_get_wcs(
    Wcs* self,
    /*@unused@*/ void* closure) {

  if (self->py_wcsprm) {
    Py_INCREF(self->py_wcsprm);
    return self->py_wcsprm;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static int
Wcs_set_wcs(
    Wcs* self,
    /*@shared@*/ PyObject* value,
    /*@unused@*/ void* closure) {

  Py_CLEAR(self->py_wcsprm);
  self->x.wcs = NULL;

  if (value != NULL && value != Py_None) {
    if (!PyObject_TypeCheck(value, &PyWcsprmType)) {
      PyErr_SetString(PyExc_TypeError,
                      "wcs must be Wcsprm object");
      return -1;
    }

    Py_INCREF(value);
    self->py_wcsprm = value;
    self->x.wcs = &(((PyWcsprm*)value)->x);
  }

  return 0;
}

static PyObject*
Wcs_get_cpdis1(
    Wcs* self,
    /*@unused@*/ void* closure) {

  if (self->py_distortion_lookup[0]) {
    Py_INCREF(self->py_distortion_lookup[0]);
    return self->py_distortion_lookup[0];
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static int
Wcs_set_cpdis1(
    Wcs* self,
    /*@shared@*/ PyObject* value,
    /*@unused@*/ void* closure) {

  Py_CLEAR(self->py_distortion_lookup[0]);
  self->x.cpdis[0] = NULL;

  if (value != NULL && value != Py_None) {
    if (!PyObject_TypeCheck(value, &PyDistLookupType)) {
      PyErr_SetString(PyExc_TypeError,
                      "cpdis1 must be DistortionLookupTable object");
      return -1;
    }

    Py_INCREF(value);
    self->py_distortion_lookup[0] = value;
    self->x.cpdis[0] = &(((PyDistLookup*)value)->x);
  }

  return 0;
}

/*@shared@*/ static PyObject*
Wcs_get_cpdis2(
    Wcs* self,
    /*@unused@*/ void* closure) {

  if (self->py_distortion_lookup[1]) {
    Py_INCREF(self->py_distortion_lookup[1]);
    return self->py_distortion_lookup[1];
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static int
Wcs_set_cpdis2(
    Wcs* self,
    /*@shared@*/ PyObject* value,
    /*@unused@*/ void* closure) {

  Py_CLEAR(self->py_distortion_lookup[1]);
  self->x.cpdis[1] = NULL;

  if (value != NULL && value != Py_None) {
    if (!PyObject_TypeCheck(value, &PyDistLookupType)) {
      PyErr_SetString(PyExc_TypeError,
                      "cpdis2 must be DistortionLookupTable object");
      return -1;
    }

    Py_INCREF(value);
    self->py_distortion_lookup[1] = value;
    self->x.cpdis[1] = &(((PyDistLookup*)value)->x);
  }

  return 0;
}

static PyObject*
Wcs_get_det2im1(
    Wcs* self,
    /*@unused@*/ void* closure) {

  if (self->py_det2im[0]) {
    Py_INCREF(self->py_det2im[0]);
    return self->py_det2im[0];
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static int
Wcs_set_det2im1(
    Wcs* self,
    /*@shared@*/ PyObject* value,
    /*@unused@*/ void* closure) {

  Py_CLEAR(self->py_det2im[0]);
  self->x.det2im[0] = NULL;

  if (value != NULL && value != Py_None) {
    if (!PyObject_TypeCheck(value, &PyDistLookupType)) {
      PyErr_SetString(PyExc_TypeError,
                      "det2im1 must be DistortionLookupTable object");
      return -1;
    }

    Py_INCREF(value);
    self->py_det2im[0] = value;
    self->x.det2im[0] = &(((PyDistLookup*)value)->x);
  }

  return 0;
}

/*@shared@*/ static PyObject*
Wcs_get_det2im2(
    Wcs* self,
    /*@unused@*/ void* closure) {

  if (self->py_det2im[1]) {
    Py_INCREF(self->py_det2im[1]);
    return self->py_det2im[1];
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static int
Wcs_set_det2im2(
    Wcs* self,
    /*@shared@*/ PyObject* value,
    /*@unused@*/ void* closure) {

  Py_CLEAR(self->py_det2im[1]);
  self->x.det2im[1] = NULL;

  if (value != NULL && value != Py_None) {
    if (!PyObject_TypeCheck(value, &PyDistLookupType)) {
      PyErr_SetString(PyExc_TypeError,
                      "det2im2 must be DistortionLookupTable object");
      return -1;
    }

    Py_INCREF(value);
    self->py_det2im[1] = value;
    self->x.det2im[1] = &(((PyDistLookup*)value)->x);
  }

  return 0;
}

/*@shared@*/ static PyObject*
Wcs_get_sip(
    Wcs* self,
    /*@unused@*/ void* closure) {

  if (self->py_sip) {
    Py_INCREF(self->py_sip);
    return self->py_sip;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static int
Wcs_set_sip(
    Wcs* self,
    /*@shared@*/ PyObject* value,
    /*@unused@*/ void* closure) {

  Py_CLEAR(self->py_sip);
  self->x.sip = NULL;

  if (value != NULL && value != Py_None) {
    if (!PyObject_TypeCheck(value, &PySipType)) {
      PyErr_SetString(PyExc_TypeError,
                      "sip must be Sip object");
      return -1;
    }

    Py_INCREF(value);
    self->py_sip = value;
    self->x.sip = &(((PySip*)value)->x);
  }

  return 0;
}

static PyObject*
_sanity_check(
    PyObject* self,
    PyObject* args,
    PyObject* kwds) {

  if (sizeof(WCSLIB_INT64) != 8) {
    Py_INCREF(Py_False);
    return Py_False;
  }

  Py_INCREF(Py_True);
  return Py_True;
}

/***************************************************************************
 * Wcs definition structures
 */

static PyGetSetDef Wcs_getset[] = {
  {"det2im1", (getter)Wcs_get_det2im1, (setter)Wcs_set_det2im1, (char *)doc_det2im1},
  {"det2im2", (getter)Wcs_get_det2im2, (setter)Wcs_set_det2im2, (char *)doc_det2im2},
  {"cpdis1", (getter)Wcs_get_cpdis1, (setter)Wcs_set_cpdis1, (char *)doc_cpdis1},
  {"cpdis2", (getter)Wcs_get_cpdis2, (setter)Wcs_set_cpdis2, (char *)doc_cpdis2},
  {"sip", (getter)Wcs_get_sip, (setter)Wcs_set_sip, (char *)doc_sip},
  {"wcs", (getter)Wcs_get_wcs, (setter)Wcs_set_wcs, (char *)doc_wcs},
  {NULL}
};

static PyMethodDef Wcs_methods[] = {
  {"_all_pix2world", (PyCFunction)Wcs_all_pix2world, METH_VARARGS|METH_KEYWORDS, doc_all_pix2world},
  {"_det2im", (PyCFunction)Wcs_det2im, METH_VARARGS|METH_KEYWORDS, doc_det2im},
  {"_p4_pix2foc", (PyCFunction)Wcs_p4_pix2foc, METH_VARARGS|METH_KEYWORDS, doc_p4_pix2foc},
  {"_pix2foc", (PyCFunction)Wcs_pix2foc, METH_VARARGS|METH_KEYWORDS, doc_pix2foc},
  {NULL}
};

static PyMethodDef module_methods[] = {
  {"_sanity_check", (PyCFunction)_sanity_check, METH_NOARGS, ""},
  {"find_all_wcs", (PyCFunction)PyWcsprm_find_all_wcs, METH_VARARGS|METH_KEYWORDS, doc_find_all_wcs},
  {"set_wtbarr_fitsio_callback", (PyCFunction)PyWcsprm_set_wtbarr_fitsio_callback, METH_VARARGS, NULL},
  {NULL}  /* Sentinel */
};

static PyTypeObject WcsType = {
  PyVarObject_HEAD_INIT(NULL, 0)
  "astropy.wcs.WCSBase",                 /*tp_name*/
  sizeof(Wcs),                /*tp_basicsize*/
  0,                            /*tp_itemsize*/
  (destructor)Wcs_dealloc,    /*tp_dealloc*/
  0,                            /*tp_print*/
  0,                            /*tp_getattr*/
  0,                            /*tp_setattr*/
  0,                            /*tp_compare*/
  0,                            /*tp_repr*/
  0,                            /*tp_as_number*/
  0,                            /*tp_as_sequence*/
  0,                            /*tp_as_mapping*/
  0,                            /*tp_hash */
  0,                            /*tp_call*/
  0,                            /*tp_str*/
  0,                            /*tp_getattro*/
  0,                            /*tp_setattro*/
  0,                            /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
  doc_Wcs,                      /* tp_doc */
  (traverseproc)Wcs_traverse, /* tp_traverse */
  (inquiry)Wcs_clear,         /* tp_clear */
  0,                            /* tp_richcompare */
  0,                            /* tp_weaklistoffset */
  0,                            /* tp_iter */
  0,                            /* tp_iternext */
  Wcs_methods,                /* tp_methods */
  0,                            /* tp_members */
  Wcs_getset,                 /* tp_getset */
  0,                            /* tp_base */
  0,                            /* tp_dict */
  0,                            /* tp_descr_get */
  0,                            /* tp_descr_set */
  0,                            /* tp_dictoffset */
  (initproc)Wcs_init,         /* tp_init */
  0,                            /* tp_alloc */
  Wcs_new,                    /* tp_new */
};


/***************************************************************************
 * Module-level
 ***************************************************************************/

int _setup_wcs_type(
    PyObject* m) {

  if (PyType_Ready(&WcsType) < 0)
    return -1;

  Py_INCREF(&WcsType);
  return PyModule_AddObject(m, "_Wcs", (PyObject *)&WcsType);
}

struct module_state {
/* The Sun compiler can't handle empty structs */
#if defined(__SUNPRO_C) || defined(_MSC_VER)
    int _dummy;
#endif
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_wcs",
    NULL,
    sizeof(struct module_state),
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__wcs(void)

{
  PyObject* m;

  wcs_errexc[0] = NULL;                         /* Success */
  wcs_errexc[1] = &PyExc_MemoryError;           /* Null wcsprm pointer passed */
  wcs_errexc[2] = &PyExc_MemoryError;           /* Memory allocation failed */
  wcs_errexc[3] = &WcsExc_SingularMatrix;       /* Linear transformation matrix is singular */
  wcs_errexc[4] = &WcsExc_InconsistentAxisTypes; /* Inconsistent or unrecognized coordinate axis types */
  wcs_errexc[5] = &PyExc_ValueError;            /* Invalid parameter value */
  wcs_errexc[6] = &WcsExc_InvalidTransform;     /* Invalid coordinate transformation parameters */
  wcs_errexc[7] = &WcsExc_InvalidTransform;     /* Ill-conditioned coordinate transformation parameters */
  wcs_errexc[8] = &WcsExc_InvalidCoordinate;    /* One or more of the pixel coordinates were invalid, */
  /* as indicated by the stat vector */
  wcs_errexc[9] = &WcsExc_InvalidCoordinate;    /* One or more of the world coordinates were invalid, */
  /* as indicated by the stat vector */
  wcs_errexc[10] = &WcsExc_InvalidCoordinate;    /* Invalid world coordinate */
  wcs_errexc[11] = &WcsExc_NoSolution;           /* no solution found in the specified interval */
  wcs_errexc[12] = &WcsExc_InvalidSubimageSpecification; /* Invalid subimage specification (no spectral axis) */
  wcs_errexc[13] = &WcsExc_NonseparableSubimageCoordinateSystem; /* Non-separable subimage coordinate system */

  m = PyModule_Create(&moduledef);

  if (m == NULL)
    return NULL;

  import_array();

  if (_setup_api(m)                 ||
      _setup_str_list_proxy_type(m) ||
      _setup_unit_list_proxy_type(m)||
      _setup_wcsprm_type(m)         ||
      _setup_auxprm_type(m)         ||
      _setup_prjprm_type(m)         ||
      _setup_celprm_type(m)         ||
      _setup_tabprm_type(m)         ||
      _setup_wtbarr_type(m)         ||
      _setup_distortion_type(m)     ||
      _setup_sip_type(m)            ||
      _setup_wcs_type(m)          ||
      _define_exceptions(m)) {
    Py_DECREF(m);
    return NULL;
  }

#ifdef HAVE_WCSLIB_VERSION
  if (PyModule_AddStringConstant(m, "__version__", wcslib_version(NULL))) {
    return NULL;
  }
#else
  if (PyModule_AddStringConstant(m, "__version__", "4.x")) {
    return NULL;
  }
#endif

  return m;
}
