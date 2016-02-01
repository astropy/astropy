/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

#include "astropy_wcs/sip_wrap.h"
#include "astropy_wcs/docstrings.h"
#include "wcs.h"

static void
PySip_dealloc(
    PySip* self) {

  sip_free(&self->x);
  Py_TYPE(self)->tp_free((PyObject*)self);
}

/*@null@*/ static PyObject *
PySip_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PySip* self;

  self = (PySip*)type->tp_alloc(type, 0);
  if (self != NULL) {
    sip_clear(&self->x);
  }
  return (PyObject*)self;
}

static int
convert_matrix(
    /*@null@*/ PyObject* pyobj,
    PyArrayObject** array,
    double** data,
    unsigned int* order) {

  if (pyobj == Py_None) {
    *array = NULL;
    *data = NULL;
    *order = 0;
    return 0;
  }

  *array = (PyArrayObject*)PyArray_ContiguousFromAny(
      pyobj, NPY_DOUBLE, 2, 2);
  if (*array == NULL) {
    return -1;
  }

  if (PyArray_DIM(*array, 0) != PyArray_DIM(*array, 1)) {
    PyErr_SetString(PyExc_ValueError,
                    "Matrix must be square.");
    return -1;
  }

  *data = (double*)PyArray_DATA(*array);
  *order = (unsigned int)PyArray_DIM(*array, 0) - 1;

  return 0;
}

static int
PySip_init(
    PySip* self,
    PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyObject*      py_a     = NULL;
  PyObject*      py_b     = NULL;
  PyObject*      py_ap    = NULL;
  PyObject*      py_bp    = NULL;
  PyObject*      py_crpix = NULL;
  PyArrayObject* a        = NULL;
  PyArrayObject* b        = NULL;
  PyArrayObject* ap       = NULL;
  PyArrayObject* bp       = NULL;
  PyArrayObject* crpix    = NULL;
  double*        a_data   = NULL;
  double*        b_data   = NULL;
  double*        ap_data  = NULL;
  double*        bp_data  = NULL;
  unsigned int   a_order  = 0;
  unsigned int   b_order  = 0;
  unsigned int   ap_order = 0;
  unsigned int   bp_order = 0;
  int            status   = -1;

  if (!PyArg_ParseTuple(args, "OOOOO:Sip.__init__",
                        &py_a, &py_b, &py_ap, &py_bp, &py_crpix)) {
    return -1;
  }

  if (convert_matrix(py_a, &a, &a_data, &a_order) ||
      convert_matrix(py_b, &b, &b_data, &b_order) ||
      convert_matrix(py_ap, &ap, &ap_data, &ap_order) ||
      convert_matrix(py_bp, &bp, &bp_data, &bp_order)) {
    goto exit;
  }

  crpix = (PyArrayObject*)PyArray_ContiguousFromAny(py_crpix, NPY_DOUBLE,
                                                    1, 1);
  if (crpix == NULL) {
    goto exit;
  }

  if (PyArray_DIM(crpix, 0) != 2) {
    PyErr_SetString(PyExc_ValueError, "CRPIX wrong length");
    goto exit;
  }

  status = sip_init(&self->x,
                    a_order, a_data,
                    b_order, b_data,
                    ap_order, ap_data,
                    bp_order, bp_data,
                    PyArray_DATA(crpix));

 exit:
  Py_XDECREF(a);
  Py_XDECREF(b);
  Py_XDECREF(ap);
  Py_XDECREF(bp);
  Py_XDECREF(crpix);

  if (status == 0) {
    return 0;
  } else if (status == -1) {
    /* Exception already set */
    return -1;
  } else {
    wcserr_to_python_exc(self->x.err);
    return -1;
  }
}

/*@null@*/ static PyObject*
PySip_pix2foc(
    PySip* self,
    PyObject* args,
    PyObject* kwds) {

  PyObject*      pixcrd_obj = NULL;
  int            origin     = 1;
  PyArrayObject* pixcrd     = NULL;
  PyArrayObject* foccrd     = NULL;
  double*        foccrd_data = NULL;
  unsigned int   nelem      = 0;
  unsigned int   i, j;
  int            status     = -1;
  const char*    keywords[] = {
    "pixcrd", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oi:pix2foc", (char **)keywords,
                                   &pixcrd_obj, &origin)) {
    return NULL;
  }

  if (self->x.a == NULL || self->x.b == NULL) {
    PyErr_SetString(
        PyExc_ValueError,
        "SIP object does not have coefficients for pix2foc transformation (A and B)");
    return NULL;
  }

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny(pixcrd_obj, NPY_DOUBLE, 2, 2);
  if (pixcrd == NULL) {
    goto exit;
  }

  if (PyArray_DIM(pixcrd, 1) != 2) {
    PyErr_SetString(PyExc_ValueError, "Pixel array must be an Nx2 array");
    goto exit;
  }

  foccrd = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(pixcrd),
                                             NPY_DOUBLE);
  if (foccrd == NULL) {
    goto exit;
  }

  Py_BEGIN_ALLOW_THREADS
  preoffset_array(pixcrd, origin);
  status = sip_pix2foc(&self->x,
                       (unsigned int)PyArray_DIM(pixcrd, 1),
                       (unsigned int)PyArray_DIM(pixcrd, 0),
                       (const double*)PyArray_DATA(pixcrd),
                       (double*)PyArray_DATA(foccrd));
  unoffset_array(pixcrd, origin);

  /* Adjust for crpix */
  foccrd_data = (double *)PyArray_DATA(foccrd);
  nelem = (unsigned int)PyArray_DIM(foccrd, 0);
  for (i = 0; i < nelem; ++i) {
    for (j = 0; j < 2; ++j) {
      foccrd_data[i*2 + j] -= self->x.crpix[j];
    }
  }
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
      wcserr_to_python_exc(self->x.err);
      return NULL;
    }
  }
}

/*@null@*/ static PyObject*
PySip_foc2pix(
    PySip* self,
    PyObject* args,
    PyObject* kwds) {

  PyObject*      foccrd_obj = NULL;
  int            origin     = 1;
  PyArrayObject* foccrd     = NULL;
  PyArrayObject* pixcrd     = NULL;
  int            status     = -1;
  double*        foccrd_data = NULL;
  unsigned int   nelem      = 0;
  unsigned int   i, j;
  const char*    keywords[] = {
    "foccrd", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oi:foc2pix", (char **)keywords,
                                   &foccrd_obj, &origin)) {
    return NULL;
  }

  if (self->x.ap == NULL || self->x.bp == NULL) {
    PyErr_SetString(
        PyExc_ValueError,
        "SIP object does not have coefficients for foc2pix transformation (AP and BP)");
    return NULL;
  }

  foccrd = (PyArrayObject*)PyArray_ContiguousFromAny(foccrd_obj, NPY_DOUBLE, 2, 2);
  if (foccrd == NULL) {
    goto exit;
  }

  if (PyArray_DIM(foccrd, 1) != 2) {
    PyErr_SetString(PyExc_ValueError, "Pixel array must be an Nx2 array");
    goto exit;
  }

  pixcrd = (PyArrayObject*)PyArray_SimpleNew(2, PyArray_DIMS(foccrd),
                                             NPY_DOUBLE);
  if (pixcrd == NULL) {
    status = 2;
    goto exit;
  }

  Py_BEGIN_ALLOW_THREADS
  preoffset_array(foccrd, origin);
  /* Adjust for crpix */
  foccrd_data = (double *)PyArray_DATA(foccrd);
  nelem = (unsigned int)PyArray_DIM(foccrd, 0);
  for (i = 0; i < nelem; ++i) {
    for (j = 0; j < 2; ++j) {
      foccrd_data[i*2 + j] += self->x.crpix[j];
    }
  }

  status = sip_foc2pix(&self->x,
                       (unsigned int)PyArray_DIM(pixcrd, 1),
                       (unsigned int)PyArray_DIM(pixcrd, 0),
                       (double*)PyArray_DATA(foccrd),
                       (double*)PyArray_DATA(pixcrd));

  /* Adjust for crpix */
  for (i = 0; i < nelem; ++i) {
    for (j = 0; j < 2; ++j) {
      foccrd_data[i*2 + j] -= self->x.crpix[j];
    }
  }
  unoffset_array(foccrd, origin);
  unoffset_array(pixcrd, origin);
  Py_END_ALLOW_THREADS

 exit:
  Py_XDECREF(foccrd);

  if (status == 0) {
    return (PyObject*)pixcrd;
  } else {
    Py_XDECREF(pixcrd);
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
PySip_get_a(
    PySip* self,
    /*@unused@*/ void* closure) {

  npy_intp dims[2];

  if (self->x.a == NULL) {
    Py_INCREF(Py_None);
    return Py_None;
  }

  dims[0] = (npy_intp)self->x.a_order + 1;
  dims[1] = (npy_intp)self->x.a_order + 1;

  return get_double_array("a", self->x.a, 2, dims, (PyObject*)self);
}

/*@null@*/ static PyObject*
PySip_get_b(
    PySip* self,
    /*@unused@*/ void* closure) {

  npy_intp dims[2];

  if (self->x.b == NULL) {
    Py_INCREF(Py_None);
    return Py_None;
  }

  dims[0] = (npy_intp)self->x.b_order + 1;
  dims[1] = (npy_intp)self->x.b_order + 1;

  return get_double_array("b", self->x.b, 2, dims, (PyObject*)self);
}

/*@null@*/ static PyObject*
PySip_get_ap(
    PySip* self,
    /*@unused@*/ void* closure) {

  npy_intp dims[2];

  if (self->x.ap == NULL) {
    Py_INCREF(Py_None);
    return Py_None;
  }

  dims[0] = (npy_intp)self->x.ap_order + 1;
  dims[1] = (npy_intp)self->x.ap_order + 1;

  return get_double_array("ap", self->x.ap, 2, dims, (PyObject*)self);
}

/*@null@*/ static PyObject*
PySip_get_bp(
    PySip* self,
    /*@unused@*/ void* closure) {

  npy_intp dims[2];

  if (self->x.bp == NULL) {
    Py_INCREF(Py_None);
    return Py_None;
  }

  dims[0] = (npy_intp)self->x.bp_order + 1;
  dims[1] = (npy_intp)self->x.bp_order + 1;

  return get_double_array("bp", self->x.bp, 2, dims, (PyObject*)self);
}

static PyObject*
PySip_get_a_order(
    PySip* self,
    /*@unused@*/ void* closure) {

  return get_int("a_order", (long int)self->x.a_order);
}

static PyObject*
PySip_get_b_order(
    PySip* self,
    /*@unused@*/ void* closure) {

  return get_int("b_order", (long int)self->x.b_order);
}

static PyObject*
PySip_get_ap_order(
    PySip* self,
    /*@unused@*/ void* closure) {

  return get_int("ap_order", (long int)self->x.ap_order);
}

static PyObject*
PySip_get_bp_order(
    PySip* self,
    /*@unused@*/ void* closure) {

  return get_int("bp_order", (long int)self->x.bp_order);
}

static PyObject*
PySip_get_crpix(
    PySip* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return get_double_array("crpix", self->x.crpix, 1, &naxis, (PyObject*)self);
}

static PyObject*
PySip___copy__(
    PySip* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PySip* copy         = NULL;

  copy = (PySip*)PySip_new(&PySipType, NULL, NULL);
  if (copy == NULL) {
    return NULL;
  }

  if (sip_init(&copy->x,
               self->x.a_order, self->x.a,
               self->x.b_order, self->x.b,
               self->x.ap_order, self->x.ap,
               self->x.bp_order, self->x.bp,
               self->x.crpix)) {
    Py_DECREF(copy);
    return NULL;
  }

  return (PyObject*)copy;
}


static PyGetSetDef PySip_getset[] = {
  {"a", (getter)PySip_get_a, NULL, (char *)doc_a},
  {"a_order", (getter)PySip_get_a_order, NULL, (char *)doc_a_order},
  {"b", (getter)PySip_get_b, NULL, (char *)doc_b},
  {"b_order", (getter)PySip_get_b_order, NULL, (char *)doc_b_order},
  {"ap", (getter)PySip_get_ap, NULL, (char *)doc_ap},
  {"ap_order", (getter)PySip_get_ap_order, NULL, (char *)doc_ap_order},
  {"bp", (getter)PySip_get_bp, NULL, (char *)doc_bp},
  {"bp_order", (getter)PySip_get_bp_order, NULL, (char *)doc_bp_order},
  {"crpix", (getter)PySip_get_crpix, NULL, (char *)doc_crpix},
  {NULL}
};

static PyMethodDef PySip_methods[] = {
  {"__copy__", (PyCFunction)PySip___copy__, METH_NOARGS, NULL},
  {"__deepcopy__", (PyCFunction)PySip___copy__, METH_O, NULL},
  {"pix2foc", (PyCFunction)PySip_pix2foc, METH_VARARGS|METH_KEYWORDS, doc_sip_pix2foc},
  {"foc2pix", (PyCFunction)PySip_foc2pix, METH_VARARGS|METH_KEYWORDS, doc_sip_foc2pix},
  {NULL}
};

PyTypeObject PySipType = {
  #if PY3K
  PyVarObject_HEAD_INIT(NULL, 0)
  #else
  PyObject_HEAD_INIT(NULL)
  0,                            /*ob_size*/
  #endif
  "astropy.wcs.Sip",            /*tp_name*/
  sizeof(PySip),                /*tp_basicsize*/
  0,                            /*tp_itemsize*/
  (destructor)PySip_dealloc,    /*tp_dealloc*/
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
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  doc_Sip,                      /* tp_doc */
  0,                            /* tp_traverse */
  0,                            /* tp_clear */
  0,                            /* tp_richcompare */
  0,                            /* tp_weaklistoffset */
  0,                            /* tp_iter */
  0,                            /* tp_iternext */
  PySip_methods,                /* tp_methods */
  0,                            /* tp_members */
  PySip_getset,                 /* tp_getset */
  0,                            /* tp_base */
  0,                            /* tp_dict */
  0,                            /* tp_descr_get */
  0,                            /* tp_descr_set */
  0,                            /* tp_dictoffset */
  (initproc)PySip_init,         /* tp_init */
  0,                            /* tp_alloc */
  PySip_new,                    /* tp_new */
};

int
_setup_sip_type(
    PyObject* m) {

  if (PyType_Ready(&PySipType) < 0)
    return -1;

  Py_INCREF(&PySipType);
  return PyModule_AddObject(m, "Sip", (PyObject *)&PySipType);
}
