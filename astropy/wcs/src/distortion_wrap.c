/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

#include "astropy_wcs/distortion_wrap.h"
#include "astropy_wcs/docstrings.h"

#include <structmember.h> /* From Python */

static int
PyDistLookup_traverse(
    PyDistLookup* self,
    visitproc visit,
    void* arg) {

  Py_VISIT(self->py_data);
  Py_VISIT((PyObject*)Py_TYPE((PyObject*)self));

  return 0;
}

static int
PyDistLookup_clear(
    PyDistLookup* self) {

  Py_CLEAR(self->py_data);

  return 0;
}

static void
PyDistLookup_dealloc(
    PyDistLookup* self) {

  PyObject_GC_UnTrack(self);
  distortion_lookup_t_free(&self->x);
  Py_XDECREF((PyObject*)self->py_data);
  PyTypeObject *tp = Py_TYPE((PyObject*)self);
  freefunc free_func = PyType_GetSlot(tp, Py_tp_free);
  free_func((PyObject*)self);
  Py_DECREF(tp);
}

/*@null@*/ static PyObject *
PyDistLookup_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyDistLookup* self;

  allocfunc alloc_func = PyType_GetSlot(type, Py_tp_alloc);
  self = (PyDistLookup*)alloc_func(type, 0);
  if (self != NULL) {
    if (distortion_lookup_t_init(&self->x)) {
      return NULL;
    }
    self->py_data = NULL;
  }
  return (PyObject*)self;
}

static int
PyDistLookup_init(
    PyDistLookup* self,
    PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyObject* py_array_obj = NULL;
  PyArrayObject* array_obj = NULL;

  if (!PyArg_ParseTuple(args, "O(dd)(dd)(dd):DistortionLookupTable.__init__",
                        &py_array_obj,
                        &(self->x.crpix[0]), &(self->x.crpix[1]),
                        &(self->x.crval[0]), &(self->x.crval[1]),
                        &(self->x.cdelt[0]), &(self->x.cdelt[1]))) {
    return -1;
  }

  array_obj = (PyArrayObject*)PyArray_ContiguousFromAny(py_array_obj, NPY_FLOAT32, 2, 2);
  if (array_obj == NULL) {
    return -1;
  }

  self->py_data = array_obj;
  self->x.naxis[0] = (unsigned int)PyArray_DIM(array_obj, 1);
  self->x.naxis[1] = (unsigned int)PyArray_DIM(array_obj, 0);
  self->x.data = (float *)PyArray_DATA(array_obj);

  return 0;
}

static PyObject*
PyDistLookup_get_cdelt(
    PyDistLookup* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return get_double_array("cdelt", self->x.cdelt, 1, &naxis, (PyObject*)self);
}

static int
PyDistLookup_set_cdelt(
    PyDistLookup* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 2;

  return set_double_array("cdelt", value, 1, &naxis, self->x.cdelt);
}

static PyObject*
PyDistLookup_get_crpix(
    PyDistLookup* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return get_double_array("crpix", self->x.crpix, 1, &naxis, (PyObject*)self);
}

static int
PyDistLookup_set_crpix(
    PyDistLookup* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 2;

  return set_double_array("crpix", value, 1, &naxis, self->x.crpix);
}

static PyObject*
PyDistLookup_get_crval(
    PyDistLookup* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return get_double_array("crval", self->x.crval, 1, &naxis, (PyObject*)self);
}

static int
PyDistLookup_set_crval(
    PyDistLookup* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 2;

  return set_double_array("crval", value, 1, &naxis, self->x.crval);
}

/*@shared@*/ static PyObject*
PyDistLookup_get_data(
    PyDistLookup* self,
    /*@unused@*/ void* closure) {

  if (self->py_data == NULL) {
    Py_INCREF(Py_None);
    return Py_None;
  } else {
    Py_INCREF((PyObject*)self->py_data);
    return (PyObject*)self->py_data;
  }
}

static int
PyDistLookup_set_data(
    PyDistLookup* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  PyArrayObject* value_array = NULL;

  if (value == NULL) {
    Py_CLEAR(self->py_data);
    self->x.data = NULL;
    return 0;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value, NPY_FLOAT32, 2, 2);

  if (value_array == NULL) {
    return -1;
  }

  Py_XDECREF((PyObject*)self->py_data);

  self->py_data = value_array;
  self->x.naxis[0] = (unsigned int)PyArray_DIM(value_array, 1);
  self->x.naxis[1] = (unsigned int)PyArray_DIM(value_array, 0);
  self->x.data = (float *)PyArray_DATA(value_array);

  return 0;
}

/*@null@*/ static PyObject*
PyDistLookup_get_offset(
    PyDistLookup* self,
    PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  double coord[NAXES];
  double result;

  if (self->x.data == NULL) {
    PyErr_SetString(PyExc_RuntimeError,
                    "No data has been set for the lookup table");
    return NULL;
  }

  if (!PyArg_ParseTuple(args, "dd:get_offset", &coord[0], &coord[1])) {
    return NULL;
  }

  result = get_distortion_offset(&self->x, coord);
  return PyFloat_FromDouble(result);
}

static PyObject*
PyDistLookup___copy__(
    PyDistLookup* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyDistLookup* copy = NULL;
  int           i    = 0;

  copy = (PyDistLookup*)PyDistLookup_new((PyTypeObject*)PyDistLookupType, NULL, NULL);
  if (copy == NULL) {
    return NULL;
  }

  for (i = 0; i < 2; ++i) {
    copy->x.naxis[i] = self->x.naxis[i];
    copy->x.crpix[i] = self->x.crpix[i];
    copy->x.crval[i] = self->x.crval[i];
    copy->x.cdelt[i] = self->x.cdelt[i];
  }

  if (self->py_data) {
    PyDistLookup_set_data(copy, (PyObject*)self->py_data, NULL);
  }

  return (PyObject*)copy;
}

static PyObject*
PyDistLookup___deepcopy__(
    PyDistLookup* self,
    PyObject* memo,
    /*@unused@*/ PyObject* kwds) {

  PyDistLookup* copy;
  PyObject*     obj_copy;
  int           i = 0;

  copy = (PyDistLookup*)PyDistLookup_new((PyTypeObject*)PyDistLookupType, NULL, NULL);
  if (copy == NULL) {
    return NULL;
  }

  for (i = 0; i < 2; ++i) {
    copy->x.naxis[i] = self->x.naxis[i];
    copy->x.crpix[i] = self->x.crpix[i];
    copy->x.crval[i] = self->x.crval[i];
    copy->x.cdelt[i] = self->x.cdelt[i];
  }

  if (self->py_data) {
    obj_copy = get_deepcopy((PyObject*)self->py_data, memo);
    if (obj_copy == NULL) {
      Py_DECREF(copy);
      return NULL;
    }
    PyDistLookup_set_data(copy, (PyObject*)obj_copy, NULL);
    Py_DECREF(obj_copy);
  }

  return (PyObject*)copy;
}


static PyGetSetDef PyDistLookup_getset[] = {
  {"cdelt", (getter)PyDistLookup_get_cdelt, (setter)PyDistLookup_set_cdelt, (char *)doc_cdelt},
  {"crpix", (getter)PyDistLookup_get_crpix, (setter)PyDistLookup_set_crpix, (char *)doc_crpix},
  {"crval", (getter)PyDistLookup_get_crval, (setter)PyDistLookup_set_crval, (char *)doc_crval},
  {"data",  (getter)PyDistLookup_get_data,  (setter)PyDistLookup_set_data,  (char *)doc_data},
  {NULL}
};

static PyMethodDef PyDistLookup_methods[] = {
  {"__copy__", (PyCFunction)PyDistLookup___copy__, METH_NOARGS, NULL},
  {"__deepcopy__", (PyCFunction)PyDistLookup___deepcopy__, METH_O, NULL},
  {"get_offset", (PyCFunction)PyDistLookup_get_offset, METH_VARARGS, doc_get_offset},
  {NULL}
};

static PyType_Spec PyDistLookupType_spec = {
  .name = "astropy.wcs.DistortionLookupTable",
  .basicsize = sizeof(PyDistLookup),
  .itemsize = 0,
  .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_IMMUTABLETYPE,
  .slots = (PyType_Slot[]){
    {Py_tp_dealloc, (destructor)PyDistLookup_dealloc},
    {Py_tp_doc, doc_DistortionLookupTable},
    {Py_tp_traverse, (traverseproc)PyDistLookup_traverse},
    {Py_tp_clear, (inquiry)PyDistLookup_clear},
    {Py_tp_methods, PyDistLookup_methods},
    {Py_tp_getset, PyDistLookup_getset},
    {Py_tp_init, (initproc)PyDistLookup_init},
    {Py_tp_new, PyDistLookup_new},
    {0, NULL},
  },
};

PyObject* PyDistLookupType = NULL;

int _setup_distortion_type(
    PyObject* m) {

  PyDistLookupType = PyType_FromSpec(&PyDistLookupType_spec);
  if (PyDistLookupType == NULL) {
    return -1;
  }

  return PyModule_AddObject(m, "DistortionLookupTable", PyDistLookupType);
}
