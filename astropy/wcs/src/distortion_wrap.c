/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

#include "astropy_wcs/distortion_wrap.h"
#include "astropy_wcs/docstrings.h"

#include <structmember.h> /* From Python */

static int
DistLookup_traverse(
    DistLookup* self,
    visitproc visit,
    void* arg) {

  Py_VISIT(self->py_data);
  Py_VISIT((PyObject*)Py_TYPE((PyObject*)self));

  return 0;
}

static int
DistLookup_clear(
    DistLookup* self) {

  Py_CLEAR(self->py_data);

  return 0;
}

static void
DistLookup_dealloc(
    DistLookup* self) {

  PyObject_GC_UnTrack(self);
  distortion_lookup_t_free(&self->x);
  Py_XDECREF((PyObject*)self->py_data);
  PyTypeObject *tp = Py_TYPE((PyObject*)self);
  freefunc free_func = PyType_GetSlot(tp, Py_tp_free);
  free_func((PyObject*)self);
  Py_DECREF(tp);
}

/*@null@*/ static PyObject *
DistLookup_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  DistLookup* self;

  allocfunc alloc_func = PyType_GetSlot(type, Py_tp_alloc);
  self = (DistLookup*)alloc_func(type, 0);
  if (self != NULL) {
    if (distortion_lookup_t_init(&self->x)) {
      return NULL;
    }
    self->py_data = NULL;
  }
  return (PyObject*)self;
}

static int
DistLookup_init(
    DistLookup* self,
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
DistLookup_get_cdelt(
    DistLookup* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return get_double_array("cdelt", self->x.cdelt, 1, &naxis, (PyObject*)self);
}

static int
DistLookup_set_cdelt(
    DistLookup* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 2;

  return set_double_array("cdelt", value, 1, &naxis, self->x.cdelt);
}

static PyObject*
DistLookup_get_crpix(
    DistLookup* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return get_double_array("crpix", self->x.crpix, 1, &naxis, (PyObject*)self);
}

static int
DistLookup_set_crpix(
    DistLookup* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 2;

  return set_double_array("crpix", value, 1, &naxis, self->x.crpix);
}

static PyObject*
DistLookup_get_crval(
    DistLookup* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 2;

  return get_double_array("crval", self->x.crval, 1, &naxis, (PyObject*)self);
}

static int
DistLookup_set_crval(
    DistLookup* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 2;

  return set_double_array("crval", value, 1, &naxis, self->x.crval);
}

/*@shared@*/ static PyObject*
DistLookup_get_data(
    DistLookup* self,
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
DistLookup_set_data(
    DistLookup* self,
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
DistLookup_get_offset(
    DistLookup* self,
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
DistLookup___copy__(
    DistLookup* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  DistLookup* copy = NULL;
  int           i    = 0;

  copy = (DistLookup*)DistLookup_new((PyTypeObject*)DistLookupType, NULL, NULL);
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
    DistLookup_set_data(copy, (PyObject*)self->py_data, NULL);
  }

  return (PyObject*)copy;
}

static PyObject*
DistLookup___deepcopy__(
    DistLookup* self,
    PyObject* memo,
    /*@unused@*/ PyObject* kwds) {

  DistLookup* copy;
  PyObject*     obj_copy;
  int           i = 0;

  copy = (DistLookup*)DistLookup_new((PyTypeObject*)DistLookupType, NULL, NULL);
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
    DistLookup_set_data(copy, (PyObject*)obj_copy, NULL);
    Py_DECREF(obj_copy);
  }

  return (PyObject*)copy;
}


static PyGetSetDef DistLookup_getset[] = {
  {"cdelt", (getter)DistLookup_get_cdelt, (setter)DistLookup_set_cdelt, (char *)doc_cdelt},
  {"crpix", (getter)DistLookup_get_crpix, (setter)DistLookup_set_crpix, (char *)doc_crpix},
  {"crval", (getter)DistLookup_get_crval, (setter)DistLookup_set_crval, (char *)doc_crval},
  {"data",  (getter)DistLookup_get_data,  (setter)DistLookup_set_data,  (char *)doc_data},
  {NULL}
};

static PyMethodDef DistLookup_methods[] = {
  {"__copy__", (PyCFunction)DistLookup___copy__, METH_NOARGS, NULL},
  {"__deepcopy__", (PyCFunction)DistLookup___deepcopy__, METH_O, NULL},
  {"get_offset", (PyCFunction)DistLookup_get_offset, METH_VARARGS, doc_get_offset},
  {NULL}
};

static PyType_Spec DistLookupType_spec = {
  .name = "astropy.wcs.DistortionLookupTable",
  .basicsize = sizeof(DistLookup),
  .itemsize = 0,
  .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_IMMUTABLETYPE,
  .slots = (PyType_Slot[]){
    {Py_tp_dealloc, (destructor)DistLookup_dealloc},
    {Py_tp_doc, doc_DistortionLookupTable},
    {Py_tp_traverse, (traverseproc)DistLookup_traverse},
    {Py_tp_clear, (inquiry)DistLookup_clear},
    {Py_tp_methods, DistLookup_methods},
    {Py_tp_getset, DistLookup_getset},
    {Py_tp_init, (initproc)DistLookup_init},
    {Py_tp_new, DistLookup_new},
    {0, NULL},
  },
};

PyObject* DistLookupType = NULL;

int _setup_distortion_type(
    PyObject* m) {

  DistLookupType = PyType_FromSpec(&DistLookupType_spec);
  if (DistLookupType == NULL) {
    return -1;
  }

  return PyModule_AddObject(m, "DistortionLookupTable", DistLookupType);
}
