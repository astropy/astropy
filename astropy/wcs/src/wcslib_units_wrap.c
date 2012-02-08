/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

#include "wcslib_units_wrap.h"

#include <wcsunits.h>

/*
 It gets to be really tedious to type long docstrings in ANSI C syntax
 (since multi-line strings literals are not valid).  Therefore, the
 docstrings are written in doc/docstrings.py, which are then converted
 by setup.py into docstrings.h, which we include here.
*/
#include "docstrings.h"

/***************************************************************************
 * PyTabprm methods
 */

static int
PyUnits_traverse(
    PyUnits* self, visitproc visit, void *arg) {

  return 0;
}

static int
PyUnits_clear(
    PyUnits* self) {

  return 0;
}

static void
PyUnits_dealloc(
    PyUnits* self) {

  PyUnits_clear(self);
  Py_TYPE(self)->tp_free((PyObject*)self);
}

PyUnits*
PyUnits_cnew(
    const char* const have,
    const char* const want,
    const double scale,
    const double offset,
    const double power) {

  PyUnits* self;
  self = (PyUnits*)(&PyUnitsType)->tp_alloc(&PyUnitsType, 0);
  if (have == NULL) {
    self->have[0] = 0;
  } else {
    strncpy(self->have, have, 80);
  }
  if (want == NULL) {
    self->want[0] = 0;
  } else {
    strncpy(self->want, want, 80);
  }
  self->scale = scale;
  self->offset = offset;
  self->power = power;
  return self;
}

static PyObject *
PyUnits_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyUnits* self;
  self = (PyUnits*)type->tp_alloc(type, 0);
  self->have[0] = 0;
  self->want[0] = 0;
  self->scale = 1.0;
  self->offset = 0.0;
  self->power = 1.0;
  return (PyObject*)self;
}

static int
PyUnits_init(
    PyUnits* self,
    PyObject* args,
    PyObject* kwds) {

  int            status     = -1;
  const char*    have;
  const char*    want;
  const char*    ctrl_str   = NULL;
  int            ctrl       = 0;
  const char*    keywords[] = {"have", "want", "translate_units", NULL};
  struct wcserr* err        = NULL;

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ss|s:UnitConverter.__init__",
                                   (char **)keywords, &have, &want,
                                   &ctrl_str)) {
    goto exit;
  }

  if (ctrl_str != NULL) {
    if (parse_unsafe_unit_conversion_spec(ctrl_str, &ctrl)) {
      goto exit;
    }
  }

  /* Copy the strings since we can't have wcslib monkeying with the
     data in a Python string */
  strncpy(self->have, have, 80);
  strncpy(self->want, want, 80);

  status = wcsutrne(ctrl, self->have, &err);
  if (status != -1 && status != 0) {
    goto exit;
  }

  status = wcsutrne(ctrl, self->want, &err);
  if (status != -1 && status != 0) {
    goto exit;
  }

  status = wcsunitse(self->have, self->want,
                     &self->scale, &self->offset, &self->power, &err);

 exit:

  if (PyErr_Occurred()) {
    return -1;
  } else if (status == 0) {
    return 0;
  } else {
    wcserr_units_to_python_exc(err);
    free(err);
    return -1;
  }
}

/*@null@*/ static PyObject*
PyUnits___str__(
    PyUnits* self) {

  #define BUF_SIZE (1 << 8)
  char buffer[BUF_SIZE];
  char scale[BUF_SIZE];
  char offset[BUF_SIZE];
  char power[BUF_SIZE];

  if (self->scale != 1.0) {
    snprintf(scale, BUF_SIZE, "*%.12g", self->scale);
  } else {
    scale[0] = 0;
  }

  if (self->offset != 0.0) {
    snprintf(offset, BUF_SIZE, " + %.12g", self->offset);
  } else {
    offset[0] = 0;
  }

  if (self->power != 1.0) {
    snprintf(power, BUF_SIZE, " ** %.12g", self->power);
  } else {
    power[0] = 0;
  }

  snprintf(buffer, 1 << 8, "<astropy.wcs.UnitConverter '%s' to '%s' (x%s%s)%s>",
           self->have, self->want, scale, offset, power);

  #if PY3K
  return PyUnicode_FromString(buffer);
  #else
  return PyString_FromString(buffer);
  #endif
}

static PyObject*
PyUnits_convert(
    PyUnits* self,
    PyObject* args,
    PyObject* kwds) {

  int            status       = 1;
  PyObject*      input        = NULL;
  PyArrayObject* input_arr    = NULL;
  PyArrayObject* output_arr   = NULL;
  PyObject*      input_iter   = NULL;
  PyObject*      output_iter  = NULL;
  double         input_val;
  double         output_val;

  if (!PyArg_ParseTuple(args, "O:UnitConverter.convert", &input)) {
    goto exit;
  }

  input_arr = (PyArrayObject*)PyArray_FromObject(
      input, NPY_DOUBLE, 0, NPY_MAXDIMS);
  if (input_arr == NULL) {
    goto exit;
  }

  output_arr = (PyArrayObject*)PyArray_SimpleNew(
      PyArray_NDIM(input_arr), PyArray_DIMS(input_arr), PyArray_DOUBLE);
  if (output_arr == NULL) {
    goto exit;
  }

  input_iter = PyArray_IterNew((PyObject*)input_arr);
  if (input_iter == NULL) {
    goto exit;
  }

  output_iter = PyArray_IterNew((PyObject*)output_arr);
  if (output_iter == NULL) {
    goto exit;
  }

  if (self->power != 1.0) {
    while (PyArray_ITER_NOTDONE(input_iter)) {
      input_val = *(double *)PyArray_ITER_DATA(input_iter);
      output_val = pow(self->scale*input_val + self->offset, self->power);
      if (errno) {
        PyErr_SetFromErrno(PyExc_ValueError);
        goto exit;
      }
      *(double *)PyArray_ITER_DATA(output_iter) = output_val;
      PyArray_ITER_NEXT(input_iter);
      PyArray_ITER_NEXT(output_iter);
    }
  } else {
    while (PyArray_ITER_NOTDONE(input_iter)) {
      input_val = *(double *)PyArray_ITER_DATA(input_iter);
      output_val = self->scale*input_val + self->offset;
      *(double *)PyArray_ITER_DATA(output_iter) = output_val;
      PyArray_ITER_NEXT(input_iter);
      PyArray_ITER_NEXT(output_iter);
    }
  }

  status = 0;

 exit:

  Py_XDECREF((PyObject*)input_arr);
  Py_XDECREF(input_iter);
  Py_XDECREF(output_iter);
  if (status) {
    Py_XDECREF((PyObject*)output_arr);
    return NULL;
  }
  return (PyObject*)output_arr;
}

/***************************************************************************
 * Member getters/setters (properties)
 */

/*@null@*/ static PyObject*
PyUnits_get_have(
    PyUnits* self,
    /*@unused@*/ void* closure) {

  #if PY3K
  return PyUnicode_FromString(self->have);
  #else
  return PyString_FromString(self->have);
  #endif
}

/*@null@*/ static PyObject*
PyUnits_get_want(
    PyUnits* self,
    /*@unused@*/ void* closure) {

  #if PY3K
  return PyUnicode_FromString(self->want);
  #else
  return PyString_FromString(self->want);
  #endif
}

/*@null@*/ static PyObject*
PyUnits_get_scale(
    PyUnits* self,
    /*@unused@*/ void* closure) {

  return get_double("scale", self->scale);
}

/*@null@*/ static PyObject*
PyUnits_get_offset(
    PyUnits* self,
    /*@unused@*/ void* closure) {

  return get_double("offset", self->offset);
}

/*@null@*/ static PyObject*
PyUnits_get_power(
    PyUnits* self,
    /*@unused@*/ void* closure) {

  return get_double("power", self->power);
}

/***************************************************************************
 * PyUnits definition structures
 */

static PyGetSetDef PyUnits_getset[] = {
  {"have", (getter)PyUnits_get_have, NULL, (char *)doc_have},
  {"want", (getter)PyUnits_get_want, NULL, (char *)doc_want},
  {"scale", (getter)PyUnits_get_scale, NULL, (char *)doc_scale},
  {"offset", (getter)PyUnits_get_offset, NULL, (char *)doc_offset},
  {"power", (getter)PyUnits_get_power, NULL, (char *)doc_power},
  {NULL}
};

static PyMethodDef PyUnits_methods[] = {
  {"convert", (PyCFunction)PyUnits_convert, METH_VARARGS, doc_convert},
  {NULL}
};

PyTypeObject PyUnitsType = {
  #if PY3K
  PyVarObject_HEAD_INIT(NULL, 0)
  #else
  PyObject_HEAD_INIT(NULL)
  0,                            /*ob_size*/
  #endif
  "astropy.wcs.UnitConverter",        /*tp_name*/
  sizeof(PyUnits),              /*tp_basicsize*/
  0,                            /*tp_itemsize*/
  (destructor)PyUnits_dealloc,  /*tp_dealloc*/
  0,                            /*tp_print*/
  0,                            /*tp_getattr*/
  0,                            /*tp_setattr*/
  0,                            /*tp_compare*/
  (reprfunc)PyUnits___str__,    /*tp_repr*/
  0,                            /*tp_as_number*/
  0,                            /*tp_as_sequence*/
  0,                            /*tp_as_mapping*/
  0,                            /*tp_hash */
  0,                            /*tp_call*/
  (reprfunc)PyUnits___str__,    /*tp_str*/
  0,                            /*tp_getattro*/
  0,                            /*tp_setattro*/
  0,                            /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  doc_UnitConverter,            /* tp_doc */
  (traverseproc)PyUnits_traverse, /* tp_traverse */
  (inquiry)PyUnits_clear,       /* tp_clear */
  0,                            /* tp_richcompare */
  0,                            /* tp_weaklistoffset */
  0,                            /* tp_iter */
  0,                            /* tp_iternext */
  PyUnits_methods,              /* tp_methods */
  0,                            /* tp_members */
  PyUnits_getset,               /* tp_getset */
  0,                            /* tp_base */
  0,                            /* tp_dict */
  0,                            /* tp_descr_get */
  0,                            /* tp_descr_set */
  0,                            /* tp_dictoffset */
  (initproc)PyUnits_init,       /* tp_init */
  0,                            /* tp_alloc */
  PyUnits_new,                  /* tp_new */
};

int
_setup_units_type(
    PyObject* m) {

  if (PyType_Ready(&PyUnitsType) < 0) {
    return -1;
  }

  Py_INCREF(&PyUnitsType);

  PyModule_AddObject(m, "UnitConverter", (PyObject *)&PyUnitsType);

  return 0;
}
