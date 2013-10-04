/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

#include "astropy_wcs/wcslib_wtbarr_wrap.h"

#include <wcs.h>

/*
 It gets to be really tedious to type long docstrings in ANSI C syntax
 (since multi-line strings literals are not valid).  Therefore, the
 docstrings are written in doc/docstrings.py, which are then converted
 by setup.py into docstrings.h, which we include here.
*/
#include "astropy_wcs/docstrings.h"

/***************************************************************************
 * Helper functions                                                        *
 ***************************************************************************/

/***************************************************************************
 * PyWtbarr methods
 */

static int
PyWtbarr_traverse(
    PyWtbarr* self, visitproc visit, void *arg) {
  int vret;

  vret = visit(self->owner, arg);
  if (vret != 0) {
    return vret;
  }

  return 0;
}

static int
PyWtbarr_clear(
    PyWtbarr* self) {
  PyObject* tmp;

  tmp = self->owner;
  self->owner = NULL;
  Py_XDECREF(tmp);

  return 0;
}

static void
PyWtbarr_dealloc(
    PyWtbarr* self) {

  PyWtbarr_clear(self);
  Py_TYPE(self)->tp_free((PyObject*)self);
}

PyWtbarr*
PyWtbarr_cnew(PyObject* wcsprm, struct wtbarr* x) {
  PyWtbarr* self;
  self = (PyWtbarr*)(&PyWtbarrType)->tp_alloc(&PyWtbarrType, 0);
  self->x = x;
  Py_INCREF(wcsprm);
  self->owner = wcsprm;
  return self;
}

/***************************************************************************
 * Member getters/setters (properties)
 */

/*@null@*/ static PyObject*
PyWtbarr_get_data(
    PyWtbarr* self,
    /*@unused@*/ void* closure) {

  int ndims;
  npy_intp dims[NPY_MAXDIMS];
  int i;

  if (is_null(self->x->arrayp)) {
    return NULL;
  }

  ndims = self->x->ndim;
  for (i = 0; i < ndims; ++i) {
    dims[i] = self->x->dimlen[i];
  }

  return get_double_array("data", *self->x->arrayp, ndims, dims, (PyObject*)self);
}

/*@null@*/ static int
PyWtbarr_set_data(
    PyWtbarr* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  int ndims;
  npy_intp dims[NPY_MAXDIMS];
  int i;

  if (is_null(self->x->arrayp)) {
    return -1;
  }

  ndims = self->x->ndim;
  for (i = 0; i < ndims; ++i) {
    dims[i] = self->x->dimlen[i];
  }

  return set_double_array("data", value, ndims, dims, *self->x->arrayp);
}

/*@null@*/ static PyObject*
PyWtbarr_get_dims(
    PyWtbarr* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t ndims = 0;

  if (is_null(self->x->dimlen)) {
    return NULL;
  }

  ndims = (Py_ssize_t)self->x->ndim;

  return get_int_array("dims", self->x->dimlen, 1, &ndims, (PyObject*)self);
}

/*@null@*/ static PyObject*
PyWtbarr_get_extlev(
    PyWtbarr* self,
    /*@unused@*/ void* closure) {

  return get_int("extlev", self->x->extlev);
}

/*@null@*/ static PyObject*
PyWtbarr_get_extnam(
    PyWtbarr* self,
    /*@unused@*/ void* closure) {

  return get_string("extnam", self->x->extnam);
}

/*@null@*/ static PyObject*
PyWtbarr_get_extver(
    PyWtbarr* self,
    /*@unused@*/ void* closure) {

  return get_int("extver", self->x->extver);
}

/*@null@*/ static PyObject*
PyWtbarr_get_i(
    PyWtbarr* self,
    /*@unused@*/ void* closure) {

  return get_int("i", self->x->i);
}

/*@null@*/ static PyObject*
PyWtbarr_get_kind(
    PyWtbarr* self,
    /*@unused@*/ void* closure) {

  char kind = (char)self->x->kind;

  #if PY3K
  return PyUnicode_FromStringAndSize(&kind, 1);
  #else
  return PyString_FromStringAndSize(&kind, 1);
  #endif
}

/*@null@*/ static PyObject*
PyWtbarr_get_m(
    PyWtbarr* self,
    /*@unused@*/ void* closure) {

  return get_int("m", self->x->m);
}

/*@null@*/ static PyObject*
PyWtbarr_get_ndim(
    PyWtbarr* self,
    /*@unused@*/ void* closure) {

  return get_int("ndim", self->x->ndim);
}

/*@null@*/ static PyObject*
PyWtbarr_get_row(
    PyWtbarr* self,
    /*@unused@*/ void* closure) {

  return get_int("row", self->x->row);
}

/*@null@*/ static PyObject*
PyWtbarr_get_ttype(
    PyWtbarr* self,
    /*@unused@*/ void* closure) {

  return get_string("ttype", self->x->ttype);
}

/***************************************************************************
 * PyWtbarr definition structures
 */

static PyGetSetDef PyWtbarr_getset[] = {
  {"data", (getter)PyWtbarr_get_data, (setter)PyWtbarr_set_data, (char *)doc_data},
  {"dims", (getter)PyWtbarr_get_dims, NULL, (char *)doc_dims},
  {"extlev", (getter)PyWtbarr_get_extlev, NULL, (char *)doc_extlev},
  {"extnam", (getter)PyWtbarr_get_extnam, NULL, (char *)doc_extnam},
  {"extver", (getter)PyWtbarr_get_extver, NULL, (char *)doc_extver},
  {"i", (getter)PyWtbarr_get_i, NULL, (char *)doc_i},
  {"kind", (getter)PyWtbarr_get_kind, NULL, (char *)doc_kind},
  {"m", (getter)PyWtbarr_get_m, NULL, (char *)doc_m},
  {"ndim", (getter)PyWtbarr_get_ndim, NULL, (char *)doc_ndim},
  {"row", (getter)PyWtbarr_get_row, NULL, (char *)doc_row},
  {"ttype", (getter)PyWtbarr_get_ttype, NULL, (char *)doc_ttype},
  {NULL}
};

static PyMethodDef PyWtbarr_methods[] = {
  {NULL}
};

PyTypeObject PyWtbarrType = {
  #if PY3K
  PyVarObject_HEAD_INIT(NULL, 0)
  #else
  PyObject_HEAD_INIT(NULL)
  0,                            /*ob_size*/
  #endif
  "astropy.wcs.Wtbarr",         /*tp_name*/
  sizeof(PyWtbarr),             /*tp_basicsize*/
  0,                            /*tp_itemsize*/
  (destructor)PyWtbarr_dealloc, /*tp_dealloc*/
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
  doc_Wtbarr,                   /* tp_doc */
  (traverseproc)PyWtbarr_traverse, /* tp_traverse */
  (inquiry)PyWtbarr_clear,         /* tp_clear */
  0,                            /* tp_richcompare */
  0,                            /* tp_weaklistoffset */
  0,                            /* tp_iter */
  0,                            /* tp_iternext */
  PyWtbarr_methods,             /* tp_methods */
  0,                            /* tp_members */
  PyWtbarr_getset,              /* tp_getset */
  0,                            /* tp_base */
  0,                            /* tp_dict */
  0,                            /* tp_descr_get */
  0,                            /* tp_descr_set */
  0,                            /* tp_dictoffset */
  0,                            /* tp_init */
  0,                            /* tp_alloc */
  0,                            /* tp_new */
};

int
_setup_wtbarr_type(
    PyObject* m) {

  if (PyType_Ready(&PyWtbarrType) < 0) {
    return -1;
  }

  Py_INCREF(&PyWtbarrType);

  PyModule_AddObject(m, "Wtbarr", (PyObject *)&PyWtbarrType);

  return 0;
}
