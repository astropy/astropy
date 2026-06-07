/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

#include "astropy_wcs/wcslib_tabprm_wrap.h"

#include <wcs.h>
#include <wcsprintf.h>
#include <tab.h>

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

static INLINE void
note_change(Tabprm* self) {
  self->x->flag = 0;
}

static int
make_fancy_dims(Tabprm* self, int* ndims, npy_intp* dims) {
  int i, M;

  M = self->x->M;
  if (M + 1 > NPY_MAXDIMS) {
    PyErr_SetString(PyExc_ValueError, "Too many dimensions");
    return -1;
  }

  *ndims = M + 1;

  for (i = 0; i < M; ++i) {
    dims[i] = self->x->K[M-1-i];
  }

  dims[M] = M;

  return 0;
}

PyObject** tab_errexc[6];

static void
wcslib_tab_to_python_exc(int status) {
  if (status > 0 && status < 6) {
    PyErr_SetString(*tab_errexc[status], tab_errmsg[status]);
  } else {
    PyErr_SetString(
        PyExc_RuntimeError,
        "Unknown error occurred.  Something is seriously wrong.");
  }
}

/***************************************************************************
 * Tabprm methods
 */

static int
Tabprm_traverse(
    Tabprm* self, visitproc visit, void *arg) {
  Py_VISIT(self->owner);
  Py_VISIT(Py_TYPE((PyObject*)self));
  return 0;
}

static int
Tabprm_clear(
    Tabprm* self) {

  Py_CLEAR(self->owner);

  return 0;
}

static void
Tabprm_dealloc(
    Tabprm* self) {

  Tabprm_clear(self);
  PyTypeObject *tp = Py_TYPE((PyObject*)self);
  freefunc free_func = PyType_GetSlot(tp, Py_tp_free);
  free_func((PyObject*)self);
  Py_DECREF(tp);
}

Tabprm*
Tabprm_cnew(PyObject* wcsprm, struct tabprm* x) {
  Tabprm* self;
  PyTypeObject* type = (PyTypeObject*)TabprmType;
  allocfunc alloc_func = PyType_GetSlot(type, Py_tp_alloc);
  self = (Tabprm*)alloc_func(type, 0);
  if (self == NULL) return NULL;
  self->x = x;
  Py_INCREF(wcsprm);
  self->owner = wcsprm;
  return self;
}

static int
Tabprm_cset(
    Tabprm* self) {

  int status = 0;

  status = tabset(self->x);

  if (status == 0) {
    return 0;
  } else {
    wcslib_tab_to_python_exc(status);
    return -1;
  }
}

/*@null@*/ static PyObject*
Tabprm_set(
    Tabprm* self) {

  if (Tabprm_cset(self)) {
    return NULL;
  }

  Py_RETURN_NONE;
}

/*@null@*/ static PyObject*
Tabprm_print_contents(
    Tabprm* self) {

  if (Tabprm_cset(self)) {
    return NULL;
  }

  /* This is not thread-safe, but since we're holding onto the GIL,
     we can assume we won't have thread conflicts */
  wcsprintf_set(NULL);
  tabprt(self->x);
  printf("%s", wcsprintf_buf());
  fflush(stdout);
  Py_RETURN_NONE;
}

/*@null@*/ static PyObject*
Tabprm___str__(
    Tabprm* self) {

  if (Tabprm_cset(self)) {
    return NULL;
  }

  /* This is not thread-safe, but since we're holding onto the GIL,
     we can assume we won't have thread conflicts */
  wcsprintf_set(NULL);

  tabprt(self->x);

  return PyUnicode_FromString(wcsprintf_buf());
}

/***************************************************************************
 * Member getters/setters (properties)
 */

/*@null@*/ static PyObject*
Tabprm_get_coord(
    Tabprm* self,
    /*@unused@*/ void* closure) {

  int ndims;
  npy_intp dims[NPY_MAXDIMS];

  if (is_null(self->x->coord)) {
    return NULL;
  }

  if (make_fancy_dims(self, &ndims, dims)) {
    return NULL;
  }

  return get_double_array("coord", self->x->coord, ndims, dims, (PyObject*)self);
}

/*@null@*/ static int
Tabprm_set_coord(
    Tabprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  int ndims;
  npy_intp dims[NPY_MAXDIMS];

  if (is_null(self->x->coord)) {
    return -1;
  }

  if (make_fancy_dims(self, &ndims, dims)) {
    return -1;
  }

  return set_double_array("coord", value, ndims, dims, self->x->coord);
}

/*@null@*/ static PyObject*
Tabprm_get_crval(
    Tabprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t M = 0;

  if (is_null(self->x->crval)) {
    return NULL;
  }

  M = (Py_ssize_t)self->x->M;

  return get_double_array("crval", self->x->crval, 1, &M, (PyObject*)self);
}

static int
Tabprm_set_crval(
    Tabprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp M = 0;

  if (is_null(self->x->crval)) {
    return -1;
  }

  M = (Py_ssize_t)self->x->M;

  note_change(self);

  return set_double_array("crval", value, 1, &M, self->x->crval);
}

/*@null@*/ static PyObject*
Tabprm_get_delta(
    Tabprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t M = 0;

  if (is_null(self->x->delta)) {
    return NULL;
  }

  M = (Py_ssize_t)self->x->M;

  return get_double_array("delta", self->x->delta, 1, &M, (PyObject*)self);
}

/*@null@*/ static PyObject*
Tabprm_get_extrema(
    Tabprm* self,
    /*@unused@*/ void* closure) {

  int ndims;
  npy_intp dims[NPY_MAXDIMS];

  if (is_null(self->x->coord)) {
    return NULL;
  }

  if (make_fancy_dims(self, &ndims, dims)) {
    return NULL;
  }

  dims[ndims-2] = 2;

  return get_double_array("extrema", self->x->extrema, ndims, dims, (PyObject*)self);
}

/*@null@*/ static PyObject*
Tabprm_get_K(
    Tabprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t M = 0;

  if (is_null(self->x->K)) {
    return NULL;
  }

  M = (Py_ssize_t)self->x->M;

  return get_int_array("K", self->x->K, 1, &M, (PyObject*)self);
}

/*@null@*/ static PyObject*
Tabprm_get_M(
    Tabprm* self,
    /*@unused@*/ void* closure) {

  return get_int("M", self->x->M);
}

/*@null@*/ static PyObject*
Tabprm_get_map(
    Tabprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t M = 0;

  if (is_null(self->x->map)) {
    return NULL;
  }

  M = (Py_ssize_t)self->x->M;

  return get_int_array("map", self->x->map, 1, &M, (PyObject*)self);
}

static int
Tabprm_set_map(
    Tabprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp M = 0;

  if (is_null(self->x->map)) {
    return -1;
  }

  M = (Py_ssize_t)self->x->M;

  note_change(self);

  return set_int_array("map", value, 1, &M, self->x->map);
}

/*@null@*/ static PyObject*
Tabprm_get_nc(
    Tabprm* self,
    /*@unused@*/ void* closure) {

  return get_int("nc", self->x->nc);
}

/*@null@*/ static PyObject*
Tabprm_get_p0(
    Tabprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t M = 0;

  if (is_null(self->x->p0)) {
    return NULL;
  }

  M = (Py_ssize_t)self->x->M;

  return get_int_array("p0", self->x->p0, 1, &M, (PyObject*)self);
}

/*@null@*/ static PyObject*
Tabprm_get_sense(
    Tabprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t M = 0;

  if (is_null(self->x->sense)) {
    return NULL;
  }

  M = (Py_ssize_t)self->x->M;

  return get_int_array("sense", self->x->sense, 1, &M, (PyObject*)self);
}

/***************************************************************************
 * Tabprm definition structures
 */

static PyGetSetDef Tabprm_getset[] = {
  {"coord", (getter)Tabprm_get_coord, (setter)Tabprm_set_coord, (char *)doc_coord},
  {"crval", (getter)Tabprm_get_crval, (setter)Tabprm_set_crval, (char *)doc_crval_tabprm},
  {"delta", (getter)Tabprm_get_delta, NULL, (char *)doc_delta},
  {"extrema", (getter)Tabprm_get_extrema, NULL, (char *)doc_extrema},
  {"K", (getter)Tabprm_get_K, NULL, (char *)doc_K},
  {"M", (getter)Tabprm_get_M, NULL, (char *)doc_M},
  {"map", (getter)Tabprm_get_map, (setter)Tabprm_set_map, (char *)doc_map},
  {"nc", (getter)Tabprm_get_nc, NULL, (char *)doc_nc},
  {"p0", (getter)Tabprm_get_p0, NULL, (char *)doc_p0},
  {"sense", (getter)Tabprm_get_sense, NULL, (char *)doc_sense},
  {NULL}
};

static PyMethodDef Tabprm_methods[] = {
  {"print_contents", (PyCFunction)Tabprm_print_contents, METH_NOARGS, doc_print_contents_tabprm},
  {"set", (PyCFunction)Tabprm_set, METH_NOARGS, doc_set_tabprm},
  {NULL}
};

static PyType_Spec TabprmType_spec = {
  .name = "astropy.wcs.Tabprm",
  .basicsize = sizeof(Tabprm),
  .itemsize = 0,
  .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_IMMUTABLETYPE,
  .slots = (PyType_Slot[]) {
    {Py_tp_dealloc, (destructor)Tabprm_dealloc},
    {Py_tp_str, (reprfunc)Tabprm___str__},
    {Py_tp_doc, doc_Tabprm},
    {Py_tp_traverse, (traverseproc)Tabprm_traverse},
    {Py_tp_clear, (inquiry)Tabprm_clear},
    {Py_tp_getset, Tabprm_getset},
    {Py_tp_methods, Tabprm_methods},
    {0, NULL},
  },
};

PyObject* TabprmType = NULL;

int
_setup_tabprm_type(
    PyObject* m) {

  TabprmType = PyType_FromSpec(&TabprmType_spec);
  if (TabprmType == NULL) {
    return -1;
  }

  PyModule_AddObject(m, "Tabprm", TabprmType);

  tab_errexc[0] = NULL;                         /* Success */
  tab_errexc[1] = &PyExc_MemoryError;           /* Null wcsprm pointer passed */
  tab_errexc[2] = &PyExc_MemoryError;           /* Memory allocation failed */
  tab_errexc[3] = &WcsExc_InvalidTabularParameters;  /* Invalid tabular parameters */
  tab_errexc[4] = &WcsExc_InvalidCoordinate; /* One or more of the x coordinates were invalid */
  tab_errexc[5] = &WcsExc_InvalidCoordinate; /* One or more of the world coordinates were invalid */

  return 0;
}
