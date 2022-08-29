#define NO_IMPORT_ARRAY

#include "astropy_wcs/wcslib_wtbarr_wrap.h"

#include <wcs.h>
#include <wcsprintf.h>
#include <tab.h>
#include <wtbarr.h>

/*
 It gets to be really tedious to type long docstrings in ANSI C syntax
 (since multi-line strings literals are not valid).  Therefore, the
 docstrings are written in doc/docstrings.py, which are then converted
 by setup.py into docstrings.h, which we include here.
*/
#include "astropy_wcs/docstrings.h"


/***************************************************************************
 * PyWtbarr methods                                                        *
 ***************************************************************************/

static PyObject*
PyWtbarr_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
  PyWtbarr* self;
  self = (PyWtbarr*)type->tp_alloc(type, 0);
  return (PyObject*)self;
}


static int
PyWtbarr_traverse(PyWtbarr* self, visitproc visit, void *arg) {
  Py_VISIT(self->owner);
  return 0;
}


static int
PyWtbarr_clear(PyWtbarr* self) {
  Py_CLEAR(self->owner);
  return 0;
}


static void PyWtbarr_dealloc(PyWtbarr* self) {
  PyWtbarr_clear(self);
  Py_TYPE(self)->tp_free((PyObject*)self);
}


PyWtbarr* PyWtbarr_cnew(PyObject* wcsprm, struct wtbarr* x) {
  PyWtbarr* self;
  self = (PyWtbarr*)(&PyWtbarrType)->tp_alloc(&PyWtbarrType, 0);
  if (self == NULL) return NULL;
  self->x = x;
  Py_INCREF(wcsprm);
  self->owner = wcsprm;
  return self;
}


static void wtbarrprt(const struct wtbarr *wtb) {
  int i, nd, ndim;

  if (wtb == 0x0) return;

  wcsprintf("     i: %d\n", wtb->i);
  wcsprintf("     m: %d\n", wtb->m);
  wcsprintf("  kind: %c\n", wtb->kind);
  wcsprintf("extnam: %s\n", wtb->extnam);
  wcsprintf("extver: %d\n", wtb->extver);
  wcsprintf("extlev: %d\n", wtb->extlev);
  wcsprintf(" ttype: %s\n", wtb->ttype);
  wcsprintf("   row: %ld\n", wtb->row);
  wcsprintf("  ndim: %d\n", wtb->ndim);
  wcsprintf("dimlen: %p\n", (void *)wtb->dimlen);

  ndim = wtb->ndim - (int)(wtb->kind == 'c');
  nd = 1 + (int) log10(ndim ? ndim : 1);
  for (i = 0; i < ndim; i++) {
    wcsprintf("        %*d:   %d\n", nd, i, wtb->dimlen[i]);
  }
  wcsprintf("arrayp: %p\n", (void *)wtb->arrayp);

  return;
}


static PyObject* PyWtbarr_print_contents(PyWtbarr* self) {
  /* This is not thread-safe, but since we're holding onto the GIL,
     we can assume we won't have thread conflicts */
  wcsprintf_set(NULL);
  wtbarrprt(self->x);
  printf("%s", wcsprintf_buf());
  fflush(stdout);
  Py_RETURN_NONE;
}


static PyObject* PyWtbarr___str__(PyWtbarr* self) {
  /* This is not thread-safe, but since we're holding onto the GIL,
     we can assume we won't have thread conflicts */
  wcsprintf_set(NULL);
  wtbarrprt(self->x);
  return PyUnicode_FromString(wcsprintf_buf());
}


/***************************************************************************
 * Member getters/setters (properties)
 */


static PyObject* PyWtbarr_get_i(PyWtbarr* self, void* closure) {
  return get_int("i", self->x->i);
}


static PyObject* PyWtbarr_get_m(PyWtbarr* self, void* closure) {
  return get_int("m", self->x->m);
}


static PyObject* PyWtbarr_get_extver(PyWtbarr* self, void* closure) {
  return get_int("extver", self->x->extver);
}


static PyObject* PyWtbarr_get_extlev(PyWtbarr* self, void* closure) {
  return get_int("extlev", self->x->extlev);
}


static PyObject* PyWtbarr_get_ndim(PyWtbarr* self, void* closure) {
  return get_int("ndim", self->x->ndim);
}


static PyObject* PyWtbarr_get_row(PyWtbarr* self, void* closure) {
  return get_int("row", self->x->row);
}


static PyObject* PyWtbarr_get_extnam(PyWtbarr* self, void* closure) {
  if (is_null(self->x->extnam)) return NULL;
  return get_string("extnam", self->x->extnam);
}


static PyObject* PyWtbarr_get_ttype(PyWtbarr* self, void* closure) {
  if (is_null(self->x->ttype)) return NULL;
  return get_string("ttype", self->x->ttype);
}


static PyObject* PyWtbarr_get_kind(PyWtbarr* self, void* closure) {
  return PyUnicode_FromFormat("%c", self->x->kind);
}


/***************************************************************************
 * PyWtbarr definition structures
 */

static PyGetSetDef PyWtbarr_getset[] = {
  {"i", (getter)PyWtbarr_get_i, NULL, (char *) doc_i},
  {"m", (getter)PyWtbarr_get_m, NULL, (char *) doc_m},
  {"kind", (getter)PyWtbarr_get_kind, NULL, (char *) doc_kind},
  {"extnam", (getter)PyWtbarr_get_extnam, NULL, (char *) doc_extnam},
  {"extver", (getter)PyWtbarr_get_extver, NULL, (char *) doc_extver},
  {"extlev", (getter)PyWtbarr_get_extlev, NULL, (char *) doc_extlev},
  {"ttype", (getter)PyWtbarr_get_ttype, NULL, (char *) doc_ttype},
  {"row", (getter)PyWtbarr_get_row, NULL, (char *) doc_row},
  {"ndim", (getter)PyWtbarr_get_ndim, NULL, (char *) doc_ndim},
/*  {"dimlen", (getter)PyWtbarr_get_dimlen, NULL, (char *) NULL}, */
/*  {"arrayp", (getter)PyWtbarr_get_arrayp, NULL, (char *) NULL}, */
  {NULL}
};


static PyMethodDef PyWtbarr_methods[] = {
  {"print_contents", (PyCFunction)PyWtbarr_print_contents, METH_NOARGS, doc_print_contents_wtbarr},
  {NULL}
};

PyTypeObject PyWtbarrType = {
  PyVarObject_HEAD_INIT(NULL, 0)
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
  (reprfunc)PyWtbarr___str__,   /*tp_str*/
  0,                            /*tp_getattro*/
  0,                            /*tp_setattro*/
  0,                            /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  doc_Wtbarr,                   /* tp_doc */
  (traverseproc)PyWtbarr_traverse, /* tp_traverse */
  (inquiry)PyWtbarr_clear,      /* tp_clear */
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
_setup_wtbarr_type(PyObject* m) {
  if (PyType_Ready(&PyWtbarrType) < 0) {
    return -1;
  }

  Py_INCREF(&PyWtbarrType);

  PyModule_AddObject(m, "Wtbarr", (PyObject *)&PyWtbarrType);

  return 0;
}
