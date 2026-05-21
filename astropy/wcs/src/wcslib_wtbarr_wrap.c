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
 * Wtbarr methods                                                        *
 ***************************************************************************/

static PyObject*
Wtbarr_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
  Wtbarr* self;
  allocfunc alloc_func = PyType_GetSlot(type, Py_tp_alloc);
  self = (Wtbarr*)alloc_func(type, 0);
  return (PyObject*)self;
}


static int
Wtbarr_traverse(Wtbarr* self, visitproc visit, void *arg) {
  Py_VISIT(self->owner);
  Py_VISIT(Py_TYPE((PyObject*)self));
  return 0;
}


static int
Wtbarr_clear(Wtbarr* self) {
  Py_CLEAR(self->owner);
  return 0;
}


static void Wtbarr_dealloc(Wtbarr* self) {
  Wtbarr_clear(self);
  PyTypeObject *tp = Py_TYPE((PyObject*)self);
  freefunc free_func = PyType_GetSlot(tp, Py_tp_free);
  free_func((PyObject*)self);
  Py_DECREF(tp);
}


Wtbarr* Wtbarr_cnew(PyObject* wcsprm, struct wtbarr* x) {
  Wtbarr* self;
  PyTypeObject* type = (PyTypeObject*)WtbarrType;
  allocfunc alloc_func = PyType_GetSlot(type, Py_tp_alloc);
  self = (Wtbarr*)alloc_func(type, 0);
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


static PyObject* Wtbarr_print_contents(Wtbarr* self) {
  /* This is not thread-safe, but since we're holding onto the GIL,
     we can assume we won't have thread conflicts */
  wcsprintf_set(NULL);
  wtbarrprt(self->x);
  printf("%s", wcsprintf_buf());
  fflush(stdout);
  Py_RETURN_NONE;
}


static PyObject* Wtbarr___str__(Wtbarr* self) {
  /* This is not thread-safe, but since we're holding onto the GIL,
     we can assume we won't have thread conflicts */
  wcsprintf_set(NULL);
  wtbarrprt(self->x);
  return PyUnicode_FromString(wcsprintf_buf());
}


/***************************************************************************
 * Member getters/setters (properties)
 */


static PyObject* Wtbarr_get_i(Wtbarr* self, void* closure) {
  return get_int("i", self->x->i);
}


static PyObject* Wtbarr_get_m(Wtbarr* self, void* closure) {
  return get_int("m", self->x->m);
}


static PyObject* Wtbarr_get_extver(Wtbarr* self, void* closure) {
  return get_int("extver", self->x->extver);
}


static PyObject* Wtbarr_get_extlev(Wtbarr* self, void* closure) {
  return get_int("extlev", self->x->extlev);
}


static PyObject* Wtbarr_get_ndim(Wtbarr* self, void* closure) {
  return get_int("ndim", self->x->ndim);
}


static PyObject* Wtbarr_get_row(Wtbarr* self, void* closure) {
  return get_int("row", self->x->row);
}


static PyObject* Wtbarr_get_extnam(Wtbarr* self, void* closure) {
  if (is_null(self->x->extnam)) return NULL;
  return get_string("extnam", self->x->extnam);
}


static PyObject* Wtbarr_get_ttype(Wtbarr* self, void* closure) {
  if (is_null(self->x->ttype)) return NULL;
  return get_string("ttype", self->x->ttype);
}


static PyObject* Wtbarr_get_kind(Wtbarr* self, void* closure) {
  return PyUnicode_FromFormat("%c", self->x->kind);
}


/***************************************************************************
 * Wtbarr definition structures
 */

static PyGetSetDef Wtbarr_getset[] = {
  {"i", (getter)Wtbarr_get_i, NULL, (char *) doc_i},
  {"m", (getter)Wtbarr_get_m, NULL, (char *) doc_m},
  {"kind", (getter)Wtbarr_get_kind, NULL, (char *) doc_kind},
  {"extnam", (getter)Wtbarr_get_extnam, NULL, (char *) doc_extnam},
  {"extver", (getter)Wtbarr_get_extver, NULL, (char *) doc_extver},
  {"extlev", (getter)Wtbarr_get_extlev, NULL, (char *) doc_extlev},
  {"ttype", (getter)Wtbarr_get_ttype, NULL, (char *) doc_ttype},
  {"row", (getter)Wtbarr_get_row, NULL, (char *) doc_row},
  {"ndim", (getter)Wtbarr_get_ndim, NULL, (char *) doc_ndim},
/*  {"dimlen", (getter)Wtbarr_get_dimlen, NULL, (char *) NULL}, */
/*  {"arrayp", (getter)Wtbarr_get_arrayp, NULL, (char *) NULL}, */
  {NULL}
};


static PyMethodDef Wtbarr_methods[] = {
  {"print_contents", (PyCFunction)Wtbarr_print_contents, METH_NOARGS, doc_print_contents_wtbarr},
  {NULL}
};

static PyType_Spec WtbarrType_spec = {
  .name = "astropy.wcs.Wtbarr",
  .basicsize = sizeof(Wtbarr),
  .itemsize = 0,
  .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_IMMUTABLETYPE,
  .slots = (PyType_Slot[]){
    {Py_tp_dealloc, (destructor)Wtbarr_dealloc},
    {Py_tp_str, (reprfunc)Wtbarr___str__},
    {Py_tp_doc, doc_Wtbarr},
    {Py_tp_traverse, (traverseproc)Wtbarr_traverse},
    {Py_tp_clear, (inquiry)Wtbarr_clear},
    {Py_tp_getset, Wtbarr_getset},
    {Py_tp_methods, Wtbarr_methods},
    // FIXME: this seems logical but this slot was not previously defined
    // maybe an error from https://github.com/astropy/astropy/pull/9641 ?
    // {Py_tp_new, (newfunc)Wtbarr_new},
    {0, NULL},
  },
};

PyObject* WtbarrType = NULL;

int
_setup_wtbarr_type(PyObject* m) {
  WtbarrType = PyType_FromSpec(&WtbarrType_spec);
  if (WtbarrType == NULL) {
    return -1;
  }

  PyModule_AddObject(m, "Wtbarr", WtbarrType);

  return 0;
}
