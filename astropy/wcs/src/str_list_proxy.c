/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

#include <stdlib.h> // malloc, free
#include "astropy_wcs/pyutil.h"

/***************************************************************************
 * List-of-strings proxy object
 ***************************************************************************/

static PyObject* PyStrListProxyType;

typedef struct {
  PyObject_HEAD
  /*@null@*/ /*@shared@*/ PyObject* pyobject;
  Py_ssize_t size;
  Py_ssize_t maxsize;
  char (*array)[72];
} PyStrListProxy;

static void
PyStrListProxy_dealloc(
    PyStrListProxy* self) {

  PyObject_GC_UnTrack(self);
  Py_XDECREF(self->pyobject);
  PyTypeObject *tp = Py_TYPE((PyObject*)self);
  freefunc free_func = PyType_GetSlot(tp, Py_tp_free);
  free_func((PyObject*)self);
  Py_DECREF(tp);
}

/*@null@*/ static PyObject *
PyStrListProxy_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyStrListProxy* self = NULL;

  allocfunc alloc_func = PyType_GetSlot(type, Py_tp_alloc);
  self = (PyStrListProxy*)alloc_func(type, 0);
  if (self != NULL) {
    self->pyobject = NULL;
  }
  return (PyObject*)self;
}

static int
PyStrListProxy_traverse(
    PyStrListProxy* self,
    visitproc visit,
    void *arg) {

  Py_VISIT(self->pyobject);
  Py_VISIT((PyObject*)Py_TYPE((PyObject*)self));
  return 0;
}

static int
PyStrListProxy_clear(
    PyStrListProxy *self) {

  Py_CLEAR(self->pyobject);

  return 0;
}

/*@null@*/ PyObject *
PyStrListProxy_New(
    /*@shared@*/ PyObject* owner,
    Py_ssize_t size,
    Py_ssize_t maxsize,
    char (*array)[72]) {

  PyStrListProxy* self = NULL;

  if (maxsize == 0) {
    maxsize = 68;
  }

  PyTypeObject* tp = (PyTypeObject*)PyStrListProxyType;
  allocfunc alloc_func = PyType_GetSlot(tp, Py_tp_alloc);
  self = (PyStrListProxy*)alloc_func(tp, 0);
  if (self == NULL) {
    return NULL;
  }

  Py_XINCREF(owner);
  self->pyobject = owner;
  self->size = size;
  self->maxsize = maxsize;
  self->array = array;
  return (PyObject*)self;
}

static Py_ssize_t
PyStrListProxy_len(
    PyStrListProxy* self) {

  return self->size;
}

/*@null@*/ static PyObject*
PyStrListProxy_getitem(
    PyStrListProxy* self,
    Py_ssize_t index) {

  if (index >= self->size || index < 0) {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    return NULL;
  }

  return get_string("string", self->array[index]);
}

static int
PyStrListProxy_setitem(
    PyStrListProxy* self,
    Py_ssize_t index,
    PyObject* arg) {

  if (index >= self->size || index < 0) {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    return -1;
  }

  return set_string("string", arg, self->array[index], self->maxsize);
}

/*@null@*/ PyObject*
str_list_proxy_repr(
    char (*array)[72],
    Py_ssize_t size,
    Py_ssize_t maxsize) {

  char*       buffer  = NULL;
  char*       wp      = NULL;
  char*       rp      = NULL;
  Py_ssize_t  i       = 0;
  Py_ssize_t  j       = 0;
  PyObject*   result  = NULL;
  /* These are in descending order, so we can exit the loop quickly.  They
     are in pairs: (char_to_escape, char_escaped) */
  const char* escapes   = "\\\\''\rr\ff\vv\nn\tt\bb\aa";
  const char* e         = NULL;
  char        next_char = '\0';

  /* Overallocating to allow for escaped characters */
  buffer = malloc((size_t)size*maxsize*2 + 2);
  if (buffer == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Could not allocate memory.");
    return NULL;
  }

  wp = buffer;
  *wp++ = '[';

  for (i = 0; i < size; ++i) {
    *wp++ = '\'';
    rp = array[i];
    for (j = 0; j < maxsize && *rp != '\0'; ++j) {
      /* Check if this character should be escaped */
      e = escapes;
      next_char = *rp++;
      do {
        if (next_char > *e) {
          break;
        } else if (next_char == *e) {
          *wp++ = '\\';
          next_char = *(++e);
          break;
        } else {
          e += 2;
        }
      } while (*e != '\0');

      *wp++ = next_char;
    }
    *wp++ = '\'';

    /* Add a comma for all but the last one */
    if (i != size - 1) {
      *wp++ = ',';
      *wp++ = ' ';
    }
  }

  *wp++ = ']';
  *wp++ = '\0';

  result = PyUnicode_FromString(buffer);
  free(buffer);
  return result;
}

/*@null@*/ static PyObject*
PyStrListProxy_repr(
    PyStrListProxy* self) {

  return str_list_proxy_repr(self->array, self->size, self->maxsize);
}

static PyType_Spec PyStrListProxyType_spec = {
  .name = "astropy.wcs.StrListProxy",
  .basicsize = sizeof(PyStrListProxy),
  .itemsize = 0,
  .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
  .slots = (PyType_Slot[]){
    {Py_tp_dealloc, (destructor)PyStrListProxy_dealloc},
    {Py_tp_repr, (reprfunc)PyStrListProxy_repr},
    {Py_sq_length, (lenfunc)PyStrListProxy_len},
    {Py_sq_item, (ssizeargfunc)PyStrListProxy_getitem},
    {Py_sq_ass_item, (ssizeobjargproc)PyStrListProxy_setitem},
    {Py_tp_str, (reprfunc)PyStrListProxy_repr},
    {Py_tp_traverse, (traverseproc)PyStrListProxy_traverse},
    {Py_tp_clear, (inquiry)PyStrListProxy_clear},
    {Py_tp_new, (newfunc)PyStrListProxy_new},
    {0, NULL},
  },
};

static PyObject* PyStrListProxyType = NULL;

int
_setup_str_list_proxy_type(
    /*@unused@*/ PyObject* m) {

  PyStrListProxyType = PyType_FromSpec(&PyStrListProxyType_spec);
  if (PyStrListProxyType == NULL) {
    return 1;
  }

  return 0;
}
