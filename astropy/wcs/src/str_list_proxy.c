/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

#include "astropy_wcs/pyutil.h"
#include <stdlib.h> // malloc, free

/***************************************************************************
 * List-of-strings proxy object
 ***************************************************************************/

static PyObject* StrListProxyType;

typedef struct {
  PyObject_HEAD
  /*@null@*/ /*@shared@*/ PyObject* pyobject;
  Py_ssize_t size;
  Py_ssize_t maxsize;
  char (*array)[72];
} StrListProxy;

static void
StrListProxy_dealloc(
    StrListProxy* self) {

  PyObject_GC_UnTrack(self);
  Py_XDECREF(self->pyobject);
  PyTypeObject *tp = Py_TYPE((PyObject*)self);
  freefunc free_func = PyType_GetSlot(tp, Py_tp_free);
  free_func((PyObject*)self);
  Py_DECREF(tp);
}

/*@null@*/ static PyObject *
StrListProxy_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  StrListProxy* self = NULL;

  allocfunc alloc_func = PyType_GetSlot(type, Py_tp_alloc);
  self = (StrListProxy*)alloc_func(type, 0);
  if (self != NULL) {
    self->pyobject = NULL;
  }
  return (PyObject*)self;
}

static int
StrListProxy_traverse(
    StrListProxy* self,
    visitproc visit,
    void *arg) {

  Py_VISIT(self->pyobject);
  Py_VISIT((PyObject*)Py_TYPE((PyObject*)self));
  return 0;
}

static int
StrListProxy_clear(
    StrListProxy *self) {

  Py_CLEAR(self->pyobject);

  return 0;
}

/*@null@*/ PyObject *
StrListProxy_New(
    /*@shared@*/ PyObject* owner,
    Py_ssize_t size,
    Py_ssize_t maxsize,
    char (*array)[72]) {

  StrListProxy* self = NULL;

  if (maxsize == 0) {
    maxsize = 68;
  }

  PyTypeObject* tp = (PyTypeObject*)StrListProxyType;
  allocfunc alloc_func = PyType_GetSlot(tp, Py_tp_alloc);
  self = (StrListProxy*)alloc_func(tp, 0);
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
StrListProxy_len(
    StrListProxy* self) {

  return self->size;
}

/*@null@*/ static PyObject*
StrListProxy_getitem(
    StrListProxy* self,
    Py_ssize_t index) {

  if (index >= self->size || index < 0) {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    return NULL;
  }

  return get_string("string", self->array[index]);
}

static int
StrListProxy_setitem(
    StrListProxy* self,
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
StrListProxy_repr(
    StrListProxy* self) {

  return str_list_proxy_repr(self->array, self->size, self->maxsize);
}

static PyType_Spec StrListProxyType_spec = {
  .name = "astropy.wcs.StrListProxy",
  .basicsize = sizeof(StrListProxy),
  .itemsize = 0,
  .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
  .slots = (PyType_Slot[]){
    {Py_tp_dealloc, (destructor)StrListProxy_dealloc},
    {Py_tp_repr, (reprfunc)StrListProxy_repr},
    {Py_sq_length, (lenfunc)StrListProxy_len},
    {Py_sq_item, (ssizeargfunc)StrListProxy_getitem},
    {Py_sq_ass_item, (ssizeobjargproc)StrListProxy_setitem},
    {Py_tp_str, (reprfunc)StrListProxy_repr},
    {Py_tp_traverse, (traverseproc)StrListProxy_traverse},
    {Py_tp_clear, (inquiry)StrListProxy_clear},
    {Py_tp_new, (newfunc)StrListProxy_new},
    {0, NULL},
  },
};

static PyObject* StrListProxyType = NULL;

int
_setup_str_list_proxy_type(
    /*@unused@*/ PyObject* m) {

  StrListProxyType = PyType_FromSpec(&StrListProxyType_spec);
  if (StrListProxyType == NULL) {
    return 1;
  }

  return 0;
}
