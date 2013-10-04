/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

#include "astropy_wcs/pyutil.h"

/***************************************************************************
 * List-of-strings proxy object
 ***************************************************************************/

static PyTypeObject PyStrListProxyType;

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

  Py_XDECREF(self->pyobject);
  Py_TYPE(self)->tp_free((PyObject*)self);
}

/*@null@*/ static PyObject *
PyStrListProxy_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyStrListProxy* self = NULL;

  self = (PyStrListProxy*)type->tp_alloc(type, 0);
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

  int vret;

  if (self->pyobject) {
    vret = visit(self->pyobject, arg);
    if (vret != 0) {
      return vret;
    }
  }

  return 0;
}

static int
PyStrListProxy_clear(
    PyStrListProxy *self) {

  PyObject *tmp;

  tmp = self->pyobject;
  self->pyobject = NULL;
  Py_XDECREF(tmp);

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

  self = (PyStrListProxy*)PyStrListProxyType.tp_alloc(&PyStrListProxyType, 0);
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

  if (index >= self->size) {
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

  if (index > self->size) {
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

  #if PY3K
  result = PyUnicode_FromString(buffer);
  #else
  result = PyString_FromString(buffer);
  #endif
  free(buffer);
  return result;
}

/*@null@*/ static PyObject*
PyStrListProxy_repr(
    PyStrListProxy* self) {

  return str_list_proxy_repr(self->array, self->size, self->maxsize);
}

static PySequenceMethods PyStrListProxy_sequence_methods = {
  (lenfunc)PyStrListProxy_len,
  NULL,
  NULL,
  (ssizeargfunc)PyStrListProxy_getitem,
  NULL,
  (ssizeobjargproc)PyStrListProxy_setitem,
  NULL,
  NULL,
  NULL,
  NULL
};

static PyTypeObject PyStrListProxyType = {
  #if PY3K
  PyVarObject_HEAD_INIT(NULL, 0)
  #else
  PyObject_HEAD_INIT(NULL)
  0,                          /*ob_size*/
  #endif
  "astropy.wcs.StrListProxy", /*tp_name*/
  sizeof(PyStrListProxy),  /*tp_basicsize*/
  0,                          /*tp_itemsize*/
  (destructor)PyStrListProxy_dealloc, /*tp_dealloc*/
  0,                          /*tp_print*/
  0,                          /*tp_getattr*/
  0,                          /*tp_setattr*/
  0,                          /*tp_compare*/
  (reprfunc)PyStrListProxy_repr, /*tp_repr*/
  0,                          /*tp_as_number*/
  &PyStrListProxy_sequence_methods, /*tp_as_sequence*/
  0,                          /*tp_as_mapping*/
  0,                          /*tp_hash */
  0,                          /*tp_call*/
  (reprfunc)PyStrListProxy_repr, /*tp_str*/
  0,                          /*tp_getattro*/
  0,                          /*tp_setattro*/
  0,                          /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
  0,                          /* tp_doc */
  (traverseproc)PyStrListProxy_traverse, /* tp_traverse */
  (inquiry)PyStrListProxy_clear, /* tp_clear */
  0,                          /* tp_richcompare */
  0,                          /* tp_weaklistoffset */
  0,                          /* tp_iter */
  0,                          /* tp_iternext */
  0,                          /* tp_methods */
  0,                          /* tp_members */
  0,                          /* tp_getset */
  0,                          /* tp_base */
  0,                          /* tp_dict */
  0,                          /* tp_descr_get */
  0,                          /* tp_descr_set */
  0,                          /* tp_dictoffset */
  0,                          /* tp_init */
  0,                          /* tp_alloc */
  PyStrListProxy_new,      /* tp_new */
};

int
_setup_str_list_proxy_type(
    /*@unused@*/ PyObject* m) {

  if (PyType_Ready(&PyStrListProxyType) < 0) {
    return 1;
  }

  return 0;
}
