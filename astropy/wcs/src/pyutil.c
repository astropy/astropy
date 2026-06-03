/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

/* util.h must be imported first */
#include "astropy_wcs/pyutil.h"
#include "astropy_wcs/docstrings.h"
#include "astropy_wcs/wcslib_wrap.h"  /* for the Wcsprm struct (write-back) */

#include <stdlib.h> // malloc, free
#include <string.h> // memcpy

#include "wcsfix.h"
#include "wcshdr.h"
#include "wcsprintf.h"
#include "wcsunits.h"

/*@null@*/ static INLINE PyObject*
_ArrayProxy_New(
    /*@shared@*/ PyObject* self,
    int nd,
    const npy_intp* dims,
    int typenum,
    const void* data,
    const int flags) {

  PyArray_Descr* type_descr = NULL;
  PyObject*      result     = NULL;

  type_descr = (PyArray_Descr*)PyArray_DescrFromType(typenum);
  if (type_descr == NULL) {
    return NULL;
  }

  result = (PyObject*)PyArray_NewFromDescr(
      &PyArray_Type,
      type_descr,
      nd, (npy_intp*)dims,
      NULL,
      (void*)data,
      NPY_ARRAY_C_CONTIGUOUS | flags,
      NULL);

  if (result == NULL) {
    return NULL;
  }
  Py_INCREF(self);
  PyArray_SetBaseObject((PyArrayObject *)result, self);
  return result;
}

/*@null@*/ PyObject*
ArrayProxy_New(
    /*@shared@*/ PyObject* self,
    int nd,
    const npy_intp* dims,
    int typenum,
    const void* data) {

  return _ArrayProxy_New(self, nd, dims, typenum, data, NPY_ARRAY_WRITEABLE);
}

/*@null@*/ PyObject*
ArrayReadOnlyProxy_New(
    /*@shared@*/ PyObject* self,
    int nd,
    const npy_intp* dims,
    int typenum,
    const void* data) {

  return _ArrayProxy_New(self, nd, dims, typenum, data, 0);
}

/****************************************************************************
 * WCSParameterArray
 *
 * The wcsprm struct is stored canonically in WCSLIB's native UNDEFINED form
 * (GH-16409).  An array-valued parameter (crpix, cdelt, cd, pc, ...) is
 * exposed to Python as a WCSParameterArray: an ndarray subclass holding a
 * *copy* of the field with UNDEFINED translated to NaN.  Item assignment is
 * written straight back into the owning wcsprm field with NaN translated to
 * UNDEFINED, so `w.wcs.cdelt[1] = np.nan` stores UNDEFINED.
 *
 * Parent/destination are carried on a tiny write-back token set as the
 * array's base object, so the subclass itself needs no extra inline fields
 * (its tp_basicsize equals ndarray's) and NumPy manages the token's lifetime.
 ****************************************************************************/

typedef struct {
  PyObject_HEAD
  PyObject*   parent;    /* owning Wcsprm, holds a reference */
  double*     dest;      /* &parent->x.<field> */
  npy_intp    nelem;
  const char* propname;  /* static field name, for the deprecation message */
} WCSParamWriteback;

static void
WCSParamWriteback_dealloc(WCSParamWriteback* self) {
  Py_XDECREF(self->parent);
  Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyTypeObject WCSParamWriteback_Type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "astropy.wcs._WCSParamWriteback",
  .tp_basicsize = sizeof(WCSParamWriteback),
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_dealloc = (destructor)WCSParamWriteback_dealloc,
};

static PyMappingMethods WCSParameterArray_as_mapping;

/* Instances carry a __dict__ (NumPy/astropy.units expect any ndarray subclass
 * to have one), stored in an extra pointer past NumPy's ndarray fields.  The
 * offset and basicsize are set in setup since ndarray's size is only known at
 * runtime.  The dict must be cleared on dealloc and traversed for GC. */
static void
WCSParameterArray_dealloc(PyObject* self) {
  PyObject** dictptr = _PyObject_GetDictPtr(self);
  if (dictptr != NULL) {
    Py_CLEAR(*dictptr);
  }
  PyArray_Type.tp_dealloc(self);
}

static int
WCSParameterArray_traverse(PyObject* self, visitproc visit, void* arg) {
  PyObject** dictptr = _PyObject_GetDictPtr(self);
  if (dictptr != NULL) {
    Py_VISIT(*dictptr);
  }
  if (PyArray_Type.tp_traverse != NULL) {
    return PyArray_Type.tp_traverse(self, visit, arg);
  }
  return 0;
}

static int
WCSParameterArray_clear(PyObject* self) {
  PyObject** dictptr = _PyObject_GetDictPtr(self);
  if (dictptr != NULL) {
    Py_CLEAR(*dictptr);
  }
  if (PyArray_Type.tp_clear != NULL) {
    return PyArray_Type.tp_clear(self);
  }
  return 0;
}

/* Expose obj.__dict__ (PyType_Ready does not add the descriptor automatically
 * for a static type that merely sets tp_dictoffset). */
static PyGetSetDef WCSParameterArray_getset[] = {
  {"__dict__", PyObject_GenericGetDict, PyObject_GenericSetDict, NULL, NULL},
  {NULL}
};

PyTypeObject WCSParameterArray_Type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "astropy.wcs.WCSParameterArray",
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
  .tp_dealloc = (destructor)WCSParameterArray_dealloc,
  .tp_traverse = WCSParameterArray_traverse,
  .tp_clear = WCSParameterArray_clear,
  .tp_getset = WCSParameterArray_getset,
  /* tp_base, tp_basicsize, tp_dictoffset and tp_as_mapping are filled in at
   * setup time because NumPy's true ndarray size is only known at runtime. */
};

static int
WCSParameterArray_ass_subscript(PyObject* self, PyObject* key, PyObject* value) {
  PyObject* base = PyArray_BASE((PyArrayObject*)self);

  /* In-place mutation of an auxiliary WCS parameter array is deprecated.
   * Only the interceptable paths (a[i]=, a[:]=) reach here; the eventual fix
   * is to return a read-only array so that every write path (including
   * np.fill_diagonal, a.flat[i]=, a+=) raises loudly instead.
   * TODO: switch these getters to read-only arrays and delete this subclass. */
  if (base != NULL && Py_TYPE(base) == &WCSParamWriteback_Type) {
    WCSParamWriteback* wb = (WCSParamWriteback*)base;
    if (PyErr_WarnFormat(
            PyExc_DeprecationWarning, 1,
            "In-place modification of wcs.wcs.%s is deprecated and will stop "
            "working in a future version; assign the whole attribute instead, "
            "e.g. wcs.wcs.%s = new_values.",
            wb->propname, wb->propname) < 0) {
      return -1;  /* warning escalated to an error */
    }
  }

  /* Perform the normal ndarray assignment first (NaN is allowed). */
  if (PyArray_Type.tp_as_mapping->mp_ass_subscript(self, key, value) != 0) {
    return -1;
  }

  /* Then sync the whole (small) buffer back into the wcsprm field,
   * translating NaN -> UNDEFINED. */
  if (base != NULL && Py_TYPE(base) == &WCSParamWriteback_Type) {
    WCSParamWriteback* wb = (WCSParamWriteback*)base;
    const double* src = (const double*)PyArray_DATA((PyArrayObject*)self);
    for (npy_intp i = 0; i < wb->nelem; ++i) {
      wb->dest[i] = npy_isnan(src[i]) ? UNDEFINED : src[i];
    }
    /* Equivalent to note_change(): force a wcsset on next use. */
    ((Wcsprm*)wb->parent)->x.flag = 0;
  }
  return 0;
}

int
_setup_wcsparameter_array_type(PyObject* m) {
  if (PyType_Ready(&WCSParamWriteback_Type) < 0) {
    return -1;
  }

  /* Inherit ndarray's mapping methods, override only assignment. */
  WCSParameterArray_as_mapping = *PyArray_Type.tp_as_mapping;
  WCSParameterArray_as_mapping.mp_ass_subscript = WCSParameterArray_ass_subscript;

  WCSParameterArray_Type.tp_base = &PyArray_Type;
  /* Reserve one extra pointer past ndarray's fields for the instance __dict__. */
  WCSParameterArray_Type.tp_basicsize =
      PyArray_Type.tp_basicsize + sizeof(PyObject*);
  WCSParameterArray_Type.tp_dictoffset = PyArray_Type.tp_basicsize;
  WCSParameterArray_Type.tp_as_mapping = &WCSParameterArray_as_mapping;
  if (PyType_Ready(&WCSParameterArray_Type) < 0) {
    return -1;
  }
  Py_INCREF(&WCSParameterArray_Type);
  if (PyModule_AddObject(m, "WCSParameterArray",
                         (PyObject*)&WCSParameterArray_Type) < 0) {
    Py_DECREF(&WCSParameterArray_Type);
    return -1;
  }
  return 0;
}

static PyObject*
new_double_array(PyTypeObject* subtype, int ndims, const npy_intp* dims,
                 double* value, npy_intp* nelem_out) {
  npy_intp nelem = 1;
  for (int i = 0; i < ndims; ++i) {
    nelem *= dims[i];
  }
  /* data == NULL here, so the flags argument is interpreted by NumPy as a
   * boolean Fortran-order flag -- pass 0 for C order (NOT
   * NPY_ARRAY_C_CONTIGUOUS, whose nonzero value would request Fortran order
   * and transpose 2-D fields like pc/cd). */
  PyObject* arr = PyArray_NewFromDescr(
      subtype, PyArray_DescrFromType(NPY_DOUBLE), ndims, (npy_intp*)dims,
      NULL, NULL, 0, NULL);
  if (arr == NULL) {
    return NULL;
  }
  double* d = (double*)PyArray_DATA((PyArrayObject*)arr);
  for (npy_intp i = 0; i < nelem; ++i) {
    d[i] = undefined(value[i]) ? (double)NPY_NAN : value[i];
  }
  if (nelem_out) {
    *nelem_out = nelem;
  }
  return arr;
}

/*@null@*/ PyObject*
WCSParameterArray_New(PyObject* owner, const char* propname, int ndims,
                      const npy_intp* dims, double* value) {
  npy_intp nelem = 0;
  PyObject* arr = new_double_array(&WCSParameterArray_Type, ndims, dims, value,
                                   &nelem);
  if (arr == NULL) {
    return NULL;
  }
  WCSParamWriteback* wb = PyObject_New(WCSParamWriteback, &WCSParamWriteback_Type);
  if (wb == NULL) {
    Py_DECREF(arr);
    return NULL;
  }
  Py_INCREF(owner);
  wb->parent = owner;
  wb->dest = value;
  wb->nelem = nelem;
  wb->propname = propname;
  if (PyArray_SetBaseObject((PyArrayObject*)arr, (PyObject*)wb) < 0) {
    Py_DECREF(wb);
    Py_DECREF(arr);
    return NULL;
  }
  return arr;
}

void
preoffset_array(
    PyArrayObject* array,
    int value) {

  npy_intp  size;
  double   *data;

  if (value == 1) {
    return;
  }

  size = PyArray_Size((PyObject*)array);
  data = (double*)PyArray_DATA(array);
  offset_c_array(data, size, (double)(1 - value));
}

void
unoffset_array(
    PyArrayObject* array,
    int value) {

  npy_intp  size;
  double   *data;

  if (value == 1) {
    return;
  }

  size = PyArray_Size((PyObject*)array);
  data = (double*)PyArray_DATA(array);
  offset_c_array(data, size, (double)-(1 - value));
}

void
copy_array_to_c_double(
    PyArrayObject* array,
    double* dest) {

  npy_intp size = 1;
  double*  data = NULL;

  size = PyArray_Size((PyObject*)array);
  data = (double*)PyArray_DATA(array);

  memcpy(dest, data, size * sizeof(double));
}

void
copy_array_to_c_int(
    PyArrayObject* array,
    int* dest) {

  npy_intp size = 1;
  int*     data = NULL;

  size = PyArray_Size((PyObject*)array);
  data = (int*)PyArray_DATA(array);

  memcpy(dest, data, size * sizeof(int));
}

int
is_null(
    /*@null@*/ void *p) {

  if (p == NULL) {
    PyErr_SetString(PyExc_AssertionError, "Underlying object is NULL.");
    return 1;
  }
  return 0;
}

/* wcslib represents undefined values using its own special constant,
   UNDEFINED, and to be consistent with NumPy/Python we expose them as NaN.

   Historically the struct was stored "canonically" in NaN form and these two
   functions flipped the whole struct to/from UNDEFINED in place around every
   wcslib call.  That in-place flip of a shared struct is not thread-safe
   (GH-16409).  The struct is now stored canonically in wcslib's native
   UNDEFINED form, and the NaN<->UNDEFINED translation happens at the Python
   attribute boundary (the getters/setters and WCSParameterArray) instead.

   These two functions are retained as no-ops only because they are part of
   astropy's exported WCS C API (see astropy_wcs_api.c).  They no longer
   mutate the struct.
*/

void
wcsprm_c2python(
    /*@null@*/ struct wcsprm* x) {
  (void)x;
}

void
wcsprm_python2c(
    /*@null@*/ struct wcsprm* x) {
  (void)x;
}

/***************************************************************************
 * Exceptions                                                              *
 ***************************************************************************/

PyObject* WcsExc_Wcs;
PyObject* WcsExc_SingularMatrix;
PyObject* WcsExc_InconsistentAxisTypes;
PyObject* WcsExc_InvalidTransform;
PyObject* WcsExc_InvalidCoordinate;
PyObject* WcsExc_NoSolution;
PyObject* WcsExc_InvalidSubimageSpecification;
PyObject* WcsExc_NonseparableSubimageCoordinateSystem;
PyObject* WcsExc_NoWcsKeywordsFound;
PyObject* WcsExc_InvalidTabularParameters;
PyObject* WcsExc_InvalidPrjParameters;

/* This is an array mapping the wcs status codes to Python exception
 * types.  The exception string is stored as part of wcslib itself in
 * wcs_errmsg.
 */
PyObject** wcs_errexc[14];

static PyObject*
_new_exception_with_doc(char *name, char *doc, PyObject *base)
{
  return PyErr_NewExceptionWithDoc(name, doc, base, NULL);
}

#define DEFINE_EXCEPTION(exc) \
  WcsExc_##exc = _new_exception_with_doc(                             \
      "astropy.wcs._wcs." #exc "Error",                                 \
      doc_##exc,                                                        \
      WcsExc_Wcs);                                                      \
  if (WcsExc_##exc == NULL) \
    return 1; \
  PyModule_AddObject(m, #exc "Error", WcsExc_##exc); \

int
_define_exceptions(
    PyObject* m) {

  WcsExc_Wcs = _new_exception_with_doc(
      "astropy.wcs._wcs.WcsError",
      doc_WcsError,
      PyExc_ValueError);
  if (WcsExc_Wcs == NULL) {
    return 1;
  }
  PyModule_AddObject(m, "WcsError", WcsExc_Wcs);

  DEFINE_EXCEPTION(SingularMatrix);
  DEFINE_EXCEPTION(InconsistentAxisTypes);
  DEFINE_EXCEPTION(InvalidTransform);
  DEFINE_EXCEPTION(InvalidCoordinate);
  DEFINE_EXCEPTION(NoSolution);
  DEFINE_EXCEPTION(InvalidSubimageSpecification);
  DEFINE_EXCEPTION(NonseparableSubimageCoordinateSystem);
  DEFINE_EXCEPTION(NoWcsKeywordsFound);
  DEFINE_EXCEPTION(InvalidTabularParameters);
  DEFINE_EXCEPTION(InvalidPrjParameters);
  return 0;
}

const char*
wcslib_get_error_message(int status) {
  return wcs_errmsg[status];
}

void
wcserr_to_python_exc(const struct wcserr *err) {
  PyObject *exc;
  if (err == NULL) {
    PyErr_SetString(PyExc_RuntimeError, "NULL error object in wcslib");
  } else {
    if (err->status > 0 && err->status <= WCS_ERRMSG_MAX) {
      exc = *wcs_errexc[err->status];
    } else {
      exc = PyExc_RuntimeError;
    }
    /* This is technically not thread-safe -- make sure we have the GIL */
    wcsprintf_set(NULL);
    wcserr_prt(err, "");
    PyErr_SetString(exc, wcsprintf_buf());
  }
}

void
wcs_to_python_exc(const struct wcsprm *wcs) {
  PyObject* exc;
  const struct wcserr *err = wcs->err;
  if (err == NULL) {
    PyErr_SetString(PyExc_RuntimeError, "NULL error object in wcslib");
  } else {
    if (err->status > 0 && err->status < WCS_ERRMSG_MAX) {
      exc = *wcs_errexc[err->status];
    } else {
      exc = PyExc_RuntimeError;
    }
    /* This is technically not thread-safe -- make sure we have the GIL */
    wcsprintf_set(NULL);
    wcsperr(wcs, "");
    PyErr_SetString(exc, wcsprintf_buf());
  }
}

void
wcserr_fix_to_python_exc(const struct wcserr *err) {
  PyObject *exc;
  if (err == NULL) {
    PyErr_SetString(PyExc_RuntimeError, "NULL error object in wcslib");
  } else {
    if (err->status > 0 && err->status <= FIXERR_NO_REF_PIX_VAL) {
      exc = PyExc_ValueError;
    } else {
      exc = PyExc_RuntimeError;
    }
    /* This is technically not thread-safe -- make sure we have the GIL */
    wcsprintf_set(NULL);
    wcserr_prt(err, "");
    PyErr_SetString(exc, wcsprintf_buf());
  }
}

void
wcshdr_err_to_python_exc(int status, const struct wcsprm *wcs) {
  /* Add error to wcslib error buffer */
  wcsperr(wcs, NULL);
  if (status > 0 && status != WCSHDRERR_PARSER) {
    PyErr_Format(
      PyExc_MemoryError,
      "Memory allocation error:\n%s",
      wcsprintf_buf()
    );
  } else {
    PyErr_Format(
      PyExc_ValueError,
      "Internal error in wcslib header parser:\n %s",
      wcsprintf_buf()
    );
  }
}


/***************************************************************************
  Property helpers
 ***************************************************************************/

#define SHAPE_STR_LEN 2048

/* Helper function to display the desired shape of an array as a
   string, eg. 2x2 */
static void
shape_to_string(
    int ndims,
    const npy_intp* dims,
    char* str /* [SHAPE_STR_LEN] */) {

  int i;
  char value[32]; /* More than large enough to hold string rep of a
                     64-bit integer (way overkill) */

  if (ndims > 3) {
    strncpy(str, "ERROR", 6);
    return;
  }

  str[0] = 0;
  for (i = 0; i < ndims; ++i) {
      snprintf(value, 32, "%d", (int)dims[i]);
    strncat(str, value, 32);
    if (i != ndims - 1) {
      strncat(str, "x", 2);
    }
  }
}

/* get_string is inlined */

int
set_string(
    const char* propname,
    PyObject* value,
    char* dest,
    Py_ssize_t maxlen) {

  char*      buffer;
  Py_ssize_t len;
  PyObject*  ascii_obj = NULL;
  int        result = -1;

  if (check_delete(propname, value)) {
    return -1;
  }

  if (PyUnicode_Check(value)) {
    ascii_obj = PyUnicode_AsASCIIString(value);
    if (ascii_obj == NULL) {
      goto end;
    }
    if (PyBytes_AsStringAndSize(ascii_obj, &buffer, &len) == -1) {
      goto end;
    }
  } else if (PyBytes_Check(value)) {
    if (PyBytes_AsStringAndSize(value, &buffer, &len) == -1) {
      goto end;
    }
  } else {
    PyErr_SetString(PyExc_TypeError, "'value' must be bytes or unicode.");
    goto end;
  }

  if (len >= maxlen) {
    PyErr_Format(
        PyExc_ValueError,
        "'%s' length must be less than %u characters.",
        propname,
        (unsigned int) maxlen);
    goto end;
  }

  strncpy(dest, buffer, (size_t)len + 1);
  result = 0;

 end:
  Py_XDECREF(ascii_obj);
  return result;
}

/* get_bool is inlined */

int
set_bool(
    const char* propname,
    PyObject* value,
    int* dest) {

  if (check_delete(propname, value)) {
    return -1;
  }

  *dest = PyObject_IsTrue(value);

  return 0;
}

/* get_int is inlined */

int
set_int(
    const char* propname,
    PyObject* value,
    int* dest) {
  long value_int;

  if (check_delete(propname, value)) {
    return -1;
  }

  value_int = PyLong_AsLong(value);
  if (value_int == -1 && PyErr_Occurred()) {
    return -1;
  }

  if ((unsigned long)value_int > 0x7fffffff) {
    PyErr_SetString(PyExc_OverflowError, "integer value too large");
    return -1;
  }

  *dest = (int)value_int;

  return 0;
}

/* get_double is inlined */

int
set_double(
    const char* propname,
    PyObject* value,
    double* dest) {

  if (check_delete(propname, value)) {
    return -1;
  }

  double v = PyFloat_AsDouble(value);

  if (PyErr_Occurred()) {
    return -1;
  }

  /* Store NaN as WCSLIB's native UNDEFINED (GH-16409). */
  *dest = npy_isnan(v) ? UNDEFINED : v;
  return 0;
}

/* get_double_array is inlined */

int
_set_double_array(
    const char* propname,
    PyObject* value,
    int ndims,
    const npy_intp* dims,
    double* dest,
    int to_undefined) {

  PyArrayObject* value_array = NULL;
  npy_int        i           = 0;
  char           shape_str[SHAPE_STR_LEN];

  if (check_delete(propname, value)) {
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value, NPY_DOUBLE,
                                                          ndims, ndims);
  if (value_array == NULL) {
    return -1;
  }

  if (dims != NULL) {
    for (i = 0; i < ndims; ++i) {
      if (PyArray_DIM(value_array, i) != dims[i]) {
        shape_to_string(ndims, dims, shape_str);
        PyErr_Format(
            PyExc_ValueError,
            "'%s' array is the wrong shape, must be %s",
            propname, shape_str);
        Py_DECREF(value_array);
        return -1;
      }
    }
  }

  if (to_undefined) {
    /* Auxiliary UNDEFINED-capable field: translate NaN -> UNDEFINED. */
    npy_intp n = PyArray_SIZE(value_array);
    const double* src = (const double*)PyArray_DATA(value_array);
    for (npy_intp j = 0; j < n; ++j) {
      dest[j] = npy_isnan(src[j]) ? UNDEFINED : src[j];
    }
  } else {
    /* Core/derived field: stored verbatim (NaN stays NaN). */
    copy_array_to_c_double(value_array, dest);
  }

  Py_DECREF(value_array);

  return 0;
}

int
set_double_array(
    const char* propname,
    PyObject* value,
    int ndims,
    const npy_intp* dims,
    double* dest) {
  return _set_double_array(propname, value, ndims, dims, dest, 0);
}

int
set_double_array_undefined(
    const char* propname,
    PyObject* value,
    int ndims,
    const npy_intp* dims,
    double* dest) {
  return _set_double_array(propname, value, ndims, dims, dest, 1);
}

int
set_int_array(
    const char* propname,
    PyObject* value,
    int ndims,
    const npy_intp* dims,
    int* dest) {
  PyArrayObject* value_array = NULL;
  npy_int        i           = 0;
  char           shape_str[SHAPE_STR_LEN];

  if (check_delete(propname, value)) {
    return -1;
  }

  value_array = (PyArrayObject*)PyArray_ContiguousFromAny(value, NPY_INT,
                                                          ndims, ndims);
  if (value_array == NULL) {
    return -1;
  }

  if (dims != NULL) {
    for (i = 0; i < ndims; ++i) {
      if (PyArray_DIM(value_array, i) != dims[i]) {
        shape_to_string(ndims, dims, shape_str);
        PyErr_Format(
            PyExc_ValueError,
            "'%s' array is the wrong shape, must be %s",
            propname, shape_str);
        Py_DECREF(value_array);
        return -1;
      }
    }
  }

  copy_array_to_c_int(value_array, dest);

  Py_DECREF(value_array);

  return 0;
}

/* get_str_list is inlined */

int
set_str_list(
    const char* propname,
    PyObject* value,
    Py_ssize_t len,
    Py_ssize_t maxlen,
    char (*dest)[72]) {

  PyObject*  str      = NULL;
  Py_ssize_t input_len;
  Py_ssize_t i        = 0;

  if (check_delete(propname, value)) {
    return -1;
  }

  if (maxlen == 0) {
    maxlen = 68;
  }

  if (!PySequence_Check(value)) {
    PyErr_Format(
        PyExc_TypeError,
        "'%s' must be a sequence of strings",
        propname);
    return -1;
  }

  if (PySequence_Size(value) != len) {
    PyErr_Format(
        PyExc_ValueError,
        "len(%s) must be %u",
        propname,
        (unsigned int)len);
    return -1;
  }

  /* We go through the list twice, once to verify that the list is
     in the correct format, and then again to do the data copy.  This
     way, we won't partially copy the contents and then throw an
     exception. */
  for (i = 0; i < len; ++i) {
    str = PySequence_GetItem(value, i);
    if (str == NULL) {
      return -1;
    }

    if (!(PyBytes_CheckExact(str) || PyUnicode_CheckExact(str))) {
      PyErr_Format(
          PyExc_TypeError,
          "'%s' must be a sequence of bytes or strings",
          propname);
      Py_DECREF(str);
      return -1;
    }

    input_len = PySequence_Size(str);
    if (input_len > maxlen) {
      PyErr_Format(
          PyExc_ValueError,
          "Each entry in '%s' must be less than %u characters",
          propname, (unsigned int)maxlen);
      Py_DECREF(str);
      return -1;
    } else if (input_len == -1) {
      Py_DECREF(str);
      return -1;
    }

    Py_DECREF(str);
  }

  for (i = 0; i < len; ++i) {
    str = PySequence_GetItem(value, i);
    if (str == NULL) {
      /* Theoretically, something has gone really wrong here, since
         we've already verified the list. */
      PyErr_Clear();
      PyErr_Format(
          PyExc_RuntimeError,
          "Input values have changed underneath us.  Something is seriously wrong.");
      return -1;
    }

    if (set_string(propname, str, dest[i], maxlen)) {
      PyErr_Clear();
      PyErr_Format(
          PyExc_RuntimeError,
          "Input values have changed underneath us.  Something is seriously wrong.");
      Py_DECREF(str);
      return -1;
    }

    Py_DECREF(str);
  }

  return 0;
}


/*@null@*/ PyObject*
get_pscards(
    /*@unused@*/ const char* propname,
    struct pscard* ps,
    int nps) {

  PyObject*  result    = NULL;
  PyObject*  subresult = NULL;
  Py_ssize_t i         = 0;

  if (nps < 0) {
    nps = 0;
  }

  result = PyList_New((Py_ssize_t)nps);
  if (result == NULL) {
    return NULL;
  }

  if (nps && ps == NULL) {
    PyErr_SetString(PyExc_MemoryError, "NULL pointer");
    return NULL;
  }

  for (i = 0; i < (Py_ssize_t)nps; ++i) {
    subresult = Py_BuildValue("iis", ps[i].i, ps[i].m, ps[i].value);
    if (subresult == NULL) {
      Py_DECREF(result);
      return NULL;
    }

    if (PyList_SetItem(result, i, subresult)) {
      Py_DECREF(result);
      return NULL;
    }
  }

  return result;
}

int
set_pscards(
    /*@unused@*/ const char* propname,
    PyObject* value,
    struct pscard** ps,
    int *nps,
    int *npsmax) {

  PyObject*   subvalue  = NULL;
  Py_ssize_t  i         = 0;
  Py_ssize_t  size      = 0;
  int         ival      = 0;
  int         mval      = 0;
  const char* strvalue  = 0;
  void*       newmem    = NULL;

  if (!PySequence_Check(value))
    return -1;
  size = PySequence_Size(value);
  if (size > 0x7fffffff) {
    /* Must be a 32-bit size */
    return -1;
  }

  if (size > (Py_ssize_t)*npsmax) {
    newmem = malloc(sizeof(struct pscard) * size);
    if (newmem == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Could not allocate memory.");
      return -1;
    }
    free(*ps);
    *ps = newmem;
    *npsmax = (int)size;
  }

  /* Verify the entire list for correct types first, so we don't have
     to undo anything copied into the canonical array. */
  for (i = 0; i < size; ++i) {
    subvalue = PySequence_GetItem(value, i);
    if (subvalue == NULL) {
      return -1;
    }
    if (!PyArg_ParseTuple(subvalue, "iis", &ival, &mval, &strvalue)) {
      Py_DECREF(subvalue);
      return -1;
    }
    Py_DECREF(subvalue);
  }

  for (i = 0; i < size; ++i) {
    subvalue = PySequence_GetItem(value, i);
    if (subvalue == NULL) {
      return -1;
    }
    if (!PyArg_ParseTuple(subvalue, "iis", &ival, &mval, &strvalue)) {
      Py_DECREF(subvalue);
      return -1;
    }
    Py_DECREF(subvalue);

    (*ps)[i].i = ival;
    (*ps)[i].m = mval;
    strncpy((*ps)[i].value, strvalue, 72);
    (*ps)[i].value[71] = '\0';
    (*nps) = (int)(i + 1);
  }

  return 0;
}

/*@null@*/ PyObject*
get_pvcards(
    /*@unused@*/ const char* propname,
    struct pvcard* pv,
    int npv) {

  PyObject*  result    = NULL;
  PyObject*  subresult = NULL;
  Py_ssize_t i         = 0;

  if (npv < 0) {
    npv = 0;
  }

  result = PyList_New((Py_ssize_t)npv);
  if (result == NULL) {
    return NULL;
  }

  if (npv && pv == NULL) {
    PyErr_SetString(PyExc_MemoryError, "NULL pointer");
    return NULL;
  }

  for (i = 0; i < (Py_ssize_t)npv; ++i) {
    subresult = Py_BuildValue("iid", pv[i].i, pv[i].m, pv[i].value);
    if (subresult == NULL) {
      Py_DECREF(result);
      return NULL;
    }

    if (PyList_SetItem(result, i, subresult)) {
      Py_DECREF(result);
      return NULL;
    }
  }

  return result;
}

int
set_pvcards(
    /*@propname@*/ const char* propname,
    PyObject* value,
    struct pvcard** pv,
    int *npv,
    int *npvmax) {

  PyObject* fastseq = NULL;
  struct pvcard* newmem = NULL;
  Py_ssize_t size;
  int ret = -1;
  int i;

  fastseq = PySequence_Fast(value, "Expected sequence type");
  if (!fastseq)
    goto done;

  size = PySequence_Size(value);
  newmem = malloc(sizeof(struct pvcard) * size);

  /* Raise exception if size is nonzero but newmem
   * could not be allocated. */
  if (size && !newmem) {
    PyErr_SetString(PyExc_MemoryError, "Could not allocate memory.");
    return -1;
  }

  PyObject* item = NULL;
  for (i = 0; i < size; ++i)
  {
    if (!PyArg_ParseTuple((item = PySequence_GetItem(value, i)), "iid",
        &newmem[i].i, &newmem[i].m, &newmem[i].value))
    {
      Py_DECREF(item);
      goto done;
    }
  }

  if (size <= (Py_ssize_t)*npvmax) {
    memcpy(*pv, newmem, sizeof(struct pvcard) * size);
  } else { /* (size > (Py_ssize_t)*npvmax) */
    free(*pv);
    *npv = (int)size;
    *pv = newmem;
    newmem = NULL;
  }
  *npv = (int)size;

  ret = 0;
done:
  Py_XDECREF(fastseq);
  free(newmem);
  return ret;
}

PyObject*
get_deepcopy(
    PyObject* obj,
    PyObject* memo) {

  if (PyObject_HasAttrString(obj, "__deepcopy__")) {
    return PyObject_CallMethod(obj, "__deepcopy__", "O", memo);
  } else {
    return PyObject_CallMethod(obj, "__copy__", "");
  }
}

/***************************************************************************
 * Miscellaneous helper functions                                          *
 ***************************************************************************/

int
parse_unsafe_unit_conversion_spec(
    const char* arg, int* ctrl) {

  const char* p = NULL;

  *ctrl = 0;

  for (p = arg; *p != '\0'; ++p) {
    switch (*p) {
    case 's':
    case 'S':
      *ctrl |= 1;
      break;
    case 'h':
    case 'H':
      *ctrl |= 2;
      break;
    case 'd':
    case 'D':
      *ctrl |= 4;
      break;
    default:
      PyErr_SetString(
          PyExc_ValueError,
          "translate_units may only contain the characters 's', 'h' or 'd'");
      return 1;
    }
  }

  return 0;
}
