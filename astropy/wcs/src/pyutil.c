/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

/* util.h must be imported first */
#include "astropy_wcs/pyutil.h"

#include "astropy_wcs/docstrings.h"

#include "wcsfix.h"
#include "wcshdr.h"
#include "wcsprintf.h"
#include "wcsunits.h"

/*@null@*/ static INLINE PyObject*
_PyArrayProxy_New(
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
PyArrayProxy_New(
    /*@shared@*/ PyObject* self,
    int nd,
    const npy_intp* dims,
    int typenum,
    const void* data) {

  return _PyArrayProxy_New(self, nd, dims, typenum, data, NPY_ARRAY_WRITEABLE);
}

/*@null@*/ PyObject*
PyArrayReadOnlyProxy_New(
    /*@shared@*/ PyObject* self,
    int nd,
    const npy_intp* dims,
    int typenum,
    const void* data) {

  return _PyArrayProxy_New(self, nd, dims, typenum, data, 0);
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
   UNDEFINED.  To be consistent with the Pythonic way of doing things,
   it's nicer to represent undefined values using NaN.  Unfortunately,
   in order to get nice mutable arrays in Python, Python must be able
   to edit the wcsprm values directly.  The solution is to store NaNs
   in the struct "canonically", but convert those NaNs to/from
   UNDEFINED around every call into a wcslib function.  It's not as
   computationally expensive as it sounds, as all these arrays are
   quite small.
*/

static INLINE void
wcsprm_fix_values(
    struct wcsprm* x,
    value_fixer_t value_fixer) {

  unsigned int naxis = (unsigned int)x->naxis;

  value_fixer(x->cd, naxis * naxis);
  value_fixer(x->cdelt, naxis);
  value_fixer(x->crder, naxis);
  value_fixer(x->crota, naxis);
  value_fixer(x->crpix, naxis);
  value_fixer(x->crval, naxis);
  value_fixer(x->csyer, naxis);
  value_fixer(&x->equinox, 1);
  value_fixer(&x->latpole, 1);
  value_fixer(&x->lonpole, 1);
  value_fixer(&x->mjdavg, 1);
  value_fixer(&x->mjdobs, 1);
  value_fixer(x->obsgeo, 3);
  value_fixer(&x->cel.phi0, 1);
  value_fixer(&x->restfrq, 1);
  value_fixer(&x->restwav, 1);
  value_fixer(&x->cel.theta0, 1);
  value_fixer(&x->velangl, 1);
  value_fixer(&x->velosys, 1);
  value_fixer(&x->zsource, 1);
}

void
wcsprm_c2python(
    /*@null@*/ struct wcsprm* x) {

  if (x != NULL) {
    wcsprm_fix_values(x, &undefined2nan);
  }
}

void
wcsprm_python2c(
    /*@null@*/ struct wcsprm* x) {

  if (x != NULL) {
    wcsprm_fix_values(x, &nan2undefined);
  }
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

/* This is an array mapping the wcs status codes to Python exception
 * types.  The exception string is stored as part of wcslib itself in
 * wcs_errmsg.
 */
PyObject** wcs_errexc[14];

static PyObject*
_new_exception_with_doc(char *name, char *doc, PyObject *base)
{
#if ((PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION >= 7) || \
     (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION >= 2))
  return PyErr_NewExceptionWithDoc(name, doc, base, NULL);
#else
  /* Python 2.6 and 3.1 don't have PyErr_NewExceptionWithDoc */
  PyObject *dict;
  PyObject *docobj;
  int result;

  dict = PyDict_New();
  if (dict == NULL) {
    return NULL;
  }

  if (doc != NULL) {
    docobj = PyUnicode_FromString(doc);
    if (docobj == NULL) {
      Py_DECREF(dict);
      return NULL;
    }

    result = PyDict_SetItemString(dict, "__doc__", docobj);
    Py_DECREF(docobj);
    if (result < 0) {
      Py_DECREF(dict);
      return NULL;
    }

    return PyErr_NewException(name, base, dict);
  }
#endif
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
wcshdr_err_to_python_exc(int status) {
  if (status > 0 && status != WCSHDRERR_PARSER) {
    PyErr_SetString(PyExc_MemoryError, "Memory allocation error");
  } else {
    PyErr_SetString(PyExc_ValueError, "Internal error in wcslib header parser");
  }
}


/***************************************************************************
  Property helpers
 ***************************************************************************/

#define SHAPE_STR_LEN 128

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
    PyErr_SetString(PyExc_TypeError, "value must be bytes or unicode");
    goto end;
  }

  if (len > maxlen) {
    PyErr_Format(
        PyExc_ValueError,
        "'%s' must be less than %u characters",
        propname,
        (unsigned int)maxlen);
    goto end;
  }

  strncpy(dest, buffer, (size_t)maxlen);

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

  #if PY3K
  value_int = PyLong_AsLong(value);
  #else
  value_int = PyInt_AsLong(value);
  #endif
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

  *dest = PyFloat_AsDouble(value);

  if (PyErr_Occurred()) {
    return -1;
  } else {
    return 0;
  }
}

/* get_double_array is inlined */

int
set_double_array(
    const char* propname,
    PyObject* value,
    int ndims,
    const npy_intp* dims,
    double* dest) {

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

  copy_array_to_c_double(value_array, dest);

  Py_DECREF(value_array);

  return 0;
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
      Py_DECREF(subresult);
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
      Py_DECREF(subresult);
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

  size = PySequence_Fast_GET_SIZE(value);
  newmem = malloc(sizeof(struct pvcard) * size);

  /* Raise exception if size is nonzero but newmem
   * could not be allocated. */
  if (size && !newmem) {
    PyErr_SetString(PyExc_MemoryError, "Could not allocate memory.");
    return -1;
  }

  for (i = 0; i < size; ++i)
  {
    if (!PyArg_ParseTuple(PySequence_Fast_GET_ITEM(value, i), "iid",
        &newmem[i].i, &newmem[i].m, &newmem[i].value))
    {
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
