/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __PYUTIL_H__
#define __PYUTIL_H__

#include "util.h"

#define PY_ARRAY_UNIQUE_SYMBOL astropy_wcs_numpy_api

#include <Python.h>

#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

#if PY_MAJOR_VERSION >= 3
#define PY3K 1
#else
#define PY3K 0
#ifndef Py_TYPE
  #define Py_TYPE(ob) (((PyObject*)(ob))->ob_type)
#endif
#endif

PyObject*
PyArrayProxy_New(
    PyObject* self,
    int nd,
    const npy_intp* dims,
    int typenum,
    const void* data);

PyObject*
PyArrayReadOnlyProxy_New(
    PyObject* self,
    int nd,
    const npy_intp* dims,
    int typenum,
    const void* data);

/*@null@*/ PyObject *
PyStrListProxy_New(
    PyObject* owner,
    Py_ssize_t size,
    Py_ssize_t maxsize,
    char (*array)[72]
    );

int
_setup_str_list_proxy_type(
    PyObject* m);

static INLINE void
offset_c_array(
    double* value,
    npy_intp size,
    double offset) {
  double* end = value + size;

  for ( ; value != end; ++value) {
    *value += offset;
  }
}

static INLINE
void nan2undefined(
    double* value,
    unsigned int nvalues) {

  double* end = value + nvalues;

  for ( ; value != end; ++value) {
    if (isnan64(*value)) {
      *value = UNDEFINED;
    }
  }
}

static INLINE
void undefined2nan(
    double* value,
    unsigned int nvalues) {

  double* end = value + nvalues;

  for ( ; value != end; ++value) {
    if (*value == UNDEFINED) {
      *value = (double)NPY_NAN;
    }
  }
}

void
preoffset_array(
    PyArrayObject* array,
    int value);

void
unoffset_array(
    PyArrayObject* array,
    int value);

void
copy_array_to_c_double(
    PyArrayObject* array,
    double* dest);

void
copy_array_to_c_int(
    PyArrayObject* array,
    int* dest);

/**
 Returns TRUE if pointer is NULL, and sets Python exception
*/
int
is_null(/*@null@*/ void *);

typedef void (*value_fixer_t)(double*, unsigned int);

void
wcsprm_c2python(
    /*@null@*/ struct wcsprm* x);

void
wcsprm_python2c(
    /*@null@*/ struct wcsprm* x);

/***************************************************************************
 * Exceptions                                                              *
 ***************************************************************************/

extern PyObject* WcsExc_SingularMatrix;
extern PyObject* WcsExc_InconsistentAxisTypes;
extern PyObject* WcsExc_InvalidTransform;
extern PyObject* WcsExc_InvalidCoordinate;
extern PyObject* WcsExc_NoSolution;
extern PyObject* WcsExc_InvalidSubimageSpecification;
extern PyObject* WcsExc_NonseparableSubimageCoordinateSystem;
extern PyObject* WcsExc_NoWcsKeywordsFound;
extern PyObject* WcsExc_InvalidTabularParameters;

/* This is an array mapping the wcs status codes to Python exception
 * types.  The exception string is stored as part of wcslib itself in
 * wcs_errmsg.
 */
extern PyObject** wcs_errexc[14];
#define WCS_ERRMSG_MAX 14
#define WCSFIX_ERRMSG_MAX 11

int
_define_exceptions(PyObject* m);

const char*
wcslib_get_error_message(int stat);

void
wcserr_to_python_exc(const struct wcserr *err);

void
wcs_to_python_exc(const struct wcsprm *wcs);

void
wcserr_fix_to_python_exc(const struct wcserr *err);

/***************************************************************************
  Property helpers
 ***************************************************************************/
static INLINE int
check_delete(
    const char* propname,
    PyObject* value) {

  if (value == NULL) {
    PyErr_Format(PyExc_TypeError, "'%s' can not be deleted", propname);
    return -1;
  }

  return 0;
}

static INLINE PyObject*
get_string(
    /*@unused@*/ const char* propname,
    const char* value) {

  #if PY3K
  return PyUnicode_FromString(value);
  #else
  return PyBytes_FromString(value);
  #endif
}

int
set_string(
    const char* propname,
    PyObject* value,
    char* dest,
    Py_ssize_t maxlen);

static INLINE PyObject*
get_bool(
    /*@unused@*/ const char* propname,
    long value) {

  return PyBool_FromLong(value);
}

int
set_bool(
    const char* propname,
    PyObject* value,
    int* dest);

static INLINE PyObject*
get_int(
    /*@unused@*/ const char* propname,
    long value) {

  #if PY3K
  return PyLong_FromLong(value);
  #else
  return PyInt_FromLong(value);
  #endif
}

int
set_int(
    const char* propname,
    PyObject* value,
    int* dest);

static INLINE PyObject*
get_double(
    const char* propname,
    double value) {

  return PyFloat_FromDouble(value);
}

int
set_double(
    const char* propname,
    PyObject* value,
    double* dest);

/*@null@*/ static INLINE PyObject*
get_double_array(
    /*@unused@*/ const char* propname,
    double* value,
    int ndims,
    const npy_intp* dims,
    /*@shared@*/ PyObject* owner) {

  return PyArrayProxy_New(owner, ndims, dims, PyArray_DOUBLE, value);
}

/*@null@*/ static INLINE PyObject*
get_double_array_readonly(
    /*@unused@*/ const char* propname,
    double* value,
    int ndims,
    const npy_intp* dims,
    /*@shared@*/ PyObject* owner) {

  return PyArrayReadOnlyProxy_New(owner, ndims, dims, PyArray_DOUBLE, value);
}

int
set_double_array(
    const char* propname,
    PyObject* value,
    int ndims,
    const npy_intp* dims,
    double* dest);

/*@null@*/ static INLINE PyObject*
get_int_array(
    /*@unused@*/ const char* propname,
    int* value,
    int ndims,
    const npy_intp* dims,
    /*@shared@*/ PyObject* owner) {

  return PyArrayProxy_New(owner, ndims, dims, PyArray_INT, value);
}

int
set_int_array(
    const char* propname,
    PyObject* value,
    int ndims,
    const npy_intp* dims,
    int* dest);

static INLINE PyObject*
get_str_list(
    /*@unused@*/ const char* propname,
    char (*array)[72],
    Py_ssize_t len,
    Py_ssize_t maxlen,
    PyObject* owner) {

  return PyStrListProxy_New(owner, len, maxlen, array);
}

int
set_str_list(
    const char* propname,
    PyObject* value,
    Py_ssize_t len,
    Py_ssize_t maxlen,
    char (*dest)[72]);

PyObject*
get_pscards(
    const char* propname,
    struct pscard* ps,
    int nps);

int
set_pscards(
    const char* propname,
    PyObject* value,
    struct pscard** ps,
    int *nps,
    int *npsmax);

PyObject*
get_pvcards(
    const char* propname,
    struct pvcard* pv,
    int npv);

int
set_pvcards(
    const char* propname,
    PyObject* value,
    struct pvcard** pv,
    int *npv,
    int *npvmax);

PyObject*
get_deepcopy(
    PyObject* obj,
    PyObject* memo);

/***************************************************************************
  Miscellaneous helper functions
 ***************************************************************************/

int
parse_unsafe_unit_conversion_spec(
    const char* arg, int* ctrl);

#endif /* __PYUTIL_H__ */
