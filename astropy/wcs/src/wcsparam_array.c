/*
 WCSParameterArray: a NumPy ndarray subclass used to expose the auxiliary
 wcsprm array parameters that WCSLIB represents with its UNDEFINED sentinel
 (obsgeo, crder, csyer, cperi, czphs, mjdref).  It holds a copy of the field
 with UNDEFINED translated to NaN, and writes element assignments back into the
 owning wcsprm field translating NaN -> UNDEFINED (GH-16409).

 It is a heap type (PyType_FromSpecWithBases) that delegates to ndarray's slots
 via PyType_GetSlot, with the write-back context carried on a PyCapsule, so the
 whole thing is compatible with the limited C API / abi3 builds.

 The core/derived array parameters (crpix, crval, cdelt, pc, cd, ...) are NOT
 exposed this way -- WCSLIB never treats them as UNDEFINED, so they are plain
 writeable views (see get_double_array in pyutil.h).
*/

#define NO_IMPORT_ARRAY

#include "astropy_wcs/pyutil.h"
#include "astropy_wcs/wcslib_wrap.h"  /* for the Wcsprm struct (write-back) */

#include <stdlib.h>  /* malloc, free */

/* Write-back context for an auxiliary WCSParameterArray.  It is carried on the
 * array's base object as a PyCapsule rather than a custom C type, because a
 * capsule is fully compatible with the limited C API / abi3 builds (a static
 * PyTypeObject is not). */
typedef struct {
  PyObject*   parent;    /* owning Wcsprm, holds a reference */
  double*     dest;      /* &parent->x.<field> */
  npy_intp    nelem;
} WCSParamWriteback;

static const char* const WCSPARAM_CAPSULE_NAME = "astropy.wcs._wcsparam_writeback";

static void
wcsparam_writeback_destructor(PyObject* capsule) {
  WCSParamWriteback* wb =
      (WCSParamWriteback*)PyCapsule_GetPointer(capsule, WCSPARAM_CAPSULE_NAME);
  if (wb != NULL) {
    Py_XDECREF(wb->parent);
    free(wb);
  }
}

/* Created at setup via PyType_FromSpecWithBases so the whole subclass works
 * under the limited API.  ndarray's own subscript slots are fetched with
 * PyType_GetSlot and delegated to. */
static PyObject*     WCSParameterArray_Type = NULL;
static binaryfunc    ndarray_subscript = NULL;
static objobjargproc ndarray_ass_subscript = NULL;

/* Read indexing decays to a plain ndarray, so a slice or element of an
 * auxiliary array is an ordinary array and downstream code (e.g.
 * astropy.units) never has to handle the subclass -- which also means the
 * subclass needs no instance __dict__. */
static PyObject*
WCSParameterArray_subscript(PyObject* self, PyObject* key) {
  PyObject* result = ndarray_subscript(self, key);
  if (result != NULL && PyArray_Check(result) &&
      Py_TYPE(result) == (PyTypeObject*)WCSParameterArray_Type) {
    PyObject* plain = PyArray_View((PyArrayObject*)result, NULL, &PyArray_Type);
    Py_DECREF(result);
    return plain;
  }
  return result;
}

static int
WCSParameterArray_ass_subscript(PyObject* self, PyObject* key, PyObject* value) {
  PyObject* base = PyArray_BASE((PyArrayObject*)self);
  WCSParamWriteback* wb = NULL;
  if (base != NULL && PyCapsule_IsValid(base, WCSPARAM_CAPSULE_NAME)) {
    wb = (WCSParamWriteback*)PyCapsule_GetPointer(base, WCSPARAM_CAPSULE_NAME);
  }

  /* Only the interceptable assignment paths (a[i]=, a[:]=) reach here; writes
   * that bypass __setitem__ (np.fill_diagonal, a.flat[i]=, a+=) are not synced
   * back.  Making every write path safe would mean returning read-only arrays
   * and steering users to whole-attribute assignment -- left to a follow-up. */

  /* Perform the normal ndarray assignment first (NaN is allowed). */
  if (ndarray_ass_subscript(self, key, value) != 0) {
    return -1;
  }

  /* Then sync the whole (small) buffer back into the wcsprm field,
   * translating NaN -> UNDEFINED. */
  if (wb != NULL) {
    const double* src = (const double*)PyArray_DATA((PyArrayObject*)self);
    for (npy_intp i = 0; i < wb->nelem; ++i) {
      wb->dest[i] = npy_isnan(src[i]) ? UNDEFINED : src[i];
    }
    /* Equivalent to note_change(): force a wcsset on next use. */
    ((Wcsprm*)wb->parent)->x.flag = 0;
  }
  return 0;
}

static PyType_Slot WCSParameterArray_slots[] = {
  {Py_mp_subscript, (void*)WCSParameterArray_subscript},
  {Py_mp_ass_subscript, (void*)WCSParameterArray_ass_subscript},
  {0, NULL}
};

static PyType_Spec WCSParameterArray_spec = {
  "astropy.wcs.WCSParameterArray",
  0,  /* basicsize == 0 inherits from the base type (ndarray) */
  0,  /* itemsize */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  WCSParameterArray_slots
};

int
_setup_wcsparameter_array_type(PyObject* m) {
  ndarray_subscript =
      (binaryfunc)PyType_GetSlot(&PyArray_Type, Py_mp_subscript);
  ndarray_ass_subscript =
      (objobjargproc)PyType_GetSlot(&PyArray_Type, Py_mp_ass_subscript);
  if (ndarray_subscript == NULL || ndarray_ass_subscript == NULL) {
    return -1;
  }

  PyObject* bases = PyTuple_Pack(1, (PyObject*)&PyArray_Type);
  if (bases == NULL) {
    return -1;
  }
  WCSParameterArray_Type =
      PyType_FromSpecWithBases(&WCSParameterArray_spec, bases);
  Py_DECREF(bases);
  if (WCSParameterArray_Type == NULL) {
    return -1;
  }
  return PyModule_AddObjectRef(m, "WCSParameterArray", WCSParameterArray_Type);
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
WCSParameterArray_New(PyObject* owner, int ndims,
                      const npy_intp* dims, double* value) {
  npy_intp nelem = 0;
  PyObject* arr = new_double_array((PyTypeObject*)WCSParameterArray_Type, ndims,
                                   dims, value, &nelem);
  if (arr == NULL) {
    return NULL;
  }

  WCSParamWriteback* wb = (WCSParamWriteback*)malloc(sizeof(WCSParamWriteback));
  if (wb == NULL) {
    Py_DECREF(arr);
    return PyErr_NoMemory();
  }
  Py_INCREF(owner);
  wb->parent = owner;
  wb->dest = value;
  wb->nelem = nelem;

  PyObject* capsule = PyCapsule_New(wb, WCSPARAM_CAPSULE_NAME,
                                    wcsparam_writeback_destructor);
  if (capsule == NULL) {
    Py_DECREF(owner);
    free(wb);
    Py_DECREF(arr);
    return NULL;
  }
  /* SetBaseObject steals the capsule reference on success. */
  if (PyArray_SetBaseObject((PyArrayObject*)arr, capsule) < 0) {
    Py_DECREF(capsule);  /* runs the destructor: DECREF owner + free wb */
    Py_DECREF(arr);
    return NULL;
  }
  return arr;
}
