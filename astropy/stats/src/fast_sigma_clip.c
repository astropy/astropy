#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>
#include "compute_bounds.h"

/* Define docstrings */
static char module_docstring[] = "Fast sigma clipping";
static char _sigma_clip_fast_docstring[] = "Compute sigma clipping";

/* Declare the C functions here. */
static PyObject *_sigma_clip_fast(PyObject *self, PyObject *args);

/* Define the methods that will be available on the module. */
static PyMethodDef module_methods[] = {{"_sigma_clip_fast", _sigma_clip_fast,
                                        METH_VARARGS,
                                        _sigma_clip_fast_docstring},
                                       {NULL, NULL, 0, NULL}};

/* This is the function that is called on import. */

#define MOD_ERROR_VAL NULL
#define MOD_SUCCESS_VAL(val) val
#define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#define MOD_DEF(ob, name, doc, methods)                                        \
  static struct PyModuleDef moduledef = {                                      \
      PyModuleDef_HEAD_INIT, name, doc, -1, methods,                           \
  };                                                                           \
  ob = PyModule_Create(&moduledef);

MOD_INIT(_fast_sigma_clip) {
  PyObject *m;
  MOD_DEF(m, "_fast_sigma_clip", module_docstring, module_methods);
  if (m == NULL)
    return MOD_ERROR_VAL;
  import_array();
  return MOD_SUCCESS_VAL(m);
}

static PyObject *_sigma_clip_fast(PyObject *self, PyObject *args) {

  long n, m;
  int i, j;
  double *buffer;
  PyObject *data_obj, *mask_obj;
  PyArrayObject *data_array, *mask_array;
  double *data;
  uint8_t *mask;
  int use_median, maxiters;
  double sigma_lower, sigma_upper;
  int iteration, count;
  double lower, upper;
  npy_intp dims[2];
  PyArrayObject *arrays[2];
  npy_uint32 op_flags[2];

  PyObject *bounds_obj;
  PyArrayObject *bounds_array;
  double *bounds;

  NpyIter *iter;
  NpyIter_IterNextFunc *iternext;
  char **dataptr;
  npy_intp *strideptr, *innersizeptr, stride, innersize;
  PyArray_Descr *dtypes[2];

  // Parse the input tuple
  if (!PyArg_ParseTuple(args, "OOiidd", &data_obj, &mask_obj, &use_median,
                        &maxiters, &sigma_lower, &sigma_upper)) {
    PyErr_SetString(PyExc_TypeError, "Error parsing input");
    return NULL;
  }

  // Interpret the input data array as a `numpy` array.
  data_array = (PyArrayObject *)PyArray_FROM_O(data_obj);
  if (data_array == NULL) {
    PyErr_SetString(PyExc_TypeError, "Couldn't parse the input data array.");
    Py_XDECREF(data_array);
    return NULL;
  }

  // Interpret the input mask array as a `numpy` array.
  mask_array = (PyArrayObject *)PyArray_FROM_O(mask_obj);
  if (mask_array == NULL) {
    PyErr_SetString(PyExc_TypeError, "Couldn't parse the input mask array.");
    Py_XDECREF(mask_array);
    return NULL;
  }

  // Get the data array shape
  n = (long)PyArray_DIM(data_array, 0);
  m = (long)PyArray_DIM(data_array, 1);

  // Numpy iterator set-up - we use this rather that iterate manually over the
  // array so that data values are automatically cast to double.

  arrays[0] = data_array;
  arrays[1] = mask_array;

  op_flags[0] = NPY_ITER_READONLY;
  op_flags[1] = NPY_ITER_READONLY;

  dtypes[0] = PyArray_DescrFromType(NPY_DOUBLE);
  dtypes[1] = PyArray_DescrFromType(NPY_UINT8);

  iter = NpyIter_MultiNew(2, arrays,
                          NPY_ITER_EXTERNAL_LOOP | NPY_ITER_BUFFERED,
                          NPY_FORTRANORDER,
                          NPY_SAFE_CASTING,
                          op_flags,
                          dtypes);

  if (iter == NULL) {
    PyErr_SetString(PyExc_RuntimeError, "Couldn't set up iterator");
    Py_DECREF(data_array);
    return NULL;
  }

  // The iternext function gets stored in a local variable
  // so it can be called repeatedly in an efficient manner.
  iternext = NpyIter_GetIterNext(iter, NULL);
  if (iternext == NULL) {
    PyErr_SetString(PyExc_RuntimeError, "Couldn't set up iterator");
    NpyIter_Deallocate(iter);
    Py_DECREF(data_array);
    return NULL;
  }

  // The location of the data pointer which the iterator may update
  dataptr = NpyIter_GetDataPtrArray(iter);

  // The location of the stride which the iterator may update
  strideptr = NpyIter_GetInnerStrideArray(iter);

  // The location of the inner loop size which the iterator may update
  innersizeptr = NpyIter_GetInnerLoopSizePtr(iter);

  // Build the temporary buffer array
  buffer = (double *)malloc(n * sizeof(double));
  if (buffer == NULL) {
    PyErr_SetString(PyExc_RuntimeError, "Couldn't build buffer array");
    Py_DECREF(data_array);
    return NULL;
  }

  // Build output bounds arrays
  dims[0] = 2;
  dims[1] = m;
  bounds_obj = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
  if (bounds_obj == NULL) {
    PyErr_SetString(PyExc_RuntimeError, "Couldn't build bounds array");
    Py_DECREF(data_array);
    Py_XDECREF(bounds_obj);
    return NULL;
  }

  bounds_array = (PyArrayObject *)bounds_obj;
  bounds = (double *)PyArray_DATA(bounds_array);

  // This function is constructed to take a 2-d array of values and assumes
  // that each 1-d array when looping over the last dimension should be
  // treated separately. We use the Numpy iterator set up above which iterates
  // in C order, so we should get each 1-d array after one another. However
  // the inner loop of the iterator may not match the size of these 1-d arrays
  // so when we update the iterator may not be in sync with the m loop.

  stride = *strideptr;
  innersize = *innersizeptr;

  for (j = 0; j < m; j++) {

    iteration = 0;

    // We copy all finite values from array into the buffer

    count = 0;
    for (i = 0; i < n; i++) {
      innersize--;
      if(innersize == 0) {
        iternext(iter);
        innersize = *innersizeptr;
      }
      if (*(uint8_t *)dataptr[1] == 0) {
        buffer[count] = *(double *)dataptr[0];
        count += 1;
      }
      dataptr[0] += strideptr[0];
      dataptr[1] += strideptr[1];
    }

    // If end == 0, no values have been copied over (this can happen
    // for example if all the values are NaN). In this case, we just
    // proceed to the next array.
    if (count == 0)
      continue;

    compute_sigma_clipped_bounds(buffer, count, use_median, maxiters,
                                 sigma_lower, sigma_upper, &lower, &upper);

    bounds[j] = lower;
    bounds[j + m] = upper;

  }

  Py_DECREF(data_array);

  return bounds_obj;
}
