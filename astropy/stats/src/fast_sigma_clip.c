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
  int i, j, istride, jstride, index;
  double *buffer;
  PyObject *data_obj, *mask_obj;
  PyArrayObject *data_array, *mask_array;
  double *data;
  uint8_t *mask;
  int use_median, maxiters;
  double sigma_lower, sigma_upper;
  int iteration, count;
  double lower, upper;

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

  // Build the temporary array and mask
  buffer = (double *)malloc(n * sizeof(double));
  if (buffer == NULL) {
    PyErr_SetString(PyExc_RuntimeError, "Couldn't build buffer array");
    Py_DECREF(data_array);
    Py_DECREF(mask_array);
    return NULL;
  }

  // At this point we should ideally use the iterator API in Numpy,
  // but for simplicity for now we assume the data is contiguous in
  // memory and will deal with the correct iteration later

  data = (double *)PyArray_DATA(data_array);
  mask = (uint8_t *)PyArray_DATA(mask_array);

  istride = PyArray_STRIDE(data_array, 0) / 8;
  jstride = PyArray_STRIDE(data_array, 1) / 8;

  // This function is constructed to take a 2-d array of values and assumes
  // that each 1-d array when looping over the last dimension should be
  // treated separately.

  for (j = 0; j < m; j++) {

    iteration = 0;

    // We copy all finite values from array into the buffer

    count = 0;
    index = j * jstride;
    for (i = 0; i < n; i++) {
      if (mask[index] == 0) {
        buffer[count] = data[index];
        count += 1;
      }
      index += istride;
    }

    // If end == 0, no values have been copied over (this can happen
    // for example if all the values are NaN). In this case, we just
    // proceed to the next array.
    if (count == 0)
      continue;

    compute_sigma_clipped_bounds(buffer, count, use_median, maxiters,
                                 sigma_lower, sigma_upper, &lower, &upper);

    // Populate the final (unsorted) mask
    index = j * jstride;
    for (i = 0; i < n; i++) {
      if (data[index] < lower || data[index] > upper) {
        mask[index] = 1;
      } else {
        mask[index] = 0;
      }
      index += istride;
   }
  }

  return mask_obj;
}
