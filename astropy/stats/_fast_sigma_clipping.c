#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

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

MOD_INIT(_fast_sigma_clipping) {
  PyObject *m;
  MOD_DEF(m, "_fast_sigma_clipping", module_docstring, module_methods);
  if (m == NULL)
    return MOD_ERROR_VAL;
  import_array();
  return MOD_SUCCESS_VAL(m);
}

static int sort_cmp(const void *a, const void *b) {
  if (*(double *)a > *(double *)b)
    return 1;
  else if (*(double *)a < *(double *)b)
    return -1;
  else
    return 0;
}

size_t locate(double *array, double target, size_t start, size_t end) {
  size_t mid;
  while (1) {
    if (end == start + 1)
      return start;
    mid = (start + end) / 2;
    if (array[mid] < target) {
      start = mid;
    } else {
      end = mid;
    }
  }
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
  double sigma_lower, sigma_upper, mean, median, std;
  int iteration, start, end, count, median_index, start_prev, end_prev;
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

  // Build the temporary array
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

  // This function is constructed to take a 2-d array of values and assumes
  // that each 1-d array when looping over the last dimension should be
  // treated separately.

  for (j = 0; j < m; j++) {

    iteration = 0;

    // We copy all finite values from array into the buffer buffer,
    // sort it, then keep track of the range of values that are
    // not rejected by the sigma clipping. Having the values in a
    // sorted array means that the sigma clipping is removing values
    // from either or both ends of the array - [start:end] gives the
    // range of values that have not been rejected. This has the
    // advantage that we don't then need to update a data or mask
    // array in each iteration - just the start and end indices.

    start = 0;
    end = 0;

    for (i = 0; i < n; i++) {
      if (data[i * m + j] == data[i * m + j]) {
        buffer[end] = data[i * m + j];
        end += 1;
      }
    }

    // If end == 0, no values have been copied over (this can happen
    // for example if all the values are NaN). In this case, we just
    // proceed to the next array.
    if (end == 0)
      continue;

    // We now sort the values in the array up to the end index.
    qsort(&buffer[0], end, 8, &sort_cmp);

    while (1) {

      count = end - start;

      // Calculate the mean and standard deviation of values so far.

      mean = 0;
      for (i = start; i < end; i++) {
        mean += buffer[i];
      }
      mean /= count;

      std = 0;
      for (i = start; i < end; i++) {
        std += pow(mean - buffer[i], 2);
      }
      std = sqrt(std / count);

      // If needed, we compute the median
      if (use_median) {

        if (count % 2 == 0) {
          median_index = start + count / 2 - 1;
          median = 0.5 * (buffer[median_index] + buffer[median_index + 1]);
        } else {
          median_index = start + (count - 1) / 2;
          median = buffer[median_index];
        }

        lower = median - sigma_lower * std;
        upper = median + sigma_upper * std;

      } else {

        lower = mean - sigma_lower * std;
        upper = mean + sigma_upper * std;
      }

      // If all array values in the [start:end] range are still inside
      // (lower, upper) then the process has converged and we can exit
      // the loop over iterations.
      if (buffer[start] > lower && buffer[end - 1] < upper)
        break;

      // We need to keep track of the previous start/end values as
      // we need the original values for both locate calls.
      start_prev = start;
      end_prev = end;

      // Update the start/end values based on the new lower/upper values
      if (buffer[start] < lower)
        start = locate(buffer, lower, start_prev, end_prev) + 1;
      if (buffer[end - 1] > upper)
        end = locate(buffer, upper, start_prev, end_prev) + 1;

      iteration += 1;

      if (maxiters != -1 && iteration >= maxiters)
        break;
    }

    // Populate the final (unsorted) mask
    for (i = 0; i < n; i++) {
      if (data[i * m + j] < lower || data[i * m + j] > upper) {
        mask[i * m + j] = 1;
      } else {
        mask[i * m + j] = 0;
      }
    }
  }

  return mask_obj;
}