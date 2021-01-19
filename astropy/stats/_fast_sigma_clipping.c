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

/*---------------------------------------------------------------------------

   Algorithm from N. Wirth's book, implementation by N. Devillard.
   This code in public domain.

   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :   use the median() macro defined below to get the median. 

                Reference:

                  Author: Wirth, Niklaus 
                   Title: Algorithms + data structures = programs 
               Publisher: Englewood Cliffs: Prentice-Hall, 1976 
    Physical description: 366 p. 
                  Series: Prentice-Hall Series in Automatic Computation 

 ---------------------------------------------------------------------------*/

#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

double kth_smallest(double a[], int n, int k)
{
    register int i,j,l,m ;
    register double x ;

    l=0 ; m=n-1 ;
    while (l<m) {
        x=a[k] ;
        i=l ;
        j=m ;
        do {
            while (a[i]<x) i++ ;
            while (x<a[j]) j-- ;
            if (i<=j) {
                ELEM_SWAP(a[i],a[j]) ;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    return a[k] ;
}

double wirth_median(double a[], int n) {
  if (n % 2 == 0) {
    return 0.5 * (kth_smallest(a,n,n/2) + kth_smallest(a,n,n/2 - 1));
  } else {
    return kth_smallest(a,n,(n-1)/2);
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
  int iteration, count;
  double lower, upper;
  int new_count;

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

  // This function is constructed to take a 2-d array of values and assumes
  // that each 1-d array when looping over the last dimension should be
  // treated separately.

  for (j = 0; j < m; j++) {

    iteration = 0;

    // We copy all finite values from array into the buffer

    count = 0;
    for (i = 0; i < n; i++) {
      if (data[i * m + j] == data[i * m + j]) {
        buffer[count] = data[i * m + j];
        count += 1;
      }
    }

    // If end == 0, no values have been copied over (this can happen
    // for example if all the values are NaN). In this case, we just
    // proceed to the next array.
    if (count == 0)
      continue;

    while (1) {

      // Calculate the mean and standard deviation of values so far.

      mean = 0;
      for (i = 0; i < count; i++) {
        mean += buffer[i];
      }
      mean /= count;

      std = 0;
      for (i = 0; i < count; i++) {
        std += pow(mean - buffer[i], 2);
      }
      std = sqrt(std / count);

      if (use_median) {

        median = wirth_median(buffer, count);

        lower = median - sigma_lower * std;
        upper = median + sigma_upper * std;

      } else {

        lower = mean - sigma_lower * std;
        upper = mean + sigma_upper * std;
      }

      // We now exclude values from the buffer using these
      // limits and shift values so that we end up with a
      // packed array of 'valid' values
      new_count = 0;
      for (i = 0; i < count; i++) {
        if (buffer[i] >= lower && buffer[i] <= upper) {
          buffer[new_count] = buffer[i];
          new_count += 1;
        }
      }

      if (new_count == count)
        break;

      iteration += 1;

      count = new_count;

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
