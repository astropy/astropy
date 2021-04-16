#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include "numpy/ufuncobject.h"
#include "compute_bounds.h"

/* Define docstrings */
static char module_docstring[] = "Fast sigma clipping";
static char _sigma_clip_fast_docstring[] = "Compute sigma clipping";

/* Declare the C functions here. */
static void _sigma_clip_fast(
    char **args, npy_intp *dimensions, npy_intp* steps, void* data);

/* Define the methods that will be available on the module. */
static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

/* This is the function that is called on import. */

#define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#define MOD_DEF(ob, name, doc, methods)                                        \
  static struct PyModuleDef moduledef = {                                      \
      PyModuleDef_HEAD_INIT, name, doc, -1, methods,                           \
      NULL, NULL, NULL, NULL                                                   \
  };                                                                           \
  ob = PyModule_Create(&moduledef);

MOD_INIT(_fast_sigma_clip) {
    PyObject *m, *d;
    PyUFuncObject *ufunc;
    static char types[8] = {
        NPY_DOUBLE, /* data array */
        NPY_BOOL, /* mask array */
        NPY_BOOL, /* use median */
        NPY_INT, /* max iter */
        NPY_DOUBLE, /* sigma low */
        NPY_DOUBLE, /* sigma high */
        NPY_DOUBLE, /* output: lower bound */
        NPY_DOUBLE /* output: upper bound */
    };
    /* In principle, can have multiple functions for multiple input types */
    static PyUFuncGenericFunction funcs[1] = { &_sigma_clip_fast };
    static void *data[1] = {NULL};
    MOD_DEF(m, "_fast_sigma_clip", module_docstring, module_methods);
    if (m == NULL) {
        goto fail;
    }
    d = PyModule_GetDict(m); /* borrowed ref. */
    if (d == NULL) {
        goto fail;
    }
    import_array();
    import_umath();

    ufunc = (PyUFuncObject *)PyUFunc_FromFuncAndDataAndSignature(
        funcs, data, types, 1, 6, 2, PyUFunc_None, "_sigma_clip_fast",
        _sigma_clip_fast_docstring, 0, "(n),(n),(),(),(),()->(),()");
    if (ufunc == NULL) {
        goto fail;
    }
    PyDict_SetItemString(d, "_sigma_clip_fast", (PyObject *)ufunc);
    Py_DECREF(ufunc);
    return m;

  fail:
    Py_XDECREF(m);
    Py_XDECREF(d);
    return NULL;
}


static void _sigma_clip_fast(
    char **args, npy_intp *dimensions, npy_intp* steps, void* data)
{
    npy_intp i_o, i;
    int count;
    /* dimensions, pointers and step sizes for outer loop */
    npy_intp n_o = *dimensions++;
    char *array = *args++;
    npy_intp s_array = *steps++;
    char *mask = *args++;
    npy_intp s_mask = *steps++;
    char *use_median = *args++;
    npy_intp s_use_median = *steps++;
    char *max_iter = *args++;
    npy_intp s_max_iter = *steps++;
    char *sigma_low = *args++;
    npy_intp s_sigma_low = *steps++;
    char *sigma_high = *args++;
    npy_intp s_sigma_high = *steps++;
    char *bound_low = *args++;
    npy_intp s_bound_low = *steps++;
    char *bound_high = *args++;
    npy_intp s_bound_high = *steps++;
    /* dimension and step sizes for inner loop */
    npy_intp n_i = dimensions[0];
    char *in_array, *in_mask;
    npy_intp is_array = *steps++;
    npy_intp is_mask = *steps++;

    double *buffer = NULL;

    buffer = (double *)PyArray_malloc(n_i * sizeof(double));
    if (buffer == NULL) {
        PyErr_NoMemory();
        return;
    }
    for (i_o = 0; i_o < n_o;
         i_o++, array += s_array, mask += s_mask,
             use_median += s_use_median, max_iter += s_max_iter,
             sigma_low += s_sigma_low, sigma_high += s_sigma_high,
             bound_low += s_bound_low, bound_high += s_bound_high) {
        /* copy to buffer */
        in_array = array;
        in_mask = mask;
        count = 0;
        for (i = 0; i < n_i; i++, in_array += is_array, in_mask += is_mask) {
            if (*(uint8_t *)in_mask == 0) {
                buffer[count] = *(double *)in_array;
                count += 1;
            }
        }
        if (count > 0) {
            compute_sigma_clipped_bounds(
                buffer, count,
                (int)(*(npy_bool *)use_median), *(int *)max_iter,
                *(double *)sigma_low, *(double *)sigma_high,
                (double *)bound_low, (double *)bound_high);
        }
        else {
            *(double *)bound_low = NPY_NAN;
            *(double *)bound_high = NPY_NAN;
        }
    }
    PyArray_free((void *)buffer);
}
