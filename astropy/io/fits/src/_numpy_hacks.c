/* Licensed under a 3-clause BSD style license - see LICENSE.rst

   This module is for functions that do tricky things with Numpy arrays and
   dtypes that are not normally supported in Numpy (but can work in limited
   cases relevant to FITS) or that otherwise require workarounds.

   Currently there is only one such function--realign_dtype--you can
   find its docstring in the source below.
*/

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>


#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif


const char realign_dtype_doc[] =
"Given a Numpy struct dtype object an a list of integer offsets, with one\n\
offset per field in the dtype, returns a new dtype where each field has the\n\
given offset.\n\n\
All offsets must be non-negative integers, but otherwise have no\n\
restrictions, and may overlap, per the usual rules for creating struct\n\
dtypes.  The new dtype will have an itemsize equal to the offset of the\n\
right-most field plus the width of that field.\n\n\
One restriction of this function is that it must not be used with object\n\
arrays--incorrect offsets may lead to invalid pointers in the arrays.\n\
However, this function is really only meant for use by astropy.io.fits\n\
and object arrays are not supported for FITS data anyhow.\n\n\
This function is used primarily to get around a shortcoming in Numpy that\n\
it is currently impossible to create dtypes with arbitrary offsets, *and*\n\
that have zero-width fields.  Both of these features are needed for full\n\
FITS support.  However, this will be fixed in a future version of Numpy\n\
at which point use of this hack can be deprecated.  See\n\
https://github.com/numpy/numpy/pull/6430";


PyObject* numpy_hacks_realign_dtype(PyObject *not_used, PyObject *args) {
    PyObject *offsets, *field_tuple, *new_field_tuple, *field_title, *tmp;
    PyArray_Descr *dtype, *field_dtype;
    int n, idx, jdx, new_elsize;
    long offset, new_offset;
    long maxoffset = 0;

    if (!PyArg_ParseTuple(args, "O!O", &PyArrayDescr_Type, &dtype, &offsets)) {
        return NULL;
    }

    if (!PyDataType_HASFIELDS(dtype)) {
        PyErr_SetString(PyExc_ValueError,
                        "given dtype does not have struct fields");
        return NULL;
    }

    n = PyObject_Length(dtype->fields);
    new_elsize = dtype->elsize;  /* Just by default */

    if (!PySequence_Check(offsets) ||
            PyObject_Length(offsets) != n) {
        PyErr_SetString(PyExc_ValueError,
                        "offsets must be a sequence of length equal "
                        "to the number of fields in the dtype");
        return NULL;
    }

    /* Replace the original dtype argument with a copy to return */
    /* Also copy the fields dict (otherwise both dtypes reference the same
       dict */
    dtype = PyArray_DescrNew(dtype);
    dtype->fields = PyDict_Copy(dtype->fields);

    for (idx = 0; idx < n; idx++) {
        field_tuple = PyObject_GetItem(
                dtype->fields, PyTuple_GET_ITEM(dtype->names, idx));
        if (PyTuple_GET_SIZE(field_tuple) == 2) {
            if (!PyArg_ParseTuple(field_tuple, "O!i", &PyArrayDescr_Type,
                                  &field_dtype, &offset)) {
                return NULL;
            }
        } else {
            /* Only other possibility is len(field_tuple) == 3 */
            if (!PyArg_ParseTuple(field_tuple, "O!iO", &PyArrayDescr_Type,
                                  &field_dtype, &offset, &field_title)) {
                return NULL;
            }
        }

        tmp = PySequence_GetItem(offsets, idx);
        if (!tmp) {
            return NULL;
        }
        new_offset = PyLong_AsLong(tmp);

        if (PyErr_Occurred()) {
            Py_DECREF(tmp);
            return NULL;
        }

        if (offset != new_offset) {
            if (new_offset > maxoffset) {
                maxoffset = new_offset;
                new_elsize = new_offset + field_dtype->elsize;
            }

            new_field_tuple = PyTuple_New(PyTuple_GET_SIZE(field_tuple));
            PyTuple_SET_ITEM(new_field_tuple, 1, tmp);

            for (jdx = 0; jdx < PyTuple_GET_SIZE(field_tuple); jdx++) {
                if (jdx == 1) {
                    continue;
                }
                tmp = PyTuple_GET_ITEM(field_tuple, jdx);
                Py_INCREF(tmp);
                PyTuple_SET_ITEM(new_field_tuple, jdx, tmp);
            }

            PyObject_SetItem(dtype->fields,
                             PyTuple_GET_ITEM(dtype->names, idx),
                             new_field_tuple);
        } else {
            Py_DECREF(tmp);
        }
    }

    dtype->elsize = new_elsize;

    return (PyObject *)dtype;
}


static PyMethodDef _numpy_hacks_methods[] =
{
   {"realign_dtype", numpy_hacks_realign_dtype, METH_VARARGS,
    realign_dtype_doc},
   {NULL, NULL}
};


const char _numpy_hacks_doc[] =
"This module is for functions that do tricky things with Numpy arrays and\n\
 dtypes that are not normally supported in Numpy (but can work in limited\n\
 cases relevant to FITS) or that otherwise require workarounds.";


#ifdef IS_PY3K
static struct PyModuleDef _numpy_hacksmodule = {
    PyModuleDef_HEAD_INIT,
    "_numpy_hacks",
    _numpy_hacks_doc,
    -1, /* No global state */
    _numpy_hacks_methods
};

PyObject *
PyInit__numpy_hacks(void)
{
    PyObject* module = PyModule_Create(&_numpy_hacksmodule);

    /* Needed to use Numpy routines */
    /* Note -- import_array() is a macro that behaves differently in Python2.x
     * vs. Python 3. See the discussion at:
     * https://groups.google.com/d/topic/astropy-dev/6_AesAsCauM/discussion
     */
    import_array();
    return module;
}
#else
PyMODINIT_FUNC init_numpy_hacks(void)
{
   PyObject* module = Py_InitModule4("_numpy_hacks", _numpy_hacks_methods,
                                     _numpy_hacks_doc,
                                     NULL, PYTHON_API_VERSION);
   import_array();
}
#endif
