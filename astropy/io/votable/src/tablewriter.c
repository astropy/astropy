/******************************************************************************
 * C extension code for vo.table.
 *
 * Everything in this file has an alternate Python implementation and
 * is included for performance reasons only.
 *
 * This contains a write_tabledata function to quickly write out a Numpy array
 * in TABLEDATA format.
 *
 ******************************************************************************/

#include <Python.h>

/******************************************************************************
 * Convenience macros and functions
 ******************************************************************************/
#undef  CLAMP
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

static Py_ssize_t
next_power_of_2(Py_ssize_t n)
{
    /* Calculate the next-highest power of two */
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;

    return n == 0 ? 2 : n;
}

/******************************************************************************
 * Python version compatibility macros
 ******************************************************************************/
#if PY_MAJOR_VERSION >= 3
#  define IS_PY3K
#endif

#ifndef Py_TYPE
#  define Py_TYPE(o) ((o)->ob_type)
#endif

/******************************************************************************
 * Write TABLEDATA
 ******************************************************************************/

#define CHAR Py_UNICODE

/*
 * Reallocate the write buffer to the requested size
 */
static int
_buffer_realloc(
        CHAR** buffer, Py_ssize_t* buffer_size, CHAR** x, Py_ssize_t req_size)
{
    Py_ssize_t  n       = req_size;
    CHAR *      new_mem = NULL;

    if (req_size < *buffer_size) {
        return 0;
    }

    /* Calculate the next-highest power of two */
    n = next_power_of_2(n);

    if (n < req_size) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory for XML text.");
        return -1;
    }

    new_mem = realloc((void *)*buffer, n * sizeof(CHAR));
    if (new_mem == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory for XML text.");
        return -1;
    }

    *x = (CHAR *)new_mem + (*x - *buffer);
    *buffer = new_mem;
    *buffer_size = n;

    return 0;
}

/*
 * Write *indent* spaces to the buffer
 */
static int
_write_indent(CHAR** buffer, Py_ssize_t* buffer_size,
              CHAR** x, Py_ssize_t indent)
{
    if (_buffer_realloc(buffer, buffer_size, x,
                        (*x - *buffer + indent))) {
        return 1;
    }

    for (; indent; --indent) {
        *(*x)++ = ' ';
    }

    return 0;
}

/*
 * Write a string into a buffer.
 */
static int
_write_string(CHAR** buffer, Py_ssize_t* buffer_size,
              CHAR** x, const CHAR* src, const Py_ssize_t len) {
    if (_buffer_realloc(buffer, buffer_size, x,
                        (*x - *buffer + len))) {
        return 1;
    }

    while (*src != (CHAR)0) {
        *(*x)++ = *src++;
    }

    return 0;
}

/*
 * Write an 8-bit ascii-encoded C string to a Unicode string.
 */
static int
_write_cstring(CHAR** buffer, Py_ssize_t* buffer_size,
               CHAR** x, const char* src, const Py_ssize_t len) {
    if (_buffer_realloc(buffer, buffer_size, x,
                        (*x - *buffer + len))) {
        return 1;
    }

    while (*src != (char)0) {
        *(*x)++ = *src++;
    }

    return 0;
}

/*
 * Write a TABLEDATA element tree to the given write method.
 *
 * The Python arguments are:
 *
 * *write_method* (callable): A Python callable that takes a unicode
 *    string and writes it to a file or buffer.
 *
 * *array* (numpy structured array): A Numpy record array containing
 *    the data
 *
 * *mask* (numpy array): A Numpy array which is True everywhere a
 *    value is missing.  Must have the same shape as *array*.
 *
 * *converters* (list of callables): A sequence of methods which
 *    convert from the native data types in the columns of *array* to
 *    a unicode string in VOTABLE XML format.  Must have the same
 *    length as the number of columns in *array*.
 *
 * *write_null_values* (boolean): When True, write null values in
 *    their entirety in the table.  When False, just write empty <TD/>
 *    elements when the data is null or missing.
 *
 * *indent* (integer): The number of spaces to indent the table.
 *
 * *buf_size* (integer): The size of the write buffer.
 *
 * Returns None.
 */
static PyObject*
write_tabledata(PyObject* self, PyObject *args, PyObject *kwds)
{
    /* Inputs */
    PyObject* write_method = NULL;
    PyObject* array = NULL;
    PyObject* mask = NULL;
    PyObject* converters = NULL;
    PyObject* py_supports_empty_values = NULL;
    Py_ssize_t indent = 0;
    Py_ssize_t buf_size = (Py_ssize_t)1 << 8;

    /* Output buffer */
    CHAR* buf = NULL;
    CHAR* x;

    Py_ssize_t nrows = 0;
    Py_ssize_t ncols = 0;
    Py_ssize_t i, j;
    int write_full;
    int all;
    PyObject* numpy_module = NULL;
    PyObject* numpy_all_method = NULL;
    PyObject* array_row = NULL;
    PyObject* mask_row = NULL;
    PyObject* array_val = NULL;
    PyObject* mask_val = NULL;
    PyObject* converter = NULL;
    PyObject* all_masked_obj = NULL;
    PyObject* str_val = NULL;
    PyObject* tmp = NULL;
    CHAR* str_tmp = NULL;
    Py_ssize_t str_len = 0;
    int* supports_empty_values = NULL;
    PyObject* result = 0;

    if (!PyArg_ParseTuple(args, "OOOOOnn:write_tabledata",
                          &write_method, &array, &mask, &converters,
                          &py_supports_empty_values, &indent, &buf_size)) {
        goto exit;
    }

    if (!PyCallable_Check(write_method)) goto exit;
    if (!PySequence_Check(array)) goto exit;
    if (!PySequence_Check(mask)) goto exit;
    if (!PyList_Check(converters)) goto exit;
    if (!PyList_Check(py_supports_empty_values)) goto exit;
    indent = CLAMP(indent, (Py_ssize_t)0, (Py_ssize_t)80);
    buf_size = CLAMP(buf_size, (Py_ssize_t)1 << 8, (Py_ssize_t)1 << 24);

    if ((numpy_module = PyImport_ImportModule("numpy")) == NULL) goto exit;
    if ((numpy_all_method = PyObject_GetAttrString(numpy_module, "all"))
        == NULL) goto exit;

    if ((nrows = PySequence_Size(array)) == -1) goto exit;
    if ((ncols = PyList_Size(converters)) == -1) goto exit;
    if (PyList_Size(py_supports_empty_values) != ncols) goto exit;

    supports_empty_values = PyMem_Malloc(sizeof(int) * ncols);
    if (!supports_empty_values) goto exit;
    for (i = 0; i < ncols; ++i) {
        supports_empty_values[i] = PyObject_IsTrue(
                PyList_GET_ITEM(py_supports_empty_values, i));
    }

    if ((buf = PyMem_Malloc((size_t)buf_size * sizeof(CHAR))) == NULL) goto exit;

    for (i = 0; i < nrows; ++i) {
        if ((array_row = PySequence_GetItem(array, i)) == NULL) goto exit;
        if ((mask_row = PySequence_GetItem(mask, i)) == NULL) goto exit;

        x = buf;
        if (_write_indent(&buf, &buf_size, &x, indent)) goto exit;
        if (_write_cstring(&buf, &buf_size, &x, " <TR>\n", 6)) goto exit;

        for (j = 0; j < ncols; ++j) {
            if ((converter = PyList_GET_ITEM(converters, j)) == NULL) goto exit;
            if ((array_val = PySequence_GetItem(array_row, j)) == NULL) goto exit;
            if ((mask_val = PySequence_GetItem(mask_row, j)) == NULL) goto exit;

            write_full = 1;
            if (mask_val == Py_True) {
                write_full = 0;
            } else if (mask_val == Py_False) {
                // pass
            } else if (supports_empty_values[j]) {
                if ((all_masked_obj =
                     PyObject_CallFunctionObjArgs(numpy_all_method, mask_val, NULL))
                    == NULL) goto exit;
                if ((all = PyObject_IsTrue(all_masked_obj)) == -1) {
                    Py_DECREF(all_masked_obj);
                    goto exit;
                }
                Py_DECREF(all_masked_obj);

                write_full = !all;
            }

            if (write_full) {
                if (_write_indent(&buf, &buf_size, &x, indent)) goto exit;

                if ((str_val =
                     PyObject_CallFunctionObjArgs(converter, array_val, mask_val, NULL))
                    == NULL) goto exit;
                if (PyBytes_Check(str_val)) {
                    tmp = PyUnicode_FromEncodedObject(str_val, "utf-8", "ignore");
                    Py_DECREF(str_val);
                    str_val = tmp;
                }
                if ((str_tmp = PyUnicode_AsUnicode(str_val)) == NULL) {
                    Py_DECREF(str_val);
                    goto exit;
                }

                str_len = PyUnicode_GetSize(str_val);
                if (str_len) {
                    if (_write_cstring(&buf, &buf_size, &x, "  <TD>", 6) ||
                        _write_string(&buf, &buf_size, &x, str_tmp, str_len) ||
                        _write_cstring(&buf, &buf_size, &x, "</TD>\n", 6)) {
                        Py_DECREF(str_val);
                        goto exit;
                    }
                } else {
                    if (_write_cstring(&buf, &buf_size, &x, "  <TD/>\n", 8)) {
                        Py_DECREF(str_val);
                        goto exit;
                    }
                }

                Py_DECREF(str_val);
            } else {
                if (_write_indent(&buf, &buf_size, &x, indent)) goto exit;
                if (_write_cstring(&buf, &buf_size, &x, "  <TD/>\n", 8)) goto exit;
            }

            Py_DECREF(array_val); array_val = NULL;
            Py_DECREF(mask_val);  mask_val = NULL;
        }

        Py_DECREF(array_row); array_row = NULL;
        Py_DECREF(mask_row);  mask_row = NULL;

        if (_write_indent(&buf, &buf_size, &x, indent)) goto exit;
        if (_write_cstring(&buf, &buf_size, &x, " </TR>\n", 7)) goto exit;

        /* NULL-terminate the string */
        *x = (CHAR)0;
        if ((tmp = PyObject_CallFunction(write_method, "u#", buf, x - buf))
            == NULL) goto exit;
        Py_DECREF(tmp);
    }

    Py_INCREF(Py_None);
    result = Py_None;

 exit:
    Py_XDECREF(numpy_module);
    Py_XDECREF(numpy_all_method);

    Py_XDECREF(array_row);
    Py_XDECREF(mask_row);
    Py_XDECREF(array_val);
    Py_XDECREF(mask_val);

    PyMem_Free(buf);
    PyMem_Free(supports_empty_values);

    return result;
}

/******************************************************************************
 * Module setup
 ******************************************************************************/

static PyMethodDef module_methods[] =
{
    {"write_tabledata", (PyCFunction)write_tabledata, METH_VARARGS,
     "Fast C method to write tabledata"},
    {NULL}  /* Sentinel */
};

struct module_state {
    void* none;
};

#ifdef IS_PY3K
static int module_traverse(PyObject* m, visitproc visit, void* arg)
{
    return 0;
}

static int module_clear(PyObject* m)
{
    return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "tablewriter",
    "Fast way to write VOTABLE TABLEDATA",
    sizeof(struct module_state),
    module_methods,
    NULL,
    module_traverse,
    module_clear,
    NULL
};

#  define INITERROR return NULL

PyMODINIT_FUNC
PyInit_tablewriter(void)
#else /* Not PY3K */
#  define INITERROR return

#  ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#    define PyMODINIT_FUNC void
#  endif

PyMODINIT_FUNC
inittablewriter(void)
#endif
{
    PyObject* m;

#ifdef IS_PY3K
    m = PyModule_Create(&moduledef);
#else
    m = Py_InitModule3("tablewriter", module_methods,
                       "Fast way to write VOTABLE TABLEDATA");
#endif

    if (m == NULL)
        INITERROR;

#ifdef IS_PY3K
    return m;
#endif
}
