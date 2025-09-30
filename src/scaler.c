#define NPY_TARGET_VERSION NPY_2_0_API_VERSION
#define PY_SSIZE_T_CLEAN
#include "numpy/arrayobject.h"
#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <numpy/ufuncobject.h>
#include <stddef.h> /* for offsetof() */

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here. */
    double factor;
    PyObject *O_factor;
    PyObject *A_factor;
    vectorcallfunc vectorcall;
} ScalerObject;

static PyUFuncObject *multiply = NULL;

static PyObject *use_contiguous_loop(PyArrayObject *arr, double *factor)
{
    static PyUFuncGenericFunction loop = NULL; // [NPY_CDOUBLE] = {NULL};
    const char type_num = PyArray_DESCR(arr)->type_num;
    if (loop == NULL) {
        // Unfortunately, new-style loops are not exposed, so get legacy one.
        int nargs = multiply->nargs;
        const char *types = multiply->types;
        int i, j;
        for (i = 0; i < multiply->ntypes; ++i) {
            for (j = 0; j < nargs; ++j) {
                if (types[j] != type_num) {
                    break;
                }
            }
            if (j == nargs) {
                loop = multiply->functions[i];
                break;
            }
            types += nargs;
        }
        if (i == multiply->ntypes) {
            PyErr_SetString(PyExc_RuntimeError, "could not find loop");
            return NULL;
        }
    }
    // do it!
    PyObject *res =
        PyArray_EMPTY(PyArray_NDIM(arr), PyArray_DIMS(arr), type_num, PyArray_ISFORTRAN(arr));
    char *data[3] = {PyArray_DATA(arr), (char *)factor, PyArray_DATA((PyArrayObject *)res)};
    npy_intp n = PyArray_SIZE(arr);
    npy_intp strides[3] = {PyArray_ITEMSIZE(arr), 0, PyArray_ITEMSIZE(arr)};
    PyUFunc_clearfperr();
    loop(data, &n, strides, NULL);
    int fpe_errors = PyUFunc_getfperr();
    if (fpe_errors) {
        if (PyUFunc_GiveFloatingpointErrors("multiply", fpe_errors) < 0) {
            Py_DECREF(res);
            return NULL;
        }
    }
    return res;
}

static PyObject *Scaler_vectorcall(
    ScalerObject *self, PyObject *const *args, size_t len_args, PyObject *kwnames
)
{
    PyObject *const obj = args[0];
    if (PyVectorcall_NARGS(len_args) != 1) {
        PyErr_Format(
            PyExc_TypeError, "scaler() takes 1 argument, not %d", PyVectorcall_NARGS(len_args)
        );
        return NULL;
    }
    // fastest paths: special-case known objects.
    if (PyFloat_CheckExact(obj)) {
        return PyFloat_FromDouble(PyFloat_AS_DOUBLE(obj) * self->factor);
    }
    else if (PyArray_CheckExact(obj)) {
        PyArrayObject *const arr = (PyArrayObject *)obj;
        if (PyArray_DESCR(arr)->type_num == NPY_FLOAT64 && PyArray_ISONESEGMENT(arr)) {
            // speed up over ufunc by using strided loop directly.
            return use_contiguous_loop(arr, &self->factor);
        }
        if (self->A_factor == NULL) {
            npy_intp const dims[0];
            self->A_factor = PyArray_SimpleNewFromData(0, dims, NPY_FLOAT64, &self->factor);
            if (self->A_factor == NULL) {
                npy_intp const dims[0];
                self->A_factor = PyArray_SimpleNewFromData(0, dims, NPY_FLOAT64, &self->factor);
                if (self->A_factor == NULL) {
                    return NULL;
                }
            }
            return Py_TYPE(obj)->tp_as_number->nb_multiply(obj, self->A_factor);
        }
    }
    else if (PyLong_CheckExact(obj)) {
        double d = PyLong_AsDouble(obj);
        if (d == -1.0 && PyErr_Occurred()) {
            return NULL;
        }
        return PyFloat_FromDouble(d * self->factor);
    }
    // For cases without special treatment, we need the float object, so
    // construct it ahead of time.
    if (self->O_factor == NULL) {
        self->O_factor = PyFloat_FromDouble(self->factor);
        if (self->O_factor == NULL) {
            return NULL;
        }
    }
    // fast path, go directly for the slot, to avoid PyNumber call.
    if (Py_TYPE(obj)->tp_as_number != NULL) {
        binaryfunc slotv = Py_TYPE(obj)->tp_as_number->nb_multiply;
        if (slotv != NULL) {
            PyObject *res = slotv(obj, self->O_factor);
            if (res != Py_NotImplemented) {
                return res;
            }
            // Fall back to slow path, which will almost certainly raise.
            Py_DECREF(res);
        }
    }
    return PyNumber_Multiply(obj, self->O_factor);
}

static PyObject *Scaler_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    if (kwds != NULL || PyTuple_Size(args) != 1) {
        PyErr_SetString(PyExc_TypeError, "Scaler takes exactly 1 positional argument.");
        return NULL;
    }
    double factor = PyFloat_AsDouble(PyTuple_GET_ITEM(args, 0));
    if (factor == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    ScalerObject *self = (ScalerObject *)type->tp_alloc(type, 0);
    if (self == NULL) {
        return NULL;
    }
    self->factor = factor;
    self->O_factor = NULL;
    self->A_factor = NULL;
    self->vectorcall = (vectorcallfunc)&Scaler_vectorcall;
    return (PyObject *)self;
}

static void Scaler_dealloc(ScalerObject *self)
{
    Py_XDECREF(self->O_factor);
    Py_XDECREF(self->A_factor);
    Py_TYPE(self)->tp_free(self);
}

static PyMemberDef Scaler_members[] = {
    {"factor",
     Py_T_DOUBLE,
     offsetof(ScalerObject, factor),
     Py_READONLY,
     "factor with which input is multiplied"},
    {NULL},
};

static PyTypeObject ScalerType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0).tp_name = "Scaler",
    .tp_doc = PyDoc_STR("Multiplies input by the given scale factor"),
    .tp_basicsize = sizeof(ScalerObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_VECTORCALL,
    .tp_new = Scaler_new,
    .tp_dealloc = (destructor)Scaler_dealloc,
    .tp_vectorcall_offset = offsetof(ScalerObject, vectorcall),
    .tp_call = &PyVectorcall_Call,
    .tp_members = Scaler_members,
};

static int scaler_module_exec(PyObject *m)
{
    if (PyType_Ready(&ScalerType) < 0) {
        return -1;
    }

    if (PyModule_AddObjectRef(m, "Scaler", (PyObject *)&ScalerType) < 0) {
        return -1;
    }

    return 0;
}

static PyModuleDef_Slot scaler_module_slots[] = {
    {Py_mod_exec, scaler_module_exec},
    // Just use this while using static types
    {Py_mod_multiple_interpreters, Py_MOD_MULTIPLE_INTERPRETERS_NOT_SUPPORTED},
    {0, NULL}
};

static PyModuleDef scaler_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "scaler",
    .m_doc = "Example module that creates an extension type.",
    .m_size = 0,
    .m_slots = scaler_module_slots,
};

PyMODINIT_FUNC PyInit_scaler(void)
{
    if (PyUFunc_ImportUFuncAPI() < 0) {
        return NULL;
    }
    if (PyArray_ImportNumPyAPI() < 0) {
        return NULL;
    }
    PyObject *mod = PyImport_ImportModule("numpy");
    if (mod != NULL) {
        multiply = (PyUFuncObject *)PyObject_GetAttrString(mod, "multiply");
        Py_DECREF(mod);
    }
    if (multiply == NULL) {
        return NULL;
    }
    return PyModuleDef_Init(&scaler_module);
}
