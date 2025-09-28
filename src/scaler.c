#define PY_SSIZE_T_CLEAN
#include "numpy/arrayobject.h"
#include <Python.h>
#include <stddef.h> /* for offsetof() */

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here. */
    double factor;
    PyObject *O_factor;
    PyObject *A_factor;
    vectorcallfunc vectorcall;
} ScalerObject;

PyObject *Scaler_vectorcall(
    ScalerObject *self, PyObject *const *args, size_t len_args, PyObject *kwnames
)
{
    PyObject *obj = args[0];
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
        return Py_TYPE(obj)->tp_as_number->nb_multiply(obj, self->A_factor);
    }
    else if (PyLong_CheckExact(obj)) {
        double d = PyLong_AsDouble(obj);
        if (d == -1.0 && PyErr_Occurred()) {
            return NULL;
        }
        return PyFloat_FromDouble(d * self->factor);
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
    ScalerObject *self;
    self = (ScalerObject *)type->tp_alloc(type, 0);
    if (self == NULL) {
        return NULL;
    }
    static char *kwlist[] = {"factor", NULL};
    int res = PyArg_ParseTupleAndKeywords(args, kwds, "d", kwlist, &self->factor);
    if (!res) {
        Py_DECREF(self);
        return NULL;
    }
    // For cases without special treatment, we need the float object, so
    // construct it ahead of time.
    self->O_factor = PyFloat_FromDouble(self->factor);
    if (self->O_factor == NULL) {
        Py_DECREF(self);
        return NULL;
    }
    PyArray_Descr *dtype = PyArray_DescrFromType(NPY_FLOAT64); // cannot fail;
    self->A_factor = PyArray_FromAny(self->O_factor, dtype, 0, 0, 0, NULL);
    if (self->A_factor == NULL) {
        Py_DECREF(self);
        Py_DECREF(self->O_factor);
        return NULL;
    }
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
     Py_T_OBJECT_EX,
     offsetof(ScalerObject, O_factor),
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
    import_array();
    return PyModuleDef_Init(&scaler_module);
}
