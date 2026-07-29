#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h> /* for offsetof() */

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here. */
    PyObject *factor;
} ScalerObject;

static PyObject *Scaler_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    ScalerObject *self;
    self = (ScalerObject *)type->tp_alloc(type, 0);
    if (self == NULL) {
        return NULL;
    }
    static char *kwlist[] = {"factor", NULL};
    double factor;
    int res = PyArg_ParseTupleAndKeywords(args, kwds, "d", kwlist, &factor);
    if (!res) {
        Py_DECREF(self);
        return NULL;
    }
    self->factor = PyFloat_FromDouble(factor);
    if (self->factor == NULL) {
        Py_DECREF(self);
        return NULL;
    }
    return (PyObject *)self;
}

static void Scaler_dealloc(ScalerObject *self)
{
    Py_XDECREF(self->factor);
    Py_TYPE(self)->tp_free(self);
}

PyObject *Scaler_call(ScalerObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *obj, *mul, *res;
    static char *kwlist[] = {"x", NULL};
    if (PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &obj) < 0) {
        return NULL;
    }
    mul = PyObject_GetAttrString(obj, "__mul__");
    if (mul == NULL) {
        return NULL;
    }
    res = PyObject_CallOneArg(mul, self->factor);
    Py_DECREF(mul);
    if (res == Py_NotImplemented) {
        // Fairly unlikely to happen, but int cannot multiply with float.
        Py_DECREF(res);
        mul = PyObject_GetAttrString(self->factor, "__rmul__");
        res = PyObject_CallOneArg(mul, obj);
        Py_DECREF(mul);
        if (res == Py_NotImplemented) { // Not sure I can hit this.
            Py_DECREF(res);
            PyErr_Format(PyExc_TypeError, "cannot scale type '%T'", obj);
            res = NULL;
        }
    }
    return res;
}

static PyMemberDef Scaler_members[] = {
    {"factor",
     Py_T_OBJECT_EX,
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
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = Scaler_new,
    .tp_dealloc = (destructor)Scaler_dealloc,
    .tp_call = (ternaryfunc)Scaler_call,
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
    return PyModuleDef_Init(&scaler_module);
}
