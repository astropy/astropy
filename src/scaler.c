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
    double factor_imag; // always set to 0, to have complex factor.
    float factor_f;
    float factor_f_imag; // always set to 0.
    PyObject *O_factor;
    PyObject *A_factor;
    PyObject *A_factor_f;
    vectorcallfunc vectorcall;
} ScalerObject;

static PyUFuncGenericFunction loops[6] = {NULL};

static inline PyObject *use_contiguous_loop(PyArrayObject *arr, char *factor_ptr)
{
    const char type_num = PyArray_TYPE(arr);
    PyUFuncGenericFunction loop = loops[(int)(type_num - NPY_FLOAT)];
    PyObject *res =
        PyArray_EMPTY(PyArray_NDIM(arr), PyArray_DIMS(arr), type_num, PyArray_ISFORTRAN(arr));
    char *data[3] = {PyArray_DATA(arr), factor_ptr, PyArray_DATA((PyArrayObject *)res)};
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
        char type_num = PyArray_TYPE(arr);
        npy_bool needs_float = type_num == NPY_FLOAT || type_num == NPY_CFLOAT;
        char *f_ptr = needs_float ? (char *)&self->factor_f : (char *)&self->factor;
        // Pass contiguous float or complex arrays directly to multiply loop,
        // bypassing ufunc setup.
        if (PyArray_ISONESEGMENT(arr) &&
            (type_num == NPY_DOUBLE || type_num == NPY_CDOUBLE || needs_float)) {
            return use_contiguous_loop(arr, f_ptr);
        }
        // If not, convert factor to array here, since that makes ufunc
        // call substantially faster.
        PyObject **A_factor = needs_float ? &self->A_factor_f : &self->A_factor;
        if (*A_factor == NULL) {
            const npy_intp dims[1] = {0};
            *A_factor =
                PyArray_SimpleNewFromData(0, dims, needs_float ? NPY_FLOAT : NPY_DOUBLE, f_ptr);
        }
        return Py_TYPE(obj)->tp_as_number->nb_multiply(obj, *A_factor);
    }
    else if (PyLong_CheckExact(obj)) {
        double d = PyLong_AsDouble(obj);
        if (d == -1.0 && PyErr_Occurred()) {
            return NULL;
        }
        return PyFloat_FromDouble(d * self->factor);
    }
    // For cases without special treatment, we need the float object;
    // construct it if not done already.
    if (self->O_factor == NULL) {
        self->O_factor = PyFloat_FromDouble(self->factor);
        if (self->O_factor == NULL) {
            return NULL;
        }
    }
    // If a known type of object, like a subclass, or something that has
    // a dtype or supports the Array API, just run multiply.
    if (PyFloat_Check(obj) || PyComplex_Check(obj) || PyLong_Check(obj) || PyArray_Check(obj) ||
        PyObject_HasAttrString(obj, "dtype") ||
        PyObject_HasAttrString(obj, "__array_namespace__")) {
        // Generally, should be possible to go directly for the slot.
        if (Py_TYPE(obj)->tp_as_number != NULL) {
            binaryfunc slotv = Py_TYPE(obj)->tp_as_number->nb_multiply;
            if (slotv != NULL) {
                PyObject *res = slotv(obj, self->O_factor);
                if (res != Py_NotImplemented) {
                    return res;
                }
                // Fall back to slow path, which will likely raise,
                // thus giving a reasonable error message.
                Py_DECREF(res);
            }
        }
        return PyNumber_Multiply(obj, self->O_factor);
    }
    // If not a known type, try converting to an array.
    PyObject *arr = PyArray_FromAny(obj, NULL, 0, 0, 0, NULL);
    if (arr == NULL) {
        return NULL;
    }
    PyObject *res = NULL;
    if (PyArray_ISNUMBER((PyArrayObject *)arr)) {
        // Re-enter to use the array path.
        res = Scaler_vectorcall(self, &arr, 1, NULL);
    }
    else {
        // If not numeric, same error message as in _condition_arg.
        PyErr_SetString(
            PyExc_ValueError,
            "Value not scalar compatible or convertible to "
            "an int, float, or complex array"
        );
    }
    Py_DECREF(arr);
    return res;
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
    self->factor_imag = 0; // allows using &self->factor for complex.
    self->factor_f = (float)factor;
    self->factor_f_imag = 0;
    self->O_factor = NULL; // following initialized only if needed in call.
    self->A_factor = NULL;
    self->A_factor_f = NULL;
    self->vectorcall = (vectorcallfunc)&Scaler_vectorcall;
    return (PyObject *)self;
}

static void Scaler_dealloc(ScalerObject *self)
{
    Py_XDECREF(self->O_factor);
    Py_XDECREF(self->A_factor);
    Py_XDECREF(self->A_factor_f);
    Py_TYPE(self)->tp_free(self);
}

static PyMemberDef Scaler_members[] = {
    {"factor",
     Py_T_DOUBLE,
     offsetof(ScalerObject, factor),
     Py_READONLY,
     "Factor with which input is multiplied."},
    {NULL},
};

static PyTypeObject ScalerType = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0).tp_name = "Scaler",
    .tp_doc = PyDoc_STR("Multiplies input by the given scale factor."),
    .tp_basicsize = sizeof(ScalerObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_VECTORCALL,
    .tp_new = Scaler_new,
    .tp_dealloc = (destructor)Scaler_dealloc,
    .tp_vectorcall_offset = offsetof(ScalerObject, vectorcall),
    .tp_call = &PyVectorcall_Call,
    .tp_members = Scaler_members,
};

static int get_multiply_loops(PyUFuncGenericFunction loops[6])
{
    // Get and cache multiplication loops for use with plain ndarray.
    PyObject *mod = PyImport_ImportModule("numpy");
    if (mod == NULL) {
        return -1;
    }
    PyUFuncObject *multiply = (PyUFuncObject *)PyObject_GetAttrString(mod, "multiply");
    Py_DECREF(mod);
    if (multiply == NULL) {
        return -1;
    }
    // Unfortunately, new-style loops are not exposed, so get legacy ones.
    for (int k = 0; k < 5; k++) {
        char type = (char)(k + (int)NPY_FLOAT);
        if (type == NPY_LONGDOUBLE) {
            continue;
        }
        const char *types = multiply->types;
        for (int i = 0; i < multiply->ntypes; i++) {
            if (types[0] == type && types[1] == type && types[2] == type) {
                loops[k] = multiply->functions[i];
                break;
            }
            types += 3;
        }
        if (loops[k] == NULL) {
            PyErr_SetString(PyExc_RuntimeError, "could not find loop");
            return -1;
        }
    }
    // Keep reference to multiply since we are using its loops
    // (not that it will ever disappear...).
    return 0;
}

static int scaler_module_exec(PyObject *m)
{
    if (PyType_Ready(&ScalerType) < 0) {
        return -1;
    }
    if (PyModule_AddObjectRef(m, "Scaler", (PyObject *)&ScalerType) < 0) {
        return -1;
    }
    if (get_multiply_loops(loops) < 0) {
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
    .m_doc = "Provides the inspectable Scaler class.",
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
    return PyModuleDef_Init(&scaler_module);
}
