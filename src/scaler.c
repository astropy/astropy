#define NPY_TARGET_VERSION NPY_2_0_API_VERSION // For PyUFunc_GiveFloatingpointErrors
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#if Py_Version < 0x03120000
#include <structmember.h> // for PyMemberDef, not in Python.h for python <= 3.11
#endif
#include <numpy/arrayobject.h>
#include <numpy/arrayscalars.h>
#include <numpy/ndarrayobject.h>
#include <numpy/ufuncobject.h>
#include <stddef.h> // for offsetof()

PyDoc_STRVAR(scaler_module_doc, "Compiled module providing the inspectable Scaler class.");
PyDoc_STRVAR(
    Scaler_doc,
    "Multiplies input by the given scale factor.\n\n"
    "Parameters\n"
    "----------\n"
    "factor : float, positional only\n"
    "    The scale factor with which, when called, an argument will be multiplied.\n"
);


typedef struct {
    PyObject_HEAD
    vectorcallfunc vectorcall;
    double factor;
    float factor_f;
    PyObject *O_factor;
    PyObject *A_factor;
    PyObject *A_factor_f;
} ScalerObject;

// Double and float multiply loops, filled in on initialization.
static PyUFuncGenericFunction loops[2] = {NULL, NULL};
// Unity scaler, which is kept unique, filled in on initialization.
static PyObject *unity_scaler = NULL;

// Multiply a contiguous array directly with the factor, using np.multiply's
// loop function.
static inline PyObject *use_contiguous_loop(PyArrayObject *arr, char *factor_ptr)
{
    const int type_num = PyArray_TYPE(arr);
    PyArrayObject *res = (PyArrayObject *)PyArray_EMPTY(
        PyArray_NDIM(arr), PyArray_DIMS(arr), type_num, PyArray_ISFORTRAN(arr)
    );
    if (res == NULL) {
        return NULL;
    }
    npy_intp n = PyArray_SIZE(arr);
    if (n == 0) {
        return (PyObject *)res; // Nothing to do.
    }
    PyUFuncGenericFunction loop;
    npy_intp strides[3];
    char *data[3] = {PyArray_DATA(arr), factor_ptr, PyArray_DATA(res)};
    if (type_num == NPY_DOUBLE || type_num == NPY_FLOAT) {
        strides[0] = PyArray_ITEMSIZE(arr);
        loop = loops[type_num - NPY_FLOAT];
    }
    else { // For complex, with real factor, we can use real loop.
        strides[0] = PyArray_ITEMSIZE(arr) / 2;
        n *= 2;
        loop = loops[type_num - NPY_CFLOAT];
    }
    strides[1] = 0;
    strides[2] = strides[0];
    PyUFunc_clearfperr();
    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS_THRESHOLDED(n);
    loop(data, &n, strides, NULL);
    NPY_END_THREADS;
    int fpe_errors = PyUFunc_getfperr();
    if (fpe_errors) {
        if (PyUFunc_GiveFloatingpointErrors("multiply", fpe_errors) < 0) {
            Py_DECREF(res);
            return NULL;
        }
    }
    return (PyObject *)res;
}

// Get the factor as a PyFloat, caching it as needed.
static inline PyObject *get_o_factor(ScalerObject *self)
{
    if (self->O_factor == NULL) {
        self->O_factor = PyFloat_FromDouble(self->factor);
        if (self->O_factor == NULL) {
            return NULL;
        }
    }
    return self->O_factor;
}

// Scale the input with the factor, using shortcuts where possible.
static PyObject *Scaler_vectorcall(
    ScalerObject *self, PyObject *const *args, size_t len_args, PyObject *kwnames
)
{
    if (kwnames != NULL && PyTuple_GET_SIZE(kwnames) > 0) {
        PyErr_Format(
            PyExc_TypeError,
            "scaler() got an unexpected keyword argument %R",
            PyTuple_GET_ITEM(kwnames, 0)
        );
        return NULL;
    }
    if (PyVectorcall_NARGS(len_args) != 1) {
        PyErr_Format(
            PyExc_TypeError,
            "scaler() takes 1 positional argument but %d were given",
            PyVectorcall_NARGS(len_args)
        );
        return NULL;
    }

    PyObject *const obj = args[0];
    // Fast paths for python double.
    if (PyFloat_CheckExact(obj)) {
        return PyFloat_FromDouble(PyFloat_AS_DOUBLE(obj) * self->factor);
    }
    // Fast path for plain ndarray, with special care for contiguous data.
    if (PyArray_CheckExact(obj)) {
        PyArrayObject *const arr = (PyArrayObject *)obj;
        const int type_num = PyArray_TYPE(arr);
        npy_bool needs_float = type_num == NPY_FLOAT || type_num == NPY_CFLOAT;
        npy_bool supported_fast_type =
            (needs_float || type_num == NPY_DOUBLE || type_num == NPY_CDOUBLE);
        char *f_ptr = needs_float ? (char *)&self->factor_f : (char *)&self->factor;
        // Pass contiguous float or complex arrays directly to multiply loop,
        // bypassing ufunc setup.
        if (supported_fast_type && PyArray_ISONESEGMENT(arr) && PyArray_ISNOTSWAPPED(arr)) {
            return use_contiguous_loop(arr, f_ptr);
        }
        // If not, convert factor to array here, since that makes ufunc
        // call substantially faster.  Use cached version if available.
        PyObject **A_factor = needs_float ? &self->A_factor_f : &self->A_factor;
        if (*A_factor == NULL) {
            const npy_intp dims[1] = {0};
            *A_factor =
                PyArray_SimpleNewFromData(0, dims, needs_float ? NPY_FLOAT : NPY_DOUBLE, f_ptr);
            if (*A_factor == NULL) {
                return NULL;
            }
        }
        return Py_TYPE(obj)->tp_as_number->nb_multiply(obj, *A_factor);
    }
// Fast paths for numpy double and float scalars; particularly useful for
// float32, where it avoids double->float conversion.
// Note: keep PyArray_VAL and ASSIGN where they are: they prevent the compiler
// from reordering the multiplication with the fperr calls.
#define FAST_PATH_NUMPY_SCALAR(obj, Type, npy_type, factor) \
    if (PyArray_IsScalar(obj, Type)) { \
        PyObject *res = PyArrayScalar_New(Type); \
        if (res == NULL) { \
            return NULL; \
        } \
        PyUFunc_clearfperr(); \
        npy_type x = PyArrayScalar_VAL(obj, Type); \
        npy_type out = x * factor; \
        PyArrayScalar_ASSIGN(res, Type, out); \
        int fpe_errors = PyUFunc_getfperr(); \
        if (fpe_errors && PyUFunc_GiveFloatingpointErrors("scalar multiply", fpe_errors) < 0) { \
            Py_CLEAR(res); \
        } \
        return res; \
    }
    FAST_PATH_NUMPY_SCALAR(obj, Double, npy_double, self->factor);
    FAST_PATH_NUMPY_SCALAR(obj, Float, npy_float, self->factor_f);
#undef FAST_PATH_NUMPY_SCALAR
    // Fast path for python integers.
    if (PyLong_CheckExact(obj)) {
        double d = PyLong_AsDouble(obj);
        if (d == -1.0 && PyErr_Occurred()) {
            return NULL;
        }
        return PyFloat_FromDouble(d * self->factor);
    }
    // For cases without special treatment, we need the float object.
    PyObject *O_factor = get_o_factor(self);
    if (O_factor == NULL) {
        return NULL;
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
                PyObject *res = slotv(obj, O_factor);
                if (res != Py_NotImplemented) {
                    return res;
                }
                // Fall back to slow path, which will likely raise,
                // thus giving a reasonable error message.
                Py_DECREF(res);
            }
        }
        return PyNumber_Multiply(obj, O_factor);
    }
    // If obj is not a known type, try converting it to an array.
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

static inline PyObject *Scaler_from_factor(PyTypeObject *type, double factor)
{
    ScalerObject *self = (ScalerObject *)type->tp_alloc(type, 0);
    if (self == NULL) {
        return NULL;
    }
    self->vectorcall = (vectorcallfunc)&Scaler_vectorcall;
    self->factor = factor;
    self->factor_f = (float)factor;
    // No need to explicitly set to zero, since zeroed by tp_alloc.
    // self->O_factor = NULL;  // initialized as needed in call.
    // self->A_factor = NULL;
    // self->A_factor_f = NULL;
    return (PyObject *)self;
}

static PyObject *Scaler_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double factor;
    Py_ssize_t nargs = PyTuple_GET_SIZE(args);
    if (nargs != 1 || kwds != NULL) {
        // Use parser to give error message.
        char *const kwlist[] = {"", NULL};
        PyArg_ParseTupleAndKeywords(args, kwds, "d:Scaler", kwlist, &factor);
        return NULL;
    }
    factor = PyFloat_AsDouble(PyTuple_GET_ITEM(args, 0));
    if (factor == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    if (factor == 1.0) { // Use singleton.
        return Py_NewRef(unity_scaler);
    }
    return Scaler_from_factor(type, factor);
}

// traverse: Visit all references from an object, including its type
// (since we're a heap type, which can get deallocated).
static int Scaler_traverse(ScalerObject *self, visitproc visit, void *arg)
{
    // Visit the type
    Py_VISIT(Py_TYPE(self));

    // Visit attributes that may hold references.
    Py_VISIT(self->O_factor);
    Py_VISIT(self->A_factor);
    Py_VISIT(self->A_factor_f);
    return 0;
}

// Clear internal references; called from finalize and dealloc.
static int Scaler_clear(ScalerObject *self)
{
    Py_CLEAR(self->O_factor);
    Py_CLEAR(self->A_factor);
    Py_CLEAR(self->A_factor_f);
    return 0;
}

static void Scaler_dealloc(ScalerObject *self)
{
    PyObject_GC_UnTrack(self);
    Scaler_clear(self);
    PyTypeObject *type = Py_TYPE(self);
    type->tp_free(self);
    Py_DECREF(type);
}

static PyObject *Scaler_repr(ScalerObject *self)
{
    PyObject *O_factor = get_o_factor(self);
    if (O_factor == NULL) {
        return NULL;
    }
    return PyUnicode_FromFormat("Scaler(%S)", O_factor);
}

static PyObject *Scaler_richcompare(PyObject *self, PyObject *other, int op)
{
    if (op != Py_EQ && op != Py_NE) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    npy_bool same =
        (Py_TYPE(other) == Py_TYPE(self) &&
         (((ScalerObject *)other)->factor == ((ScalerObject *)self)->factor));
    if ((op == Py_EQ) ^ same) {
        Py_RETURN_FALSE;
    }
    else {
        Py_RETURN_TRUE;
    }
}

static Py_hash_t Scaler_hash(PyObject *self)
{
    PyObject *O_factor = get_o_factor((ScalerObject *)self);
    if (O_factor == NULL) {
        return -1;
    }
    return (PyObject_Hash((PyObject *)Py_TYPE(self)) ^ PyObject_Hash(O_factor));
}

static PyObject *Scaler___reduce__(ScalerObject *self)
{
    return Py_BuildValue("(O(d))", Py_TYPE(self), self->factor);
}

static PyMemberDef Scaler_members[] = {
#if Py_Version < 0x03120000
    {"factor",
     T_DOUBLE,
     offsetof(ScalerObject, factor),
     READONLY,
     "Factor with which input is multiplied."},
    {"__vectorcalloffset__", T_PYSSIZET, offsetof(ScalerObject, vectorcall), READONLY},
#else
    {"factor",
     Py_T_DOUBLE,
     offsetof(ScalerObject, factor),
     Py_READONLY,
     "Factor with which input is multiplied."},
    {"__vectorcalloffset__", Py_T_PYSSIZET, offsetof(ScalerObject, vectorcall), Py_READONLY},
#endif
    {NULL},
};

static PyMethodDef Scaler_methods[] = {
    {"__reduce__", (PyCFunction)Scaler___reduce__, METH_NOARGS, NULL},
    {NULL, NULL, 0, NULL},
};

static PyType_Slot Scaler_slots[] = {
    {Py_tp_doc, Scaler_doc},
    {Py_tp_new, (newfunc)Scaler_new},
    {Py_tp_traverse, (traverseproc)Scaler_traverse},
    {Py_tp_clear, (inquiry)Scaler_clear},
    {Py_tp_dealloc, (destructor)Scaler_dealloc},
    {Py_tp_repr, (reprfunc)Scaler_repr},
    {Py_tp_richcompare, (richcmpfunc)Scaler_richcompare},
    {Py_tp_hash, (hashfunc)Scaler_hash},
    {Py_tp_call, &PyVectorcall_Call},
    {Py_tp_members, Scaler_members},
    {Py_tp_methods, Scaler_methods},
    {0, NULL},
};

static PyType_Spec Scaler_spec = {
    .name = "scaler.Scaler",
    .basicsize = sizeof(ScalerObject),
    .flags =
        Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_VECTORCALL | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_BASETYPE,
    .slots = Scaler_slots,
};

// Get and cache multiplication loops for use with plain ndarray.
// Only gets double and float loops, since those can be used for
// complex too. We ignore long double for our fast path.
static int get_multiply_loops(PyUFuncGenericFunction loops[2])
{
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
    for (int k = 0; k < 2; k++) {
        char type = (char)(k + NPY_FLOAT);
        const char *types = multiply->types;
        for (int i = 0; i < multiply->ntypes; i++) {
            if (types[0] == type && types[1] == type && types[2] == type) {
                loops[k] = multiply->functions[i];
                break;
            }
            types += 3;
        }
        if (loops[k] == NULL) { // Should never happen.
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
    PyTypeObject *Scaler = (PyTypeObject *)PyType_FromModuleAndSpec(m, &Scaler_spec, NULL);
    if (Scaler == NULL) {
        return -1;
    }
    if (PyModule_AddType(m, Scaler) < 0) {
        return -1;
    }
    unity_scaler = Scaler_from_factor(Scaler, 1.0);
    if (unity_scaler == NULL) {
        return -1;
    }
    if (get_multiply_loops(loops) < 0) {
        return -1;
    }
    return 0;
}

static PyModuleDef_Slot scaler_module_slots[] = {
    {Py_mod_exec, scaler_module_exec},
#if Py_Version >= 0x03120000
    {Py_mod_multiple_interpreters, Py_MOD_MULTIPLE_INTERPRETERS_NOT_SUPPORTED},
#endif
    {0, NULL}
};

static PyModuleDef scaler_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "scaler",
    .m_doc = scaler_module_doc,
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
