// This extension is adapted from the positional_defaults PyPI package
// https://pypi.org/project/positional-defaults/ version 2023.4.19
// MIT. see licenses/POSITIONAL_DEFAULTS.rst


#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <structmember.h>

#define SINCE_CHAR_SIZE 32
#define NAMES_CHAR_SIZE 128
#define MSG_SIZE 512

typedef struct {
    PyObject_HEAD
    PyObject* dict;
    PyObject* wrapped;
    PyObject* names;
    PyObject* since;
} DeprKwsObject;



static void
depr_kws_wrap_dealloc(DeprKwsObject* self)
{
    Py_XDECREF(self->wrapped);
    Py_XDECREF(self->names);
    Py_XDECREF(self->since);
    Py_TYPE(self)->tp_free((PyObject*)self);
}


static PyObject*
depr_kws_wrap_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    DeprKwsObject* self = (DeprKwsObject*)type->tp_alloc(type, 0);

    if (self != NULL) {
        self->names = PyTuple_New(0);
        if (self->names == NULL) {
            Py_DECREF(self);
            return NULL;
        }

        Py_INCREF(Py_None);
        self->wrapped = Py_None;

        Py_INCREF(Py_None);
        self->since = Py_None;
    }

    return (PyObject*)self;
}


static int
depr_kws_wrap_init(DeprKwsObject* self, PyObject* args, PyObject* kwds)
{
    static char *kwlist[] = {"wrapped", "names", "since", NULL};
    Py_ssize_t i, n_names;
    PyObject* wrapped, *names, *since, *tmp;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOO:wrap", kwlist,
                                     &wrapped, &names, &since))
        return -1;

    if (!PyTuple_Check(names)) {
        PyErr_SetString(PyExc_TypeError, "names must be a tuple");
        return -1;
    }

    n_names = PyTuple_GET_SIZE(names);

    for (i = 0; i < n_names; ++i) {
        PyObject* name = PyTuple_GET_ITEM(names, i);
        if (!PyUnicode_Check(name)) {
            PyErr_Format(PyExc_TypeError, "names[%zd] must be a string", i);
            return -1;
        }
    }

    if (!PyUnicode_Check(since)) {
        PyErr_Format(PyExc_TypeError, "since must be a string", i);
        return -1;
    }

    tmp = self->wrapped;
    Py_INCREF(wrapped);
    self->wrapped = wrapped;
    Py_XDECREF(tmp);

    tmp = self->names;
    Py_INCREF(names);
    self->names = names;
    Py_XDECREF(tmp);

    tmp = self->since;
    Py_INCREF(since);
    self->since = since;
    Py_XDECREF(tmp);
    return 0;
}


static PyMemberDef depr_kws_wrap_members[] = {
    {"__dict__", T_OBJECT, offsetof(DeprKwsObject, dict), READONLY},
    {"wrapped", T_OBJECT, offsetof(DeprKwsObject, wrapped), READONLY},
    {"names", T_OBJECT, offsetof(DeprKwsObject, names), READONLY},
    {"since", T_OBJECT, offsetof(DeprKwsObject, since), READONLY},
    {NULL}
};


static PyObject*
depr_kws_wrap_call(DeprKwsObject* self, PyObject* args, PyObject* kwds) {
    // step 0: return early whenever possible
    if (self->wrapped == NULL)
        Py_RETURN_NONE;

    if (kwds == NULL)
        return PyObject_Call(self->wrapped, args, kwds);

    // step 1: detect any deprecated keyword arguments, return if none.
    Py_ssize_t n_names = PyTuple_GET_SIZE(self->names);
    PyObject *deprecated_kwargs = PyList_New(n_names);
    Py_INCREF(deprecated_kwargs);
    PyObject *name = NULL;
    Py_ssize_t i = 0;
    int has_kw = -2;

    Py_ssize_t n_depr = 0;
    for (i=0 ; i < n_names ; ++i) {
        name = PyTuple_GET_ITEM(self->names, i);
        has_kw = PyDict_Contains(kwds, name);
        if (has_kw) {
            PyList_SET_ITEM(deprecated_kwargs, n_depr, name);
            ++n_depr;
        }
    }

    if (n_depr == 0)
        return PyObject_Call(self->wrapped, args, kwds);

    // step 2: create and emit warning message
    char names_char[NAMES_CHAR_SIZE];
    char *s, *arguments, *respectively, *pronoun;

    PyObject *names_unicode;
    if (n_depr > 1) {
        names_unicode = PyObject_Str(PyList_GetSlice(deprecated_kwargs, 0, n_depr));
        s = "s";
        arguments = " arguments";
        respectively = ", respectively";
        pronoun = "them";
    } else {
        names_unicode = PyObject_Repr(PyList_GET_ITEM(deprecated_kwargs, 0));
        s = arguments = respectively = "";
        pronoun = "it";
    }
    const char* names_utf8 = PyUnicode_AsUTF8(names_unicode);
    snprintf(names_char, NAMES_CHAR_SIZE, "%s", names_utf8);

    PyObject *since_unicode = PyObject_Str(self->since);
    const char* since_utf8 = PyUnicode_AsUTF8(since_unicode);
    char since_char[SINCE_CHAR_SIZE];
    snprintf(since_char, SINCE_CHAR_SIZE, "%s", since_utf8);

    char msg[MSG_SIZE];
    snprintf(
        msg,
        MSG_SIZE,
        "Passing %s%s as keyword%s "
        "is deprecated since version %s "
        "and will stop working in a future release. "
        "Pass %s positionally to suppress this warning.",
        names_char, arguments, s, since_char, pronoun
    );
    const char* msg_ptr = msg;

    int status = PyErr_WarnEx(PyExc_FutureWarning, msg_ptr, 2);
    if (status == -1) {
        // avoid leaking memory if Warning is promoted to Exception
        Py_DECREF(deprecated_kwargs);
    }

    return PyObject_Call(self->wrapped, args, kwds);
}


static PyObject*
depr_kws_wrap_get(PyObject* self, PyObject* obj, PyObject* type) {
    if (obj == Py_None || obj == NULL) {
        Py_INCREF(self);
        return self;
    }
    return PyMethod_New(self, obj);
}


static PyTypeObject DeprKwsWrap = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "signature_deprecations.wrap",
    .tp_doc = PyDoc_STR("wrap a function with deprecated keyword arguments"),
    .tp_basicsize = sizeof(DeprKwsObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_dictoffset = offsetof(DeprKwsObject, dict),
    .tp_new = depr_kws_wrap_new,
    .tp_init = (initproc)depr_kws_wrap_init,
    .tp_dealloc = (destructor)depr_kws_wrap_dealloc,
    .tp_members = depr_kws_wrap_members,
    .tp_call = (ternaryfunc)depr_kws_wrap_call,
    .tp_descr_get = depr_kws_wrap_get,
};


static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "signature_deprecations",
    .m_doc = PyDoc_STR("fast decorators to mark signature details as deprecated"),
    .m_size = -1,
};


PyMODINIT_FUNC
PyInit_signature_deprecations(void) {
    PyObject* m;

    if (PyType_Ready(&DeprKwsWrap) < 0)
        return NULL;

    m = PyModule_Create(&module);
    if (m == NULL)
        return NULL;

    Py_INCREF(&DeprKwsWrap);
    if (PyModule_AddObject(m, "_depr_kws_wrap", (PyObject*)&DeprKwsWrap) < 0) {
        Py_DECREF(&DeprKwsWrap);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
