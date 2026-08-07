// Licensed under a 3-clause BSD style license - see LICENSE.rst
#ifndef _SCALER_LIMITED_API_WORKAROUNDS_H
#define _SCALER_LIMITED_API_WORKAROUNDS_H

// Once we're at 3.13, move SCALER_TP_FLAGS to its use in scaler.c,
// and remove this whole include file.
#if !defined(Py_LIMITED_API) || Py_LIMITED_API + 0 >= 0x030D0000
#define SCALER_TP_FLAGS \
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_VECTORCALL | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_BASETYPE
#else // Limited API 3.12 or 3.11
#if Py_LIMITED_API + 0 < 0x030C0000
// Some vectorcall stuff was not yet part of the limited API in 3.11
// Includes a simple substitution for PyVectorcall_Call that is specific
// to our own vectorcall implementation (hence the forward definition).
#define Py_TPFLAGS_HAVE_VECTORCALL (1UL << 11)
#define PY_VECTORCALL_ARGUMENTS_OFFSET (_Py_STATIC_CAST(size_t, 1) << (8 * sizeof(size_t) - 1))
#define PyVectorcall_NARGS(nargsf) (Py_ssize_t)(nargsf & ~PY_VECTORCALL_ARGUMENTS_OFFSET)
#define vectorcallfunc void *
static PyObject *Scaler_vectorcall(
    PyObject *self, PyObject *const *args, size_t len_args, PyObject *kwnames
);
static PyObject *PyVectorcall_Call(PyObject *self, PyObject *tuple, PyObject *dict)
{
    if (dict != NULL && PyDict_Size(dict) > 0) {
        Py_ssize_t pos = 0;
        PyObject *key, *value;
        PyDict_Next(dict, &pos, &key, &value);
        PyErr_Format(PyExc_TypeError, "scaler() got an unexpected keyword argument %R", key);
        return NULL;
    }
    if (PyTuple_Size(tuple) != 1) {
        PyErr_Format(
            PyExc_TypeError,
            "scaler() takes 1 positional argument but %d were given",
            PyTuple_Size(tuple)
        );
        return NULL;
    }
    PyObject *args[1];
    args[0] = PyTuple_GetItem(tuple, 0);
    if (args[0] == NULL) {
        return NULL;
    }
    return Scaler_vectorcall(self, args, 1, NULL);
}
#endif // Py_LIMITED_API + 0 < 0x030C0000
// PyType_GetModuleByDef is part of the limited API only since 3.13.
// Its definition is rather long so we use a nearly good enough substitute,
// the only downside of which is that one cannot subclass (hence the adjustment
// to the flags). Not that important really, since it should not be needed
// in normal use, but we might as well allow it once we can.
#define PyType_GetModuleByDef(type, unused) PyType_GetModule(type)
#define SCALER_TP_FLAGS Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_VECTORCALL | Py_TPFLAGS_HAVE_GC
#endif // defined(Py_LIMITED_API) && Py_LIMITED_API + 0 >= 0x030D0000

// Other changes needed for python 3.11, with our without limited API.
#if PY_VERSION_HEX < 0x030C0000
#include <structmember.h> // For PyMemberDef, included by default only for >=3.12
#define Py_T_DOUBLE T_DOUBLE
#define Py_READONLY READONLY
#define Py_T_PYSSIZET T_PYSSIZET
#endif // PY_VERSION_HEX < 0x030C0000

#endif // _SCALER_LIMITED_API_WORKAROUNDS_H
