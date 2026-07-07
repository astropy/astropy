/*
 Exercises the public astropy.wcs C API exposed through the AstropyWcs_API
 function-pointer table (discoverable from astropy.wcs.get_include()).  This is
 the same entry point a downstream package such as drizzlepac uses.  The module
 provides one Python-callable function per member, or group of members, of the
 table that astropy still supports, plus one that calls a deprecated member.
*/

#include "Python.h"
#include "astropy_wcs_api.h"
#include "astropy_wcs.h"
#include "wcslib_wrap.h"

/* Slot 0: the version returned through the table. */
static PyObject*
get_c_version(PyObject* self, PyObject* args) {
    return PyLong_FromLong((long)AstropyWcs_GetCVersion());
}

/* Slot 23: wcslib_get_error_message. */
static PyObject*
get_error_message(PyObject* self, PyObject* arg) {
    long code = PyLong_AsLong(arg);
    const char* msg;
    if (code == -1 && PyErr_Occurred()) {
        return NULL;
    }
    msg = wcslib_get_error_message((int)code);
    return PyUnicode_FromString(msg == NULL ? "" : msg);
}

/* Slots 1, 2, 20, 21: bracket wcsp2s / wcss2p with the wcsprm_python2c /
   wcsprm_c2python pair (exactly as drizzlepac does) and round-trip a single
   two-dimensional coordinate through pixel -> world -> pixel. */
static PyObject*
roundtrip(PyObject* self, PyObject* args) {
    PyObject* wcsprm_obj;
    double x, y;
    struct wcsprm* w;
    double pixcrd[2], imgcrd[2], world[2], pix2[2];
    double phi[1], theta[1];
    int stat[1];
    int s1, s2;

    if (!PyArg_ParseTuple(args, "Odd", &wcsprm_obj, &x, &y)) {
        return NULL;
    }
    w = &((Wcsprm*)wcsprm_obj)->x;

    pixcrd[0] = x;
    pixcrd[1] = y;

    wcsprm_python2c(w);
    s1 = wcsp2s(w, 1, 2, pixcrd, imgcrd, phi, theta, world, stat);
    s2 = wcss2p(w, 1, 2, world, phi, theta, imgcrd, pix2, stat);
    wcsprm_c2python(w);

    if (s1 != 0 || s2 != 0) {
        PyErr_Format(PyExc_RuntimeError,
                     "wcslib transform failed (p2s=%d, s2p=%d)", s1, s2);
        return NULL;
    }
    return Py_BuildValue("(dd)", pix2[0], pix2[1]);
}

/* Slot 22: wcsprt. */
static PyObject*
print_wcsprm(PyObject* self, PyObject* arg) {
    struct wcsprm* w = &((Wcsprm*)arg)->x;
    return PyLong_FromLong((long)wcsprt(w));
}

/* Slot 18: pipeline_all_pixel2world.  astropy applies preoffset_array to the
   pixel coordinates before calling it, which adds (1 - origin) so that the
   pipeline always sees the one-based (FITS) convention wcslib uses; we do the
   same here. */
static PyObject*
all_pix2world(PyObject* self, PyObject* args) {
    PyObject* wcs_obj;
    double x, y;
    int origin;
    pipeline_t* p;
    double pixcrd[2], world[2];
    int status;

    if (!PyArg_ParseTuple(args, "Oddi", &wcs_obj, &x, &y, &origin)) {
        return NULL;
    }
    p = &((Wcs*)wcs_obj)->x;

    pixcrd[0] = x + (1.0 - (double)origin);
    pixcrd[1] = y + (1.0 - (double)origin);

    status = pipeline_all_pixel2world(p, 1, 2, pixcrd, world);
    if (status != 0 && status != 8) {
        PyErr_Format(PyExc_RuntimeError,
                     "pipeline_all_pixel2world failed (status=%d)", status);
        return NULL;
    }
    return Py_BuildValue("(dd)", world[0], world[1]);
}

/* A deprecated member of the table, mirroring the call the original
   wcsapi_test.c (removed in #12489) made. */
static PyObject*
call_deprecated(PyObject* self, PyObject* args) {
    pipeline_t p;
    pipeline_clear(&p);
    Py_RETURN_NONE;
}

static PyMethodDef module_methods[] = {
    {"get_c_version", get_c_version, METH_NOARGS, ""},
    {"get_error_message", get_error_message, METH_O, ""},
    {"roundtrip", roundtrip, METH_VARARGS, ""},
    {"print_wcsprm", print_wcsprm, METH_O, ""},
    {"all_pix2world", all_pix2world, METH_VARARGS, ""},
    {"call_deprecated", call_deprecated, METH_NOARGS, ""},
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "wcsapi_test",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit_wcsapi_test(void) {
    PyObject* m = PyModule_Create(&moduledef);
    if (m == NULL) {
        return NULL;
    }

    import_astropy_wcs();

    if (PyErr_Occurred()) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
