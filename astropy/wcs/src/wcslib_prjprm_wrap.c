#define NO_IMPORT_ARRAY
#include <math.h>
#include <float.h>

#include "astropy_wcs/wcslib_celprm_wrap.h"
#include "astropy_wcs/wcslib_prjprm_wrap.h"

#include <wcs.h>
#include <wcsprintf.h>
#include <prj.h>
#include <numpy/npy_math.h>
#include "astropy_wcs/docstrings.h"


PyObject** prj_errexc[5];


static int is_dbl_equal(double x1, double x2)
{
    double ax1 = fabs(x1);
    double ax2 = fabs(x2);
    double minx = (ax1 < ax2) ? ax1 : ax2;
    double diff = fabs(x1 - x2);
    return (diff <= (2.0 * DBL_EPSILON * minx) || diff < DBL_MIN);
}


static int wcslib_prj_to_python_exc(int status)
{
    if (status > 0 && status < 5) {
        PyErr_SetString(*prj_errexc[status], prj_errmsg[status]);
    } else if (status > 5) {
        PyErr_SetString(
            PyExc_RuntimeError,
            "Unknown WCSLIB prjprm-related error occurred.");
    }
    return status;
}


static int is_readonly(PyPrjprm* self)
{
    if (self != NULL && self->owner != NULL &&
        ((PyCelprm*)self->owner)->owner != NULL) {
        PyErr_SetString(
            PyExc_AttributeError,
            "Attribute 'prj' of 'astropy.wcs.Wcsprm.cel' objects is read-only.");
        return 1;
    } else {
        return 0;
    }
}


static int is_prj_null(PyPrjprm* self)
{
    if (self->x == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Underlying 'prjprm' object is NULL.");
        return 1;
    } else {
        return 0;
    }
}


/***************************************************************************
 * PyPrjprm methods                                                        *
 ***************************************************************************/

static PyObject* PyPrjprm_new(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
    PyPrjprm* self;
    self = (PyPrjprm*)type->tp_alloc(type, 0);
    if (self == NULL) return NULL;
    self->owner = NULL;
    self->x = NULL;
    self->prefcount = NULL;
    if ((self->x = calloc(1, sizeof(struct prjprm))) == 0x0) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory.");
        return NULL;
    }
    if ((self->prefcount = (int*) malloc(sizeof(int))) == 0x0) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory.");
        free(self->x);
        return NULL;
    }
    if (wcslib_prj_to_python_exc(prjini(self->x)))
    {
        free(self->x);
        free(self->prefcount);
        return NULL;
    }
    *(self->prefcount) = 1;
    return (PyObject*)self;
}


static int PyPrjprm_traverse(PyPrjprm* self, visitproc visit, void *arg)
{
    Py_VISIT(self->owner);
    return 0;
}


static int PyPrjprm_clear(PyPrjprm* self)
{
    Py_CLEAR(self->owner);
    return 0;
}


static void PyPrjprm_dealloc(PyPrjprm* self)
{
    PyPrjprm_clear(self);
    if (self->prefcount && (--(*self->prefcount)) == 0) {
        wcslib_prj_to_python_exc(prjfree(self->x));
        free(self->x);
        free(self->prefcount);
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}


PyPrjprm* PyPrjprm_cnew(PyObject* celprm_obj, struct prjprm* x, int* prefcount)
{
    PyPrjprm* self;
    self = (PyPrjprm*)(&PyPrjprmType)->tp_alloc(&PyPrjprmType, 0);
    if (self == NULL) return NULL;
    self->x = x;
    Py_XINCREF(celprm_obj);
    self->owner = celprm_obj;
    self->prefcount = prefcount;
    if (prefcount) (*prefcount)++;
    return self;
}


static PyObject* PyPrjprm_copy(PyPrjprm* self)
{
    PyPrjprm* copy = NULL;
    copy = PyPrjprm_cnew(self->owner, self->x, self->prefcount);
    if (copy == NULL) return NULL;
    return (PyObject*)copy;
}


static PyObject* PyPrjprm_deepcopy(PyPrjprm* self)
{
    PyPrjprm* copy = (PyPrjprm*) PyPrjprm_new(&PyPrjprmType, NULL, NULL);
    if (copy == NULL) return NULL;

    memcpy(copy->x, self->x, sizeof(struct prjprm));
    copy->x->err = NULL;
    return (PyObject*)copy;
}


static PyObject* PyPrjprm___str__(PyPrjprm* self)
{
    wcsprintf_set(NULL);
    if (wcslib_prj_to_python_exc(prjprt(self->x))) {
        return NULL;
    }
    return PyUnicode_FromString(wcsprintf_buf());
}


static int PyPrjprm_cset(PyPrjprm* self)
{
    if (wcslib_prj_to_python_exc(prjset(self->x))) {
        return -1;
    }
    return 0;
}


static PyObject* PyPrjprm_set(PyPrjprm* self)
{
    if (is_readonly(self) || PyPrjprm_cset(self)) return NULL;
    Py_RETURN_NONE;
}


static PyObject* _prj_eval(PyPrjprm* self, int (*prjfn)(PRJX2S_ARGS),
                           PyObject* x1_in, PyObject* x2_in)
{
    Py_ssize_t i, ndim;
    npy_intp *x1_dims, *x2_dims;
    Py_ssize_t     nelem      = 1;
    PyArrayObject* x1         = NULL;
    PyArrayObject* x2         = NULL;
    PyArrayObject* prj_x1     = NULL;
    PyArrayObject* prj_x2     = NULL;
    PyArrayObject* stat       = NULL;
    PyObject*      result     = NULL;
    int            status     = -1;

    // TODO: This assumes the same shape for the input arrays.
    //       Instead, we should broadcast.

    x1 = (PyArrayObject *) PyArray_ContiguousFromObject(x1_in, NPY_DOUBLE, 1, NPY_MAXDIMS);
    if (x1 == NULL) {
      goto exit;
    }
    x2 = (PyArrayObject *) PyArray_ContiguousFromObject(x2_in, NPY_DOUBLE, 1, NPY_MAXDIMS);
    if (x2 == NULL) {
      goto exit;
    }

    ndim = PyArray_NDIM(x1);

    if (ndim != PyArray_NDIM(x2)) {
        PyErr_SetString(PyExc_ValueError, "Input array dimensions do not match.");
        goto exit;
    }

    x1_dims = PyArray_DIMS(x1);
    x2_dims = PyArray_DIMS(x2);
    for (i = 0; i < ndim; i++) {
        if (x1_dims[i] != x2_dims[i]) {
            PyErr_SetString(PyExc_ValueError, "Input array dimensions do not match.");
            goto exit;
        }
        nelem *= x1_dims[i];
    }

    prj_x1 = (PyArrayObject*)PyArray_SimpleNew(ndim, x1_dims, NPY_DOUBLE);
    if (prj_x1 == NULL) {
      goto exit;
    }

    prj_x2 = (PyArrayObject*)PyArray_SimpleNew(ndim, x1_dims, NPY_DOUBLE);
    if (prj_x2 == NULL) {
      goto exit;
    }

    stat = (PyArrayObject*)PyArray_SimpleNew(ndim, x1_dims, NPY_INT);
    if (stat == NULL) {
      goto exit;
    }

    Py_BEGIN_ALLOW_THREADS
    status = prjfn(
        self->x,
        nelem, 0, 1, 1,
        (double*)PyArray_DATA(x1),
        (double*)PyArray_DATA(x2),
        (double*)PyArray_DATA(prj_x1),
        (double*)PyArray_DATA(prj_x2),
        (int*)PyArray_DATA(stat));
    Py_END_ALLOW_THREADS

    switch (status) {
    case 3:
    case 4:
        for (i = 0; i < nelem; ++i) {
            if (((int *)PyArray_DATA(stat))[i]) {
                ((double *)PyArray_DATA(prj_x1))[i] = NPY_NAN;
                ((double *)PyArray_DATA(prj_x2))[i] = NPY_NAN;
            }
        }
    case 0:
        result = Py_BuildValue("(OO)", prj_x1, prj_x2);
        break;
    default:
        wcslib_prj_to_python_exc(status);
        break;
    }

    exit:
        Py_XDECREF(x1);
        Py_XDECREF(x2);
        Py_XDECREF(prj_x1);
        Py_XDECREF(prj_x2);
        Py_XDECREF(stat);

    return result;
}


static PyObject* PyPrjprm_prjx2s(PyPrjprm* self, PyObject* args, PyObject* kwds)
{
    PyObject* x = NULL;
    PyObject* y = NULL;
    const char* keywords[] = { "x", "y", NULL };

    if (is_prj_null(self)) return NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO:prjx2s",
        (char **)keywords, &x, &y)) {
        return NULL;
    }

    if (self->x->prjx2s == NULL || self->x->flag == 0) {
        if (is_readonly(self)) {
            PyErr_SetString(
                PyExc_AttributeError,
                "Attribute 'prj' of 'astropy.wcs.Wcsprm.cel' objects is "
                "read-only and cannot be automatically set.");
            return NULL;
        } else if (PyPrjprm_cset(self)) {
            return NULL;
        }
    }

    return _prj_eval(self, self->x->prjx2s, x, y);
}


static PyObject* PyPrjprm_prjs2x(PyPrjprm* self, PyObject* args, PyObject* kwds)
{
    PyObject* phi = NULL;
    PyObject* theta = NULL;
    const char* keywords[] = { "phi", "theta", NULL };

    if (is_prj_null(self)) return NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO:prjs2x",
        (char **)keywords, &phi, &theta)) {
        return NULL;
    }

    if (self->x->prjs2x == NULL || self->x->flag == 0) {
        if (is_readonly(self)) {
            PyErr_SetString(
                PyExc_AttributeError,
                "Attribute 'prj' of 'astropy.wcs.Wcsprm.cel' objects is "
                "read-only and cannot be automatically set.");
            return NULL;
        } else if (PyPrjprm_cset(self)) {
            return NULL;
        }
    }

    return _prj_eval(self, self->x->prjs2x, phi, theta);
}


/***************************************************************************
 * Member getters/setters (properties)
 */

static PyObject* PyPrjprm_get_flag(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return get_int("flag", self->x->flag);
    }
}


static PyObject* PyPrjprm_get_code(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return get_string("code", self->x->code);
    }
}


static int PyPrjprm_set_code(PyPrjprm* self, PyObject* value, void* closure)
{
    char code[4];
    int code_len;

    if (is_prj_null(self) || is_readonly(self)) {
        return -1;
    } else if (value == Py_None) {
        if (strcmp("   ", self->x->code)) {
            strcpy(self->x->code, "   ");
            self->x->flag = 0;
            if (self->owner) ((PyCelprm*)self->owner)->x->flag = 0;
        }
    } else {
        if (set_string("code", value, code, 4)) return -1;
        code_len = strlen(code);
        if (code_len != 3) {
            PyErr_Format(PyExc_ValueError,
                "'code' must be exactly a three character string. "
                "Provided 'code' ('%s') is %d characters long.",
                code, code_len);
            return -1;
        }
        if (strcmp(code, self->x->code)) {
            strncpy(self->x->code, code, 4);
            self->x->code[3] = '\0';  /* just to be safe */
            self->x->flag = 0;
            if (self->owner) ((PyCelprm*)self->owner)->x->flag = 0;
        }
    }
    return 0;
}


static PyObject* PyPrjprm_get_r0(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else if (self->x->r0 == UNDEFINED) {
        Py_RETURN_NONE;
    } else {
        return get_double("r0", self->x->r0);
    }
}


static int PyPrjprm_set_r0(PyPrjprm* self, PyObject* value, void* closure)
{
    int result;
    double r0;
    if (is_prj_null(self) || is_readonly(self)) {
        return -1;
    } else if (value == Py_None) {
        if (self->x->r0 != UNDEFINED) {
            self->x->r0 = UNDEFINED;
            self->x->flag = 0;
            if (self->owner) ((PyCelprm*)self->owner)->x->flag = 0;
        }
    } else {
        result = set_double("r0", value, &r0);
        if (result) return result;
        if (r0 != self->x->r0) {
            self->x->r0 = r0;
            self->x->flag = 0;
            if (self->owner) ((PyCelprm*)self->owner)->x->flag = 0;
        }
    }
    return 0;
}


static PyObject* PyPrjprm_get_phi0(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else if (self->x->phi0 == UNDEFINED) {
        Py_RETURN_NONE;
    } else {
        return get_double("phi0", self->x->phi0);
    }
}


static int PyPrjprm_set_phi0(PyPrjprm* self, PyObject* value, void* closure)
{
    int result;
    double phi0;
    if (is_prj_null(self) || is_readonly(self)) {
        return -1;
    } else if (value == Py_None) {
        if (self->x->phi0 != UNDEFINED) {
            self->x->phi0 = UNDEFINED;
            self->x->flag = 0;
            if (self->owner) ((PyCelprm*)self->owner)->x->flag = 0;
        }
    } else {
        result = set_double("phi0", value, &phi0);
        if (result) return result;
        if (phi0 != self->x->phi0) {
            self->x->phi0 = phi0;
            self->x->flag = 0;
            if (self->owner) ((PyCelprm*)self->owner)->x->flag = 0;
        }
    }
    return 0;
}


static PyObject* PyPrjprm_get_theta0(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else if (self->x->theta0 == UNDEFINED) {
        Py_RETURN_NONE;
    } else {
        return get_double("theta0", self->x->theta0);
    }
}


static int PyPrjprm_set_theta0(PyPrjprm* self, PyObject* value, void* closure)
{
    int result;
    double theta0;
    if (is_prj_null(self) || is_readonly(self)) {
        return -1;
    } else if (value == Py_None) {
        if (self->x->theta0 != UNDEFINED) {
            self->x->theta0 = UNDEFINED;
            self->x->flag = 0;
            if (self->owner) ((PyCelprm*)self->owner)->x->flag = 0;
        }
    } else {
        result = set_double("theta0", value, &theta0);
        if (result) return result;
        if (theta0 != self->x->theta0) {
            self->x->theta0 = theta0;
            self->x->flag = 0;
            if (self->owner) ((PyCelprm*)self->owner)->x->flag = 0;
        }
    }
    return 0;
}


static PyObject* PyPrjprm_get_pv(PyPrjprm* self, void* closure)
{
    int k;
    Py_ssize_t size = PVN;
    double *pv;
    PyObject* pv_pyobj;
    PyArrayObject* pv_array;
    if (is_prj_null(self)) return NULL;

    pv_pyobj = PyArray_SimpleNew(1, &size, NPY_DOUBLE);
    pv_array = (PyArrayObject*) pv_pyobj;
    if (pv_array == NULL) return NULL;
    pv = (double*) PyArray_DATA(pv_array);

    for (k = 0; k < PVN; k++) {
        if (self->x->pv[k] == UNDEFINED) {
            pv[k] = (double) NPY_NAN;
        } else {
            pv[k] = self->x->pv[k];
        }
    }

    return pv_pyobj;
}


static int PyPrjprm_set_pv(PyPrjprm* self, PyObject* value, void* closure)
{
    int k, modified;
    npy_intp size;
    double *data;
    PyArrayObject* value_array = NULL;
    int skip[PVN];

    if (is_prj_null(self) || is_readonly(self)) return -1;

    if (value == Py_None) {
        /* If pv is set to None - reset pv to prjini values: */
        self->x->pv[0] = 0.0;
        for (k = 1; k < 4; self->x->pv[k++] = UNDEFINED);
        for (k = 4; k < PVN; self->x->pv[k++] = 0.0);
        self->x->flag = 0;
        if (self->owner) ((PyCelprm*)self->owner)->x->flag = 0;
        return 0;
    }

    value_array = (PyArrayObject*) PyArray_ContiguousFromAny(value, NPY_DOUBLE, 1, 1);
    if (!value_array) return -1;

    size = PyArray_SIZE(value_array);

    if (size < 1) {
        Py_DECREF(value_array);
        PyErr_SetString(PyExc_ValueError,
            "PV must be a non-empty 1-dimentional list of values or None.");
        return -1;
    }

    if (size > PVN) {
        Py_DECREF(value_array);
        PyErr_Format(PyExc_RuntimeError, "Number of PV values cannot exceed %d.", PVN);
        return -1;
    }

    if (PyList_Check(value)) {
        for (k = 0; k < size; k++) {
            skip[k] = (PyList_GetItem(value, k) == Py_None);
        }
    } else if (PyTuple_Check(value)) {
        for (k = 0; k < size; k++) {
            skip[k] = (PyTuple_GetItem(value, k) == Py_None);
        }
    } else {
        for (k = 0; k < size; k++) skip[k] = 0;
    }

    data = (double*) PyArray_DATA(value_array);

    modified = 0;
    for (k = 0; k < size; k++) {
        if (skip[k]) continue;
        if (is_dbl_equal(self->x->pv[k], data[k])) {
            /* update PV but do not flag it as modified since values are
               essentially the same.
            */
            self->x->pv[k] = data[k];
        } else if (npy_isnan(data[k])) {
            self->x->pv[k] = UNDEFINED;
            modified = 1;
        } else {
            self->x->pv[k] = data[k];
            modified = 1;
        }
    }
    Py_DECREF(value_array);

    if (modified) {
        self->x->flag = 0;
        if (self->owner) ((PyCelprm*)self->owner)->x->flag = 0;
    }
    return 0;
}


static PyObject* PyPrjprm_get_pvi(PyPrjprm* self, PyObject* args, PyObject* kwds)
{
    int idx;
    PyObject* index = NULL;
    PyObject* value = NULL;
    const char* keywords[] = { "index", NULL };

    if (is_prj_null(self)) return NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O:get_pvi",
        (char **)keywords, &index)) {
        return NULL;
    }

    if (!PyLong_Check(index)) {
        PyErr_SetString(PyExc_TypeError,
        "PV index must be an integer number.");
    }

    idx = PyLong_AsLong(index);
    if (idx == -1 && PyErr_Occurred()) {
        return NULL;
    }

    if (idx < 0 || idx >= PVN) {
        PyErr_Format(PyExc_ValueError,
            "PV index must be an integer number between 0 and %d.", PVN - 1);
        return NULL;
    }

    if (self->x->pv[idx] == UNDEFINED) {
        return PyFloat_FromDouble((double) NPY_NAN);
    } else {
        return PyFloat_FromDouble(self->x->pv[idx]);
    }
}


static PyObject* PyPrjprm_set_pvi(PyPrjprm* self, PyObject* args, PyObject* kwds)
{
    int idx, size;
    double data;
    PyObject* scalar= NULL;
    PyObject* index = NULL;
    PyObject* value = NULL;
    PyObject* flt_value = NULL;
    PyObject* value_array_pyobj = NULL;
    PyArrayObject* value_array = NULL;
    const char* keywords[] = { "index", "value", NULL };
    PyArray_Descr* dbl_descr = PyArray_DescrNewFromType(NPY_DOUBLE);

    if (is_prj_null(self) || is_readonly(self)) return NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO:set_pvi",
        (char **)keywords, &index, &value)) {
        return NULL;
    }

    if (!PyLong_Check(index)) {
        PyErr_SetString(PyExc_TypeError,
        "PV index must be an integer number.");
    }

    idx = PyLong_AsLong(index);
    if (idx == -1 && PyErr_Occurred()) {
        return NULL;
    }

    if (idx < 0 || idx >= PVN) {
        PyErr_Format(PyExc_ValueError,
            "PV index must be an integer number between 0 and %d.", PVN - 1);
        return NULL;
    }

    if (value == Py_None) {
        /* If pv is set to None - reset pv to prjini values: */
        self->x->pv[idx] = (idx > 0 && idx < 4) ? UNDEFINED : 0.0;
        self->x->flag = 0;
        if (self->owner) ((PyCelprm*)self->owner)->x->flag = 0;
        Py_RETURN_NONE;
    }

    if (PyFloat_Check(value) || PyLong_Check(value)) {
        data = PyFloat_AsDouble(value);
        if (data == -1.0 && PyErr_Occurred()) {
            return NULL;
        }

    } else if (PyUnicode_Check(value)) {
        flt_value = PyFloat_FromString(value);
        if (!flt_value) return NULL;
        data = PyFloat_AsDouble(flt_value);
        Py_DECREF(flt_value);
        if (data == -1.0 && PyErr_Occurred()) {
            return NULL;
        }

    } else {
        if (PyArray_Converter(value, &value_array_pyobj) == NPY_FAIL) {
            return NULL;
        }
        value_array = (PyArrayObject*) value_array_pyobj;

        size = PyArray_SIZE(value_array);
        if (size != 1) {
            Py_DECREF(value_array);
            PyErr_SetString(PyExc_ValueError,
                "PV value must be a scalar-like object or None.");
            return NULL;
        }

        scalar = PyArray_ToScalar(PyArray_DATA(value_array), value_array);
        Py_DECREF(value_array);
        if (!scalar) {
            Py_DECREF(scalar);
            PyErr_SetString(PyExc_TypeError, "Unable to convert value to scalar.");
        }

        PyArray_CastScalarToCtype(scalar, &data, dbl_descr);
        Py_DECREF(scalar);
        if (PyErr_Occurred()) {
            return NULL;
        }
    }

    data = (isnan(data)) ? UNDEFINED : data;

    if (!is_dbl_equal(self->x->pv[idx], data)) {
        self->x->flag = 0;
        if (self->owner) ((PyCelprm*)self->owner)->x->flag = 0;
    }
    self->x->pv[idx] = data;

    Py_RETURN_NONE;
}


static PyObject* PyPrjprm_get_bounds(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return get_int("bounds", self->x->bounds);
    }
}


static int PyPrjprm_set_bounds(PyPrjprm* self, PyObject* value, void* closure)
{
    if (is_prj_null(self) || is_readonly(self)) {
        return -1;
    } else if (value == Py_None) {
        self->x->bounds = 0;
        return 0;
    } else {
        return set_int("bounds", value, &self->x->bounds);
    }
}


static PyObject* PyPrjprm_get_w(PyPrjprm* self, void* closure)
{
    Py_ssize_t size = 10;
    int k;
    double *w;
    PyArrayObject* w_array;

    if (is_prj_null(self)) return NULL;

    w_array = (PyArrayObject*) PyArray_SimpleNew(1, &size, NPY_DOUBLE);
    if (w_array == NULL) return NULL;
    w = (double*) PyArray_DATA(w_array);

    for (k = 0; k < size; k++) {
        if (self->x->w[k] == UNDEFINED) {
            w[k] = (double) NPY_NAN;
        } else {
            w[k] = self->x->w[k];
        }
    }

    return (PyObject*) w_array;
}


static PyObject* PyPrjprm_get_name(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return get_string("name", self->x->name);
    }
}


static PyObject* PyPrjprm_get_category(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return get_int("category", self->x->category);
    }
}


static PyObject* PyPrjprm_get_pvrange(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return get_int("pvrange", self->x->pvrange);
    }
}


static PyObject* PyPrjprm_get_simplezen(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return PyBool_FromLong(self->x->simplezen);
    }
}


static PyObject* PyPrjprm_get_equiareal(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return PyBool_FromLong(self->x->equiareal);
    }
}


static PyObject* PyPrjprm_get_conformal(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return PyBool_FromLong(self->x->conformal);
    }
}


static PyObject* PyPrjprm_get_global_projection(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return PyBool_FromLong(self->x->global);
    }
}


static PyObject* PyPrjprm_get_divergent(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return PyBool_FromLong(self->x->divergent);
    }
}


static PyObject* PyPrjprm_get_x0(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return get_double("x0", self->x->x0);
    }
}


static PyObject* PyPrjprm_get_y0(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return get_double("y0", self->x->y0);
    }
}


static PyObject* PyPrjprm_get_m(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return get_int("m", self->x->m);
    }
}


static PyObject* PyPrjprm_get_n(PyPrjprm* self, void* closure)
{
    if (is_prj_null(self)) {
        return NULL;
    } else {
        return get_int("n", self->x->n);
    }
}


/***************************************************************************
 * PyPrjprm definition structures
 */

static PyGetSetDef PyPrjprm_getset[] = {
    {"r0", (getter)PyPrjprm_get_r0, (setter)PyPrjprm_set_r0, (char *)doc_prjprm_r0},
    {"phi0", (getter)PyPrjprm_get_phi0, (setter)PyPrjprm_set_phi0, (char *)doc_prjprm_phi0},
    {"theta0", (getter)PyPrjprm_get_theta0, (setter)PyPrjprm_set_theta0, (char *)doc_prjprm_theta0},
    {"pv", (getter)PyPrjprm_get_pv, (setter)PyPrjprm_set_pv, (char *)doc_prjprm_pv},
    {"w", (getter)PyPrjprm_get_w, NULL, (char *)doc_prjprm_w},
    {"name", (getter)PyPrjprm_get_name, NULL, (char *)doc_prjprm_name},
    {"code", (getter)PyPrjprm_get_code, (setter)PyPrjprm_set_code, (char *)doc_prjprm_code},
    {"bounds", (getter)PyPrjprm_get_bounds, (setter)PyPrjprm_set_bounds, (char *)doc_prjprm_bounds},
    {"category", (getter)PyPrjprm_get_category, NULL, (char *)doc_prjprm_category},
    {"pvrange", (getter)PyPrjprm_get_pvrange, NULL, (char *)doc_prjprm_pvrange},
    {"simplezen", (getter)PyPrjprm_get_simplezen, NULL, (char *)doc_prjprm_simplezen},
    {"equiareal", (getter)PyPrjprm_get_equiareal, NULL, (char *)doc_prjprm_equiareal},
    {"conformal", (getter)PyPrjprm_get_conformal, NULL, (char *)doc_prjprm_conformal},
    {"global_projection", (getter)PyPrjprm_get_global_projection, NULL, (char *)doc_prjprm_global_projection},
    {"divergent", (getter)PyPrjprm_get_divergent, NULL, (char *)doc_prjprm_divergent},
    {"x0", (getter)PyPrjprm_get_x0, NULL, (char *)doc_prjprm_x0},
    {"y0", (getter)PyPrjprm_get_y0, NULL, (char *)doc_prjprm_y0},
    {"m", (getter)PyPrjprm_get_m, NULL, (char *)doc_prjprm_m},
    {"n", (getter)PyPrjprm_get_n, NULL, (char *)doc_prjprm_n},
    {"_flag", (getter)PyPrjprm_get_flag, NULL, ""},
    {NULL}
};


static PyMethodDef PyPrjprm_methods[] = {
    {"set", (PyCFunction)PyPrjprm_set, METH_NOARGS, (char*)doc_prjprm_set},
    {"prjx2s", (PyCFunction)PyPrjprm_prjx2s, METH_VARARGS|METH_KEYWORDS, (char*)doc_prjprm_prjx2s},
    {"prjs2x", (PyCFunction)PyPrjprm_prjs2x, METH_VARARGS|METH_KEYWORDS, (char*)doc_prjprm_prjs2x},
    {"set_pvi", (PyCFunction)PyPrjprm_set_pvi, METH_VARARGS|METH_KEYWORDS, (char*)doc_prjprm_pvi},
    {"get_pvi", (PyCFunction)PyPrjprm_get_pvi, METH_VARARGS|METH_KEYWORDS, (char*)doc_prjprm_pvi},
    {"__copy__", (PyCFunction)PyPrjprm_copy, METH_NOARGS, ""},
    {"__deepcopy__", (PyCFunction)PyPrjprm_deepcopy, METH_O, ""},
    {NULL}
};


PyTypeObject PyPrjprmType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "astropy.wcs.Prjprm",         /*tp_name*/
    sizeof(PyPrjprm),             /*tp_basicsize*/
    0,                            /*tp_itemsize*/
    (destructor)PyPrjprm_dealloc, /*tp_dealloc*/
    0,                            /*tp_print*/
    0,                            /*tp_getattr*/
    0,                            /*tp_setattr*/
    0,                            /*tp_compare*/
    0,                            /*tp_repr*/
    0,                            /*tp_as_number*/
    0,                            /*tp_as_sequence*/
    0,                            /*tp_as_mapping*/
    0,                            /*tp_hash */
    0,                            /*tp_call*/
    (reprfunc)PyPrjprm___str__,   /*tp_str*/
    0,                            /*tp_getattro*/
    0,                            /*tp_setattro*/
    0,                            /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    doc_Prjprm,                   /* tp_doc */
    (traverseproc)PyPrjprm_traverse, /* tp_traverse */
    (inquiry)PyPrjprm_clear,      /* tp_clear */
    0,                            /* tp_richcompare */
    0,                            /* tp_weaklistoffset */
    0,                            /* tp_iter */
    0,                            /* tp_iternext */
    PyPrjprm_methods,             /* tp_methods */
    0,                            /* tp_members */
    PyPrjprm_getset,              /* tp_getset */
    0,                            /* tp_base */
    0,                            /* tp_dict */
    0,                            /* tp_descr_get */
    0,                            /* tp_descr_set */
    0,                            /* tp_dictoffset */
    0,                            /* tp_init */
    0,                            /* tp_alloc */
    PyPrjprm_new,                 /* tp_new */
};


int _setup_prjprm_type(PyObject* m)
{
    if (PyType_Ready(&PyPrjprmType) < 0) return -1;
    Py_INCREF(&PyPrjprmType);
    PyModule_AddObject(m, "Prjprm", (PyObject *)&PyPrjprmType);

    prj_errexc[0] = NULL;                         /* Success */
    prj_errexc[1] = &PyExc_MemoryError;           /* Null prjprm pointer passed */
    prj_errexc[2] = &WcsExc_InvalidPrjParameters; /* Invalid projection parameters */
    prj_errexc[3] = &WcsExc_InvalidCoordinate;    /* One or more of the (x,y) coordinates were invalid */
    prj_errexc[4] = &WcsExc_InvalidCoordinate;    /* One or more of the (lng,lat) coordinates were invalid */

    return 0;
}
