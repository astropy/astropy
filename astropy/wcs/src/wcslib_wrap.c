/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

#include "astropy_wcs/wcslib_wrap.h"
#include "astropy_wcs/wcslib_auxprm_wrap.h"
#include "astropy_wcs/wcslib_prjprm_wrap.h"
#include "astropy_wcs/wcslib_celprm_wrap.h"
#include "astropy_wcs/wcslib_tabprm_wrap.h"
#include "astropy_wcs/wcslib_wtbarr_wrap.h"
#include "astropy_wcs/wcslib_units_wrap.h"
#include "astropy_wcs/unit_list_proxy.h"
#include <structmember.h> /* from Python */

#include <stdlib.h> // calloc, malloc, free
#include <string.h> // memset, strlen, strncpy

#include <wcs.h>
#include <wcsfix.h>
#include <wcshdr.h>
#include <wcsmath.h>
#include <wcsprintf.h>
#include <wcsunits.h>
#include <cel.h>
#include <prj.h>
#include <tab.h>
#include <wtbarr.h>
#include <stdio.h>

#include "astropy_wcs/isnan.h"
#include "astropy_wcs/distortion.h"

/*
 It gets to be really tedious to type long docstrings in ANSI C syntax
 (since multi-line strings literals are not valid).  Therefore, the
 docstrings are written in doc/docstrings.py, which are then converted
 by setup.py into docstrings.h, which we include here.
*/
#include "astropy_wcs/docstrings.h"

/***************************************************************************
 * Helper functions                                                        *
 ***************************************************************************/

enum e_altlin {
  has_pc = 1,
  has_cd = 2,
  has_crota = 4
};

static int
is_valid_alt_key(
    const char* key) {

  if (key[1] != '\0' ||
      !(key[0] == ' ' ||
        (key[0] >= 'A' && key[0] <= 'Z'))) {
    PyErr_SetString(PyExc_ValueError, "key must be ' ' or 'A'-'Z'");
    return 0;
  }

  return 1;
}

static int
convert_rejections_to_warnings() {
  char buf[1024];
  const char *src;
  char *dst;
  int last_was_space;
  PyObject *wcs_module = NULL;
  PyObject *FITSFixedWarning = NULL;
  int status = -1;
  char delimiter;

#ifdef HAVE_WCSLIB_VERSION
  delimiter = ',';
#else
  delimiter = ':';
#endif

  if (wcsprintf_buf()[0] == 0) {
    return 0;
  }

  wcs_module = PyImport_ImportModule("astropy.wcs");
  if (wcs_module == NULL) {
    goto exit;
  }

  FITSFixedWarning = PyObject_GetAttrString(
      wcs_module, "FITSFixedWarning");
  if (FITSFixedWarning == NULL) {
    goto exit;
  }

  src = wcsprintf_buf();
  while (*src != 0) {
    dst = buf;

    /* Read the first line, removing any repeated spaces */
    last_was_space = 0;
    for (; *src != 0; ++src) {
      if (*src == ' ') {
        if (!last_was_space) {
          *(dst++) = *src;
          last_was_space = 1;
        }
      } else if (*src == '\n') {
        ++src;
        break;
      } else {
        *(dst++) = *src;
        last_was_space = 0;
      }
    }

    *(dst++) = '\n';

    /* For the second line, remove everything up to and including the
       first colon */
    for (; *src != 0; ++src) {
      if (*src == delimiter) {
        ++src;
        break;
      }
    }

    /* Read to the end of the second line, removing any repeated
       spaces */
    last_was_space = 1;
    for (; *src != 0; ++src) {
      if (*src == ' ') {
        if (!last_was_space) {
          *(dst++) = *src;
          last_was_space = 1;
        }
      } else if (*src == '\n') {
        ++src;
        break;
      } else {
        *(dst++) = *src;
        last_was_space = 0;
      }
    }

    /* NULL terminate the string */
    *dst = 0;

    /* Raise the warning.  Depending on the user's configuration, this
       may raise an exception, and PyErr_WarnEx returns -1. */
    if (PyErr_WarnEx(FITSFixedWarning, buf, 1)) {
      goto exit;
    }
  }

  status = 0;

 exit:

  Py_XDECREF(wcs_module);
  Py_XDECREF(FITSFixedWarning);

  return status;
}

PyObject* get_scaled_double_array_readonly(
    /*@unused@*/ const char* propname,
    double* value,
    double *scale,
    const npy_intp naxis,
    const npy_intp* dims,
    /*@shared@*/ PyObject* owner) {

  /*
  This function is equivalent to get_double_array_readonly except that values
  can be scaled before constructing the final array. The input parameters are as
  in get_double_array_readonly except for the addition of scale which should be
  an array of scalefactors with length dims[0] - the number of coordinate axes.
  For now this supports only 1D or 2D arrays, as we do not need support for
  higher dimensions.
  */

  PyArrayObject* array;

  array = (PyArrayObject*)PyArray_SimpleNew(naxis, dims, NPY_DOUBLE);
  if (array == NULL) {
    return NULL;
  }

  double *array_data = (double *)PyArray_DATA(array);

  if (naxis == 1) {
    for (npy_intp i = 0; i < dims[0]; ++i) {
      array_data[i] = value[i] / scale[i];
    }
  } else if (naxis == 2) {
    for (npy_intp i = 0; i < dims[0]; ++i) {
      for (npy_intp j = 0; j < dims[1]; ++j) {
        array_data[i * dims[1] + j] = value[i * dims[1] + j] / scale[i];
      }
    }
  } else {
    return NULL;
  }

  // Make the array read-only
  PyArray_CLEARFLAGS(array, NPY_ARRAY_WRITEABLE);

  return (PyObject*)array;

}

int cunit_empty(const char cunit[][72], int naxis) {
  /* This function returns 0 if all the cunit values are empty, and 1 otherwise */
  for (int i = 0; i < naxis; ++i) {
      if (strcmp(cunit[i], "") != 0) {
        return 0;
      }
  }
  return 1;
}

/***************************************************************************
 * wtbarr-related global variables and functions                           *
 ***************************************************************************/

static PyObject *get_wtbarr_data = NULL;


void _set_wtbarr_callback(PyObject* callback) {
  Py_XINCREF(callback);         /* Add a reference to new callback */
  Py_XDECREF(get_wtbarr_data);  /* Dispose of previous callback */
  get_wtbarr_data = callback;   /* Remember new callback */
}


int _update_wtbarr_from_hdulist(PyObject *hdulist, struct wtbarr *wtb) {
  PyArrayObject *arrayp=NULL;
  PyObject *result=NULL;
  int i, naxis, nelem, naxes[NPY_MAXDIMS];
  npy_intp *npy_naxes;
  npy_double *appayp_data;

  if (hdulist == NULL || hdulist == Py_None) {
    PyErr_SetString(PyExc_ValueError,
                    "HDUList is required to retrieve -TAB coordinates "
                    "and/or indices.");
    return 0;
  }

  if (wtb->ndim < 1) {
    PyErr_SetString(PyExc_ValueError, "Number of dimensions should be positive.");
    return 0;
  }

  result = PyObject_CallFunction(get_wtbarr_data, "(OsiiCsli)", hdulist,
      wtb->extnam, wtb->extver, wtb->extlev, wtb->kind, wtb->ttype, wtb->row,
      wtb->ndim);

  if (result == NULL) return 0;

  arrayp = (PyArrayObject *)PyArray_FromAny(result,
      PyArray_DescrFromType(NPY_DOUBLE), 0, 0, NPY_ARRAY_CARRAY, NULL);

  Py_DECREF(result);

  if (arrayp == NULL) {
    PyErr_SetString(PyExc_TypeError, "Unable to convert wtbarr callback "
                    "result to a numpy.ndarray.");
    return 0;
  }

  if (!PyArray_Check((PyObject*)arrayp)) {
    PyErr_SetString(PyExc_TypeError,
                    "wtbarr callback must return a numpy.ndarray type "
                    "coordinate or index array.");
    Py_DECREF(arrayp);
    return 0;
  }

  naxis = PyArray_NDIM(arrayp);

  if (naxis == 0) {
    PyErr_SetString(PyExc_ValueError, "-TAB coordinate or index arrays "
                    "cannot be 0-dimensional.");
    Py_DECREF(arrayp);
    return 0;
  }

  npy_naxes = PyArray_DIMS(arrayp);
  for (i = 0; i < naxis; i++) {
    naxes[i] = (int) npy_naxes[i];
  }

  if (naxis != wtb->ndim) {
    if (wtb->kind == 'c' && wtb->ndim == 2 && naxis == 1) {
      /* Allow TDIMn to be omitted for degenerate coordinate arrays. */
      naxis = 2;
      naxes[1] = 1;
    } else {
      PyErr_Format(PyExc_ValueError,
          "An array with an unexpected number of axes was "
          "received from the callback. Expected %d but got %d.",
          wtb->ndim, (int) naxis);
      Py_DECREF(arrayp);
      return 0;
    }
  }

  if (wtb->kind == 'c') {
    /* Coordinate array; calculate the array size. */
    nelem = naxes[naxis-1];
    for (i = 0; i < naxis-1; i++) {
      *(wtb->dimlen + i) = naxes[naxis-2-i];
      nelem *= naxes[i];
    }
  } else {
    /* Index vector; check length. */
    if ((nelem = naxes[naxis-1]) != *(wtb->dimlen)) {
      /* N.B. coordinate array precedes the index vectors. */
      PyErr_Format(PyExc_ValueError,
          "An index array with an unexpected number of dimensions was "
          "received from the callback. Expected %d but got %d.",
          *(wtb->dimlen), (int) nelem);
      Py_DECREF(arrayp);
      return 0;
    }
  }

  /* Allocate memory for the array. */
  if (!((*wtb->arrayp) = calloc((size_t)nelem, sizeof(double)))) {
    PyErr_SetString(PyExc_MemoryError, "Out of memory: can't allocate "
                    "coordinate or index array.");
    Py_DECREF(arrayp);
    return 0;
  }

  /* Read the array from the table. */
  appayp_data = (npy_double*)PyArray_DATA(arrayp);
  for (i = 0; i < nelem; i++) {
    (*wtb->arrayp)[i] = (double)appayp_data[i];
  }

  Py_DECREF(arrayp);
  return 1;
}


/***************************************************************************
 * Wcsprm methods
 */

int
Wcsprm_cset(Wcsprm* self, const int convert);

static INLINE void
note_change(Wcsprm* self) {
  self->x.flag = 0;
}

static void
Wcsprm_dealloc(
    Wcsprm* self) {

  wcsfree(&self->x);
  PyTypeObject *tp = Py_TYPE((PyObject*)self);
  freefunc free_func = PyType_GetSlot(tp, Py_tp_free);
  free_func((PyObject*)self);
  Py_DECREF(tp);
}

static Wcsprm*
Wcsprm_cnew(void) {
  Wcsprm* self;
  PyTypeObject* type = (PyTypeObject*)WcsprmType;
  allocfunc alloc_func = PyType_GetSlot(type, Py_tp_alloc);
  self = (Wcsprm*)alloc_func(type, 0);
  return self;
}

static PyObject *
Wcsprm_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  Wcsprm* self;
  allocfunc alloc_func = PyType_GetSlot(type, Py_tp_alloc);
  self = (Wcsprm*)alloc_func(type, 0);
  return (PyObject*)self;
}

static int
Wcsprm_init(
    Wcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  int            status;
  PyObject*      header_obj    = NULL;
  PyObject*      hdulist       = NULL;
  char *         header        = NULL;
  Py_ssize_t     header_length = 0;
  Py_ssize_t     nkeyrec       = 0;
  const char *   key           = " ";
  PyObject*      relax_obj     = NULL;
  int            relax         = 0;
  int            naxis         = -1;
  int            keysel        = -1;
  PyObject*      colsel        = Py_None;
  PyArrayObject* colsel_array  = NULL;
  int*           colsel_data  = NULL;
  int*           colsel_ints   = NULL;
  int            warnings      = 1;
  int            nreject       = 0;
  int            nwcs          = 0;
  struct wcsprm* wcs           = NULL;
  int            i, j;
  int            preserve_units = 0;
  const char*    keywords[]    = {"header", "key", "relax", "naxis", "keysel",
                                  "colsel", "warnings", "hdulist", "preserve_units", NULL};

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "|OsOiiOiOp:WCSBase.__init__",
          (char **)keywords, &header_obj, &key, &relax_obj, &naxis, &keysel,
          &colsel, &warnings, &hdulist, &preserve_units)) {
    return -1;
  }

  // The user can optionally choose to preserve the original units rather
  // than use the units WCSLIB converts to (e.g. arcsec -> deg)
  self->preserve_units = preserve_units;

  if (header_obj == NULL || header_obj == Py_None) {
    if (keysel > 0) {
      PyErr_SetString(
          PyExc_ValueError,
          "If no header is provided, keysel may not be provided either.");
      return -1;
    }

    if (colsel != Py_None) {
      PyErr_SetString(
          PyExc_ValueError,
          "If no header is provided, colsel may not be provided either.");
      return -1;
    }

    /* Default number of axes is 2 */
    if (naxis < 0) {
        naxis = 2;
    }

    if (naxis < 1 || naxis > 15) {
      PyErr_SetString(
          PyExc_ValueError,
          "naxis must be in range 1-15");
      return -1;
    }

    self->x.flag = -1;
    status = wcsini(1, naxis, &self->x);

    if (status != 0) {
      PyErr_SetString(
          PyExc_MemoryError,
          self->x.err->msg);
      return -1;
    }

    self->x.alt[0] = key[0];

    if (Wcsprm_cset(self, 0)) {
      return -1;
    }
    wcsprm_c2python(&self->x);

    return 0;
  } else { /* header != NULL */
    if (PyBytes_AsStringAndSize(header_obj, &header, &header_length)) {
      return -1;
    }

    if (relax_obj == Py_True) {
      relax = WCSHDR_all;
    } else if (relax_obj == NULL || relax_obj == Py_False) {
      relax = WCSHDR_none;
    } else {
      relax = (int)PyLong_AsLong(relax_obj);
      if (relax == -1) {
        PyErr_SetString(
            PyExc_ValueError,
            "relax must be True, False or an integer.");
        return -1;
      }
    }

    if (!is_valid_alt_key(key)) {
      return -1;
    }

    if (naxis >= 0) {
      PyErr_SetString(
          PyExc_ValueError,
          "naxis may not be provided if a header is provided.");
      return -1;
    }

    nkeyrec = header_length / 80;
    if (nkeyrec > 0x7fffffff) {
      PyErr_SetString(
          PyExc_MemoryError,
          "header is too long");
      return -1;
    }

    if (colsel != Py_None) {
      colsel_array = (PyArrayObject*) PyArray_ContiguousFromAny(
        colsel, NPY_INT, 1, 1);
      if (colsel_array == NULL) {
        return -1;
      }

      colsel_ints = malloc(sizeof(int) * (PyArray_DIM(colsel_array, 0) + 1));
      if (colsel_ints == NULL) {
        Py_DECREF(colsel_array);
        PyErr_SetString(
            PyExc_MemoryError,
            "Memory allocation error.");
        return -1;
      }

      colsel_ints[0] = (int)PyArray_DIM(colsel_array, 0);
      colsel_data = (int *)PyArray_DATA(colsel_array);
      for (i = 0; i < colsel_ints[0]; ++i) {
        colsel_ints[i+1] = colsel_data[i];
      }

      Py_DECREF(colsel_array);
    }

    wcsprintf_set(NULL);

    /* Call the header parser twice, the first time to get warnings
       out about "rejected" keywords (which we can then send to Python
       as warnings), and the second time to get a corrected wcsprm
       object. */

    if (keysel < 0) {
      status = wcspih(
          header,
          (int)nkeyrec,
          WCSHDR_reject,
          2,
          &nreject,
          &nwcs,
          &wcs);
    } else {
      status = wcsbth(
          header,
          (int)nkeyrec,
          WCSHDR_reject,
          2,
          keysel,
          colsel_ints,
          &nreject,
          &nwcs,
          &wcs);
    }

    if (status != 0) {
      free(colsel_ints);
      wcshdr_err_to_python_exc(status, wcs);
      return -1;
    }

    wcsvfree(&nwcs, &wcs);

    if (warnings && convert_rejections_to_warnings()) {
      free(colsel_ints);
      return -1;
    }

    if (keysel < 0) {
      status = wcspih(
          header,
          (int)nkeyrec,
          relax,
          0,
          &nreject,
          &nwcs,
          &wcs);
    } else {
      status = wcsbth(
          header,
          (int)nkeyrec,
          relax,
          0,
          keysel,
          colsel_ints,
          &nreject,
          &nwcs,
          &wcs);
    }

    free(colsel_ints);

    if (status != 0) {
      wcshdr_err_to_python_exc(status, wcs);
      return -1;
    }

    if (nwcs == 0) {
      wcsvfree(&nwcs, &wcs);
      PyErr_SetString(
          WcsExc_NoWcsKeywordsFound,
          "No WCS keywords found in the given header");
      return -1;
    }

    /* Find the desired WCS */
    for (i = 0; i < nwcs; ++i) {
      if (wcs[i].alt[0] == key[0]) {
        break;
      }
    }

    if (i >= nwcs) {
      wcsvfree(&nwcs, &wcs);
      PyErr_Format(
          PyExc_KeyError,
          "No WCS with key '%s' was found in the given header",
          key);
      return -1;
    }

    if (wcscopy(1, wcs + i, &self->x) != 0) {
      wcsvfree(&nwcs, &wcs);
      PyErr_SetString(
          PyExc_MemoryError,
          self->x.err->msg);
      return -1;
    }

    if (self->x.ntab) {
      wcstab(&self->x);
      for (j = 0; j < self->x.nwtb; j++) {
        if (!_update_wtbarr_from_hdulist(hdulist, &(self->x.wtb[j]))) {
          wcsfree(&self->x);
          return -1;
        }
      }
    }

    note_change(self);
    wcsprm_c2python(&self->x);
    wcsvfree(&nwcs, &wcs);
    return 0;
  }
}

/*@null@*/ static PyObject*
Wcsprm_bounds_check(
    Wcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  unsigned char pix2sky    = 1;
  unsigned char sky2pix    = 1;
  int           bounds     = 0;
  const char*   keywords[] = {"pix2world", "world2pix", NULL};

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "|bb:bounds_check", (char **)keywords,
          &pix2sky, &sky2pix)) {
    return NULL;
  }

  if (pix2sky) {
      bounds |= 2|4;
  }

  if (sky2pix) {
      bounds |= 1;
  }

  wcsprm_python2c(&self->x);
  wcsbchk(&self->x, bounds);

  Py_RETURN_NONE;
}


/*@null@*/ static PyObject*
Wcsprm_copy(
    Wcsprm* self) {

  Wcsprm*     copy = NULL;
  int           status;

  copy = Wcsprm_cnew();
  if (copy == NULL) {
    return NULL;
  }

  wcsini(0, self->x.naxis, &copy->x);

  wcsprm_python2c(&self->x);
  status = wcscopy(1, &self->x, &copy->x);
  wcsprm_c2python(&self->x);

  // Also copy over any information related to preserving units

  copy->preserve_units = self->preserve_units;

  if (self->original_cunit != NULL) {
    copy->original_cunit = malloc(copy->x.naxis * sizeof(*copy->original_cunit));
    if (self->original_cunit == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    for (int i = 0; i < copy->x.naxis; ++i) {
        strncpy(copy->original_cunit[i], self->original_cunit[i], 72);
    }
  }

  if (self->unit_scaling != NULL) {
      copy->unit_scaling = malloc(copy->x.naxis * sizeof(double));
      if (copy->unit_scaling == NULL) {
          PyErr_NoMemory();
          return NULL;
      }
    for (int i = 0; i < copy->x.naxis; ++i) {
        copy->unit_scaling[i] = self->unit_scaling[i];
    }
  }


  if (status == 0) {
    if (Wcsprm_cset(copy, 0)) {
      Py_XDECREF((PyObject*)copy);
      return NULL;
    }

    wcsprm_c2python(&copy->x);
    return (PyObject*)copy;
  } else {
    Py_XDECREF((PyObject*)copy);
    wcs_to_python_exc(&(self->x));
    return NULL;
  }
}

static Wcsprm* Wcsprm_copy_with_patched_units(Wcsprm* source) {

    // This function returns a copy of a Wcsprm with units patched
    // to match the original units (before WCSLIB changed them). The
    // returned object should not be used to do any kind of transformations
    // and is only for use in e.g. converting to a header, or printing
    // contents.

    int original_flag;
    Wcsprm* copy = (Wcsprm*)Wcsprm_copy(source);

    // We make sure wcsset is happy
    wcsset(&copy->x);

    // We preserve the value of the flag which indicates that wcsset
    // has been called and there have been no further modifications
    original_flag = copy->x.flag;

    // We now override the unit, cdelt and crval values
    for (int i = 0; i < copy->x.naxis; ++i) {
      strncpy(copy->x.cunit[i], source->original_cunit[i], 72);
      copy->x.cdelt[i] /= source->unit_scaling[i];
      copy->x.crval[i] /= source->unit_scaling[i];
    }

    // We restore the flag to trick WCSLIB into thinking it no longer
    // needs to change the units.
    copy->x.flag = original_flag;

    return copy;

}

PyObject*
Wcsprm_find_all_wcs(
    PyObject* __,
    PyObject* args,
    PyObject* kwds) {

  PyObject*      header_obj    = NULL;
  char *         header        = NULL;
  Py_ssize_t     header_length = 0;
  Py_ssize_t     nkeyrec       = 0;
  PyObject*      relax_obj     = NULL;
  int            relax         = 0;
  int            keysel        = 0;
  int            warnings      = 1;
  int            nreject       = 0;
  int            nwcs          = 0;
  struct wcsprm* wcs           = NULL;
  PyObject*      result        = NULL;
  Wcsprm*      subresult     = NULL;
  int            i             = 0;
  const char*    keywords[]    = {"header", "relax", "keysel", "warnings", NULL};
  int            status        = -1;

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "O|Oii:find_all_wcs",
          (char **)keywords, &header_obj, &relax_obj, &keysel, &warnings)) {
    return NULL;
  }

  if (PyBytes_AsStringAndSize(header_obj, &header, &header_length)) {
    return NULL;
  }

  nkeyrec = header_length / 80;
  if (nkeyrec > 0x7fffffff) {
    PyErr_SetString(
        PyExc_MemoryError,
        "header is too long");
    return NULL;
  }

  if (relax_obj == Py_True) {
    relax = WCSHDR_all;
  } else if (relax_obj == NULL || relax_obj == Py_False) {
    relax = WCSHDR_none;
  } else {
    relax = (int)PyLong_AsLong(relax_obj);
    if (relax == -1) {
      PyErr_SetString(
          PyExc_ValueError,
          "relax must be True, False or an integer.");
      return NULL;
    }
  }

  /* Call the header parser twice, the first time to get warnings
     out about "rejected" keywords (which we can then send to Python
     as warnings), and the second time to get a corrected wcsprm
     object. */

  Py_BEGIN_ALLOW_THREADS
  if (keysel < 0) {
    status = wcspih(
        header,
        (int)nkeyrec,
        WCSHDR_reject,
        2,
        &nreject,
        &nwcs,
        &wcs);
  } else {
    status = wcsbth(
        header,
        (int)nkeyrec,
        WCSHDR_reject,
        2,
        keysel,
        NULL,
        &nreject,
        &nwcs,
        &wcs);
  }
  Py_END_ALLOW_THREADS

  if (status != 0) {
    wcshdr_err_to_python_exc(status, wcs);
    return NULL;
  }

  wcsvfree(&nwcs, &wcs);

  if (warnings && convert_rejections_to_warnings()) {
    return NULL;
  }

  Py_BEGIN_ALLOW_THREADS
  if (keysel < 0) {
    status = wcspih(
        header,
        (int)nkeyrec,
        relax,
        0,
        &nreject,
        &nwcs,
        &wcs);
  } else {
    status = wcsbth(
        header,
        (int)nkeyrec,
        relax,
        0,
        keysel,
        NULL,
        &nreject,
        &nwcs,
        &wcs);
  }
  Py_END_ALLOW_THREADS

  if (status != 0) {
    wcshdr_err_to_python_exc(status, wcs);
    return NULL;
  }

  result = PyList_New(nwcs);
  if (result == NULL) {
    wcsvfree(&nwcs, &wcs);
    return NULL;
  }

  for (i = 0; i < nwcs; ++i) {
    subresult = Wcsprm_cnew();
    if (wcscopy(1, wcs + i, &subresult->x) != 0) {
      Py_DECREF(result);
      wcsvfree(&nwcs, &wcs);
      PyErr_SetString(
          PyExc_MemoryError,
          "Could not initialize wcsprm object");
      return NULL;
    }

    if (PyList_SetItem(result, i, (PyObject *)subresult) == -1) {
      Py_DECREF(result);
      wcsvfree(&nwcs, &wcs);
      return NULL;
    }

    subresult->x.flag = 0;
    wcsprm_c2python(&subresult->x);
  }

  wcsvfree(&nwcs, &wcs);
  return result;
}

static PyObject*
Wcsprm_cdfix(
    Wcsprm* self) {

  int status = 0;

  wcsprm_python2c(&self->x);
  status = cdfix(&self->x);
  wcsprm_c2python(&self->x);

  if (status == -1 || status == 0) {
    return PyLong_FromLong((long)status);
  } else {
    wcserr_fix_to_python_exc(self->x.err);
    return NULL;
  }
}

static PyObject*
Wcsprm_celfix(
    Wcsprm* self) {

  int status = 0;

  wcsprm_python2c(&self->x);
  status = celfix(&self->x);
  wcsprm_c2python(&self->x);

  if (status == -1 || status == 0) {
    return PyLong_FromLong((long)status);
  } else {
    wcserr_fix_to_python_exc(self->x.err);
    return NULL;
  }
}

static PyObject *
Wcsprm_compare(
    Wcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  int cmp = 0;
  Wcsprm *other;
  double tolerance = 0.0;
  int equal;
  int status;

  const char* keywords[] = {"other", "cmp", "tolerance", NULL};

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "O!|id:compare", (char **)keywords,
          (PyTypeObject*)WcsprmType, &other, &cmp, &tolerance)) {
    return NULL;
  }


  wcsprm_python2c(&self->x);
  wcsprm_python2c(&other->x);
  status = wcscompare(cmp, tolerance, &self->x, &other->x, &equal);
  wcsprm_c2python(&self->x);
  wcsprm_c2python(&other->x);

  if (status) {
    wcserr_fix_to_python_exc(self->x.err);
    return NULL;
  } else {
    if (equal) {
      Py_RETURN_TRUE;
    } else {
      Py_RETURN_FALSE;
    }
  }
}

/*@null@*/ static PyObject*
Wcsprm_cylfix(
    Wcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  PyObject*      naxis_obj   = NULL;
  PyArrayObject* naxis_array = NULL;
  int*           naxis       = NULL;
  int            status      = 0;
  const char*    keywords[]  = {"naxis", NULL};

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "|O:cylfix", (char **)keywords,
          &naxis_obj)) {
    return NULL;
  }

  if (naxis_obj != NULL && naxis_obj != Py_None) {
    naxis_array = (PyArrayObject*)PyArray_ContiguousFromAny(
        naxis_obj, NPY_INT, 1, 1);
    if (naxis_array == NULL) {
      return NULL;
    }
    if (PyArray_DIM(naxis_array, 0) != self->x.naxis) {
      PyErr_Format(
          PyExc_ValueError,
          "naxis must be same length as the number of axes of "
          "the Wcsprm object (%d).",
          self->x.naxis);
      Py_DECREF(naxis_array);
      return NULL;
    }
    naxis = (int*)PyArray_DATA(naxis_array);
  }

  wcsprm_python2c(&self->x);
  status = cylfix(naxis, &self->x);
  wcsprm_c2python(&self->x);

  Py_XDECREF((PyObject*)naxis_array);

  if (status == -1 || status == 0) {
    return PyLong_FromLong((long)status);
  } else {
    wcserr_fix_to_python_exc(self->x.err);
    return NULL;
  }
}

static PyObject*
Wcsprm_datfix(
    Wcsprm* self) {

  int status = 0;

  wcsprm_python2c(&self->x);
  status = datfix(&self->x);
  wcsprm_c2python(&self->x);

  if (status == -1 || status == 0) {
    return PyLong_FromLong((long)status);
  } else {
    wcserr_fix_to_python_exc(self->x.err);
    return NULL;
  }
}

int initialize_preserve_units(Wcsprm* self) {

  // If the user has requested to preserve units, we keep track of what CUNIT
  // was before and after fixing so that we can store the scaling factor if
  // needed. However, we don't do this yet if the units have not been set.
  if (self->preserve_units == 1 && !cunit_empty(self->x.cunit, self->x.naxis)) {

    if (self->original_cunit == NULL) {
      self->original_cunit = malloc(self->x.naxis * sizeof(*self->original_cunit));
      if (self->original_cunit == NULL) {
          PyErr_NoMemory();
          return -1;
      }
      for (int i = 0; i < self->x.naxis; ++i) {
          strncpy(self->original_cunit[i], self->x.cunit[i], 72);
      }
    }

  }

  return 0;

}

int check_unit_changes(Wcsprm* self) {

  double         scale, offset, power;
  int            status;

  // Check which units have changed, and store scaling if changed.
  if (self->original_cunit != NULL) {

    if (self->unit_scaling == NULL) {
      self->unit_scaling = malloc(self->x.naxis * sizeof(double));
      if (self->unit_scaling == NULL) {
          PyErr_NoMemory();
          return -1;
      }
    }

    for (int i = 0; i < self->x.naxis; ++i) {
      if (strcmp(self->original_cunit[i], self->x.cunit[i]) != 0) {
        status = wcsunits(self->original_cunit[i], self->x.cunit[i], &scale, &offset, &power);
        if (status == 0) {
          // We don't support offset and power currently because it should not be needed,
          // and would introduce additional complexity. However, if needed this can be
          // supported in future.
          if (offset != 0 || power != 1)
                PyErr_Format(
                  PyExc_ValueError,
                 "Preserving original units with non-trivial offset and power is not supported "
                );
          self->unit_scaling[i] = scale;
        }
      } else {
        self->unit_scaling[i] = 1.0;
      }

    }

    // We now check if all the scales are 1. If so, then we actually deallocate the array
    // again to not waste time doing unit conversions subsequently. It might seem wasteful
    // to allocate the array and then deallocate it, but the other option would be to
    // allocate a temporary array for the scales in this function and then only if all the
    // values are not 1 allocate self->unit_scaling, but then we'd be doing two allocations.

    int scaling_needed = 0;
    for (int i = 0; i < self->x.naxis; ++i) {
        if (self->unit_scaling[i] != 1.0) {
            scaling_needed = 1;
            break;
        }
    }
    if (!scaling_needed) {
      free(self->unit_scaling);
      self->unit_scaling = NULL;
    }

  }
  return 0;
}


/*@null@*/ static PyObject*
Wcsprm_fix(
    Wcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  const char*    translate_units = NULL;
  int            ctrl            = 0;
  PyObject*      naxis_obj       = NULL;
  PyArrayObject* naxis_array     = NULL;
  int*           naxis           = NULL;
  int            stat[NWCSFIX];
  struct wcserr  err[NWCSFIX];
  PyObject*      subresult;
  PyObject*      result;
  int            i               = 0;
  int            msg_index       = 0;
  const char*    message;

  struct message_map_entry {
    const char* name;
    const int index;
  };
  const struct message_map_entry message_map[NWCSFIX] = {
    {"cdfix", CDFIX},
    {"datfix", DATFIX},
#if (NWCSFIX > 6)
    {"obsfix", OBSFIX},
#endif
    {"unitfix", UNITFIX},
    {"celfix", CELFIX},
    {"spcfix", SPCFIX},
    {"cylfix", CYLFIX}
  };
  const char* keywords[] = {"translate_units", "naxis", NULL};

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "|sO:fix", (char **)keywords,
          &translate_units, &naxis_obj)) {
    return NULL;
  }

  if (translate_units != NULL) {
    if (parse_unsafe_unit_conversion_spec(translate_units, &ctrl)) {
      return NULL;
    }
  }

  if (naxis_obj != NULL && naxis_obj != Py_None) {
    naxis_array = (PyArrayObject*)PyArray_ContiguousFromAny(
        naxis_obj, NPY_INT, 1, 1);
    if (naxis_array == NULL) {
      return NULL;
    }
    if (PyArray_DIM(naxis_array, 0) != self->x.naxis) {
      PyErr_Format(
          PyExc_ValueError,
          "naxis must be same length as the number of axes of "
          "the Wcprm object (%d).",
          self->x.naxis);
      Py_DECREF(naxis_array);
      return NULL;
    }
    naxis = (int*)PyArray_DATA(naxis_array);
  }

  memset(err, 0, sizeof(struct wcserr) * NWCSFIX);

  initialize_preserve_units(self);

  wcsprm_python2c(&self->x);
  wcsfixi(ctrl, naxis, &self->x, stat, err);
  wcsprm_c2python(&self->x);

  check_unit_changes(self);

  /* We're done with this already, so deref now so we don't have to remember
     later */
  Py_XDECREF((PyObject*)naxis_array);

  result = PyDict_New();
  if (result == NULL) {
    return NULL;
  }

  for (i = 0; i < NWCSFIX; ++i) {
    msg_index = stat[message_map[i].index];
    message = err[message_map[i].index].msg;
    if (message == NULL || message[0] == 0) {
      if (msg_index == FIXERR_SUCCESS) {
        message = "Success";
      } else {
        message = "No change";
      }
    }
    subresult = PyUnicode_FromString(message);
    if (subresult == NULL ||
        PyDict_SetItemString(result, message_map[i].name, subresult)) {
      Py_XDECREF(subresult);
      Py_XDECREF(result);
      return NULL;
    }
    Py_XDECREF(subresult);
  }

  return result;
}

/*@null@*/ static PyObject*
Wcsprm_get_cdelt_func(
    Wcsprm* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  Py_ssize_t naxis = 0;
  PyObject *result;

  if (is_null(self->x.cdelt)) {
    return NULL;
  }

  if (Wcsprm_cset(self, 1)) {
    return NULL;
  }

  naxis = self->x.naxis;

  if (self->unit_scaling != NULL) {
    result = get_scaled_double_array_readonly("cdelt", self->x.cdelt, self->unit_scaling, 1, &naxis, (PyObject*)self);
  } else {
    result = get_double_array_readonly("cdelt", self->x.cdelt, 1, &naxis, (PyObject*)self);
  }

  return result;
}

/*@null@*/ static PyObject*
Wcsprm_get_pc_func(
    Wcsprm* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  npy_intp dims[2];

  if (is_null(self->x.pc)) {
    return NULL;
  }

  if (Wcsprm_cset(self, 1)) {
    return NULL;
  }

  dims[0] = self->x.naxis;
  dims[1] = self->x.naxis;

  return get_double_array_readonly("pc", self->x.pc, 2, dims, (PyObject*)self);
}

/*@null@*/ static PyObject*
Wcsprm_get_ps(
    Wcsprm* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  return get_pscards("ps", self->x.ps, self->x.nps);
}

/*@null@*/ static PyObject*
Wcsprm_get_pv(
    Wcsprm* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  return get_pvcards("pv", self->x.pv, self->x.npv);
}

static PyObject*
Wcsprm_has_cdi_ja(
    Wcsprm* self) {

  int result = 0;

  result = self->x.altlin & has_cd;

  return PyBool_FromLong(result);
}

static PyObject*
Wcsprm_has_crotaia(
    Wcsprm* self) {

  int result = 0;

  result = self->x.altlin & has_crota;

  return PyBool_FromLong(result);
}

static PyObject*
Wcsprm_has_pci_ja(
    Wcsprm* self) {

  int result = 0;

  result = (self->x.altlin == 0 || self->x.altlin & has_pc);

  return PyBool_FromLong(result);
}

static PyObject*
Wcsprm_is_unity(
    Wcsprm* self) {

  if (Wcsprm_cset(self, 1)) {
    return NULL;
  }

  return PyBool_FromLong(self->x.lin.unity);
}

/*@null@*/ static PyObject*
Wcsprm_mix(
    Wcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  int            mixpix     = 0;
  int            mixcel     = 0;
  double         vspan[2]   = {0, 0};
  double         vstep      = 0;
  int            viter      = 0;
  Py_ssize_t     naxis      = 0;
  PyObject*      world_obj  = NULL;
  PyObject*      pixcrd_obj = NULL;
  int            origin     = 1;
  PyArrayObject* world      = NULL;
  PyArrayObject* phi        = NULL;
  PyArrayObject* theta      = NULL;
  PyArrayObject* imgcrd     = NULL;
  PyArrayObject* pixcrd     = NULL;
  int            status     = -1;
  PyObject*      result     = NULL;
  const char*    keywords[] = {
    "mixpix", "mixcel", "vspan", "vstep", "viter", "world", "pixcrd", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(
        args, kwds, "ii(dd)diOOi:mix", (char **)keywords,
        &mixpix, &mixcel, &vspan[0], &vspan[1], &vstep, &viter, &world_obj,
        &pixcrd_obj, &origin)) {
    return NULL;
  }

  if (viter < 5 || viter > 10) {
    PyErr_SetString(
        PyExc_ValueError,
        "viter must be in the range 5 - 10");
    goto exit;
  }

  world = (PyArrayObject*)PyArray_ContiguousFromAny
    (world_obj, NPY_DOUBLE, 1, 1);
  if (world == NULL) {
    PyErr_SetString(
        PyExc_TypeError,
        "Argument 6 (world) must be a 1-dimensional numpy array");
    goto exit;
  }
  if ((int)PyArray_DIM(world, 0) != self->x.naxis) {
    PyErr_Format(
        PyExc_TypeError,
        "Argument 6 (world) must be the same length as the number "
        "of axes (%d)",
        self->x.naxis);
    goto exit;
  }

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny
    (pixcrd_obj, NPY_DOUBLE, 1, 1);
  if (pixcrd == NULL) {
    PyErr_SetString(
        PyExc_TypeError,
        "Argument 7 (pixcrd) must be a 1-dimensional numpy array");
    goto exit;
  }
  if ((int)PyArray_DIM(pixcrd, 0) != self->x.naxis) {
    PyErr_Format(
        PyExc_TypeError,
        "Argument 7 (pixcrd) must be the same length as the "
        "number of axes (%d)",
        self->x.naxis);
    goto exit;
  }

  if (mixpix < 1 || mixpix > self->x.naxis) {
    PyErr_SetString(
        PyExc_ValueError,
        "Argument 1 (mixpix) must specify a pixel coordinate "
        "axis number");
    goto exit;
  }

  if (mixcel < 1 || mixcel > 2) {
    PyErr_SetString(
        PyExc_ValueError,
        "Argument 2 (mixcel) must specify a celestial coordinate "
        "axis number (1 for latitude, 2 for longitude)");
    goto exit;
  }

  /* Now we allocate a bunch of numpy arrays to store the
   * results in.
   */
  naxis = (Py_ssize_t)self->x.naxis;
  phi = (PyArrayObject*)PyArray_SimpleNew
    (1, &naxis, NPY_DOUBLE);
  if (phi == NULL) {
    goto exit;
  }

  theta = (PyArrayObject*)PyArray_SimpleNew
    (1, &naxis, NPY_DOUBLE);
  if (theta == NULL) {
    goto exit;
  }

  imgcrd = (PyArrayObject*)PyArray_SimpleNew
    (1, &naxis, NPY_DOUBLE);
  if (imgcrd == NULL) {
    goto exit;
  }

  /* Force a call to wcsset via Wcsprm_cset before entering the parallel
   * region (see the matching note in Wcs_all_pix2world).  This ensures
   * wcs->flag == WCSSET so that wcsmix below does not invoke wcsset
   * itself, allowing the wcsprm_python2c / wcsprm_c2python round-trip
   * to be dropped. */
  if (Wcsprm_cset(self, 1)) {
    goto exit;
  }

  /* Convert pixel coordinates to 1-based. */
  Py_BEGIN_ALLOW_THREADS
  preoffset_array(pixcrd, origin);
  status = wcsmix(
      &self->x,
      mixpix,
      mixcel,
      vspan,
      vstep,
      viter,
      (double*)PyArray_DATA(world),
      (double*)PyArray_DATA(phi),
      (double*)PyArray_DATA(theta),
      (double*)PyArray_DATA(imgcrd),
      (double*)PyArray_DATA(pixcrd));
  unoffset_array(pixcrd, origin);
  unoffset_array(imgcrd, origin);
  Py_END_ALLOW_THREADS

  if (status == 0) {
    result = PyDict_New();
    if (result == NULL ||
        PyDict_SetItemString(result, "imgcrd", (PyObject*)imgcrd) ||
        PyDict_SetItemString(result, "phi", (PyObject*)phi) ||
        PyDict_SetItemString(result, "theta", (PyObject*)theta) ||
        PyDict_SetItemString(result, "world", (PyObject*)world)) {
      goto exit;
    }
  }

 exit:
  Py_XDECREF((PyObject*)world);
  Py_XDECREF((PyObject*)phi);
  Py_XDECREF((PyObject*)theta);
  Py_XDECREF((PyObject*)imgcrd);
  Py_XDECREF((PyObject*)pixcrd);

  if (status == 0) {
    return result;
  } else {
    Py_XDECREF(result);
    if (status == -1) {
      /* The error message has already been set */
      return NULL;
    } else {
      wcs_to_python_exc(&(self->x));
      return NULL;
    }
  }
}

/*@null@*/ static PyObject*
Wcsprm_p2s(
    Wcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  int            naxis      = 2;
  int            ncoord     = 0;
  int            nelem      = 0;
  PyObject*      pixcrd_obj = NULL;
  int            origin     = 1;
  PyArrayObject* pixcrd     = NULL;
  PyArrayObject* imgcrd     = NULL;
  PyArrayObject* phi        = NULL;
  PyArrayObject* theta      = NULL;
  PyArrayObject* world      = NULL;
  PyArrayObject* stat       = NULL;
  PyObject*      result     = NULL;
  int            status     = 0;
  const char*    keywords[] = {
    "pixcrd", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "Oi:p2s", (char **)keywords,
          &pixcrd_obj, &origin)) {
    return NULL;
  }

  naxis = self->x.naxis;

  pixcrd = (PyArrayObject*)PyArray_ContiguousFromAny
    (pixcrd_obj, NPY_DOUBLE, 2, 2);
  if (pixcrd == NULL) {
    return NULL;
  }

  ncoord = PyArray_DIM(pixcrd, 0);
  nelem = PyArray_DIM(pixcrd, 1);

  if (nelem < naxis) {
    PyErr_Format(
      PyExc_RuntimeError,
      "Input array must be 2-dimensional, where the second dimension >= %d",
      naxis);
    goto exit;
  }

  /* Now we allocate a bunch of numpy arrays to store the results in.
   */
  imgcrd = (PyArrayObject*)PyArray_SimpleNew(
      2, PyArray_DIMS(pixcrd), NPY_DOUBLE);
  if (imgcrd == NULL) {
    goto exit;
  }

  phi = (PyArrayObject*)PyArray_SimpleNew(
      1, PyArray_DIMS(pixcrd), NPY_DOUBLE);
  if (phi == NULL) {
    goto exit;
  }

  theta = (PyArrayObject*)PyArray_SimpleNew(
      1, PyArray_DIMS(pixcrd), NPY_DOUBLE);
  if (theta == NULL) {
    goto exit;
  }

  world = (PyArrayObject*)PyArray_SimpleNew(
      2, PyArray_DIMS(pixcrd), NPY_DOUBLE);
  if (world == NULL) {
    goto exit;
  }

  stat = (PyArrayObject*)PyArray_SimpleNew(
      1, PyArray_DIMS(pixcrd), NPY_INT);
  if (stat == NULL) {
    goto exit;
  }

  // Force a call to wcsset via Wcsprm_cset before entering the parallel
  // region (see the matching note in Wcs_all_pix2world).  This ensures
  // wcs->flag == WCSSET so that wcsp2s below does not invoke wcsset
  // itself -- with wcsset moved out of the parallel region, wcsp2s is
  // read-only on the wcsprm struct, and the wcsprm_python2c /
  // wcsprm_c2python round-trip can be dropped.
  if (Wcsprm_cset(self, 1)) {
    return NULL;
  }

  /* Make the call.
   *
   * After wcsset has run (now hoisted into Wcsprm_cset above), wcsp2s is
   * read-only on the wcsprm struct: it consumes the precomputed sub-structs
   * wcs->lin / wcs->cel / wcs->spc plus wcs->crval[i], and never re-reads
   * the raw arrays (cd, cdelt, crpix, crota, obsgeo, mjdobs, ...) that the
   * wcsprm_python2c / wcsprm_c2python pair was rewriting NaN <-> UNDEFINED
   * in place.  The round-trip can therefore be dropped here without
   * changing the transform's output, and concurrent threads no longer
   * race on those raw arrays.
   *
   * We must likewise avoid the historical in-place 0->1-based offset of the
   * caller's pixcrd array: that array may be shared between threads, and the
   * transient mutation corrupts other threads' results (GH-16244).  When an
   * offset is needed (origin == 0) we copy the input in bounded-size blocks
   * into a scratch buffer and offset there, so the extra memory is O(block)
   * rather than O(ncoord).  When origin == 1 no offset is needed, pixscratch
   * stays NULL and pixcrd is read directly in a single pass.
   */
  double delta = (origin == 1) ? 0.0 : (double)(1 - origin);

  /* Size the scratch tile at NPY_BUFSIZE bytes, numpy's default buffer size,
   * so it stays L1-resident (also friendly to threads sharing a cache).  We
   * pick the number of coordinates per block so that block * nelem doubles fit
   * in that budget, which keeps the scratch footprint independent of nelem.
   * Throughput is flat over a wide range, so the exact value is not critical. */
  const npy_intp block_doubles = NPY_BUFSIZE / (npy_intp)sizeof(double);
  npy_intp block = ncoord;
  double* pixscratch = NULL;
  if (delta != 0.0 && ncoord > 0) {
    block = block_doubles / nelem;
    if (block < 1) {
      block = 1;
    }
    if (block > ncoord) {
      block = ncoord;
    }
    pixscratch = (double*)malloc((size_t)block * nelem * sizeof(double));
    if (pixscratch == NULL) {
      PyErr_NoMemory();
      status = -1;
      goto exit;
    }
  }

  const double* pix_data = (const double*)PyArray_DATA(pixcrd);
  double* img_data   = (double*)PyArray_DATA(imgcrd);
  double* phi_data   = (double*)PyArray_DATA(phi);
  double* theta_data = (double*)PyArray_DATA(theta);
  double* world_data = (double*)PyArray_DATA(world);
  int*    stat_data  = (int*)PyArray_DATA(stat);
  int     any_invalid = 0;

  Py_BEGIN_ALLOW_THREADS
  for (npy_intp c = 0; c < ncoord; c += block) {
    int n = (int)((ncoord - c < block) ? (ncoord - c) : block);
    const double* src = pix_data + c * nelem;
    if (pixscratch != NULL) {
      /* fused copy + offset: one read of the input, one write of scratch */
      for (npy_intp j = 0; j < (npy_intp)n * nelem; ++j) {
        pixscratch[j] = src[j] + delta;
      }
      src = pixscratch;
    }
    status = wcsp2s(
        &self->x, n, nelem, src,
        img_data + c * nelem,
        phi_data + c,
        theta_data + c,
        world_data + c * nelem,
        stat_data + c);
    if (status == 8) {
      any_invalid = 1;
      set_invalid_to_nan(n, nelem, img_data + c * nelem, stat_data + c);
      set_invalid_to_nan(n, 1, phi_data + c, stat_data + c);
      set_invalid_to_nan(n, 1, theta_data + c, stat_data + c);
      set_invalid_to_nan(n, nelem, world_data + c * nelem, stat_data + c);
      status = 0;
    } else if (status != 0) {
      break;
    }
  }
  if (status == 0 && any_invalid) {
    status = 8;
  }
  /* unoffset only the freshly-allocated imgcrd output; the caller's pixcrd is
   * never modified. */
  unoffset_array(imgcrd, origin);
  Py_END_ALLOW_THREADS

  if (pixscratch != NULL) {
    free(pixscratch);
  }

  if (status == 0 || status == 8) {
    // Since the conversion succeeded, if user has requested to preserve units,
    // we convert the world coordinates to the original units
    if (self->unit_scaling != NULL) {
      double *world_data = (double *)PyArray_DATA(world);
      for (npy_intp i = 0; i < nelem; ++i) {
        for (npy_intp j = 0; j < ncoord; ++j) {
            world_data[j * nelem + i] /= self->unit_scaling[i];
        }
      }
    }
    result = PyDict_New();
    if (result == NULL ||
        PyDict_SetItemString(result, "imgcrd", (PyObject*)imgcrd) ||
        PyDict_SetItemString(result, "phi", (PyObject*)phi) ||
        PyDict_SetItemString(result, "theta", (PyObject*)theta) ||
        PyDict_SetItemString(result, "world", (PyObject*)world) ||
        PyDict_SetItemString(result, "stat", (PyObject*)stat)) {
      goto exit;
    }
  }

 exit:
  Py_XDECREF((PyObject*)pixcrd);
  Py_XDECREF((PyObject*)imgcrd);
  Py_XDECREF((PyObject*)phi);
  Py_XDECREF((PyObject*)theta);
  Py_XDECREF((PyObject*)world);
  Py_XDECREF((PyObject*)stat);

  if (status == 0 || status == 8) {
    return result;
  } else {
    Py_XDECREF(result);
    if (status == -1) {
      /* Exception already set */
      return NULL;
    } else {
      wcs_to_python_exc(&(self->x));
      return NULL;
    }
  }
}

/*@null@*/ static PyObject*
Wcsprm_s2p(
    Wcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  int            naxis     = 2;
  int            ncoord    = 0;
  int            nelem     = 0;
  PyObject*      world_obj = NULL;
  int            origin    = 1;
  PyArrayObject* world     = NULL;
  PyArrayObject* phi       = NULL;
  PyArrayObject* theta     = NULL;
  PyArrayObject* imgcrd    = NULL;
  PyArrayObject* pixcrd    = NULL;
  PyArrayObject* stat      = NULL;
  PyObject*      result    = NULL;
  int            status    = -1;
  double*        world_copy = NULL;
  const char*    keywords[] = {
    "world", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "Oi:s2p", (char **)keywords,
          &world_obj, &origin)) {
    return NULL;
  }

  naxis = self->x.naxis;

  world = (PyArrayObject*)PyArray_ContiguousFromAny(
      world_obj, NPY_DOUBLE, 2, 2);
  if (world == NULL) {
    return NULL;
  }

  ncoord = (int)PyArray_DIM(world, 0);
  nelem = (int)PyArray_DIM(world, 1);

  if (nelem < naxis) {
    PyErr_Format(
      PyExc_RuntimeError,
      "Input array must be 2-dimensional, where the second dimension >= %d",
      naxis);
    goto exit;
  }

  // Force a call to wcsset via Wcsprm_cset before entering the parallel
  // region (see the matching note in Wcs_all_pix2world).  This ensures
  // wcs->flag == WCSSET so that wcss2p below does not invoke wcsset
  // itself, allowing the wcsprm_python2c / wcsprm_c2python round-trip to
  // be dropped.

  if (Wcsprm_cset(self, 1)) {
    return NULL;
  }

  // Convert world units to native WCSLIB units if needed. Here we use a check that
  // original_cunit is allocated, because the above line does not guarantee that
  // original_cunit will be allocated (for instance if the units are not set)

  double *world_data = (double *)PyArray_DATA(world);

  if (self->unit_scaling != NULL) {

    world_copy = malloc(sizeof(double) * nelem * ncoord);
    if (!world_copy) {
      return PyErr_NoMemory();
      return NULL;
    }

    for (npy_intp i = 0; i < nelem; ++i) {
      for (npy_intp j = 0; j < ncoord; ++j) {
          world_copy[j * nelem + i] = world_data[j * nelem + i] * self->unit_scaling[i];
      }
    }

    world_data = world_copy;

  }


  /* Now we allocate a bunch of numpy arrays to store the
   * results in.
   */
  phi = (PyArrayObject*)PyArray_SimpleNew(
      1, PyArray_DIMS(world), NPY_DOUBLE);
  if (phi == NULL) {
    goto exit;
  }

  theta = (PyArrayObject*)PyArray_SimpleNew(
      1, PyArray_DIMS(world), NPY_DOUBLE);
  if (phi == NULL) {
    goto exit;
  }

  imgcrd = (PyArrayObject*)PyArray_SimpleNew(
      2, PyArray_DIMS(world), NPY_DOUBLE);
  if (theta == NULL) {
    goto exit;
  }

  pixcrd = (PyArrayObject*)PyArray_SimpleNew(
      2, PyArray_DIMS(world), NPY_DOUBLE);
  if (pixcrd == NULL) {
    goto exit;
  }

  stat = (PyArrayObject*)PyArray_SimpleNew(
      1, PyArray_DIMS(world), NPY_INT);
  if (stat == NULL) {
    goto exit;
  }

  /* Make the call.  See the comment in Wcsprm_p2s for the rationale --
   * wcss2p is read-only on the wcsprm struct after wcsset has run, so
   * the wcsprm_python2c / wcsprm_c2python round-trip can be dropped and
   * concurrent transforms no longer race on the raw arrays
   *. */
  Py_BEGIN_ALLOW_THREADS
  /* preoffset_array(world, origin); */
  status = wcss2p(
      &self->x,
      ncoord,
      nelem,
      world_data,
      (double*)PyArray_DATA(phi),
      (double*)PyArray_DATA(theta),
      (double*)PyArray_DATA(imgcrd),
      (double*)PyArray_DATA(pixcrd),
      (int*)PyArray_DATA(stat));
  /* unoffset_array(world, origin); */
  unoffset_array(pixcrd, origin);
  unoffset_array(imgcrd, origin);
  if (status == 9) {
    set_invalid_to_nan(
        ncoord, 1, (double*)PyArray_DATA(phi), (int*)PyArray_DATA(stat));
    set_invalid_to_nan(
        ncoord, 1, (double*)PyArray_DATA(theta), (int*)PyArray_DATA(stat));
    set_invalid_to_nan(
        ncoord, nelem, (double*)PyArray_DATA(imgcrd), (int*)PyArray_DATA(stat));
    set_invalid_to_nan(
        ncoord, nelem, (double*)PyArray_DATA(pixcrd), (int*)PyArray_DATA(stat));
  }
  Py_END_ALLOW_THREADS

  if (status == 0 || status == 9) {
    result = PyDict_New();
    if (result == NULL ||
        PyDict_SetItemString(result, "phi", (PyObject*)phi) ||
        PyDict_SetItemString(result, "theta", (PyObject*)theta) ||
        PyDict_SetItemString(result, "imgcrd", (PyObject*)imgcrd) ||
        PyDict_SetItemString(result, "pixcrd", (PyObject*)pixcrd) ||
        PyDict_SetItemString(result, "stat", (PyObject*)stat)) {
      goto exit;
    }
  }

 exit:

  if (world_copy!=NULL) {
    free(world_copy);
  }

  Py_XDECREF((PyObject*)pixcrd);
  Py_XDECREF((PyObject*)imgcrd);
  Py_XDECREF((PyObject*)phi);
  Py_XDECREF((PyObject*)theta);
  Py_XDECREF((PyObject*)world);
  Py_XDECREF((PyObject*)stat);

  if (status == 0 || status == 9) {
    return result;
  } else {
    Py_XDECREF(result);
    if (status == -1) {
      /* Exception already set */
      return NULL;
    } else {
      wcs_to_python_exc(&(self->x));
      return NULL;
    }
  }
}

int
Wcsprm_cset(
    Wcsprm* self,
    const int convert) {

  int status = 0;

  // We want to avoid calling wcsset whenever possible as it is not thread-safe. We use wcsenq
  // to see if the checksum of the wcsprm elements has changed since wcsset was last called.
  if (wcsenq(&self->x, WCSENQ_CHK)) {
    return 0;
  }

#ifdef Py_GIL_DISABLED
  // On a free-threaded build the GIL no longer serialises concurrent callers
  // of Wcsprm_cset, so the wcsenq guard alone cannot prevent two threads from
  // running wcsset simultaneously on the same struct (which is not thread
  // safe).  A single process-wide PyMutex around the wcsset path is enough:
  // the fast path (wcsenq-succeeds) above is lock-free, so contention only
  // happens the first time a WCS is set or after a user mutation, which is
  // rare.  On a GIL build this macro is undefined and the lock is compiled
  // out -- the GIL itself serialises us since Wcsprm_cset never releases it.
  static PyMutex wcsset_mutex;
  PyMutex_Lock(&wcsset_mutex);
  // Double-checked locking: re-check under the lock so that if two threads
  // both missed the fast path above, only the first one actually runs
  // wcsset and the second returns immediately.
  if (wcsenq(&self->x, WCSENQ_CHK)) {
    PyMutex_Unlock(&wcsset_mutex);
    return 0;
  }
#endif

  initialize_preserve_units(self);

  if (convert) wcsprm_python2c(&self->x);
  status = wcsset(&self->x);
  if (convert) wcsprm_c2python(&self->x);

  check_unit_changes(self);

#ifdef Py_GIL_DISABLED
  PyMutex_Unlock(&wcsset_mutex);
#endif

  if (status == 0) {
    return 0;
  } else {
    wcs_to_python_exc(&(self->x));
    return 1;
  }
}

/*@null@*/ static PyObject*
Wcsprm_set(
    Wcsprm* self) {

  if (Wcsprm_cset(self, 1)) {
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

/*@null@*/ static PyObject*
Wcsprm_set_ps(
    Wcsprm* self,
    PyObject* arg,
    /*@unused@*/ PyObject* kwds) {

  if (is_null(self->x.ps)) {
    return NULL;
  }

  if (set_pscards("ps", arg, &self->x.ps, &self->x.nps, &self->x.npsmax)) {
    self->x.m_ps = self->x.ps;
    return NULL;
  }
  self->x.m_ps = self->x.ps;

  note_change(self);

  Py_INCREF(Py_None);
  return Py_None;
}

/*@null@*/ static PyObject*
Wcsprm_set_pv(
    Wcsprm* self,
    PyObject* arg,
    /*@unused@*/ PyObject* kwds) {

  if (is_null(self->x.pv)) {
    return NULL;
  } else if (set_pvcards("pv", arg, &self->x.pv, &self->x.npv, &self->x.npvmax)) {
    return NULL;
  } else {
    self->x.m_pv = self->x.pv;
    note_change(self);
    Py_INCREF(Py_None);
    return Py_None;
  }
}

/* TODO: This is convenient for debugging for now -- but it's not very
 * Pythonic.  It should probably be hooked into __str__ or something.
 */
/*@null@*/ static PyObject*
Wcsprm_print_contents(
    Wcsprm* self) {

  /* This is not thread-safe, but since we're holding onto the GIL,
     we can assume we won't have thread conflicts */
  wcsprintf_set(NULL);

  wcsprm_python2c(&self->x);
  if (Wcsprm_cset(self, 0)) {
    wcsprm_c2python(&self->x);
    return NULL;
  }

if (self->unit_scaling != NULL) {
    Wcsprm* copy = Wcsprm_copy_with_patched_units(self);
    wcsprt(&copy->x);
    Wcsprm_dealloc(copy);
  } else {
    wcsprt(&self->x);
    wcsprm_c2python(&self->x);
  }

  printf("%s", wcsprintf_buf());
  fflush(stdout);

  Py_INCREF(Py_None);
  return Py_None;
}

/*@null@*/ static PyObject*
Wcsprm_spcfix(
    Wcsprm* self) {

  int status = 0;

  wcsprm_python2c(&self->x);
  status = spcfix(&self->x);
  wcsprm_c2python(&self->x);

  if (status == -1 || status == 0) {
    return PyLong_FromLong((long)status);
  } else {
    wcserr_fix_to_python_exc(self->x.err);
    return NULL;
  }
}

/*@null@*/ static PyObject*
Wcsprm_sptr(
    Wcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  int         i          = -1;
  const char* py_ctype   = NULL;
  char        ctype[9];
  int         status     = 0;
  const char* keywords[] = {"ctype", "i", NULL};

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "s|i:sptr", (char **)keywords,
          &py_ctype, &i)) {
    return NULL;
  }

  if (strlen(py_ctype) > 8) {
    PyErr_SetString(
        PyExc_ValueError,
        "ctype string has more than 8 characters.");
  }

  strncpy(ctype, py_ctype, 9);

  wcsprm_python2c(&self->x);
  status = wcssptr(&self->x, &i, ctype);
  wcsprm_c2python(&self->x);

  if (status == 0) {
    Py_INCREF(Py_None);
    return Py_None;
  } else {
    wcs_to_python_exc(&(self->x));
    return NULL;
  }
}

/*@null@*/ static PyObject*
Wcsprm___str__(
    Wcsprm* self) {

  /* This is not thread-safe, but since we're holding onto the GIL,
     we can assume we won't have thread conflicts */
  wcsprintf_set(NULL);

  wcsprm_python2c(&self->x);
  if (Wcsprm_cset(self, 0)) {
    wcsprm_c2python(&self->x);
    return NULL;
  }

  if (self->unit_scaling != NULL) {
    Wcsprm* copy = Wcsprm_copy_with_patched_units(self);
    wcsprt(&copy->x);
    Wcsprm_dealloc(copy);
  } else {
    wcsprt(&self->x);
    wcsprm_c2python(&self->x);
  }

  return PyUnicode_FromString(wcsprintf_buf());
}

PyObject *Wcsprm_richcompare(PyObject *a, PyObject *b, int op) {
  int equal;
  int status;

  struct wcsprm *ax;
  struct wcsprm *bx;

  if ((op == Py_EQ || op == Py_NE) &&
      PyObject_TypeCheck(b, (PyTypeObject*)WcsprmType)) {
    ax = &((Wcsprm *)a)->x;
    bx = &((Wcsprm *)b)->x;

    wcsprm_python2c(ax);
    wcsprm_python2c(bx);
    status = wcscompare(
        WCSCOMPARE_ANCILLARY, 0.0,
        ax, bx, &equal);
    wcsprm_c2python(ax);
    wcsprm_c2python(bx);

    if (status == 0) {
      if (op == Py_NE) {
        equal = !equal;
      }
      if (equal) {
        Py_RETURN_TRUE;
      } else {
        Py_RETURN_FALSE;
      }
    } else {
      wcs_to_python_exc(&(((Wcsprm *)a)->x));
      return NULL;
    }
  }

  Py_INCREF(Py_NotImplemented);
  return Py_NotImplemented;
}

/*@null@*/ static PyObject*
Wcsprm_sub(
    Wcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  int        i            = -1;
  Py_ssize_t tmp          = 0;
  PyObject*  py_axes      = NULL;
  Wcsprm*  py_dest_wcs  = NULL;
  PyObject*  element      = NULL;
  PyObject*  element_utf8 = NULL;
  char*      element_str  = NULL;
  int        element_val  = 0;
  int        nsub         = 0;
  int*       axes         = NULL;
  int        status       = -1;
  int        wcslib_ver[] = {0, 0, 0};
  const char* keywords[]  = {"axes", NULL};

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "|O:sub", (char **)keywords,
          &py_axes)) {
    goto exit;
  }

  wcslib_version(wcslib_ver);

  if (py_axes == NULL || py_axes == Py_None) {
    /* leave all variables as is */
  } else if (PyList_Check(py_axes) || PyTuple_Check(py_axes)) {
    tmp = PySequence_Size(py_axes);
    if (tmp == -1) {
      goto exit;
    }
    nsub = (int)tmp;

    axes = malloc(nsub * sizeof(int) * 2);
    if (axes == NULL) {
      PyErr_SetString(PyExc_MemoryError, "Out of memory");
      goto exit;
    }

    for (i = 0; i < nsub; ++i) {
      element = PySequence_GetItem(py_axes, i);
      if (element == NULL) {
        goto exit;
      }

      if (PyUnicode_Check(element) || PyBytes_Check(element)) {
        if (PyUnicode_Check(element)) {
          element_utf8 = PyUnicode_AsUTF8String(element);
          if (element_utf8 == NULL) {
            goto exit;
          }

          element_str = PyBytes_AsString(element_utf8);
        } else if (PyBytes_Check(element)) {
          element_str = PyBytes_AsString(element);
        }

        if (strncmp(element_str, "longitude", 10) == 0) {
          element_val = WCSSUB_LONGITUDE;
        } else if (strncmp(element_str, "latitude", 9) == 0) {
          element_val = WCSSUB_LATITUDE;
        } else if (strncmp(element_str, "cubeface", 9) == 0) {
          element_val = WCSSUB_CUBEFACE;
        } else if (strncmp(element_str, "spectral", 9) == 0) {
          element_val = WCSSUB_SPECTRAL;
        } else if (strncmp(element_str, "stokes", 7) == 0) {
          element_val = WCSSUB_STOKES;
#if defined(WCSSUB_TIME)
        } else if ((wcslib_ver[0] > 7 || (wcslib_ver[0] == 7 && wcslib_ver[1] >= 8)) &&
                   strncmp(element_str, "temporal", 9) == 0) {
          element_val = WCSSUB_TIME;
#else
        } else if (strncmp(element_str, "temporal", 9) == 0) {
          PyErr_Format(
            PyExc_NotImplementedError,
            "Support for 'temporal' axis requires WCSLIB version 7.8 or greater "\
            "but linked WCSLIB version is %s", wcslib_version(NULL)
          );
          goto exit;
#endif
        } else if (strncmp(element_str, "celestial", 10) == 0) {
          element_val = WCSSUB_CELESTIAL;
        } else {
          PyErr_SetString(
            PyExc_ValueError,
#if defined(WCSSUB_TIME)
            "string values for axis sequence must be one of 'latitude', 'longitude', 'cubeface', 'spectral', 'stokes', 'temporal', or 'celestial'"
#else
            "string values for axis sequence must be one of 'latitude', 'longitude', 'cubeface', 'spectral', 'stokes', or 'celestial'"
#endif
            );
          goto exit;
        }
        Py_CLEAR(element_utf8);
      } else if (PyLong_Check(element)) {
        tmp = (Py_ssize_t)PyLong_AsSsize_t(element);
        if (tmp == -1 && PyErr_Occurred()) {
          goto exit;
        }
        element_val = (int)tmp;
      } else {
        PyErr_SetString(
          PyExc_TypeError,
          "axes sequence must contain either strings or ints");
        goto exit;
      }

      axes[i] = element_val;

      Py_CLEAR(element);
    }
  } else if (PyLong_Check(py_axes)) {
    tmp = (Py_ssize_t)PyLong_AsSsize_t(py_axes);
    if (tmp == -1 && PyErr_Occurred()) {
      goto exit;
    }
    nsub = (int)tmp;

    if (nsub < 0 || nsub > self->x.naxis) {
      PyErr_Format(
        PyExc_ValueError,
        "If axes is an int, it must be in the range 0-self.naxis (%d)",
        self->x.naxis);
      goto exit;
    }
  } else {
    PyErr_SetString(
      PyExc_TypeError,
      "axes must None, a sequence or an integer");
    goto exit;
  }

  py_dest_wcs = (Wcsprm*)Wcsprm_cnew();
  py_dest_wcs->x.flag = -1;
  status = wcsini(0, nsub, &py_dest_wcs->x);
  if (status != 0) {
    goto exit;
  }

  wcsprm_python2c(&self->x);
  status = wcssub(1, &self->x, &nsub, axes, &py_dest_wcs->x);
  wcsprm_c2python(&self->x);

  // At this point, we need to copy over any relevant preserve_units
  // information, but we need to be careful since axes may have been removed or
  // re-ordered

  py_dest_wcs->preserve_units = self->preserve_units;

  // We check whether original_cunit is allocated to know whether we need to copy
  // anything over - if it is not allocated then nothing needs to be done even
  // if preserve_units is set.

  if (self->original_cunit != NULL) {

  initialize_preserve_units(py_dest_wcs);

  // Now override original_unit with the actual original units from the original WCS

  for (int i = 0; i < py_dest_wcs->x.naxis; ++i) {

    // The axis variable has now been modified by wcssub to give the mapping
    // from new axis to original axis, with the axis numbers being 1-based.
    // If axis is zero, this means that the axis is a new one and we should
    // not do anything.

    if (axes[i] > 0) {
          strncpy(py_dest_wcs->original_cunit[i], self->original_cunit[axes[i] - 1], 72);
    }

  }

    check_unit_changes(py_dest_wcs);

}
  if (Wcsprm_cset(py_dest_wcs, 0)) {
    status = -1;
    goto exit;
  }
  wcsprm_c2python(&py_dest_wcs->x);

  if (status != 0) {
    goto exit;
  }

 exit:
  free(axes);
  Py_XDECREF(element);
  Py_XDECREF(element_utf8);

  if (status == 0) {
    return (PyObject*)py_dest_wcs;
  } else if (status == -1) {
    Py_XDECREF((PyObject*)py_dest_wcs);
    /* Exception already set */
    return NULL;
  } else {
    wcs_to_python_exc(&(py_dest_wcs->x));
    Py_XDECREF((PyObject*)py_dest_wcs);
    return NULL;
  }
}

/*@null@*/ static PyObject*
Wcsprm_to_header(
    Wcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  PyObject* relax_obj    = NULL;
  int       relax        = 0;
  int       nkeyrec      = 0;
  char*     header       = NULL;
  int       status       = -1;
  PyObject* result       = NULL;
  const char* keywords[] = {"relax", NULL};

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "|O:to_header",
          (char **)keywords, &relax_obj)) {
    goto exit;
  }

  if (relax_obj == Py_True) {
    relax = WCSHDO_all;
  } else if (relax_obj == NULL || relax_obj == Py_False) {
    relax = WCSHDO_safe;
  } else {
    relax = (int)PyLong_AsLong(relax_obj);
    if (relax == -1) {
      PyErr_SetString(
          PyExc_ValueError,
          "relax must be True, False or an integer.");
      return NULL;
    }
  }

  // If the user has requested to preserve the original units, then we
  // could try and edit the header returned by wcshdo - however, this
  // might be tricky and not robust to future WCSLIB changes. Instead,
  // we make a copy of the Wcsprm object, change the units on that and
  // prevent WCSLIB fixing the units, then convert to a header.
  if (self->unit_scaling != NULL) {

    Wcsprm* copy = Wcsprm_copy_with_patched_units(self);

    wcsprm_python2c(&copy->x);
    status = wcshdo(relax, &copy->x, &nkeyrec, &header);

    if (status != 0) {
      wcs_to_python_exc(&(copy->x));
      wcsfree(&copy->x);
      goto exit;
    }

    wcsfree(&copy->x);

  } else {

    wcsprm_python2c(&self->x);
    status = wcshdo(relax, &self->x, &nkeyrec, &header);
    wcsprm_c2python(&self->x);

    if (status != 0) {
      wcs_to_python_exc(&(self->x));
      goto exit;
    }

  }

  /* Just return the raw header string.  astropy.io.fits on the Python side will
     help to parse and use this information. */
  result = PyUnicode_FromStringAndSize(header, (Py_ssize_t)nkeyrec * 80);

 exit:
  free(header);
  return result;
}

/*@null@*/ static PyObject*
Wcsprm_unitfix(
    Wcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  const char* translate_units = NULL;
  int         ctrl            = 0;
  int         status          = 0;
  const char* keywords[]      = {"translate_units", NULL};

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "|s:unitfix", (char **)keywords,
          &translate_units)) {
    return NULL;
  }

  if (translate_units != NULL) {
    if (parse_unsafe_unit_conversion_spec(translate_units, &ctrl)) {
      return NULL;
    }
  }

  status = unitfix(ctrl, &self->x);

  if (status == -1 || status == 0) {
    return PyLong_FromLong((long)status);
  } else {
    wcserr_fix_to_python_exc(self->x.err);
    return NULL;
  }
}


/***************************************************************************
 * Member getters/setters (properties)
 */
/*@null@*/ static PyObject*
Wcsprm_get_alt(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.alt)) {
    return NULL;
  }

  /* Force a null-termination of this single-character string */
  self->x.alt[1] = '\0';
  return get_string("alt", self->x.alt);
}

static int
Wcsprm_set_alt(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  char value_string[2];

  if (is_null(self->x.alt)) {
    return -1;
  }

  if (value == NULL) { /* deletion */
    self->x.alt[0] = ' ';
    self->x.alt[1] = '\0';
    note_change(self);
    return 0;
  }

  if (set_string("alt", value, value_string, 2)) {
    return -1;
  }

  if (!is_valid_alt_key(value_string)) {
    return -1;
  }

  strncpy(self->x.alt, value_string, 2);

  return 0;
}

/*@null@*/ static PyObject*
Wcsprm_get_axis_types(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.types)) {
    return NULL;
  }

  if (Wcsprm_cset(self, 1)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_int_array("axis_types", self->x.types, 1, &naxis, (PyObject*)self);
}

static PyObject*
Wcsprm_get_bepoch(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("bepoch", self->x.bepoch);
}

static int
Wcsprm_set_bepoch(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.bepoch = (double)NPY_NAN;
    return 0;
  }

  return set_double("bepoch", value, &self->x.bepoch);
}


/*@null@*/ static PyObject*
Wcsprm_get_cd(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  npy_intp dims[2];
  PyObject* result;

  if (is_null(self->x.cd)) {
    return NULL;
  }

  if ((self->x.altlin & has_cd) == 0) {
    PyErr_SetString(PyExc_AttributeError, "No cd is present.");
    return NULL;
  }

  dims[0] = self->x.naxis;
  dims[1] = self->x.naxis;

  if (self->unit_scaling != NULL) {
    result = get_scaled_double_array_readonly("cd", self->x.cd, self->unit_scaling, 2, dims, (PyObject*)self);
  } else {
    result = get_double_array("cd", self->x.cd, 2, dims, (PyObject*)self);
  }

  return result;
}

static int
Wcsprm_set_cd(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp dims[2];

  if (is_null(self->x.cd)) {
    return -1;
  }


  if (value == NULL) {
    self->x.altlin &= ~has_cd;
    note_change(self);
    return 0;
  }

  dims[0] = self->x.naxis;
  dims[1] = self->x.naxis;

  if (set_double_array("cd", value, 2, dims, self->x.cd)) {
    return -1;
  }

  if (self->unit_scaling != NULL) {
    for (npy_intp i = 0; i < dims[0]; ++i) {
      for (npy_intp j = 0; j < dims[1]; ++j) {
        self->x.cd[i * dims[1] + j] *= self->unit_scaling[i];
      }
    }
  }

  self->x.altlin |= has_cd;

  note_change(self);

  return 0;
}

 /*@null@*/ static PyObject*
Wcsprm_get_cdelt(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;
  PyObject* result;

  if (is_null(self->x.cdelt)) {
    return NULL;
  }

  naxis = self->x.naxis;

  if (self->x.altlin & has_cd) {
    PyErr_WarnEx(NULL, "cdelt will be ignored since cd is present", 1);
  }

  if (self->unit_scaling != NULL) {
    result = get_scaled_double_array_readonly("cdelt", self->x.cdelt, self->unit_scaling, 1, &naxis, (PyObject*)self);
  } else {
    result = get_double_array("cdelt", self->x.cdelt, 1, &naxis, (PyObject*)self);
  }

  return result;


}

/*@null@*/ static int
Wcsprm_set_cdelt(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp dims;
  int status;

  if (is_null(self->x.cdelt)) {
    return -1;
  }

  dims = (npy_int)self->x.naxis;

  if (self->x.altlin & has_cd) {
    PyErr_WarnEx(NULL, "cdelt will be ignored since cd is present", 1);
  }

  note_change(self);

  status = set_double_array("cdelt", value, 1, &dims, self->x.cdelt);

  if (status == 0 && self->original_cunit != NULL) {
    for (npy_intp i = 0; i < dims; ++i) {
      self->x.cdelt[i] *= self->unit_scaling[i];
    }
  }

  return status;
}

static PyObject*
Wcsprm_get_cel_offset(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_bool("cel_offset", self->x.cel.offset);
}

static int
Wcsprm_set_cel_offset(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  return set_bool("cel_offset", value, &self->x.cel.offset);
}


/*@null@*/ static PyObject*
Wcsprm_get_cname(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.cname)) {
    return NULL;
  }

  return get_str_list("cname", self->x.cname, (Py_ssize_t)self->x.naxis, 68, (PyObject*)self);
}

/*@null@*/ static int
Wcsprm_set_cname(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {
  if (is_null(self->x.cname)) {
    return -1;
  }

  return set_str_list("cname", value, (Py_ssize_t)self->x.naxis, 0, self->x.cname);
}

/*@null@*/ static PyObject*
Wcsprm_get_colax(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.colax)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_int_array("colax", self->x.colax, 1, &naxis, (PyObject*)self);
}

static int
Wcsprm_set_colax(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 0;

  if (is_null(self->x.colax)) {
    return -1;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return set_int_array("colax", value, 1, &naxis, self->x.colax);
}

static PyObject*
Wcsprm_get_colnum(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("colnum", self->x.colnum);
}

static int
Wcsprm_set_colnum(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  return set_int("colnum", value, &self->x.colnum);
}

/*@null@*/ static PyObject*
Wcsprm_get_crder(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.crder)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("crder", self->x.crder, 1, &naxis, (PyObject*)self);
}

static int
Wcsprm_set_crder(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 0;

  if (is_null(self->x.crder)) {
    return -1;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return set_double_array("crder", value, 1, &naxis, self->x.crder);
}

/*@null@*/ static PyObject*
Wcsprm_get_crota(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.crota)) {
    return NULL;
  }

  if ((self->x.altlin & has_crota) == 0) {
    PyErr_SetString(PyExc_AttributeError, "No crota is present.");
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("crota", self->x.crota, 1, &naxis, (PyObject*)self);
}

static int
Wcsprm_set_crota(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 0;

  if (is_null(self->x.crota)) {
    return -1;
  }

  if (value == NULL) { /* Deletion */
    self->x.altlin &= ~has_crota;
    note_change(self);
    return 0;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  if (set_double_array("crota", value, 1, &naxis, self->x.crota)) {
    return -1;
  }

  self->x.altlin |= has_crota;

  note_change(self);

  return 0;
}

/*@null@*/ static PyObject*
Wcsprm_get_crpix(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.crpix)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("crpix", self->x.crpix, 1, &naxis, (PyObject*)self);
}

static int
Wcsprm_set_crpix(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis = 0;

  if (is_null(self->x.crpix)) {
    return -1;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  note_change(self);

  return set_double_array("crpix", value, 1, &naxis, self->x.crpix);
}

/*@null@*/ static PyObject*
Wcsprm_get_crval(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;
  PyObject* result;

  if (is_null(self->x.crval)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  if (self->unit_scaling != NULL) {
    result = get_scaled_double_array_readonly("cdelt", self->x.crval, self->unit_scaling, 1, &naxis, (PyObject*)self);
  } else {
    result = get_double_array("crval", self->x.crval, 1, &naxis, (PyObject*)self);
  }

  return result;

}

static int
Wcsprm_set_crval(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis;
  int status;

  if (is_null(self->x.crval)) {
    return -1;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  note_change(self);

  status = set_double_array("crval", value, 1, &naxis, self->x.crval);

  if (status == 0 && self->original_cunit != NULL) {
    for (npy_intp i = 0; i < naxis; ++i) {
      self->x.crval[i] *= self->unit_scaling[i];
    }
  }

  return status;

}

/*@null@*/ static PyObject*
Wcsprm_get_csyer(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis;

  if (is_null(self->x.csyer)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("csyer", self->x.csyer, 1, &naxis, (PyObject*)self);
}

static int
Wcsprm_set_csyer(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis;

  if (is_null(self->x.csyer)) {
    return -1;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return set_double_array("csyer", value, 1, &naxis, self->x.csyer);
}

/*@null@*/ static PyObject*
Wcsprm_get_ctype(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ctype)) {
    return NULL;
  }

  return get_str_list("ctype", self->x.ctype, self->x.naxis, 68, (PyObject*)self);
}

static int
Wcsprm_set_ctype(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ctype)) {
    return -1;
  }

  note_change(self);

  return set_str_list("ctype", value, (Py_ssize_t)self->x.naxis, 0, self->x.ctype);
}

static PyObject*
Wcsprm_get_cubeface(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("cubeface", self->x.cubeface);
}

static int
Wcsprm_set_cubeface(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  return set_int("cubeface", value, &self->x.cubeface);
}

/*@null@*/ static PyObject*
Wcsprm_get_cunit(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.cunit)) {
    return NULL;
  }

  if (self->original_cunit != NULL) {
    return get_unit_list(
      "cunit", self->original_cunit, (Py_ssize_t)self->x.naxis, (PyObject*)self, 1);
  } else {
    return get_unit_list(
      "cunit", self->x.cunit, (Py_ssize_t)self->x.naxis, (PyObject*)self, 0);
  }
}

static int
Wcsprm_set_cunit(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  double scale;

  if (is_null(self->x.cunit)) {
    return -1;
  }

  note_change(self);

  // Note that if the user explicitly changes the units we should now preserve
  // those new units. The way the WCS class works normally when units are not
  // preserved is that changing units on an existing WCS is not just a unit
  // conversion but truly overriding units. The simplest and safest way to achieve
  // this here is to undo any changes that were made by the WCS class related to
  // units (scaling cdelt, crval, and optionally cd), set the new units on the
  // WCS class, and let it figure out if the new units warrant re-instating the
  // scaling.

  if (self->unit_scaling != NULL) {

    for (int i = 0; i < self->x.naxis; ++i) {

      scale = self->unit_scaling[i];

      if (scale != 1.) {

        self->x.cdelt[i] /= scale;
        self->x.crval[i] /= scale;

        if (self->x.cd) {
            for (int j = 0; j < self->x.naxis; j++) {
                *(self->x.cd + i*self->x.naxis + j) /= scale;
            }
        }

      }

    }

    free(self->unit_scaling);
    self->unit_scaling = NULL;

  }

  if (self->original_cunit != NULL) {
    free(self->original_cunit);
    self->original_cunit = NULL;
  }

  return set_unit_list(
      (PyObject *)self, "cunit", value, (Py_ssize_t)self->x.naxis, self->x.cunit);

}

/*@null@*/ static PyObject*
Wcsprm_get_czphs(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis;

  if (is_null(self->x.czphs)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("czphs", self->x.czphs, 1, &naxis, (PyObject*)self);
}

static int
Wcsprm_set_czphs(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis;

  if (is_null(self->x.czphs)) {
    return -1;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return set_double_array("czphs", value, 1, &naxis, self->x.czphs);
}

/*@null@*/ static PyObject*
Wcsprm_get_cperi(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis;

  if (is_null(self->x.cperi)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("cperi", self->x.cperi, 1, &naxis, (PyObject*)self);
}

static int
Wcsprm_set_cperi(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis;

  if (is_null(self->x.cperi)) {
    return -1;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return set_double_array("cperi", value, 1, &naxis, self->x.cperi);
}

/*@null@*/ static PyObject*
Wcsprm_get_dateavg(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateavg)) {
    return NULL;
  }

  return get_string("dateavg", self->x.dateavg);
}

static int
Wcsprm_set_dateavg(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateavg)) {
    return -1;
  }

  /* TODO: Verify that this looks like a date string */

  return set_string("dateavg", value, self->x.dateavg, 72);
}

/*@null@*/ static PyObject*
Wcsprm_get_datebeg(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.datebeg)) {
    return NULL;
  }

  return get_string("datebeg", self->x.datebeg);
}

static int
Wcsprm_set_datebeg(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.datebeg)) {
    return -1;
  }

  return set_string("datebeg", value, self->x.datebeg, 72);
}

/*@null@*/ static PyObject*
Wcsprm_get_dateend(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateend)) {
    return NULL;
  }

  return get_string("dateend", self->x.dateend);
}

static int
Wcsprm_set_dateend(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateend)) {
    return -1;
  }

  return set_string("dateend", value, self->x.dateend, 72);
}


/*@null@*/ static PyObject*
Wcsprm_get_dateobs(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateobs)) {
    return NULL;
  }

  return get_string("dateobs", self->x.dateobs);
}

static int
Wcsprm_set_dateobs(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateobs)) {
    return -1;
  }

  return set_string("dateobs", value, self->x.dateobs, 72);
}

/*@null@*/ static PyObject*
Wcsprm_get_dateref(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateref)) {
    return NULL;
  }

  return get_string("dateref", self->x.dateref);
}

static int
Wcsprm_set_dateref(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateref)) {
    return -1;
  }

  return set_string("dateref", value, self->x.dateref, 72);
}

static PyObject*
Wcsprm_get_equinox(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("equinox", self->x.equinox);
}

static int
Wcsprm_set_equinox(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.equinox = (double)NPY_NAN;
    return 0;
  }

  return set_double("equinox", value, &self->x.equinox);
}

/*@null@*/ static PyObject*
Wcsprm_get_imgpix_matrix(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  npy_intp dims[2];

  if (is_null(self->x.lin.imgpix)) {
    return NULL;
  }

  if (Wcsprm_cset(self, 1)) {
    return NULL;
  }

  dims[0] = self->x.naxis;
  dims[1] = self->x.naxis;

  return get_double_array("imgpix_matrix", self->x.lin.imgpix, 2, dims,
                          (PyObject*)self);
}

static PyObject*
Wcsprm_get_jepoch(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("jepoch", self->x.jepoch);
}

static int
Wcsprm_set_jepoch(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  if (value == NULL) {
    self->x.jepoch = (double)NPY_NAN;
    return 0;
  }

  return set_double("jepoch", value, &self->x.jepoch);
}


static PyObject*
Wcsprm_get_lat(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (Wcsprm_cset(self, 1)) {
    return NULL;
  }

  return get_int("lat", self->x.lat);
}

static PyObject*
Wcsprm_get_latpole(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("latpole", self->x.latpole);
}

static int
Wcsprm_set_latpole(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  if (value == NULL) {
    self->x.latpole = 90.0;
    return 0;
  }

  return set_double("latpole", value, &self->x.latpole);
}

/*@null@*/ static PyObject*
Wcsprm_get_lattyp(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.lattyp)) {
    return NULL;
  }

  if (Wcsprm_cset(self, 1)) {
    return NULL;
  }

  return get_string("lattyp", self->x.lattyp);
}

static PyObject*
Wcsprm_get_lng(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (Wcsprm_cset(self, 1)) {
    return NULL;
  }

  return get_int("lng", self->x.lng);
}

/*@null@*/ static PyObject*
Wcsprm_get_lngtyp(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.lngtyp)) {
    return NULL;
  }

  if (Wcsprm_cset(self, 1)) {
    return NULL;
  }

  return get_string("lngtyp", self->x.lngtyp);
}

static PyObject*
Wcsprm_get_lonpole(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("lonpole", self->x.lonpole);
}

static int
Wcsprm_set_lonpole(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  if (value == NULL) {
    self->x.lonpole = (double)NPY_NAN;
    return 0;
  }

  return set_double("lonpole", value, &self->x.lonpole);
}

static PyObject*
Wcsprm_get_mjdavg(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("mjdavg", self->x.mjdavg);
}

static int
Wcsprm_set_mjdavg(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.mjdavg = (double)NPY_NAN;
    return 0;
  }

  return set_double("mjdavg", value, &self->x.mjdavg);
}

static PyObject*
Wcsprm_get_mjdbeg(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("mjdbeg", self->x.mjdbeg);
}

static int
Wcsprm_set_mjdbeg(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.mjdbeg = (double)NPY_NAN;
    return 0;
  }

  return set_double("mjdbeg", value, &self->x.mjdbeg);
}

static PyObject*
Wcsprm_get_mjdend(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("mjdend", self->x.mjdend);
}

static int
Wcsprm_set_mjdend(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.mjdend = (double)NPY_NAN;
    return 0;
  }

  return set_double("mjdend", value, &self->x.mjdend);
}

static PyObject*
Wcsprm_get_mjdobs(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("mjdobs", self->x.mjdobs);
}

static int
Wcsprm_set_mjdobs(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  if (value == NULL) {
    self->x.mjdobs = (double)NPY_NAN;
    return 0;
  }

  return set_double("mjdobs", value, &self->x.mjdobs);
}

static PyObject*
Wcsprm_get_mjdref(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  npy_intp size = 2;

  return get_double_array("mjdref", self->x.mjdref, 1, &size, (PyObject*)self);
}

static int
Wcsprm_set_mjdref(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp size = 2;

  if (value == NULL) {
    self->x.mjdref[0] = NPY_NAN;
    self->x.mjdref[1] = NPY_NAN;
    return 0;
  }
  return set_double_array("mjdref", value, 1, &size, self->x.mjdref);
}


/*@null@*/ static PyObject*
Wcsprm_get_timesys(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.timesys)) {
    return NULL;
  }

  return get_string("timesys", self->x.timesys);
}

static int
Wcsprm_set_timesys(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.timesys)) {
    return -1;
  }

  return set_string("timesys", value, self->x.timesys, 72);
}

/*@null@*/ static PyObject*
Wcsprm_get_trefpos(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.trefpos)) {
    return NULL;
  }

  return get_string("trefpos", self->x.trefpos);
}

static int
Wcsprm_set_trefpos(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.trefpos)) {
    return -1;
  }

  return set_string("trefpos", value, self->x.trefpos, 72);
}

/*@null@*/ static PyObject*
Wcsprm_get_trefdir(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.trefdir)) {
    return NULL;
  }

  return get_string("trefdir", self->x.trefdir);
}

static int
Wcsprm_set_trefdir(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.trefdir)) {
    return -1;
  }

  return set_string("trefdir", value, self->x.trefdir, 72);
}

/*@null@*/ static PyObject*
Wcsprm_get_timeunit(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.timeunit)) {
    return NULL;
  }

  return get_string("timeunit", self->x.timeunit);
}

static int
Wcsprm_set_timeunit(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.timeunit)) {
    return -1;
  }

  return set_string("timeunit", value, self->x.timeunit, 72);
}

/*@null@*/ static PyObject*
Wcsprm_get_plephem(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.plephem)) {
    return NULL;
  }

  return get_string("plephem", self->x.plephem);
}

static int
Wcsprm_set_plephem(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.plephem)) {
    return -1;
  }

  return set_string("plephem", value, self->x.plephem, 72);
}

static PyObject*
Wcsprm_get_tstart(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("tstart", self->x.tstart);
}

static int
Wcsprm_set_tstart(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.tstart = (double)NPY_NAN;
    return 0;
  }

  return set_double("tstart", value, &self->x.tstart);
}

static PyObject*
Wcsprm_get_tstop(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("tstop", self->x.tstop);
}

static int
Wcsprm_set_tstop(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.tstop = (double)NPY_NAN;
    return 0;
  }

  return set_double("tstop", value, &self->x.tstop);
}

static PyObject*
Wcsprm_get_telapse(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("telapse", self->x.telapse);
}

static int
Wcsprm_set_telapse(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.telapse = (double)NPY_NAN;
    return 0;
  }

  return set_double("telapse", value, &self->x.telapse);
}

static PyObject*
Wcsprm_get_timeoffs(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("timeoffs", self->x.timeoffs);
}

static int
Wcsprm_set_timeoffs(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.timeoffs = (double)NPY_NAN;
    return 0;
  }

  return set_double("timeoffs", value, &self->x.timeoffs);
}

static PyObject*
Wcsprm_get_timsyer(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("timsyer", self->x.timsyer);
}

static int
Wcsprm_set_timsyer(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.timsyer = (double)NPY_NAN;
    return 0;
  }

  return set_double("timsyer", value, &self->x.timsyer);
}

static PyObject*
Wcsprm_get_timrder(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("timrder", self->x.timrder);
}

static int
Wcsprm_set_timrder(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.timrder = (double)NPY_NAN;
    return 0;
  }

  return set_double("timrder", value, &self->x.timrder);
}

static PyObject*
Wcsprm_get_timedel(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("timedel", self->x.timedel);
}

static int
Wcsprm_set_timedel(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.timedel = (double)NPY_NAN;
    return 0;
  }

  return set_double("timedel", value, &self->x.timedel);
}

static PyObject*
Wcsprm_get_timepixr(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("timepixr", self->x.timepixr);
}

static int
Wcsprm_set_timepixr(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.timepixr = (double)NPY_NAN;
    return 0;
  }

  return set_double("timepixr", value, &self->x.timepixr);
}

/*@null@*/ static PyObject*
Wcsprm_get_obsorbit(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.obsorbit)) {
    return NULL;
  }

  return get_string("obsorbit", self->x.obsorbit);
}

static int
Wcsprm_set_obsorbit(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.obsorbit)) {
    return -1;
  }

  return set_string("obsorbit", value, self->x.obsorbit, 72);
}

static PyObject*
Wcsprm_get_xposure(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("xposure", self->x.xposure);
}

static int
Wcsprm_set_xposure(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.xposure = (double)NPY_NAN;
    return 0;
  }

  return set_double("xposure", value, &self->x.xposure);
}

/*@null@*/ static PyObject*
Wcsprm_get_name(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.wcsname)) {
    return NULL;
  }

  return get_string("name", self->x.wcsname);
}

static int
Wcsprm_set_name(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.wcsname)) {
    return -1;
  }

  return set_string("name", value, self->x.wcsname, 72);
}

static PyObject*
Wcsprm_get_naxis(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("naxis", self->x.naxis);
}

/*@null@*/ static PyObject*
Wcsprm_get_obsgeo(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t size = 6;

  if (is_null(self->x.obsgeo)) {
    return NULL;
  }

  return get_double_array("obsgeo", self->x.obsgeo, 1, &size, (PyObject*)self);
}

static int
Wcsprm_set_obsgeo(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp size = 6;

  if (is_null(self->x.obsgeo)) {
    return -1;
  }

  if (value == NULL) {
    self->x.obsgeo[0] = NPY_NAN;
    self->x.obsgeo[1] = NPY_NAN;
    self->x.obsgeo[2] = NPY_NAN;
    self->x.obsgeo[3] = NPY_NAN;
    self->x.obsgeo[4] = NPY_NAN;
    self->x.obsgeo[5] = NPY_NAN;
    return 0;
  }

  return set_double_array("obsgeo", value, 1, &size, self->x.obsgeo);
}

/*@null@*/ static PyObject*
Wcsprm_get_pc(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  npy_intp dims[2];

  if (is_null(self->x.pc)) {
    return NULL;
  }

  if (self->x.altlin != 0 && (self->x.altlin & has_pc) == 0) {
    PyErr_SetString(PyExc_AttributeError, "No pc is present.");
    return NULL;
  }

  dims[0] = self->x.naxis;
  dims[1] = self->x.naxis;

  return get_double_array("pc", self->x.pc, 2, dims, (PyObject*)self);
}

static int
Wcsprm_set_pc(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp dims[2];
  int i, j, naxis;
  double* pc;

  if (is_null(self->x.pc)) {
    return -1;
  }

  note_change(self);

  if (value == NULL) { /* deletion */
    self->x.altlin &= ~has_pc;

    /* If this results in deleting all flags, pc is still the default,
       so we should set the pc matrix itself to default values. */
    naxis = self->x.naxis;
    pc = self->x.pc;
    for (i = 0; i < naxis; i++) {
      for (j = 0; j < naxis; j++) {
        if (j == i) {
          *pc = 1.0;
        } else {
          *pc = 0.0;
        }
        pc++;
      }
    }

    note_change(self);

    return 0;
  }

  dims[0] = self->x.naxis;
  dims[1] = self->x.naxis;

  if (set_double_array("pc", value, 2, dims, self->x.pc)) {
    return -1;
  }

  self->x.altlin |= has_pc;

  note_change(self);

  return 0;
}

static PyObject*
Wcsprm_get_phi0(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("phi0", self->x.cel.phi0);
}

static int
Wcsprm_set_phi0(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  if (value == NULL) {
    self->x.cel.phi0 = (double)NPY_NAN;
    return 0;
  }

  return set_double("phi0", value, &(self->x.cel.phi0));
}

static PyObject*
Wcsprm_get_piximg_matrix(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  npy_intp dims[2];

  if (is_null(self->x.lin.piximg)) {
    return NULL;
  }

  if (Wcsprm_cset(self, 1)) {
    return NULL;
  }

  dims[0] = self->x.naxis;
  dims[1] = self->x.naxis;

  return get_double_array("piximg_matrix", self->x.lin.piximg, 2, dims,
                          (PyObject*)self);
}

static PyObject*
Wcsprm_get_radesys(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.radesys)) {
    return NULL;
  }

  return get_string("radesys", self->x.radesys);
}

static int
Wcsprm_set_radesys(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.radesys)) {
    return -1;
  }

  return set_string("radesys", value, self->x.radesys, 72);
}

static PyObject*
Wcsprm_get_restfrq(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("restfrq", self->x.restfrq);
}

static int
Wcsprm_set_restfrq(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.restfrq = (double)NPY_NAN;
    return 0;
  }

  note_change(self);

  return set_double("restfrq", value, &self->x.restfrq);
}

static PyObject*
Wcsprm_get_restwav(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("restwav", self->x.restwav);
}

static int
Wcsprm_set_restwav(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.restwav = (double)NPY_NAN;
    return 0;
  }

  note_change(self);

  return set_double("restwav", value, &self->x.restwav);
}

static PyObject*
Wcsprm_get_spec(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("spec", self->x.spec);
}

/*@null@*/ static PyObject*
Wcsprm_get_specsys(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.specsys)) {
    return NULL;
  }

  return get_string("specsys", self->x.specsys);
}

static int
Wcsprm_set_specsys(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.specsys)) {
    return -1;
  }

  return set_string("specsys", value, self->x.specsys, 72);
}

/*@null@*/ static PyObject*
Wcsprm_get_ssysobs(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ssysobs)) {
    return NULL;
  }

  return get_string("ssysobs", self->x.ssysobs);
}

static int
Wcsprm_set_ssysobs(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ssysobs)) {
    return -1;
  }

  note_change(self);

  return set_string("ssysobs", value, self->x.ssysobs, 72);
}

/*@null@*/ static PyObject*
Wcsprm_get_ssyssrc(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ssyssrc)) {
    return NULL;
  }

  return get_string("ssyssrc", self->x.ssyssrc);
}

static int
Wcsprm_set_ssyssrc(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ssyssrc)) {
    return -1;
  }

  return set_string("ssyssrc", value, self->x.ssyssrc, 72);
}

static PyObject*
Wcsprm_get_tab(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  PyObject* result;
  PyObject* subresult;
  int i, ntab;

  ntab = self->x.ntab;

  result = PyList_New(ntab);
  if (result == NULL) {
    return NULL;
  }

  for (i = 0; i < ntab; ++i) {
    subresult = (PyObject *)Tabprm_cnew((PyObject *)self, &(self->x.tab[i]));
    if (subresult == NULL) {
      Py_DECREF(result);
      return NULL;
    }

    if (PyList_SetItem(result, i, subresult) == -1) {
      Py_DECREF(result);
      return NULL;
    }
  }

  return result;
}

static PyObject*
Wcsprm_get_theta0(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("theta0", self->x.cel.theta0);
}

static int
Wcsprm_set_theta0(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  if (value == NULL) {
    self->x.cel.theta0 = (double)NPY_NAN;
    return 0;
  }

  return set_double("theta0", value, &self->x.cel.theta0);
}

static PyObject*
Wcsprm_get_velangl(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("velangl", self->x.velangl);
}

static int
Wcsprm_set_velangl(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.velangl = (double)NPY_NAN;
    return 0;
  }

  return set_double("velangl", value, &self->x.velangl);
}

static PyObject*
Wcsprm_get_velosys(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("velosys", self->x.velosys);
}

static int
Wcsprm_set_velosys(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.velosys = (double)NPY_NAN;
    return 0;
  }

  return set_double("velosys", value, &self->x.velosys);
}

static PyObject*
Wcsprm_get_velref(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("velref", self->x.velref);
}

static int
Wcsprm_set_velref(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.velref = 0;
    return 0;
  }

  return set_int("velref", value, &self->x.velref);
}


static PyObject* Wcsprm_get_wtb(Wcsprm* self, void* closure) {
  PyObject* list;
  PyObject* elem;
  int i, nwtb;

  nwtb = self->x.nwtb;

  list = PyList_New(nwtb);
  if (list == NULL) return NULL;

  for (i = 0; i < nwtb; ++i) {
    elem = (PyObject *)Wtbarr_cnew((PyObject *)self, &(self->x.wtb[i]));
    if (elem == NULL) {
      Py_DECREF(list);
      return NULL;
    }

    PyList_SetItem(list, i, elem);
  }

  return list;
}


static PyObject*
Wcsprm_get_zsource(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("zsource", self->x.zsource);
}

static int
Wcsprm_set_zsource(
    Wcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.zsource = (double)NPY_NAN;
    return 0;
  }

  return set_double("zsource", value, &self->x.zsource);
}


 /*@null@*/ static PyObject*
Wcsprm_get_unit_scaling(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;
  PyObject* result;

  if (self->unit_scaling == NULL) {
    return Py_None;
  }

  naxis = self->x.naxis;

  return get_double_array("unit_scaling", self->unit_scaling, 1, &naxis, (PyObject*)self);

}

static PyObject*
Wcsprm_get_aux(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  PyObject* result;

    // If wcsprm.aux is not initialized, we should do so here so that users can
    // set auxiliary parameters on an empty WCS.

    if (self->x.aux == 0x0) {
      wcsauxi(1, &self->x);
    }

  result = (PyObject *)Auxprm_cnew((PyObject *)self, self->x.aux);

  return result;
}


static PyObject*
Wcsprm_get_cel(
    Wcsprm* self,
    /*@unused@*/ void* closure) {

  return (PyObject *)Celprm_cnew((PyObject *)self, &(self->x.cel), NULL);
}

/***************************************************************************
 * Wcsprm definition structures
 */

static PyGetSetDef Wcsprm_getset[] = {
  {"alt", (getter)Wcsprm_get_alt, (setter)Wcsprm_set_alt, (char *)doc_alt},
  {"aux", (getter)Wcsprm_get_aux, NULL, (char *)doc_aux},
  {"cel", (getter)Wcsprm_get_cel, NULL, (char *)doc_cel},
  {"axis_types", (getter)Wcsprm_get_axis_types, NULL, (char *)doc_axis_types},
  {"bepoch", (getter)Wcsprm_get_bepoch, (setter)Wcsprm_set_bepoch, (char *)doc_bepoch},
  {"cd", (getter)Wcsprm_get_cd, (setter)Wcsprm_set_cd, (char *)doc_cd},
  {"cdelt", (getter)Wcsprm_get_cdelt, (setter)Wcsprm_set_cdelt, (char *)doc_cdelt},
  {"cel_offset", (getter)Wcsprm_get_cel_offset, (setter)Wcsprm_set_cel_offset, (char *)doc_cel_offset},
  {"cname", (getter)Wcsprm_get_cname, (setter)Wcsprm_set_cname, (char *)doc_cname},
  {"colax", (getter)Wcsprm_get_colax, (setter)Wcsprm_set_colax, (char *)doc_colax},
  {"colnum", (getter)Wcsprm_get_colnum, (setter)Wcsprm_set_colnum, (char *)doc_colnum},
  {"crder", (getter)Wcsprm_get_crder, (setter)Wcsprm_set_crder, (char *)doc_crder},
  {"crota", (getter)Wcsprm_get_crota, (setter)Wcsprm_set_crota, (char *)doc_crota},
  {"crpix", (getter)Wcsprm_get_crpix, (setter)Wcsprm_set_crpix, (char *)doc_crpix},
  {"crval", (getter)Wcsprm_get_crval, (setter)Wcsprm_set_crval, (char *)doc_crval},
  {"csyer", (getter)Wcsprm_get_csyer, (setter)Wcsprm_set_csyer, (char *)doc_csyer},
  {"ctype", (getter)Wcsprm_get_ctype, (setter)Wcsprm_set_ctype, (char *)doc_ctype},
  {"cubeface", (getter)Wcsprm_get_cubeface, (setter)Wcsprm_set_cubeface, (char *)doc_cubeface},
  {"cunit", (getter)Wcsprm_get_cunit, (setter)Wcsprm_set_cunit, (char *)doc_cunit},
  {"czphs", (getter)Wcsprm_get_czphs, (setter)Wcsprm_set_czphs, (char *)doc_czphs},
  {"cperi", (getter)Wcsprm_get_cperi, (setter)Wcsprm_set_cperi, (char *)doc_cperi},
  {"dateavg", (getter)Wcsprm_get_dateavg, (setter)Wcsprm_set_dateavg, (char *)doc_dateavg},
  {"datebeg", (getter)Wcsprm_get_datebeg, (setter)Wcsprm_set_datebeg, (char *)doc_datebeg},
  {"dateend", (getter)Wcsprm_get_dateend, (setter)Wcsprm_set_dateend, (char *)doc_dateend},
  {"dateobs", (getter)Wcsprm_get_dateobs, (setter)Wcsprm_set_dateobs, (char *)doc_dateobs},
  {"dateref", (getter)Wcsprm_get_dateref, (setter)Wcsprm_set_dateref, (char *)doc_dateref},
  {"equinox", (getter)Wcsprm_get_equinox, (setter)Wcsprm_set_equinox, (char *)doc_equinox},
  {"imgpix_matrix", (getter)Wcsprm_get_imgpix_matrix, NULL, (char *)doc_imgpix_matrix},
  {"jepoch", (getter)Wcsprm_get_jepoch, (setter)Wcsprm_set_jepoch, (char *)doc_jepoch},
  {"lat", (getter)Wcsprm_get_lat, NULL, (char *)doc_lat},
  {"latpole", (getter)Wcsprm_get_latpole, (setter)Wcsprm_set_latpole, (char *)doc_latpole},
  {"lattyp", (getter)Wcsprm_get_lattyp, NULL, (char *)doc_lattyp},
  {"lng", (getter)Wcsprm_get_lng, NULL, (char *)doc_lng},
  {"lngtyp", (getter)Wcsprm_get_lngtyp, NULL, (char *)doc_lngtyp},
  {"lonpole", (getter)Wcsprm_get_lonpole, (setter)Wcsprm_set_lonpole, (char *)doc_lonpole},
  {"mjdavg", (getter)Wcsprm_get_mjdavg, (setter)Wcsprm_set_mjdavg, (char *)doc_mjdavg},
  {"mjdbeg", (getter)Wcsprm_get_mjdbeg, (setter)Wcsprm_set_mjdbeg, (char *)doc_mjdbeg},
  {"mjdend", (getter)Wcsprm_get_mjdend, (setter)Wcsprm_set_mjdend, (char *)doc_mjdend},
  {"mjdobs", (getter)Wcsprm_get_mjdobs, (setter)Wcsprm_set_mjdobs, (char *)doc_mjdobs},
  {"mjdref", (getter)Wcsprm_get_mjdref, (setter)Wcsprm_set_mjdref, (char *)doc_mjdref},
  {"name", (getter)Wcsprm_get_name, (setter)Wcsprm_set_name, (char *)doc_name},
  {"naxis", (getter)Wcsprm_get_naxis, NULL, (char *)doc_naxis},
  {"obsgeo", (getter)Wcsprm_get_obsgeo, (setter)Wcsprm_set_obsgeo, (char *)doc_obsgeo},
  {"obsorbit", (getter)Wcsprm_get_obsorbit, (setter)Wcsprm_set_obsorbit, (char *)doc_obsorbit},
  {"pc", (getter)Wcsprm_get_pc, (setter)Wcsprm_set_pc, (char *)doc_pc},
  {"phi0", (getter)Wcsprm_get_phi0, (setter)Wcsprm_set_phi0, (char *)doc_phi0},
  {"piximg_matrix", (getter)Wcsprm_get_piximg_matrix, NULL, (char *)doc_piximg_matrix},
  {"plephem", (getter)Wcsprm_get_plephem, (setter)Wcsprm_set_plephem, (char *) doc_plephem},
  {"radesys", (getter)Wcsprm_get_radesys, (setter)Wcsprm_set_radesys, (char *)doc_radesys},
  {"restfrq", (getter)Wcsprm_get_restfrq, (setter)Wcsprm_set_restfrq, (char *)doc_restfrq},
  {"restwav", (getter)Wcsprm_get_restwav, (setter)Wcsprm_set_restwav, (char *)doc_restwav},
  {"spec", (getter)Wcsprm_get_spec, NULL, (char *)doc_spec},
  {"specsys", (getter)Wcsprm_get_specsys, (setter)Wcsprm_set_specsys, (char *)doc_specsys},
  {"ssysobs", (getter)Wcsprm_get_ssysobs, (setter)Wcsprm_set_ssysobs, (char *)doc_ssysobs},
  {"ssyssrc", (getter)Wcsprm_get_ssyssrc, (setter)Wcsprm_set_ssyssrc, (char *)doc_ssyssrc},
  {"tab", (getter)Wcsprm_get_tab, NULL, (char *)doc_tab},
  {"theta0", (getter)Wcsprm_get_theta0, (setter)Wcsprm_set_theta0, (char *)doc_theta0},
  {"timesys", (getter)Wcsprm_get_timesys, (setter)Wcsprm_set_timesys, (char *) doc_timesys},
  {"trefpos", (getter)Wcsprm_get_trefpos, (setter)Wcsprm_set_trefpos, (char *) doc_trefpos},
  {"trefdir", (getter)Wcsprm_get_trefdir, (setter)Wcsprm_set_trefdir, (char *) doc_trefdir},
  {"tstart", (getter)Wcsprm_get_tstart, (setter)Wcsprm_set_tstart, (char *) doc_tstart},
  {"tstop", (getter)Wcsprm_get_tstop, (setter)Wcsprm_set_tstop, (char *) doc_tstop},
  {"telapse", (getter)Wcsprm_get_telapse, (setter)Wcsprm_set_telapse, (char *) doc_telapse},
  {"timeoffs", (getter)Wcsprm_get_timeoffs, (setter)Wcsprm_set_timeoffs, (char *) doc_timeoffs},
  {"timsyer", (getter)Wcsprm_get_timsyer, (setter)Wcsprm_set_timsyer, (char *) doc_timsyer},
  {"timrder", (getter)Wcsprm_get_timrder, (setter)Wcsprm_set_timrder, (char *) doc_timrder},
  {"timedel", (getter)Wcsprm_get_timedel, (setter)Wcsprm_set_timedel, (char *) doc_timedel},
  {"timepixr", (getter)Wcsprm_get_timepixr, (setter)Wcsprm_set_timepixr, (char *) doc_timepixr},
  {"timeunit", (getter)Wcsprm_get_timeunit, (setter)Wcsprm_set_timeunit, (char *) doc_timeunit},
  {"velangl", (getter)Wcsprm_get_velangl, (setter)Wcsprm_set_velangl, (char *)doc_velangl},
  {"velosys", (getter)Wcsprm_get_velosys, (setter)Wcsprm_set_velosys, (char *)doc_velosys},
  {"velref", (getter)Wcsprm_get_velref, (setter)Wcsprm_set_velref, (char *)doc_velref},
  {"xposure", (getter)Wcsprm_get_xposure, (setter)Wcsprm_set_xposure, (char *)doc_xposure},
  {"wtb", (getter)Wcsprm_get_wtb, NULL, (char *) doc_wtb},
  {"zsource", (getter)Wcsprm_get_zsource, (setter)Wcsprm_set_zsource, (char *)doc_zsource},
  {"_unit_scaling", (getter)Wcsprm_get_unit_scaling, NULL, NULL},  // For debugging
  {NULL}
};

static PyMethodDef Wcsprm_methods[] = {
  {"bounds_check", (PyCFunction)Wcsprm_bounds_check, METH_VARARGS|METH_KEYWORDS, doc_bounds_check},
  {"cdfix", (PyCFunction)Wcsprm_cdfix, METH_NOARGS, doc_cdfix},
  {"celfix", (PyCFunction)Wcsprm_celfix, METH_NOARGS, doc_celfix},
  {"compare", (PyCFunction)Wcsprm_compare, METH_VARARGS|METH_KEYWORDS, doc_compare},
  {"__copy__", (PyCFunction)Wcsprm_copy, METH_NOARGS, doc_copy},
  {"cylfix", (PyCFunction)Wcsprm_cylfix, METH_VARARGS|METH_KEYWORDS, doc_cylfix},
  {"datfix", (PyCFunction)Wcsprm_datfix, METH_NOARGS, doc_datfix},
  {"__deepcopy__", (PyCFunction)Wcsprm_copy, METH_O, doc_copy},
  {"fix", (PyCFunction)Wcsprm_fix, METH_VARARGS|METH_KEYWORDS, doc_fix},
  {"get_cdelt", (PyCFunction)Wcsprm_get_cdelt_func, METH_NOARGS, doc_get_cdelt},
  {"get_pc", (PyCFunction)Wcsprm_get_pc_func, METH_NOARGS, doc_get_pc},
  {"get_ps", (PyCFunction)Wcsprm_get_ps, METH_NOARGS, doc_get_ps},
  {"get_pv", (PyCFunction)Wcsprm_get_pv, METH_NOARGS, doc_get_pv},
  {"has_cd", (PyCFunction)Wcsprm_has_cdi_ja, METH_NOARGS, doc_has_cd},
  {"has_cdi_ja", (PyCFunction)Wcsprm_has_cdi_ja, METH_NOARGS, doc_has_cdi_ja},
  {"has_crota", (PyCFunction)Wcsprm_has_crotaia, METH_NOARGS, doc_has_crota},
  {"has_crotaia", (PyCFunction)Wcsprm_has_crotaia, METH_NOARGS, doc_has_crotaia},
  {"has_pc", (PyCFunction)Wcsprm_has_pci_ja, METH_NOARGS, doc_has_pc},
  {"has_pci_ja", (PyCFunction)Wcsprm_has_pci_ja, METH_NOARGS, doc_has_pci_ja},
  {"is_unity", (PyCFunction)Wcsprm_is_unity, METH_NOARGS, doc_is_unity},
  {"mix", (PyCFunction)Wcsprm_mix, METH_VARARGS|METH_KEYWORDS, doc_mix},
  {"p2s", (PyCFunction)Wcsprm_p2s, METH_VARARGS|METH_KEYWORDS, doc_p2s},
  {"print_contents", (PyCFunction)Wcsprm_print_contents, METH_NOARGS, doc_print_contents},
  {"s2p", (PyCFunction)Wcsprm_s2p, METH_VARARGS|METH_KEYWORDS, doc_s2p},
  {"set", (PyCFunction)Wcsprm_set, METH_NOARGS, doc_set},
  {"set_ps", (PyCFunction)Wcsprm_set_ps, METH_O, doc_set_ps},
  {"set_pv", (PyCFunction)Wcsprm_set_pv, METH_O, doc_set_pv},
  {"spcfix", (PyCFunction)Wcsprm_spcfix, METH_NOARGS, doc_spcfix},
  {"sptr", (PyCFunction)Wcsprm_sptr, METH_VARARGS|METH_KEYWORDS, doc_sptr},
  {"sub", (PyCFunction)Wcsprm_sub, METH_VARARGS|METH_KEYWORDS, doc_sub},
  {"to_header", (PyCFunction)Wcsprm_to_header, METH_VARARGS|METH_KEYWORDS, doc_to_header},
  {"unitfix", (PyCFunction)Wcsprm_unitfix, METH_VARARGS|METH_KEYWORDS, doc_unitfix},
  {NULL}
};

static PyType_Spec Wcsprm_spec = {
  .name = "astropy.wcs.Wcsprm",
  .basicsize = sizeof(Wcsprm),
  .itemsize = 0,
  .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_IMMUTABLETYPE,
  .slots = (PyType_Slot[]) {
    {Py_tp_dealloc, (destructor)Wcsprm_dealloc},
    {Py_tp_repr, (reprfunc)Wcsprm___str__},
    {Py_tp_str, (reprfunc)Wcsprm___str__},
    {Py_tp_doc, doc_Wcsprm},
    {Py_tp_richcompare, Wcsprm_richcompare},
    {Py_tp_methods, Wcsprm_methods},
    {Py_tp_getset, Wcsprm_getset},
    {Py_tp_init, (initproc)Wcsprm_init},
    {Py_tp_new, Wcsprm_new},
    {0, NULL},
  },
};

PyObject* WcsprmType = NULL;

#define CONSTANT(a) PyModule_AddIntConstant(m, #a, a)
#define CONSTANT2(n, v) PyModule_AddIntConstant(m, n, v)

#define XSTRINGIFY(s) STRINGIFY(s)
#define STRINGIFY(s) #s

int add_prj_codes(PyObject* module)
{
    int k;
    PyObject* code;
    PyObject* list = PyList_New(prj_ncode);
    if (list == NULL) {
        return -1;
    }

    for (k = 0; k < prj_ncode; k++) {
        code = PyUnicode_FromString(prj_codes[k]);
        if (PyList_SetItem(list, k, code)) {
            Py_DECREF(list);
            return -1;
        }
    }

    if (PyModule_AddObject(module, "PRJ_CODES", list)) {
        Py_DECREF(list);
        return -1;
    }
    return 0;
}

int
_setup_wcsprm_type(
    PyObject* m) {
  WcsprmType = PyType_FromSpec(&Wcsprm_spec);

  if (WcsprmType == NULL) {
    return -1;
  }

  wcsprintf_set(NULL);
  wcserr_enable(1);

  return (
    PyModule_AddObject(m, "Wcsprm", WcsprmType) ||
    CONSTANT(WCSSUB_LONGITUDE) ||
    CONSTANT(WCSSUB_LATITUDE)  ||
    CONSTANT(WCSSUB_CUBEFACE)  ||
    CONSTANT(WCSSUB_SPECTRAL)  ||
    CONSTANT(WCSSUB_STOKES)    ||
#if defined(WCSSUB_TIME)
    CONSTANT(WCSSUB_TIME)      ||
#endif
    CONSTANT(WCSSUB_CELESTIAL) ||
    CONSTANT(WCSHDR_IMGHEAD)   ||
    CONSTANT(WCSHDR_BIMGARR)   ||
    CONSTANT(WCSHDR_PIXLIST)   ||
    CONSTANT(WCSHDR_none)      ||
    CONSTANT(WCSHDR_all)       ||
    CONSTANT(WCSHDR_reject)    ||
#ifdef WCSHDR_strict
    CONSTANT(WCSHDR_strict)    ||
#endif
    CONSTANT(WCSHDR_CROTAia)   ||
    CONSTANT(WCSHDR_EPOCHa)    ||
    CONSTANT(WCSHDR_VELREFa)   ||
    CONSTANT(WCSHDR_CD00i00j)  ||
    CONSTANT(WCSHDR_PC00i00j)  ||
    CONSTANT(WCSHDR_PROJPn)    ||
#ifdef WCSHDR_CD0i_0ja
    CONSTANT(WCSHDR_CD0i_0ja)  ||
#endif
#ifdef WCSHDR_PC0i_0ja
    CONSTANT(WCSHDR_PC0i_0ja)  ||
#endif
#ifdef WCSHDR_PV0i_0ma
    CONSTANT(WCSHDR_PV0i_0ma)  ||
#endif
#ifdef WCSHDR_PS0i_0ma
    CONSTANT(WCSHDR_PS0i_0ma)  ||
#endif
    CONSTANT(WCSHDR_RADECSYS)  ||
    CONSTANT(WCSHDR_VSOURCE)   ||
    CONSTANT(WCSHDR_DOBSn)     ||
    CONSTANT(WCSHDR_LONGKEY)   ||
    CONSTANT(WCSHDR_CNAMn)     ||
    CONSTANT(WCSHDR_AUXIMG)    ||
    CONSTANT(WCSHDR_ALLIMG)    ||
    CONSTANT(WCSHDO_none)      ||
    CONSTANT(WCSHDO_all)       ||
    CONSTANT(WCSHDO_safe)      ||
    CONSTANT(WCSHDO_DOBSn)     ||
    CONSTANT(WCSHDO_TPCn_ka)   ||
    CONSTANT(WCSHDO_PVn_ma)    ||
    CONSTANT(WCSHDO_CRPXna)    ||
    CONSTANT(WCSHDO_CNAMna)    ||
    CONSTANT(WCSHDO_WCSNna)    ||
    CONSTANT(WCSHDO_P12)       ||
    CONSTANT(WCSHDO_P13)       ||
    CONSTANT(WCSHDO_P14)       ||
    CONSTANT(WCSHDO_P15)       ||
    CONSTANT(WCSHDO_P16)       ||
    CONSTANT(WCSHDO_P17)       ||
    CONSTANT(WCSHDO_EFMT)      ||
    CONSTANT(WCSCOMPARE_ANCILLARY) ||
    CONSTANT(WCSCOMPARE_TILING) ||
    CONSTANT(WCSCOMPARE_CRPIX)  ||
    CONSTANT2("PRJ_PVN", PVN)       ||
    add_prj_codes(m) ||
    CONSTANT2("PRJ_ZENITHAL", ZENITHAL)                   ||
    CONSTANT2("PRJ_CYLINDRICAL", CYLINDRICAL)             ||
    CONSTANT2("PRJ_PSEUDOCYLINDRICAL", PSEUDOCYLINDRICAL) ||
    CONSTANT2("PRJ_CONVENTIONAL", CONVENTIONAL)           ||
    CONSTANT2("PRJ_CONIC", CONIC)                         ||
    CONSTANT2("PRJ_POLYCONIC", POLYCONIC)                 ||
    CONSTANT2("PRJ_QUADCUBE", QUADCUBE)                   ||
    CONSTANT2("PRJ_HEALPIX", HEALPIX));
}
