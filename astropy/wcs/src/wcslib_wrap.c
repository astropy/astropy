/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#define NO_IMPORT_ARRAY

#include "astropy_wcs/wcslib_wrap.h"
#include "astropy_wcs/wcslib_tabprm_wrap.h"
#include "astropy_wcs/wcslib_wtbarr_wrap.h"
#include "astropy_wcs/wcslib_units_wrap.h"
#include "astropy_wcs/unit_list_proxy.h"
#include <structmember.h> /* from Python */

#include <wcs.h>
#include <wcsfix.h>
#include <wcshdr.h>
#include <wcsmath.h>
#include <wcsprintf.h>
#include <wcsunits.h>
#include <tab.h>

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

/***************************************************************************
 * PyWcsprm methods
 */

static int
PyWcsprm_cset(PyWcsprm* self, const int convert);

static INLINE void
note_change(PyWcsprm* self) {
  self->x.flag = 0;
}

static void
PyWcsprm_dealloc(
    PyWcsprm* self) {

  wcsfree(&self->x);
  Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyWcsprm*
PyWcsprm_cnew(void) {
  PyWcsprm* self;
  self = (PyWcsprm*)(&PyWcsprmType)->tp_alloc(&PyWcsprmType, 0);
  return self;
}

static PyObject *
PyWcsprm_new(
    PyTypeObject* type,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  PyWcsprm* self;
  self = (PyWcsprm*)type->tp_alloc(type, 0);
  return (PyObject*)self;
}

static int
PyWcsprm_init(
    PyWcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  int            status;
  PyObject*      header_obj    = NULL;
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
  int*           colsel_ints   = NULL;
  int            warnings      = 1;
  int            nreject       = 0;
  int            nwcs          = 0;
  struct wcsprm* wcs           = NULL;
  int            i             = 0;
  const char*    keywords[]    = {"header", "key", "relax", "naxis", "keysel",
                                  "colsel", "warnings", NULL};

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "|OsOiiOi:WCSBase.__init__",
          (char **)keywords, &header_obj, &key, &relax_obj, &naxis, &keysel,
          &colsel, &warnings)) {
    return -1;
  }

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

    note_change(self);
    self->x.flag = -1;
    status = wcsini(1, naxis, &self->x);

    if (status != 0) {
      PyErr_SetString(
          PyExc_MemoryError,
          self->x.err->msg);
      return -1;
    }

    self->x.alt[0] = key[0];

    if (PyWcsprm_cset(self, 0)) {
      return -1;
    }
    wcsprm_c2python(&self->x);

    return 0;
  } else { /* header != NULL */
    #if PY3K
    if (PyBytes_AsStringAndSize(header_obj, &header, &header_length)) {
    #else
    if (PyString_AsStringAndSize(header_obj, &header, &header_length)) {
    #endif
      return -1;
    }

    if (relax_obj == Py_True) {
      relax = WCSHDR_all;
    } else if (relax_obj == NULL || relax_obj == Py_False) {
      relax = WCSHDR_none;
    } else {
      #if PY3K
      relax = (int)PyLong_AsLong(relax_obj);
      #else
      relax = (int)PyInt_AsLong(relax_obj);
      #endif
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
        colsel, 1, 1, PyArray_INT);
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
      for (i = 0; i < colsel_ints[0]; ++i) {
        colsel_ints[i+1] = colsel_array->data[i];
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
      PyErr_SetString(
          PyExc_MemoryError,
          "Memory allocation error.");
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
      PyErr_SetString(
          PyExc_MemoryError,
          "Memory allocation error.");
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

    note_change(self);
    wcsprm_c2python(&self->x);
    wcsvfree(&nwcs, &wcs);
    return 0;
  }
}

/*@null@*/ static PyObject*
PyWcsprm_bounds_check(
    PyWcsprm* self,
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

  wcsbchk(&self->x, bounds);

  Py_RETURN_NONE;
}


/*@null@*/ static PyObject*
PyWcsprm_copy(
    PyWcsprm* self) {

  PyWcsprm*      copy      = NULL;
  int            status;

  copy = PyWcsprm_cnew();
  if (copy == NULL) {
    return NULL;
  }

  wcsini(0, self->x.naxis, &copy->x);

  wcsprm_python2c(&self->x);
  status = wcscopy(1, &self->x, &copy->x);
  wcsprm_c2python(&self->x);

  if (status == 0) {
    if (PyWcsprm_cset(copy, 0)) {
      Py_XDECREF(copy);
      return NULL;
    }
    wcsprm_c2python(&copy->x);
    return (PyObject*)copy;
  } else {
    Py_XDECREF(copy);
    wcs_to_python_exc(&(self->x));
    return NULL;
  }
}

PyObject*
PyWcsprm_find_all_wcs(
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
  PyWcsprm*      subresult     = NULL;
  int            i             = 0;
  const char*    keywords[]    = {"header", "relax", "keysel", "warnings", NULL};
  int            status        = -1;

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "O|Oii:find_all_wcs",
          (char **)keywords, &header_obj, &relax_obj, &keysel, &warnings)) {
    return NULL;
  }

  #if PY3K
  if (PyBytes_AsStringAndSize(header_obj, &header, &header_length)) {
  #else
  if (PyString_AsStringAndSize(header_obj, &header, &header_length)) {
  #endif
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
    #if PY3K
    relax = (int)PyLong_AsLong(relax_obj);
    #else
    relax = (int)PyInt_AsLong(relax_obj);
    #endif
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
    PyErr_SetString(
        PyExc_MemoryError,
        "Memory allocation error.");
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
    PyErr_SetString(
        PyExc_MemoryError,
        "Memory allocation error.");
    return NULL;
  }

  result = PyList_New(nwcs);
  if (result == NULL) {
    wcsvfree(&nwcs, &wcs);
    return NULL;
  }

  for (i = 0; i < nwcs; ++i) {
    subresult = PyWcsprm_cnew();
    if (wcscopy(1, wcs + i, &subresult->x) != 0) {
      Py_DECREF(result);
      wcsvfree(&nwcs, &wcs);
      PyErr_SetString(
          PyExc_MemoryError,
          "Could not initialize wcsprm object");
      return NULL;
    }

    if (PyList_SetItem(result, i, (PyObject *)subresult) == -1) {
      Py_DECREF(subresult);
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
PyWcsprm_cdfix(
    PyWcsprm* self) {

  int status = 0;

  wcsprm_python2c(&self->x);
  status = cdfix(&self->x);
  wcsprm_c2python(&self->x);

  if (status == -1 || status == 0) {
    #if PY3K
    return PyLong_FromLong((long)status);
    #else
    return PyInt_FromLong((long)status);
    #endif
  } else {
    wcserr_fix_to_python_exc(self->x.err);
    return NULL;
  }
}

static PyObject*
PyWcsprm_celfix(
    PyWcsprm* self) {

  int status = 0;

  wcsprm_python2c(&self->x);
  status = celfix(&self->x);
  wcsprm_c2python(&self->x);

  if (status == -1 || status == 0) {
    #if PY3K
    return PyLong_FromLong((long)status);
    #else
    return PyInt_FromLong((long)status);
    #endif
  } else {
    wcserr_fix_to_python_exc(self->x.err);
    return NULL;
  }
}

static PyObject *
PyWcsprm_compare(
    PyWcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  int cmp = 0;
  PyWcsprm *other;
  double tolerance = 0.0;
  int equal;
  int status;

  const char* keywords[] = {"other", "cmp", "tolerance", NULL};

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "O!|id:compare", (char **)keywords,
          &PyWcsprmType, &other, &cmp, &tolerance)) {
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
PyWcsprm_cylfix(
    PyWcsprm* self,
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
        naxis_obj, 1, 1, PyArray_INT);
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

  Py_XDECREF(naxis_array);

  if (status == -1 || status == 0) {
    #if PY3K
    return PyLong_FromLong((long)status);
    #else
    return PyInt_FromLong((long)status);
    #endif
  } else {
    wcserr_fix_to_python_exc(self->x.err);
    return NULL;
  }
}

static PyObject*
PyWcsprm_datfix(
    PyWcsprm* self) {

  int status = 0;

  wcsprm_python2c(&self->x);
  status = datfix(&self->x);
  wcsprm_c2python(&self->x);

  if (status == -1 || status == 0) {
    #if PY3K
    return PyLong_FromLong((long)status);
    #else
    return PyInt_FromLong((long)status);
    #endif
  } else {
    wcserr_fix_to_python_exc(self->x.err);
    return NULL;
  }
}

/*@null@*/ static PyObject*
PyWcsprm_fix(
    PyWcsprm* self,
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
        naxis_obj, 1, 1, PyArray_INT);
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

  wcsprm_python2c(&self->x);
  wcsfixi(ctrl, naxis, &self->x, stat, err);
  wcsprm_c2python(&self->x);

  /* We're done with this already, so deref now so we don't have to remember
     later */
  Py_XDECREF(naxis_array);

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
    #if PY3K
    subresult = PyUnicode_FromString(message);
    #else
    subresult = PyString_FromString(message);
    #endif
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
PyWcsprm_get_cdelt_func(
    PyWcsprm* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.cdelt)) {
    return NULL;
  }

  if (PyWcsprm_cset(self, 1)) {
    return NULL;
  }

  naxis = self->x.naxis;

  return get_double_array_readonly("cdelt", self->x.cdelt, 1, &naxis, (PyObject*)self);
}

/*@null@*/ static PyObject*
PyWcsprm_get_pc_func(
    PyWcsprm* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  npy_intp dims[2];

  if (is_null(self->x.pc)) {
    return NULL;
  }

  if (PyWcsprm_cset(self, 1)) {
    return NULL;
  }

  dims[0] = self->x.naxis;
  dims[1] = self->x.naxis;

  return get_double_array_readonly("pc", self->x.pc, 2, dims, (PyObject*)self);
}

/*@null@*/ static PyObject*
PyWcsprm_get_ps(
    PyWcsprm* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  return get_pscards("ps", self->x.ps, self->x.nps);
}

/*@null@*/ static PyObject*
PyWcsprm_get_pv(
    PyWcsprm* self,
    /*@unused@*/ PyObject* args,
    /*@unused@*/ PyObject* kwds) {

  return get_pvcards("pv", self->x.pv, self->x.npv);
}

static PyObject*
PyWcsprm_has_cdi_ja(
    PyWcsprm* self) {

  int result = 0;

  result = self->x.altlin & has_cd;

  return PyBool_FromLong(result);
}

static PyObject*
PyWcsprm_has_crotaia(
    PyWcsprm* self) {

  int result = 0;

  result = self->x.altlin & has_crota;

  return PyBool_FromLong(result);
}

static PyObject*
PyWcsprm_has_pci_ja(
    PyWcsprm* self) {

  int result = 0;

  result = (self->x.altlin == 0 || self->x.altlin & has_pc);

  return PyBool_FromLong(result);
}

static PyObject*
PyWcsprm_is_unity(
    PyWcsprm* self) {

  if (PyWcsprm_cset(self, 1)) {
    return NULL;
  }

  return PyBool_FromLong(self->x.lin.unity);
}

/*@null@*/ static PyObject*
PyWcsprm_mix(
    PyWcsprm* self,
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
    (world_obj, PyArray_DOUBLE, 1, 1);
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
    (pixcrd_obj, PyArray_DOUBLE, 1, 1);
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
    (1, &naxis, PyArray_DOUBLE);
  if (phi == NULL) {
    goto exit;
  }

  theta = (PyArrayObject*)PyArray_SimpleNew
    (1, &naxis, PyArray_DOUBLE);
  if (theta == NULL) {
    goto exit;
  }

  imgcrd = (PyArrayObject*)PyArray_SimpleNew
    (1, &naxis, PyArray_DOUBLE);
  if (imgcrd == NULL) {
    goto exit;
  }

  /* Convert pixel coordinates to 1-based */
  Py_BEGIN_ALLOW_THREADS
  preoffset_array(pixcrd, origin);
  wcsprm_python2c(&self->x);
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
  wcsprm_c2python(&self->x);
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
  Py_XDECREF(world);
  Py_XDECREF(phi);
  Py_XDECREF(theta);
  Py_XDECREF(imgcrd);
  Py_XDECREF(pixcrd);

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
PyWcsprm_p2s(
    PyWcsprm* self,
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
    (pixcrd_obj, PyArray_DOUBLE, 2, 2);
  if (pixcrd == NULL) {
    return NULL;
  }

  if (PyArray_DIM(pixcrd, 1) < naxis) {
    PyErr_Format(
      PyExc_RuntimeError,
      "Input array must be 2-dimensional, where the second dimension >= %d",
      naxis);
    goto exit;
  }

  /* Now we allocate a bunch of numpy arrays to store the results in.
   */
  imgcrd = (PyArrayObject*)PyArray_SimpleNew(
      2, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (imgcrd == NULL) {
    goto exit;
  }

  phi = (PyArrayObject*)PyArray_SimpleNew(
      1, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (phi == NULL) {
    goto exit;
  }

  theta = (PyArrayObject*)PyArray_SimpleNew(
      1, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (theta == NULL) {
    goto exit;
  }

  world = (PyArrayObject*)PyArray_SimpleNew(
      2, PyArray_DIMS(pixcrd), PyArray_DOUBLE);
  if (world == NULL) {
    goto exit;
  }

  stat = (PyArrayObject*)PyArray_SimpleNew(
      1, PyArray_DIMS(pixcrd), PyArray_INT);
  if (stat == NULL) {
    goto exit;
  }

  /* Make the call */
  Py_BEGIN_ALLOW_THREADS
  ncoord = PyArray_DIM(pixcrd, 0);
  nelem = PyArray_DIM(pixcrd, 1);
  preoffset_array(pixcrd, origin);
  wcsprm_python2c(&self->x);
  status = wcsp2s(
      &self->x,
      ncoord,
      nelem,
      (double*)PyArray_DATA(pixcrd),
      (double*)PyArray_DATA(imgcrd),
      (double*)PyArray_DATA(phi),
      (double*)PyArray_DATA(theta),
      (double*)PyArray_DATA(world),
      (int*)PyArray_DATA(stat));
  wcsprm_c2python(&self->x);
  unoffset_array(pixcrd, origin);
  /* unoffset_array(world, origin); */
  unoffset_array(imgcrd, origin);
  if (status == 8) {
    set_invalid_to_nan(
        ncoord, nelem, (double*)PyArray_DATA(imgcrd), (int*)PyArray_DATA(stat));
    set_invalid_to_nan(
        ncoord, 1, (double*)PyArray_DATA(phi), (int*)PyArray_DATA(stat));
    set_invalid_to_nan(
        ncoord, 1, (double*)PyArray_DATA(theta), (int*)PyArray_DATA(stat));
    set_invalid_to_nan(
        ncoord, nelem, (double*)PyArray_DATA(world), (int*)PyArray_DATA(stat));
  }
  Py_END_ALLOW_THREADS

  if (status == 0 || status == 8) {
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
  Py_XDECREF(pixcrd);
  Py_XDECREF(imgcrd);
  Py_XDECREF(phi);
  Py_XDECREF(theta);
  Py_XDECREF(world);
  Py_XDECREF(stat);

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
PyWcsprm_s2p(
    PyWcsprm* self,
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
  const char*    keywords[] = {
    "world", "origin", NULL };

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "Oi:s2p", (char **)keywords,
          &world_obj, &origin)) {
    return NULL;
  }

  naxis = self->x.naxis;

  world = (PyArrayObject*)PyArray_ContiguousFromAny(
      world_obj, PyArray_DOUBLE, 2, 2);
  if (world == NULL) {
    return NULL;
  }

  if (PyArray_DIM(world, 1) < naxis) {
    PyErr_Format(
      PyExc_RuntimeError,
      "Input array must be 2-dimensional, where the second dimension >= %d",
      naxis);
    goto exit;
  }

  /* Now we allocate a bunch of numpy arrays to store the
   * results in.
   */
  phi = (PyArrayObject*)PyArray_SimpleNew(
      1, PyArray_DIMS(world), PyArray_DOUBLE);
  if (phi == NULL) {
    goto exit;
  }

  theta = (PyArrayObject*)PyArray_SimpleNew(
      1, PyArray_DIMS(world), PyArray_DOUBLE);
  if (phi == NULL) {
    goto exit;
  }

  imgcrd = (PyArrayObject*)PyArray_SimpleNew(
      2, PyArray_DIMS(world), PyArray_DOUBLE);
  if (theta == NULL) {
    goto exit;
  }

  pixcrd = (PyArrayObject*)PyArray_SimpleNew(
      2, PyArray_DIMS(world), PyArray_DOUBLE);
  if (pixcrd == NULL) {
    goto exit;
  }

  stat = (PyArrayObject*)PyArray_SimpleNew(
      1, PyArray_DIMS(world), PyArray_INT);
  if (stat == NULL) {
    goto exit;
  }

  /* Make the call */
  Py_BEGIN_ALLOW_THREADS
  ncoord = (int)PyArray_DIM(world, 0);
  nelem = (int)PyArray_DIM(world, 1);
  /* preoffset_array(world, origin); */
  wcsprm_python2c(&self->x);
  status = wcss2p(
      &self->x,
      ncoord,
      nelem,
      (double*)PyArray_DATA(world),
      (double*)PyArray_DATA(phi),
      (double*)PyArray_DATA(theta),
      (double*)PyArray_DATA(imgcrd),
      (double*)PyArray_DATA(pixcrd),
      (int*)PyArray_DATA(stat));
  wcsprm_c2python(&self->x);
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
  Py_XDECREF(pixcrd);
  Py_XDECREF(imgcrd);
  Py_XDECREF(phi);
  Py_XDECREF(theta);
  Py_XDECREF(world);
  Py_XDECREF(stat);

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

static int
PyWcsprm_cset(
    PyWcsprm* self,
    const int convert) {

  int status = 0;

  if (convert) wcsprm_python2c(&self->x);
  status = wcsset(&self->x);
  if (convert) wcsprm_c2python(&self->x);

  if (status == 0) {
    return 0;
  } else {
    wcs_to_python_exc(&(self->x));
    return 1;
  }
}

/*@null@*/ static PyObject*
PyWcsprm_set(
    PyWcsprm* self) {

  if (PyWcsprm_cset(self, 1)) {
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

/*@null@*/ static PyObject*
PyWcsprm_set_ps(
    PyWcsprm* self,
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
PyWcsprm_set_pv(
    PyWcsprm* self,
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
PyWcsprm_print_contents(
    PyWcsprm* self) {

  /* This is not thread-safe, but since we're holding onto the GIL,
     we can assume we won't have thread conflicts */
  wcsprintf_set(NULL);

  wcsprm_python2c(&self->x);
  if (PyWcsprm_cset(self, 0)) {
    wcsprm_c2python(&self->x);
    return NULL;
  }
  wcsprt(&self->x);
  wcsprm_c2python(&self->x);

  printf("%s", wcsprintf_buf());

  Py_INCREF(Py_None);
  return Py_None;
}

/*@null@*/ static PyObject*
PyWcsprm_spcfix(
    PyWcsprm* self) {

  int status = 0;

  wcsprm_python2c(&self->x);
  status = spcfix(&self->x);
  wcsprm_c2python(&self->x);

  if (status == -1 || status == 0) {
    #if PY3K
    return PyLong_FromLong((long)status);
    #else
    return PyInt_FromLong((long)status);
    #endif
  } else {
    wcserr_fix_to_python_exc(self->x.err);
    return NULL;
  }
}

/*@null@*/ static PyObject*
PyWcsprm_sptr(
    PyWcsprm* self,
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
PyWcsprm___str__(
    PyWcsprm* self) {

  /* This is not thread-safe, but since we're holding onto the GIL,
     we can assume we won't have thread conflicts */
  wcsprintf_set(NULL);

  wcsprm_python2c(&self->x);
  if (PyWcsprm_cset(self, 0)) {
    wcsprm_c2python(&self->x);
    return NULL;
  }
  wcsprt(&self->x);
  wcsprm_c2python(&self->x);

  #if PY3K
  return PyUnicode_FromString(wcsprintf_buf());
  #else
  return PyString_FromString(wcsprintf_buf());
  #endif
}

PyObject *PyWcsprm_richcompare(PyObject *a, PyObject *b, int op) {
  int equal;
  int status;

  struct wcsprm *ax;
  struct wcsprm *bx;

  if ((op == Py_EQ || op == Py_NE) &&
      PyObject_TypeCheck(b, &PyWcsprmType)) {
    ax = &((PyWcsprm *)a)->x;
    bx = &((PyWcsprm *)b)->x;

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
      wcs_to_python_exc(&(((PyWcsprm *)a)->x));
      return NULL;
    }
  }

  Py_INCREF(Py_NotImplemented);
  return Py_NotImplemented;
}

/*@null@*/ static PyObject*
PyWcsprm_sub(
    PyWcsprm* self,
    PyObject* args,
    PyObject* kwds) {

  int        i              = -1;
  Py_ssize_t tmp            = 0;
  PyObject*  py_axes        = NULL;
  PyWcsprm*  py_dest_wcs    = NULL;
  PyObject*  element        = NULL;
  PyObject*  element_utf8   = NULL;
  char*      element_str    = NULL;
  int        element_val    = 0;
  int        nsub           = 0;
  int*       axes           = NULL;
  int        status         = -1;
  const char*    keywords[] = {"axes", NULL};

  if (!PyArg_ParseTupleAndKeywords(
          args, kwds, "|O:sub", (char **)keywords,
          &py_axes)) {
    goto exit;
  }

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
          Py_DECREF(element_utf8); element_utf8 = NULL;
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
        } else if (strncmp(element_str, "celestial", 10) == 0) {
          element_val = WCSSUB_CELESTIAL;
        } else {
          PyErr_SetString(
            PyExc_ValueError,
            "string values for axis sequence must be one of 'latitude', 'longitude', 'cubeface', 'spectral', 'stokes', or 'celestial'");
          goto exit;
        }
      #if PY3K
      } else if (PyLong_Check(element)) {
        tmp = (Py_ssize_t)PyLong_AsSsize_t(element);
      #else
      } else if (PyInt_Check(element)) {
        tmp = (Py_ssize_t)PyInt_AsLong(element);
      #endif
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

      Py_DECREF(element);
      element = NULL;
    }
  #if PY3K
  } else if (PyLong_Check(py_axes)) {
    tmp = (Py_ssize_t)PyLong_AsSsize_t(py_axes);
  #else
  } else if (PyInt_Check(py_axes)) {
    tmp = (Py_ssize_t)PyInt_AsLong(py_axes);
  #endif
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

  py_dest_wcs = (PyWcsprm*)PyWcsprm_cnew();
  py_dest_wcs->x.flag = -1;
  status = wcsini(0, nsub, &py_dest_wcs->x);
  if (status != 0) {
    goto exit;
  }

  wcsprm_python2c(&self->x);
  status = wcssub(1, &self->x, &nsub, axes, &py_dest_wcs->x);
  wcsprm_c2python(&self->x);
  if (PyWcsprm_cset(py_dest_wcs, 0)) {
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
  #if PY3K
  Py_XDECREF(element_utf8);
  #endif

  if (status == 0) {
    return (PyObject*)py_dest_wcs;
  } else if (status == -1) {
    Py_XDECREF(py_dest_wcs);
    /* Exception already set */
    return NULL;
  } else {
    wcs_to_python_exc(&(py_dest_wcs->x));
    Py_XDECREF(py_dest_wcs);
    return NULL;
  }
}

/*@null@*/ static PyObject*
PyWcsprm_to_header(
    PyWcsprm* self,
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
    #if PY3K
    relax = (int)PyLong_AsLong(relax_obj);
    #else
    relax = (int)PyInt_AsLong(relax_obj);
    #endif
    if (relax == -1) {
      PyErr_SetString(
          PyExc_ValueError,
          "relax must be True, False or an integer.");
      return NULL;
    }
  }

  wcsprm_python2c(&self->x);
  status = wcshdo(relax, &self->x, &nkeyrec, &header);
  wcsprm_c2python(&self->x);

  if (status != 0) {
    wcs_to_python_exc(&(self->x));
    goto exit;
  }

  /* Just return the raw header string.  PyFITS on the Python side will help
     to parse and use this information. */
  #if PY3K
  result = PyUnicode_FromStringAndSize(header, (Py_ssize_t)nkeyrec * 80);
  #else
  result = PyString_FromStringAndSize(header, (Py_ssize_t)nkeyrec * 80);
  #endif

 exit:
  free(header);
  return result;
}

/*@null@*/ static PyObject*
PyWcsprm_unitfix(
    PyWcsprm* self,
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
    #if PY3K
    return PyLong_FromLong((long)status);
    #else
    return PyInt_FromLong((long)status);
    #endif
  } else {
    wcserr_fix_to_python_exc(self->x.err);
    return NULL;
  }
}


/***************************************************************************
 * Member getters/setters (properties)
 */
/*@null@*/ static PyObject*
PyWcsprm_get_alt(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.alt)) {
    return NULL;
  }

  /* Force a null-termination of this single-character string */
  self->x.alt[1] = '\0';
  return get_string("alt", self->x.alt);
}

static int
PyWcsprm_set_alt(
    PyWcsprm* self,
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
PyWcsprm_get_axis_types(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.types)) {
    return NULL;
  }

  if (PyWcsprm_cset(self, 1)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_int_array("axis_types", self->x.types, 1, &naxis, (PyObject*)self);
}

/*@null@*/ static PyObject*
PyWcsprm_get_cd(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  npy_intp dims[2];

  if (is_null(self->x.cd)) {
    return NULL;
  }

  if ((self->x.altlin & has_cd) == 0) {
    PyErr_SetString(PyExc_AttributeError, "No cd is present.");
    return NULL;
  }

  dims[0] = self->x.naxis;
  dims[1] = self->x.naxis;

  return get_double_array("cd", self->x.cd, 2, dims, (PyObject*)self);
}

static int
PyWcsprm_set_cd(
    PyWcsprm* self,
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

  self->x.altlin |= has_cd;

  note_change(self);

  return 0;
}

 /*@null@*/ static PyObject*
PyWcsprm_get_cdelt(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.cdelt)) {
    return NULL;
  }

  naxis = self->x.naxis;

  if (self->x.altlin & has_cd) {
    PyErr_WarnEx(NULL, "cdelt will be ignored since cd is present", 1);
  }

  return get_double_array("cdelt", self->x.cdelt, 1, &naxis, (PyObject*)self);
}

/*@null@*/ static int
PyWcsprm_set_cdelt(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp dims;

  if (is_null(self->x.cdelt)) {
    return -1;
  }

  dims = (npy_int)self->x.naxis;

  if (self->x.altlin & has_cd) {
    PyErr_WarnEx(NULL, "cdelt will be ignored since cd is present", 1);
  }

  note_change(self);

  return set_double_array("cdelt", value, 1, &dims, self->x.cdelt);
}

static PyObject*
PyWcsprm_get_cel_offset(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_bool("cel_offset", self->x.cel.offset);
}

static int
PyWcsprm_set_cel_offset(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  return set_bool("cel_offset", value, &self->x.cel.offset);
}


/*@null@*/ static PyObject*
PyWcsprm_get_cname(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.cname)) {
    return NULL;
  }

  return get_str_list("cname", self->x.cname, (Py_ssize_t)self->x.naxis, 68, (PyObject*)self);
}

/*@null@*/ static int
PyWcsprm_set_cname(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {
  if (is_null(self->x.cname)) {
    return -1;
  }

  return set_str_list("cname", value, (Py_ssize_t)self->x.naxis, 0, self->x.cname);
}

/*@null@*/ static PyObject*
PyWcsprm_get_colax(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.colax)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_int_array("colax", self->x.colax, 1, &naxis, (PyObject*)self);
}

static int
PyWcsprm_set_colax(
    PyWcsprm* self,
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
PyWcsprm_get_colnum(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("colnum", self->x.colnum);
}

static int
PyWcsprm_set_colnum(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  return set_int("colnum", value, &self->x.colnum);
}

/*@null@*/ static PyObject*
PyWcsprm_get_crder(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.crder)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("crder", self->x.crder, 1, &naxis, (PyObject*)self);
}

static int
PyWcsprm_set_crder(
    PyWcsprm* self,
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
PyWcsprm_get_crota(
    PyWcsprm* self,
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
PyWcsprm_set_crota(
    PyWcsprm* self,
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
PyWcsprm_get_crpix(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.crpix)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("crpix", self->x.crpix, 1, &naxis, (PyObject*)self);
}

static int
PyWcsprm_set_crpix(
    PyWcsprm* self,
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
PyWcsprm_get_crval(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis = 0;

  if (is_null(self->x.crval)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("crval", self->x.crval, 1, &naxis, (PyObject*)self);
}

static int
PyWcsprm_set_crval(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp naxis;

  if (is_null(self->x.crval)) {
    return -1;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  note_change(self);

  return set_double_array("crval", value, 1, &naxis, self->x.crval);
}

/*@null@*/ static PyObject*
PyWcsprm_get_csyer(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t naxis;

  if (is_null(self->x.csyer)) {
    return NULL;
  }

  naxis = (Py_ssize_t)self->x.naxis;

  return get_double_array("csyer", self->x.csyer, 1, &naxis, (PyObject*)self);
}

static int
PyWcsprm_set_csyer(
    PyWcsprm* self,
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
PyWcsprm_get_ctype(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ctype)) {
    return NULL;
  }

  return get_str_list("ctype", self->x.ctype, self->x.naxis, 68, (PyObject*)self);
}

static int
PyWcsprm_set_ctype(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ctype)) {
    return -1;
  }

  note_change(self);

  return set_str_list("ctype", value, (Py_ssize_t)self->x.naxis, 0, self->x.ctype);
}

static PyObject*
PyWcsprm_get_cubeface(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("cubeface", self->x.cubeface);
}

static int
PyWcsprm_set_cubeface(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  return set_int("cubeface", value, &self->x.cubeface);
}

/*@null@*/ static PyObject*
PyWcsprm_get_cunit(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.cunit)) {
    return NULL;
  }

  return get_unit_list(
    "cunit", self->x.cunit, (Py_ssize_t)self->x.naxis, (PyObject*)self);
}

static int
PyWcsprm_set_cunit(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.cunit)) {
    return -1;
  }

  note_change(self);

  return set_unit_list(
    (PyObject *)self, "cunit", value, (Py_ssize_t)self->x.naxis, self->x.cunit);
}

/*@null@*/ static PyObject*
PyWcsprm_get_dateavg(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateavg)) {
    return NULL;
  }

  return get_string("dateavg", self->x.dateavg);
}

static int
PyWcsprm_set_dateavg(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateavg)) {
    return -1;
  }

  /* TODO: Verify that this looks like a date string */

  return set_string("dateavg", value, self->x.dateavg, 72);
}

/*@null@*/ static PyObject*
PyWcsprm_get_dateobs(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateobs)) {
    return NULL;
  }

  return get_string("dateobs", self->x.dateobs);
}

static int
PyWcsprm_set_dateobs(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.dateobs)) {
    return -1;
  }

  return set_string("dateobs", value, self->x.dateobs, 72);
}

static PyObject*
PyWcsprm_get_equinox(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("equinox", self->x.equinox);
}

static int
PyWcsprm_set_equinox(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.equinox = (double)NPY_NAN;
    return 0;
  }

  return set_double("equinox", value, &self->x.equinox);
}

/*@null@*/ static PyObject*
PyWcsprm_get_imgpix_matrix(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  npy_intp dims[2];

  if (is_null(self->x.lin.imgpix)) {
    return NULL;
  }

  if (PyWcsprm_cset(self, 1)) {
    return NULL;
  }

  dims[0] = self->x.naxis;
  dims[1] = self->x.naxis;

  return get_double_array("imgpix_matrix", self->x.lin.imgpix, 2, dims,
                          (PyObject*)self);
}

static PyObject*
PyWcsprm_get_lat(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (PyWcsprm_cset(self, 1)) {
    return NULL;
  }

  return get_int("lat", self->x.lat);
}

static PyObject*
PyWcsprm_get_latpole(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("latpole", self->x.latpole);
}

static int
PyWcsprm_set_latpole(
    PyWcsprm* self,
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
PyWcsprm_get_lattyp(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.lattyp)) {
    return NULL;
  }

  if (PyWcsprm_cset(self, 1)) {
    return NULL;
  }

  return get_string("lattyp", self->x.lattyp);
}

static PyObject*
PyWcsprm_get_lng(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (PyWcsprm_cset(self, 1)) {
    return NULL;
  }

  return get_int("lng", self->x.lng);
}

/*@null@*/ static PyObject*
PyWcsprm_get_lngtyp(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.lngtyp)) {
    return NULL;
  }

  if (PyWcsprm_cset(self, 1)) {
    return NULL;
  }

  return get_string("lngtyp", self->x.lngtyp);
}

static PyObject*
PyWcsprm_get_lonpole(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("lonpole", self->x.lonpole);
}

static int
PyWcsprm_set_lonpole(
    PyWcsprm* self,
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
PyWcsprm_get_mjdavg(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("mjdavg", self->x.mjdavg);
}

static int
PyWcsprm_set_mjdavg(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) {
    self->x.mjdavg = (double)NPY_NAN;
    return 0;
  }

  return set_double("mjdavg", value, &self->x.mjdavg);
}

static PyObject*
PyWcsprm_get_mjdobs(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("mjdobs", self->x.mjdobs);
}

static int
PyWcsprm_set_mjdobs(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  note_change(self);

  if (value == NULL) {
    self->x.mjdobs = (double)NPY_NAN;
    return 0;
  }

  return set_double("mjdobs", value, &self->x.mjdobs);
}

/*@null@*/ static PyObject*
PyWcsprm_get_name(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.wcsname)) {
    return NULL;
  }

  return get_string("name", self->x.wcsname);
}

static int
PyWcsprm_set_name(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.wcsname)) {
    return -1;
  }

  return set_string("name", value, self->x.wcsname, 72);
}

static PyObject*
PyWcsprm_get_naxis(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("naxis", self->x.naxis);
}

/*@null@*/ static PyObject*
PyWcsprm_get_obsgeo(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  Py_ssize_t size = 3;

  if (is_null(self->x.obsgeo)) {
    return NULL;
  }

  return get_double_array("obsgeo", self->x.obsgeo, 1, &size, (PyObject*)self);
}

static int
PyWcsprm_set_obsgeo(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  npy_intp size = 3;

  if (is_null(self->x.obsgeo)) {
    return -1;
  }

  if (value == NULL) {
    self->x.obsgeo[0] = NPY_NAN;
    self->x.obsgeo[1] = NPY_NAN;
    self->x.obsgeo[2] = NPY_NAN;
    return 0;
  }

  return set_double_array("obsgeo", value, 1, &size, self->x.obsgeo);
}

/*@null@*/ static PyObject*
PyWcsprm_get_pc(
    PyWcsprm* self,
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
PyWcsprm_set_pc(
    PyWcsprm* self,
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
PyWcsprm_get_phi0(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("phi0", self->x.cel.phi0);
}

static int
PyWcsprm_set_phi0(
    PyWcsprm* self,
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
PyWcsprm_get_piximg_matrix(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  npy_intp dims[2];

  if (is_null(self->x.lin.piximg)) {
    return NULL;
  }

  if (PyWcsprm_cset(self, 1)) {
    return NULL;
  }

  dims[0] = self->x.naxis;
  dims[1] = self->x.naxis;

  return get_double_array("piximg_matrix", self->x.lin.piximg, 2, dims,
                          (PyObject*)self);
}

static PyObject*
PyWcsprm_get_radesys(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.radesys)) {
    return NULL;
  }

  return get_string("radesys", self->x.radesys);
}

static int
PyWcsprm_set_radesys(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.radesys)) {
    return -1;
  }

  return set_string("radesys", value, self->x.radesys, 72);
}

static PyObject*
PyWcsprm_get_restfrq(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("restfrq", self->x.restfrq);
}

static int
PyWcsprm_set_restfrq(
    PyWcsprm* self,
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
PyWcsprm_get_restwav(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("restwav", self->x.restwav);
}

static int
PyWcsprm_set_restwav(
    PyWcsprm* self,
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
PyWcsprm_get_spec(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("spec", self->x.spec);
}

/*@null@*/ static PyObject*
PyWcsprm_get_specsys(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.specsys)) {
    return NULL;
  }

  return get_string("specsys", self->x.specsys);
}

static int
PyWcsprm_set_specsys(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.specsys)) {
    return -1;
  }

  return set_string("specsys", value, self->x.specsys, 72);
}

/*@null@*/ static PyObject*
PyWcsprm_get_ssysobs(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ssysobs)) {
    return NULL;
  }

  return get_string("ssysobs", self->x.ssysobs);
}

static int
PyWcsprm_set_ssysobs(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ssysobs)) {
    return -1;
  }

  note_change(self);

  return set_string("ssysobs", value, self->x.ssysobs, 72);
}

/*@null@*/ static PyObject*
PyWcsprm_get_ssyssrc(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ssyssrc)) {
    return NULL;
  }

  return get_string("ssyssrc", self->x.ssyssrc);
}

static int
PyWcsprm_set_ssyssrc(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (is_null(self->x.ssyssrc)) {
    return -1;
  }

  return set_string("ssyssrc", value, self->x.ssyssrc, 72);
}

static PyObject*
PyWcsprm_get_tab(
    PyWcsprm* self,
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
    subresult = (PyObject *)PyTabprm_cnew((PyObject *)self, &(self->x.tab[i]));
    if (subresult == NULL) {
      Py_DECREF(result);
      return NULL;
    }

    if (PyList_SetItem(result, i, subresult) == -1) {
      Py_DECREF(subresult);
      Py_DECREF(result);
      return NULL;
    }
  }

  return result;
}

static PyObject*
PyWcsprm_get_theta0(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("theta0", self->x.cel.theta0);
}

static int
PyWcsprm_set_theta0(
    PyWcsprm* self,
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
PyWcsprm_get_velangl(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("velangl", self->x.velangl);
}

static int
PyWcsprm_set_velangl(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.velangl = (double)NPY_NAN;
    return 0;
  }

  return set_double("velangl", value, &self->x.velangl);
}

static PyObject*
PyWcsprm_get_velosys(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("velosys", self->x.velosys);
}

static int
PyWcsprm_set_velosys(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.velosys = (double)NPY_NAN;
    return 0;
  }

  return set_double("velosys", value, &self->x.velosys);
}

static PyObject*
PyWcsprm_get_velref(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_int("velref", self->x.velref);
}

static int
PyWcsprm_set_velref(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.velref = 0;
    return 0;
  }

  return set_int("velref", value, &self->x.velref);
}

/* static PyObject* */
/* PyWcsprm_get_wtb( */
/*     PyWcsprm* self, */
/*     /\*@unused@*\/ void* closure) { */

/*   PyObject* result; */
/*   PyObject* subresult; */
/*   int i, nwtb; */

/*   nwtb = self->x.nwtb; */

/*   result = PyList_New(nwtb); */
/*   if (result == NULL) { */
/*     return NULL; */
/*   } */

/*   for (i = 0; i < nwtb; ++i) { */
/*     subresult = (PyObject *)PyWtbarr_cnew((PyObject *)self, &(self->x.wtb[i])); */
/*     if (subresult == NULL) { */
/*       Py_DECREF(result); */
/*       return NULL; */
/*     } */

/*     if (PyList_SetItem(result, i, subresult) == -1) { */
/*       Py_DECREF(subresult); */
/*       Py_DECREF(result); */
/*       return NULL; */
/*     } */
/*   } */

/*   return result; */
/* } */

static PyObject*
PyWcsprm_get_zsource(
    PyWcsprm* self,
    /*@unused@*/ void* closure) {

  return get_double("zsource", self->x.zsource);
}

static int
PyWcsprm_set_zsource(
    PyWcsprm* self,
    PyObject* value,
    /*@unused@*/ void* closure) {

  if (value == NULL) { /* deletion */
    self->x.zsource = (double)NPY_NAN;
    return 0;
  }

  return set_double("zsource", value, &self->x.zsource);
}

/***************************************************************************
 * PyWcsprm definition structures
 */

static PyGetSetDef PyWcsprm_getset[] = {
  {"alt", (getter)PyWcsprm_get_alt, (setter)PyWcsprm_set_alt, (char *)doc_alt},
  {"axis_types", (getter)PyWcsprm_get_axis_types, NULL, (char *)doc_axis_types},
  {"cd", (getter)PyWcsprm_get_cd, (setter)PyWcsprm_set_cd, (char *)doc_cd},
  {"cdelt", (getter)PyWcsprm_get_cdelt, (setter)PyWcsprm_set_cdelt, (char *)doc_cdelt},
  {"cel_offset", (getter)PyWcsprm_get_cel_offset, (setter)PyWcsprm_set_cel_offset, (char *)doc_cel_offset},
  {"cname", (getter)PyWcsprm_get_cname, (setter)PyWcsprm_set_cname, (char *)doc_cname},
  {"colax", (getter)PyWcsprm_get_colax, (setter)PyWcsprm_set_colax, (char *)doc_colax},
  {"colnum", (getter)PyWcsprm_get_colnum, (setter)PyWcsprm_set_colnum, (char *)doc_colnum},
  {"crder", (getter)PyWcsprm_get_crder, (setter)PyWcsprm_set_crder, (char *)doc_crder},
  {"crota", (getter)PyWcsprm_get_crota, (setter)PyWcsprm_set_crota, (char *)doc_crota},
  {"crpix", (getter)PyWcsprm_get_crpix, (setter)PyWcsprm_set_crpix, (char *)doc_crpix},
  {"crval", (getter)PyWcsprm_get_crval, (setter)PyWcsprm_set_crval, (char *)doc_crval},
  {"csyer", (getter)PyWcsprm_get_csyer, (setter)PyWcsprm_set_csyer, (char *)doc_csyer},
  {"ctype", (getter)PyWcsprm_get_ctype, (setter)PyWcsprm_set_ctype, (char *)doc_ctype},
  {"cubeface", (getter)PyWcsprm_get_cubeface, (setter)PyWcsprm_set_cubeface, (char *)doc_cubeface},
  {"cunit", (getter)PyWcsprm_get_cunit, (setter)PyWcsprm_set_cunit, (char *)doc_cunit},
  {"dateavg", (getter)PyWcsprm_get_dateavg, (setter)PyWcsprm_set_dateavg, (char *)doc_dateavg},
  {"dateobs", (getter)PyWcsprm_get_dateobs, (setter)PyWcsprm_set_dateobs, (char *)doc_dateobs},
  {"equinox", (getter)PyWcsprm_get_equinox, (setter)PyWcsprm_set_equinox, (char *)doc_equinox},
  {"imgpix_matrix", (getter)PyWcsprm_get_imgpix_matrix, NULL, (char *)doc_imgpix_matrix},
  {"lat", (getter)PyWcsprm_get_lat, NULL, (char *)doc_lat},
  {"latpole", (getter)PyWcsprm_get_latpole, (setter)PyWcsprm_set_latpole, (char *)doc_latpole},
  {"lattyp", (getter)PyWcsprm_get_lattyp, NULL, (char *)doc_lattyp},
  {"lng", (getter)PyWcsprm_get_lng, NULL, (char *)doc_lng},
  {"lngtyp", (getter)PyWcsprm_get_lngtyp, NULL, (char *)doc_lngtyp},
  {"lonpole", (getter)PyWcsprm_get_lonpole, (setter)PyWcsprm_set_lonpole, (char *)doc_lonpole},
  {"mjdavg", (getter)PyWcsprm_get_mjdavg, (setter)PyWcsprm_set_mjdavg, (char *)doc_mjdavg},
  {"mjdobs", (getter)PyWcsprm_get_mjdobs, (setter)PyWcsprm_set_mjdobs, (char *)doc_mjdobs},
  {"name", (getter)PyWcsprm_get_name, (setter)PyWcsprm_set_name, (char *)doc_name},
  {"naxis", (getter)PyWcsprm_get_naxis, NULL, (char *)doc_naxis},
  {"obsgeo", (getter)PyWcsprm_get_obsgeo, (setter)PyWcsprm_set_obsgeo, (char *)doc_obsgeo},
  {"pc", (getter)PyWcsprm_get_pc, (setter)PyWcsprm_set_pc, (char *)doc_pc},
  {"phi0", (getter)PyWcsprm_get_phi0, (setter)PyWcsprm_set_phi0, (char *)doc_phi0},
  {"piximg_matrix", (getter)PyWcsprm_get_piximg_matrix, NULL, (char *)doc_piximg_matrix},
  {"radesys", (getter)PyWcsprm_get_radesys, (setter)PyWcsprm_set_radesys, (char *)doc_radesys},
  {"restfrq", (getter)PyWcsprm_get_restfrq, (setter)PyWcsprm_set_restfrq, (char *)doc_restfrq},
  {"restwav", (getter)PyWcsprm_get_restwav, (setter)PyWcsprm_set_restwav, (char *)doc_restwav},
  {"spec", (getter)PyWcsprm_get_spec, NULL, (char *)doc_spec},
  {"specsys", (getter)PyWcsprm_get_specsys, (setter)PyWcsprm_set_specsys, (char *)doc_specsys},
  {"ssysobs", (getter)PyWcsprm_get_ssysobs, (setter)PyWcsprm_set_ssysobs, (char *)doc_ssysobs},
  {"ssyssrc", (getter)PyWcsprm_get_ssyssrc, (setter)PyWcsprm_set_ssyssrc, (char *)doc_ssyssrc},
  {"tab", (getter)PyWcsprm_get_tab, NULL, (char *)doc_tab},
  {"theta0", (getter)PyWcsprm_get_theta0, (setter)PyWcsprm_set_theta0, (char *)doc_theta0},
  {"velangl", (getter)PyWcsprm_get_velangl, (setter)PyWcsprm_set_velangl, (char *)doc_velangl},
  {"velosys", (getter)PyWcsprm_get_velosys, (setter)PyWcsprm_set_velosys, (char *)doc_velosys},
  {"velref", (getter)PyWcsprm_get_velref, (setter)PyWcsprm_set_velref, (char *)doc_velref},
  /* {"wtb", (getter)PyWcsprm_get_wtb, NULL, (char *)doc_tab}, */
  {"zsource", (getter)PyWcsprm_get_zsource, (setter)PyWcsprm_set_zsource, (char *)doc_zsource},
  {NULL}
};

static PyMethodDef PyWcsprm_methods[] = {
  {"bounds_check", (PyCFunction)PyWcsprm_bounds_check, METH_VARARGS|METH_KEYWORDS, doc_bounds_check},
  {"cdfix", (PyCFunction)PyWcsprm_cdfix, METH_NOARGS, doc_cdfix},
  {"celfix", (PyCFunction)PyWcsprm_celfix, METH_NOARGS, doc_celfix},
  {"compare", (PyCFunction)PyWcsprm_compare, METH_VARARGS|METH_KEYWORDS, doc_compare},
  {"__copy__", (PyCFunction)PyWcsprm_copy, METH_NOARGS, doc_copy},
  {"cylfix", (PyCFunction)PyWcsprm_cylfix, METH_VARARGS|METH_KEYWORDS, doc_cylfix},
  {"datfix", (PyCFunction)PyWcsprm_datfix, METH_NOARGS, doc_datfix},
  {"__deepcopy__", (PyCFunction)PyWcsprm_copy, METH_O, doc_copy},
  {"fix", (PyCFunction)PyWcsprm_fix, METH_VARARGS|METH_KEYWORDS, doc_fix},
  {"get_cdelt", (PyCFunction)PyWcsprm_get_cdelt_func, METH_NOARGS, doc_get_cdelt},
  {"get_pc", (PyCFunction)PyWcsprm_get_pc_func, METH_NOARGS, doc_get_pc},
  {"get_ps", (PyCFunction)PyWcsprm_get_ps, METH_NOARGS, doc_get_ps},
  {"get_pv", (PyCFunction)PyWcsprm_get_pv, METH_NOARGS, doc_get_pv},
  {"has_cd", (PyCFunction)PyWcsprm_has_cdi_ja, METH_NOARGS, doc_has_cd},
  {"has_cdi_ja", (PyCFunction)PyWcsprm_has_cdi_ja, METH_NOARGS, doc_has_cdi_ja},
  {"has_crota", (PyCFunction)PyWcsprm_has_crotaia, METH_NOARGS, doc_has_crota},
  {"has_crotaia", (PyCFunction)PyWcsprm_has_crotaia, METH_NOARGS, doc_has_crotaia},
  {"has_pc", (PyCFunction)PyWcsprm_has_pci_ja, METH_NOARGS, doc_has_pc},
  {"has_pci_ja", (PyCFunction)PyWcsprm_has_pci_ja, METH_NOARGS, doc_has_pci_ja},
  {"is_unity", (PyCFunction)PyWcsprm_is_unity, METH_NOARGS, doc_is_unity},
  {"mix", (PyCFunction)PyWcsprm_mix, METH_VARARGS|METH_KEYWORDS, doc_mix},
  {"p2s", (PyCFunction)PyWcsprm_p2s, METH_VARARGS|METH_KEYWORDS, doc_p2s},
  {"print_contents", (PyCFunction)PyWcsprm_print_contents, METH_NOARGS, doc_print_contents},
  {"s2p", (PyCFunction)PyWcsprm_s2p, METH_VARARGS|METH_KEYWORDS, doc_s2p},
  {"set", (PyCFunction)PyWcsprm_set, METH_NOARGS, doc_set},
  {"set_ps", (PyCFunction)PyWcsprm_set_ps, METH_O, doc_set_ps},
  {"set_pv", (PyCFunction)PyWcsprm_set_pv, METH_O, doc_set_pv},
  {"spcfix", (PyCFunction)PyWcsprm_spcfix, METH_NOARGS, doc_spcfix},
  {"sptr", (PyCFunction)PyWcsprm_sptr, METH_VARARGS|METH_KEYWORDS, doc_sptr},
  {"sub", (PyCFunction)PyWcsprm_sub, METH_VARARGS|METH_KEYWORDS, doc_sub},
  {"to_header", (PyCFunction)PyWcsprm_to_header, METH_VARARGS|METH_KEYWORDS, doc_to_header},
  {"unitfix", (PyCFunction)PyWcsprm_unitfix, METH_VARARGS|METH_KEYWORDS, doc_unitfix},
  {NULL}
};

PyTypeObject PyWcsprmType = {
  #if PY3K
  PyVarObject_HEAD_INIT(NULL, 0)
  #else
  PyObject_HEAD_INIT(NULL)
  0,                            /*ob_size*/
  #endif
  "astropy.wcs.Wcsprm",              /*tp_name*/
  sizeof(PyWcsprm),             /*tp_basicsize*/
  0,                            /*tp_itemsize*/
  (destructor)PyWcsprm_dealloc, /*tp_dealloc*/
  0,                            /*tp_print*/
  0,                            /*tp_getattr*/
  0,                            /*tp_setattr*/
  0,                            /*tp_compare*/
  (reprfunc)PyWcsprm___str__,   /*tp_repr*/
  0,                            /*tp_as_number*/
  0,                            /*tp_as_sequence*/
  0,                            /*tp_as_mapping*/
  0,                            /*tp_hash */
  0,                            /*tp_call*/
  (reprfunc)PyWcsprm___str__,   /*tp_str*/
  0,                            /*tp_getattro*/
  0,                            /*tp_setattro*/
  0,                            /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  doc_Wcsprm,                   /* tp_doc */
  0,                            /* tp_traverse */
  0,                            /* tp_clear */
  PyWcsprm_richcompare,         /* tp_richcompare */
  0,                            /* tp_weaklistoffset */
  0,                            /* tp_iter */
  0,                            /* tp_iternext */
  PyWcsprm_methods,             /* tp_methods */
  0,                            /* tp_members */
  PyWcsprm_getset,              /* tp_getset */
  0,                            /* tp_base */
  0,                            /* tp_dict */
  0,                            /* tp_descr_get */
  0,                            /* tp_descr_set */
  0,                            /* tp_dictoffset */
  (initproc)PyWcsprm_init,      /* tp_init */
  0,                            /* tp_alloc */
  PyWcsprm_new,                 /* tp_new */
};

#define CONSTANT(a) PyModule_AddIntConstant(m, #a, a)

#define XSTRINGIFY(s) STRINGIFY(s)
#define STRINGIFY(s) #s

int
_setup_wcsprm_type(
    PyObject* m) {

  if (PyType_Ready(&PyWcsprmType) < 0) {
    return -1;
  }

  Py_INCREF(&PyWcsprmType);

  wcsprintf_set(NULL);
  wcserr_enable(1);

  return (
    PyModule_AddObject(m, "Wcsprm", (PyObject *)&PyWcsprmType) ||
    CONSTANT(WCSSUB_LONGITUDE) ||
    CONSTANT(WCSSUB_LATITUDE)  ||
    CONSTANT(WCSSUB_CUBEFACE)  ||
    CONSTANT(WCSSUB_SPECTRAL)  ||
    CONSTANT(WCSSUB_STOKES)    ||
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
    CONSTANT(WCSCOMPARE_ANCILLARY) ||
    CONSTANT(WCSCOMPARE_TILING) ||
    CONSTANT(WCSCOMPARE_CRPIX));
}
