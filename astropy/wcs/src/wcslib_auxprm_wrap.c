#define NO_IMPORT_ARRAY

#include "astropy_wcs/wcslib_auxprm_wrap.h"

#include <wcs.h>
#include <wcsprintf.h>
#include <tab.h>

/*
 It gets to be really tedious to type long docstrings in ANSI C syntax
 (since multi-line strings literals are not valid).  Therefore, the
 docstrings are written in doc/docstrings.py, which are then converted
 by setup.py into docstrings.h, which we include here.
*/
#include "astropy_wcs/docstrings.h"


/***************************************************************************
 * PyAuxprm methods                                                        *
 ***************************************************************************/

static PyObject*
PyAuxprm_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
  PyAuxprm* self;

  allocfunc alloc_func = PyType_GetSlot(type, Py_tp_alloc);
  self = (PyAuxprm*)alloc_func(type, 0);
  return (PyObject*)self;
}


static int
PyAuxprm_traverse(PyAuxprm* self, visitproc visit, void *arg) {
  Py_VISIT(self->owner);
  Py_VISIT((PyObject*)Py_TYPE((PyObject*)self));
  return 0;
}


static int
PyAuxprm_clear(PyAuxprm* self) {
  Py_CLEAR(self->owner);
  return 0;
}


static void PyAuxprm_dealloc(PyAuxprm* self) {
  PyAuxprm_clear(self);
  PyTypeObject *tp = Py_TYPE((PyObject*)self);
  freefunc free_func = PyType_GetSlot(tp, Py_tp_free);
  free_func((PyObject*)self);
  Py_DECREF(tp);
}


PyAuxprm* PyAuxprm_cnew(PyObject* wcsprm, struct auxprm* x) {
  PyAuxprm* self;
  PyTypeObject* type = (PyTypeObject*)PyAuxprmType;
  allocfunc alloc_func = PyType_GetSlot(type, Py_tp_alloc);
  self = (PyAuxprm*)alloc_func(type, 0);
  if (self == NULL) return NULL;
  self->x = x;
  Py_INCREF(wcsprm);
  self->owner = wcsprm;
  return self;
}


static void auxprmprt(const struct auxprm *aux) {

  if (aux == 0x0) return;

  wcsprintf("rsun_ref:");
  if (aux->rsun_ref != UNDEFINED) wcsprintf(" %f", aux->rsun_ref);
  wcsprintf("\ndsun_obs:");
  if (aux->dsun_obs != UNDEFINED) wcsprintf(" %f", aux->dsun_obs);
  wcsprintf("\ncrln_obs:");
  if (aux->crln_obs != UNDEFINED) wcsprintf(" %f", aux->crln_obs);
  wcsprintf("\nhgln_obs:");
  if (aux->hgln_obs != UNDEFINED) wcsprintf(" %f", aux->hgln_obs);
  wcsprintf("\nhglt_obs:");
  if (aux->hglt_obs != UNDEFINED) wcsprintf(" %f", aux->hglt_obs);
  wcsprintf("\na_radius:");
  if (aux->a_radius != UNDEFINED) wcsprintf(" %f", aux->a_radius);
  wcsprintf("\nb_radius:");
  if (aux->b_radius != UNDEFINED) wcsprintf(" %f", aux->b_radius);
  wcsprintf("\nc_radius:");
  if (aux->c_radius != UNDEFINED) wcsprintf(" %f", aux->c_radius);
  wcsprintf("\nbdis_obs:");
  if (aux->bdis_obs != UNDEFINED) wcsprintf(" %f", aux->bdis_obs);
  wcsprintf("\nblon_obs:");
  if (aux->blon_obs != UNDEFINED) wcsprintf(" %f", aux->blon_obs);
  wcsprintf("\nblat_obs:");
  if (aux->blat_obs != UNDEFINED) wcsprintf(" %f", aux->blat_obs);
  return;
}


static PyObject* PyAuxprm___str__(PyAuxprm* self) {
  /* This is not thread-safe, but since we're holding onto the GIL,
     we can assume we won't have thread conflicts */
  wcsprintf_set(NULL);
  auxprmprt(self->x);
  return PyUnicode_FromString(wcsprintf_buf());
}


/***************************************************************************
 * Member getters/setters (properties)
 */

static PyObject* PyAuxprm_get_rsun_ref(PyAuxprm* self, void* closure) {
  if(self->x == NULL || self->x->rsun_ref == UNDEFINED) {
    Py_RETURN_NONE;
  } else {
    return get_double("rsun_ref", self->x->rsun_ref);
  }
}

static int PyAuxprm_set_rsun_ref(PyAuxprm* self, PyObject* value, void* closure) {
  if(self->x == NULL) {
    return -1;
  } else if (value == Py_None) {
    self->x->rsun_ref = UNDEFINED;
    return 0;
  } else {
    return set_double("rsun_ref", value, &self->x->rsun_ref);
  }
}

static PyObject* PyAuxprm_get_dsun_obs(PyAuxprm* self, void* closure) {
  if(self->x == NULL || self->x->dsun_obs == UNDEFINED) {
    Py_RETURN_NONE;
  } else {
    return get_double("dsun_obs", self->x->dsun_obs);
  }
}

static int PyAuxprm_set_dsun_obs(PyAuxprm* self, PyObject* value, void* closure) {
  if(self->x == NULL) {
    return -1;
  } else if (value == Py_None) {
    self->x->dsun_obs = UNDEFINED;
    return 0;
  } else {
    return set_double("dsun_obs", value, &self->x->dsun_obs);
  }
}

static PyObject* PyAuxprm_get_crln_obs(PyAuxprm* self, void* closure) {
  if(self->x == NULL || self->x->crln_obs == UNDEFINED) {
    Py_RETURN_NONE;
  } else {
    return get_double("crln_obs", self->x->crln_obs);
  }
}

static int PyAuxprm_set_crln_obs(PyAuxprm* self, PyObject* value, void* closure) {
  if(self->x == NULL) {
    return -1;
  } else if (value == Py_None) {
    self->x->crln_obs = UNDEFINED;
    return 0;
  } else {
    return set_double("crln_obs", value, &self->x->crln_obs);
  }
}

static PyObject* PyAuxprm_get_hgln_obs(PyAuxprm* self, void* closure) {
  if(self->x == NULL || self->x->hgln_obs == UNDEFINED) {
    Py_RETURN_NONE;
  } else {
    return get_double("hgln_obs", self->x->hgln_obs);
  }
}

static int PyAuxprm_set_hgln_obs(PyAuxprm* self, PyObject* value, void* closure) {
  if(self->x == NULL) {
    return -1;
  } else if (value == Py_None) {
    self->x->hgln_obs = UNDEFINED;
    return 0;
  } else {
    return set_double("hgln_obs", value, &self->x->hgln_obs);
  }
}

static PyObject* PyAuxprm_get_hglt_obs(PyAuxprm* self, void* closure) {
  if(self->x == NULL || self->x->hglt_obs == UNDEFINED) {
    Py_RETURN_NONE;
  } else {
    return get_double("hglt_obs", self->x->hglt_obs);
  }
}

static int PyAuxprm_set_hglt_obs(PyAuxprm* self, PyObject* value, void* closure) {
  if(self->x == NULL) {
    return -1;
  } else if (value == Py_None) {
    self->x->hglt_obs = UNDEFINED;
    return 0;
  } else {
    return set_double("hglt_obs", value, &self->x->hglt_obs);
  }
}

static PyObject* PyAuxprm_get_a_radius(PyAuxprm* self, void* closure) {
  if(self->x == NULL || self->x->a_radius == UNDEFINED) {
    Py_RETURN_NONE;
  } else {
    return get_double("a_radius", self->x->a_radius);
  }
}

static int PyAuxprm_set_a_radius(PyAuxprm* self, PyObject* value, void* closure) {
  if(self->x == NULL) {
    return -1;
  } else if (value == Py_None) {
    self->x->a_radius = UNDEFINED;
    return 0;
  } else {
    return set_double("a_radius", value, &self->x->a_radius);
  }
}

static PyObject* PyAuxprm_get_b_radius(PyAuxprm* self, void* closure) {
  if(self->x == NULL || self->x->b_radius == UNDEFINED) {
    Py_RETURN_NONE;
  } else {
    return get_double("b_radius", self->x->b_radius);
  }
}

static int PyAuxprm_set_b_radius(PyAuxprm* self, PyObject* value, void* closure) {
  if(self->x == NULL) {
    return -1;
  } else if (value == Py_None) {
    self->x->b_radius = UNDEFINED;
    return 0;
  } else {
    return set_double("b_radius", value, &self->x->b_radius);
  }
}

static PyObject* PyAuxprm_get_c_radius(PyAuxprm* self, void* closure) {
  if(self->x == NULL || self->x->c_radius == UNDEFINED) {
    Py_RETURN_NONE;
  } else {
    return get_double("c_radius", self->x->c_radius);
  }
}

static int PyAuxprm_set_c_radius(PyAuxprm* self, PyObject* value, void* closure) {
  if(self->x == NULL) {
    return -1;
  } else if (value == Py_None) {
    self->x->c_radius = UNDEFINED;
    return 0;
  } else {
    return set_double("c_radius", value, &self->x->c_radius);
  }
}

static PyObject* PyAuxprm_get_bdis_obs(PyAuxprm* self, void* closure) {
  if(self->x == NULL || self->x->bdis_obs == UNDEFINED) {
    Py_RETURN_NONE;
  } else {
    return get_double("bdis_obs", self->x->bdis_obs);
  }
}

static int PyAuxprm_set_bdis_obs(PyAuxprm* self, PyObject* value, void* closure) {
  if(self->x == NULL) {
    return -1;
  } else if (value == Py_None) {
    self->x->bdis_obs = UNDEFINED;
    return 0;
  } else {
    return set_double("bdis_obs", value, &self->x->bdis_obs);
  }
}

static PyObject* PyAuxprm_get_blon_obs(PyAuxprm* self, void* closure) {
  if(self->x == NULL || self->x->blon_obs == UNDEFINED) {
    Py_RETURN_NONE;
  } else {
    return get_double("blon_obs", self->x->blon_obs);
  }
}

static int PyAuxprm_set_blon_obs(PyAuxprm* self, PyObject* value, void* closure) {
  if(self->x == NULL) {
    return -1;
  } else if (value == Py_None) {
    self->x->blon_obs = UNDEFINED;
    return 0;
  } else {
    return set_double("blon_obs", value, &self->x->blon_obs);
  }
}

static PyObject* PyAuxprm_get_blat_obs(PyAuxprm* self, void* closure) {
  if(self->x == NULL || self->x->blat_obs == UNDEFINED) {
    Py_RETURN_NONE;
  } else {
    return get_double("blat_obs", self->x->blat_obs);
  }
}

static int PyAuxprm_set_blat_obs(PyAuxprm* self, PyObject* value, void* closure) {
  if(self->x == NULL) {
    return -1;
  } else if (value == Py_None) {
    self->x->blat_obs = UNDEFINED;
    return 0;
  } else {
    return set_double("blat_obs", value, &self->x->blat_obs);
  }
}


/***************************************************************************
 * PyAuxprm definition structures
 */

static PyGetSetDef PyAuxprm_getset[] = {
  {"rsun_ref", (getter)PyAuxprm_get_rsun_ref, (setter)PyAuxprm_set_rsun_ref, (char *)doc_rsun_ref},
  {"dsun_obs", (getter)PyAuxprm_get_dsun_obs, (setter)PyAuxprm_set_dsun_obs, (char *)doc_dsun_obs},
  {"crln_obs", (getter)PyAuxprm_get_crln_obs, (setter)PyAuxprm_set_crln_obs, (char *)doc_crln_obs},
  {"hgln_obs", (getter)PyAuxprm_get_hgln_obs, (setter)PyAuxprm_set_hgln_obs, (char *)doc_hgln_obs},
  {"hglt_obs", (getter)PyAuxprm_get_hglt_obs, (setter)PyAuxprm_set_hglt_obs, (char *)doc_hglt_obs},
  {"a_radius", (getter)PyAuxprm_get_a_radius, (setter)PyAuxprm_set_a_radius, (char *)doc_a_radius},
  {"b_radius", (getter)PyAuxprm_get_b_radius, (setter)PyAuxprm_set_b_radius, (char *)doc_b_radius},
  {"c_radius", (getter)PyAuxprm_get_c_radius, (setter)PyAuxprm_set_c_radius, (char *)doc_c_radius},
  {"bdis_obs", (getter)PyAuxprm_get_bdis_obs, (setter)PyAuxprm_set_bdis_obs, (char *)doc_bdis_obs},
  {"blon_obs", (getter)PyAuxprm_get_blon_obs, (setter)PyAuxprm_set_blon_obs, (char *)doc_blon_obs},
  {"blat_obs", (getter)PyAuxprm_get_blat_obs, (setter)PyAuxprm_set_blat_obs, (char *)doc_blat_obs},
  {NULL}
};

PyType_Spec PyAuxprmType_spec = {
  .name = "astropy.wcs.Auxprm",
  .basicsize = sizeof(PyAuxprm),
  .itemsize = 0,
  .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_IMMUTABLETYPE,
  .slots = (PyType_Slot[]) {
    {Py_tp_dealloc, (destructor)PyAuxprm_dealloc},
    {Py_tp_str, (reprfunc)PyAuxprm___str__},
    {Py_tp_doc, doc_Auxprm},
    {Py_tp_traverse, (traverseproc)PyAuxprm_traverse},
    {Py_tp_clear, (inquiry)PyAuxprm_clear},
    {Py_tp_getset, PyAuxprm_getset},
    // FIXME: this seems logical but this slot wasn't previously set
    // maybe a mistake from https://github.com/astropy/astropy/pull/10333 ?
    // {Py_tp_new, (void*)PyAuxprm_new},
    {0, NULL}
  },
};

PyObject* PyAuxprmType = NULL;

int
_setup_auxprm_type(PyObject* m) {
  PyAuxprmType = PyType_FromSpec(&PyAuxprmType_spec);
  if (PyAuxprmType == NULL) {
    return -1;
  }

  PyModule_AddObject(m, "Auxprm", PyAuxprmType);

  return 0;
}
