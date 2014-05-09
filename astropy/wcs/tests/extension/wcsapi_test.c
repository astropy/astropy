#include "Python.h"
#include "astropy_wcs_api.h"

static PyObject*
test(
    PyObject* self,
    PyObject* args,
    PyObject* kwds) {

    pipeline_t p;

    pipeline_clear(&p);

    return Py_None;
}

static PyMethodDef module_methods[] = {
  {"test", (PyCFunction)test, METH_NOARGS, ""},
  {NULL}  /* Sentinel */
};

struct module_state {
/* The Sun compiler can't handle empty structs */
#if defined(__SUNPRO_C) || defined(_MSC_VER)
    int _dummy;
#endif
};

#if PY3K
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "wcsapi_test",
        NULL,
        sizeof(struct module_state),
        module_methods,
        NULL,
        NULL,
        NULL,
        NULL
    };

    #define INITERROR return NULL

    PyMODINIT_FUNC
    PyInit_wcsapi_test(void)

#else
    #define INITERROR return

    PyMODINIT_FUNC
    initwcsapi_test(void)
#endif

{
  PyObject* m;

#if PY3K
  m = PyModule_Create(&moduledef);
#else
  m = Py_InitModule3("wcsapi_test", module_methods, NULL);
#endif

  if (m == NULL) {
      printf("HERE\n");
    INITERROR;
  }

  import_astropy_wcs();

#if PY3K
  if (PyErr_Occurred())
      return NULL;
  else
      return m;
#endif
}
