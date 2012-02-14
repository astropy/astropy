#include <Python.h>

/***************************************************************************
 * Macros for determining the compiler version.
 *
 * These are borrowed from boost, and majorly abridged to include only
 * the compilers we care about.
 ***************************************************************************/

#ifndef PY3K
#if PY_MAJOR_VERSION >= 3
#define PY3K 1
#else
#define PY3K 0
#endif
#endif


#define STRINGIZE(X) DO_STRINGIZE(X)
#define DO_STRINGIZE(X) #X

#if defined __clang__
/*  Clang C++ emulates GCC, so it has to appear early. */
#    define COMPILER "Clang version " __clang_version__

#elif defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC) || defined(__ECC)
/* Intel */
#    if defined(__INTEL_COMPILER)
#        define INTEL_VERSION __INTEL_COMPILER
#    elif defined(__ICL)
#        define INTEL_VERSION __ICL
#    elif defined(__ICC)
#        define INTEL_VERSION __ICC
#    elif defined(__ECC)
#        define INTEL_VERSION __ECC
#    endif
#    define COMPILER "Intel C compiler version " STRINGIZE(INTEL_VERSION)

#elif defined(__GNUC__)
/* gcc */
#    define COMPILER "GCC version " __VERSION__

#elif defined(__SUNPRO_CC)
/* Sun Workshop Compiler */
#    define COMPILER "Sun compiler version " STRINGIZE(__SUNPRO_CC)

#elif defined(_MSC_VER)
/* Microsoft Visual C/C++
   Must be last since other compilers define _MSC_VER for compatibility as well */
#    if _MSC_VER < 1200
#        define COMPILER_VERSION 5.0
#    elif _MSC_VER < 1300
#        define COMPILER_VERSION 6.0
#    elif _MSC_VER == 1300
#        define COMPILER_VERSION 7.0
#    elif _MSC_VER == 1310
#        define COMPILER_VERSION 7.1
#    elif _MSC_VER == 1400
#        define COMPILER_VERSION 8.0
#    elif _MSC_VER == 1500
#        define COMPILER_VERSION 9.0
#    elif _MSC_VER == 1600
#        define COMPILER_VERSION 10.0
#    else
#        define COMPILER_VERSION _MSC_VER
#    endif
#    define COMPILER "Microsoft Visual C++ version " STRINGIZE(COMPILER_VERSION)

#else
/* Fallback */
#    define COMPILER "Unknown compiler"

#endif


/***************************************************************************
 * Module-level
 ***************************************************************************/

struct module_state {
/* The Sun compiler can't handle empty structs */
#if defined(__SUNPRO_C) || defined(_MSC_VER)
    int _dummy;
#endif
};

#if PY3K
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_compiler",
        NULL,
        sizeof(struct module_state),
        NULL,
        NULL,
        NULL,
        NULL,
        NULL
    };

    #define INITERROR return NULL

    PyMODINIT_FUNC
    PyInit__compiler(void)

#else
    #define INITERROR return

    PyMODINIT_FUNC
    init_compiler(void)
#endif

{
  PyObject* m;

#if PY3K
  m = PyModule_Create(&moduledef);
#else
  m = Py_InitModule3("_compiler", NULL, NULL);
#endif

  if (m == NULL)
    INITERROR;

  PyModule_AddStringConstant(m, "compiler", COMPILER);

#if PY3K
  return m;
#endif
}
