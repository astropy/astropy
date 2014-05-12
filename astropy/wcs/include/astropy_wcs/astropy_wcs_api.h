#ifndef ASTROPY_WCS_API_H
#define ASTROPY_WCS_API_H

#include "wcsconfig.h"
#include "pyutil.h"
#include "distortion.h"
#include "pipeline.h"
#include "sip.h"
#include "wcs.h"
#include "wcsprintf.h"

/*
HOW TO UPDATE THE PUBLIC API

This code uses a table of function pointers to dynamically expose the
public API to other code that wants to use astropy.wcs from C.

Each function should be:

  1) Declared, as usual for C, in a .h file

  2) Defined in a .c file that is compiled as part of the _wcs.so file

  3) Have a macro that maps the function name to a position in the
     function table.  That macro should go in this file
     (astropy_wcs_api.h)

  4) An entry in the function table, which lives in astropy_wcs_api.c

Every time the function signatures change, or functions are added or
removed from the table, the value of REVISION should be incremented.
This allows for a rudimentary version check upon dynamic linking to
the astropy._wcs module.
 */

#define REVISION 4

#ifdef ASTROPY_WCS_BUILD

int _setup_api(PyObject* m);

#else

#if defined(NO_IMPORT_ASTROPY_WCS_API)
extern void** AstropyWcs_API;
#else
void** AstropyWcs_API;
#endif /* defined(NO_IMPORT_ASTROPY_PYWCS_API) */

/* Function macros that delegate to a function pointer in the AstropyWcs_API table */
#define AstropyWcs_GetCVersion (*(int (*)(void)) AstropyWcs_API[0])
#define wcsprm_python2c (*(void (*)(struct wcsprm*)) AstropyWcs_API[1])
#define wcsprm_c2python (*(void (*)(struct wcsprm*)) AstropyWcs_API[2])
#define distortion_lookup_t_init (*(int (*)(distortion_lookup_t* lookup)) AstropyWcs_API[3])
#define distortion_lookup_t_free (*(void (*)(distortion_lookup_t* lookup)) AstropyWcs_API[4])
#define get_distortion_offset (*(double (*)(const distortion_lookup_t*, const double* const)) AstropyWcs_API[5])
#define p4_pix2foc (*(int (*)(const unsigned int, const distortion_lookup_t**, const unsigned int, const double *, double *)) AstropyWcs_API[6])
#define p4_pix2deltas (*(int (*)(const unsigned int, const distortion_lookup_t**, const unsigned int, const double *, double *)) AstropyWcs_API[7])
#define sip_clear (*(void (*)(sip_t*) AstropyWcs_API[8]))
#define sip_init (*(int (*)(sip_t*, unsigned int, double*, unsigned int, double*, unsigned int, double*, unsigned int, double*, double*)) AstropyWcs_API[9])
#define sip_free (*(void (*)(sip_t*) AstropyWcs_API[10]))
#define sip_pix2foc (*(int (*)(sip_t*, unsigned int, unsigned int, double*, double*)) AstropyWcs_API[11])
#define sip_pix2deltas (*(int (*)(sip_t*, unsigned int, unsigned int, double*, double*)) AstropyWcs_API[12])
#define sip_foc2pix (*(int (*)(sip_t*, unsigned int, unsigned int, double*, double*)) AstropyWcs_API[13])
#define sip_foc2deltas (*(int (*)(sip_t*, unsigned int, unsigned int, double*, double*)) AstropyWcs_API[14])
#define pipeline_clear (*(void (*)(pipeline_t*)) AstropyWcs_API[15])
#define pipeline_init (*(void (*)(pipeline_t*, sip_t*, distortion_lookup_t**, struct wcsprm*)) AstropyWcs_API[16])
#define pipeline_free (*(void (*)(pipeline_t*)) AstropyWcs_API[17])
#define pipeline_all_pixel2world (*(int (*)(pipeline_t*, unsigned int, unsigned int, double*, double*)) AstropyWcs_API[18])
#define pipeline_pix2foc (*(int (*)(pipeline_t*, unsigned int, unsigned int, double*, double*)) AstropyWcs_API[19])
#define wcsp2s (*(int (*)(struct wcsprm *, int, int, const double[], double[], double[], double[], double[], int[])) AstropyWcs_API[20])
#define wcss2p (*(int (*)(struct wcsprm *, int, int, const double[], double[], double[], double[], double[], int[])) AstropyWcs_API[21])
#define wcsprt (*(int (*)(struct wcsprm *)) AstropyWcs_API[22])
#define wcslib_get_error_message (*(const char* (*)(int)) AstropyWcs_API[23])
#define wcsprintf_buf (*(const char * (*)()) AstropyWcs_API[24])

#ifndef NO_IMPORT_ASTROPY_WCS_API
int
import_astropy_wcs(void) {
  PyObject *wcs_module   = NULL;
  PyObject *c_api        = NULL;
  int       status       = -1;

  wcs_module = PyImport_ImportModule("astropy.wcs._wcs");
  if (wcs_module == NULL) goto exit;

  c_api = PyObject_GetAttrString(wcs_module, "_ASTROPY_WCS_API");
  if (c_api == NULL) goto exit;

  #if PY_VERSION_HEX >= 0x03020000
    AstropyWcs_API = (void **)PyCapsule_GetPointer(c_api, "_wcs._ASTROPY_WCS_API");
    if (AstropyWcs_API == NULL)
        goto exit;
  #else
    if (PyCObject_Check(c_api)) {
      AstropyWcs_API = (void **)PyCObject_AsVoidPtr(c_api);
    } else {
      goto exit;
    }
  #endif

  /* Perform runtime check of C API version */
  if (REVISION != AstropyWcs_GetCVersion()) {
    PyErr_Format(
                 PyExc_ImportError, "module compiled against "        \
                 "ABI version '%x' but this version of astropy.wcs is '%x'", \
                 (int)REVISION, (int)AstropyWcs_GetCVersion());
    return -1;
  }

 exit:
  Py_XDECREF(wcs_module);
  Py_XDECREF(c_api);

  return status;
}

#endif /* !defined(NO_IMPORT_ASTROPY_WCS_API) */

#endif /* ASTROPY_WCS_BUILD */

#endif /* ASTROPY_WCS_API_H */
