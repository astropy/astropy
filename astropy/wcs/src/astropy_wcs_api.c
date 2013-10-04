#define NO_IMPORT_ARRAY

#include "astropy_wcs/astropy_wcs_api.h"

int
AstropyWcs_GetCVersion(void) {
  return REVISION;
}

void* AstropyWcs_API[] = {
  /*  0 */ (void *)AstropyWcs_GetCVersion,
  /* pyutil.h */
  /*  1 */ (void *)wcsprm_python2c,
  /*  2 */ (void *)wcsprm_c2python,
  /* distortion.h */
  /*  3 */ (void *)distortion_lookup_t_init,
  /*  4 */ (void *)distortion_lookup_t_free,
  /*  5 */ (void *)get_distortion_offset,
  /*  6 */ (void *)p4_pix2foc,
  /*  7 */ (void *)p4_pix2deltas,
  /* sip.h */
  /*  8 */ (void *)sip_clear,
  /*  9 */ (void *)sip_init,
  /* 10 */ (void *)sip_free,
  /* 11 */ (void *)sip_pix2foc,
  /* 12 */ (void *)sip_pix2deltas,
  /* 13 */ (void *)sip_foc2pix,
  /* 14 */ (void *)sip_foc2deltas,
  /* pipeline.h */
  /* 15 */ (void *)pipeline_clear,
  /* 16 */ (void *)pipeline_init,
  /* 17 */ (void *)pipeline_free,
  /* 18 */ (void *)pipeline_all_pixel2world,
  /* 19 */ (void *)pipeline_pix2foc,
  /* wcs.h */
  /* 20 */ (void *)wcsp2s,
  /* 21 */ (void *)wcss2p,
  /* 22 */ (void *)wcsprt,
  /* new for api version 2 */
  /* 23 */ (void *)wcslib_get_error_message,
  /* new for api version 3 */
  /* 24 */ (void *)wcsprintf_buf
};

int _setup_api(PyObject *m) {
  PyObject* c_api;

  #if PY_VERSION_HEX >= 0x03020000
    c_api = PyCapsule_New((void *)AstropyWcs_API, "_wcs._ASTROPY_WCS_API", NULL);
  #else
    c_api = PyCObject_FromVoidPtr((void *)AstropyWcs_API, NULL);
  #endif
  PyModule_AddObject(m, "_ASTROPY_WCS_API", c_api);

  return 0;
}
