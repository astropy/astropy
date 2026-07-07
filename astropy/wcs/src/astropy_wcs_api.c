#define NO_IMPORT_ARRAY

#include "astropy_wcs/astropy_wcs_api.h"

int
AstropyWcs_GetCVersion(void) {
  return REVISION;
}

/*
 * Deprecated members of the public astropy.wcs C API.
 *
 * Of the symbols exported through the AstropyWcs_API function-pointer table
 * (and discoverable by downstream C/Cython code via astropy.wcs.get_include()),
 * only six are known to be used by a downstream package (drizzlepac):
 * wcsprm_python2c, wcsprm_c2python, pipeline_all_pixel2world, wcss2p, wcsprt
 * and wcslib_get_error_message.  wcsp2s is also kept, for symmetry with wcss2p.
 *
 * The remaining members are deprecated.  Rather than remove them (which would
 * break the ABI and the table layout that the REVISION check relies on), each
 * deprecated slot points at a thin wrapper that emits a DeprecationWarning via
 * numpy's DEPRECATE macro before forwarding to the real function.  The wrappers
 * deliberately wrap only the *exported* table entry, not the underlying
 * functions, which astropy.wcs itself still calls internally by name and which
 * must not warn.
 *
 * Note: DEPRECATE expands to PyErr_WarnEx, so the GIL must be held when a
 * deprecated entry is called through the table.  If the warning is escalated to
 * an exception (for example under -W error) the exception is left set and
 * surfaces at the next Python API boundary.
 */

#ifndef DEPRECATE
#define DEPRECATE(msg) PyErr_WarnEx(PyExc_DeprecationWarning, msg, 1)
#endif

#define ASTROPY_WCS_DEPRECATE(name)                                     \
  DEPRECATE(name " is a deprecated part of the public astropy.wcs C "   \
                 "API and will be removed in a future version of astropy")

/* distortion.h */

static int
deprecated_distortion_lookup_t_init(distortion_lookup_t* lookup) {
  ASTROPY_WCS_DEPRECATE("distortion_lookup_t_init");
  return distortion_lookup_t_init(lookup);
}

static void
deprecated_distortion_lookup_t_free(distortion_lookup_t* lookup) {
  ASTROPY_WCS_DEPRECATE("distortion_lookup_t_free");
  distortion_lookup_t_free(lookup);
}

static double
deprecated_get_distortion_offset(
    const distortion_lookup_t * const lookup,
    const double * const img) {
  ASTROPY_WCS_DEPRECATE("get_distortion_offset");
  return get_distortion_offset(lookup, img);
}

static int
deprecated_p4_pix2foc(
    const unsigned int naxes,
    const distortion_lookup_t** lookups,
    const unsigned int nelem,
    const double* pix,
    double* foc) {
  ASTROPY_WCS_DEPRECATE("p4_pix2foc");
  return p4_pix2foc(naxes, lookups, nelem, pix, foc);
}

static int
deprecated_p4_pix2deltas(
    const unsigned int naxes,
    const distortion_lookup_t** lookups,
    const unsigned int nelem,
    const double* pix,
    double* foc) {
  ASTROPY_WCS_DEPRECATE("p4_pix2deltas");
  return p4_pix2deltas(naxes, lookups, nelem, pix, foc);
}

/* sip.h */

static void
deprecated_sip_clear(sip_t* sip) {
  ASTROPY_WCS_DEPRECATE("sip_clear");
  sip_clear(sip);
}

static int
deprecated_sip_init(
    sip_t* sip,
    const unsigned int a_order, const double* a,
    const unsigned int b_order, const double* b,
    const unsigned int ap_order, const double* ap,
    const unsigned int bp_order, const double* bp,
    const double* crpix) {
  ASTROPY_WCS_DEPRECATE("sip_init");
  return sip_init(sip, a_order, a, b_order, b, ap_order, ap, bp_order, bp,
                  crpix);
}

static void
deprecated_sip_free(sip_t* sip) {
  ASTROPY_WCS_DEPRECATE("sip_free");
  sip_free(sip);
}

static int
deprecated_sip_pix2foc(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* pix,
    double* foc) {
  ASTROPY_WCS_DEPRECATE("sip_pix2foc");
  return sip_pix2foc(sip, naxes, nelem, pix, foc);
}

static int
deprecated_sip_pix2deltas(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* pix,
    double* foc) {
  ASTROPY_WCS_DEPRECATE("sip_pix2deltas");
  return sip_pix2deltas(sip, naxes, nelem, pix, foc);
}

static int
deprecated_sip_foc2pix(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* foc,
    double* pix) {
  ASTROPY_WCS_DEPRECATE("sip_foc2pix");
  return sip_foc2pix(sip, naxes, nelem, foc, pix);
}

static int
deprecated_sip_foc2deltas(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* foc,
    double* deltas) {
  ASTROPY_WCS_DEPRECATE("sip_foc2deltas");
  return sip_foc2deltas(sip, naxes, nelem, foc, deltas);
}

/* pipeline.h */

static void
deprecated_pipeline_clear(pipeline_t* pipeline) {
  ASTROPY_WCS_DEPRECATE("pipeline_clear");
  pipeline_clear(pipeline);
}

static void
deprecated_pipeline_init(
    pipeline_t* pipeline,
    distortion_lookup_t** det2im,
    sip_t* sip,
    distortion_lookup_t** cpdis,
    struct wcsprm* wcs) {
  ASTROPY_WCS_DEPRECATE("pipeline_init");
  pipeline_init(pipeline, det2im, sip, cpdis, wcs);
}

static void
deprecated_pipeline_free(pipeline_t* pipeline) {
  ASTROPY_WCS_DEPRECATE("pipeline_free");
  pipeline_free(pipeline);
}

static int
deprecated_pipeline_pix2foc(
    pipeline_t* pipeline,
    const unsigned int ncoord,
    const unsigned int nelem,
    const double* const pixcrd,
    double* foc) {
  ASTROPY_WCS_DEPRECATE("pipeline_pix2foc");
  return pipeline_pix2foc(pipeline, ncoord, nelem, pixcrd, foc);
}

/* wcsprintf.h */

static const char *
deprecated_wcsprintf_buf(void) {
  ASTROPY_WCS_DEPRECATE("wcsprintf_buf");
  return wcsprintf_buf();
}

void* AstropyWcs_API[] = {
  /*  0 */ (void *)AstropyWcs_GetCVersion,
  /* pyutil.h */
  /*  1 */ (void *)wcsprm_python2c,
  /*  2 */ (void *)wcsprm_c2python,
  /* distortion.h */
  /*  3 */ (void *)deprecated_distortion_lookup_t_init,
  /*  4 */ (void *)deprecated_distortion_lookup_t_free,
  /*  5 */ (void *)deprecated_get_distortion_offset,
  /*  6 */ (void *)deprecated_p4_pix2foc,
  /*  7 */ (void *)deprecated_p4_pix2deltas,
  /* sip.h */
  /*  8 */ (void *)deprecated_sip_clear,
  /*  9 */ (void *)deprecated_sip_init,
  /* 10 */ (void *)deprecated_sip_free,
  /* 11 */ (void *)deprecated_sip_pix2foc,
  /* 12 */ (void *)deprecated_sip_pix2deltas,
  /* 13 */ (void *)deprecated_sip_foc2pix,
  /* 14 */ (void *)deprecated_sip_foc2deltas,
  /* pipeline.h */
  /* 15 */ (void *)deprecated_pipeline_clear,
  /* 16 */ (void *)deprecated_pipeline_init,
  /* 17 */ (void *)deprecated_pipeline_free,
  /* 18 */ (void *)pipeline_all_pixel2world,
  /* 19 */ (void *)deprecated_pipeline_pix2foc,
  /* wcs.h */
  /* 20 */ (void *)wcsp2s,
  /* 21 */ (void *)wcss2p,
  /* 22 */ (void *)wcsprt,
  /* new for api version 2 */
  /* 23 */ (void *)wcslib_get_error_message,
  /* new for api version 3 */
  /* 24 */ (void *)deprecated_wcsprintf_buf
};

int _setup_api(PyObject *m) {
  PyObject* c_api;

  c_api = PyCapsule_New((void *)AstropyWcs_API, "_wcs._ASTROPY_WCS_API", NULL);
  PyModule_AddObject(m, "_ASTROPY_WCS_API", c_api);

  return 0;
}
