#ifndef CONVOLVE_INCLUDE
#define CONVOLVE_INCLUDE

#include <stddef.h>

// Forcibly disable OpenMP support at the src level
#undef _OPENMP

#if defined(_MSC_VER)

#define FORCE_INLINE  __forceinline
#define NEVER_INLINE  __declspec(noinline)

// Other compilers (including GCC & Clang)
#else

#define FORCE_INLINE inline __attribute__((always_inline))
#define NEVER_INLINE __attribute__((noinline))

#endif

// MSVC implements OpenMP 2.0 which mandates singed integers for its parallel loops
#if defined(_MSC_VER)
typedef signed omp_iter_var;
#else
typedef size_t omp_iter_var;
#endif

// MSVC exports
#if defined(_MSC_VER)
#define LIB_CONVOLVE_EXPORT __declspec(dllexport)
#else
#define LIB_CONVOLVE_EXPORT // nothing
#endif

// Distutils on Windows will automatically exports ``PyInit_lib_convolve``,
// create dummy to prevent linker complaining about missing symbol.
#if defined(_MSC_VER)
void PyInit__convolve(void);
#endif

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define NO_IMPORT_ARRAY
#include "numpy/ndarrayobject.h"
#define DTYPE npy_float64


LIB_CONVOLVE_EXPORT void convolveNd_c(DTYPE * const result,
        const DTYPE * const f,
        const unsigned n_dim,
        const size_t * const image_shape,
        const DTYPE * const g,
        const size_t * const kernel_shape,
        const bool nan_interpolate,
        const bool embed_result_within_padded_region,
        const unsigned n_threads);

// 1D
void convolve1d_c(DTYPE * const result,
        const DTYPE * const f, const size_t nx,
        const DTYPE * const g, const size_t nkx,
        const bool nan_interpolate,
        const bool embed_result_within_padded_region,
        const unsigned n_threads);
FORCE_INLINE void convolve1d(DTYPE * const result,
        const DTYPE * const f, const size_t nx,
        const DTYPE * const g, const size_t nkx,
        const bool nan_interpolate,
        const bool embed_result_within_padded_region,
        const unsigned n_threads);

// 2D
void convolve2d_c(DTYPE * const result,
        const DTYPE * const f, const size_t nx, const size_t ny,
        const DTYPE * const g, const size_t nkx, const size_t nky,
        const bool nan_interpolate,
        const bool embed_result_within_padded_region,
        const unsigned n_threads);
FORCE_INLINE void convolve2d(DTYPE * const result,
        const DTYPE * const f, const size_t nx, const size_t ny,
        const DTYPE * const g, const size_t nkx, const size_t nky,
        const bool nan_interpolate,
        const bool embed_result_within_padded_region,
        const unsigned n_threads);

// 3D
void convolve3d_c(DTYPE * const result,
        const DTYPE * const f, const size_t nx, const size_t ny, const size_t nz,
        const DTYPE * const g, const size_t nkx, const size_t nky, const size_t nkz,
        const bool nan_interpolate,
        const bool embed_result_within_padded_region,
        const unsigned n_threads);
FORCE_INLINE void convolve3d(DTYPE * const result,
        const DTYPE * const f, const size_t nx, const size_t ny, const size_t nz,
        const DTYPE * const g, const size_t nkx, const size_t nky, const size_t nkz,
        const bool nan_interpolate,
        const bool embed_result_within_padded_region,
        const unsigned n_threads);


#endif
