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
void PyInit_lib_convolve(void);
#endif

#include "numpy/ndarrayobject.h"
#define DTYPE npy_float64

#endif
