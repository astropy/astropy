// Licensed under a 3-clause BSD style license - see LICENSE.rst

/*----------------------------- WARNING! -----------------------------
 * The C functions below are NOT designed to be called externally to
 * the Python function astropy/astropy/convolution/convolve.py.
 * They do NOT include any of the required correct usage checking.
 *
 *------------------------------- NOTES ------------------------------
 *
 * The simplest implementation of convolution does not deal with any boundary
 * treatment, and pixels within half a kernel width of the edge of the image are
 * set to zero. In cases where a boundary mode is set, we pad the input array in
 * the Python code. In the 1D case, this means that the input array to the C
 * code has a size nx + nkx where nx is the original array size and nkx is the
 * size of the kernel. If we also padded the results array, then we could use
 * the exact same C code for the convolution, provided that the results array
 * was 'unpadded' in the Python code after the C code.
 *
 * However, to avoid needlessly padding the results array, we instead adjust the
 * index when accessing the results array - for example in the 1D case we shift
 * the index in the results array compared to the input array by half the kernel
 * size. This is done via the 'result_index' variable, and this behavior is
 * triggered by the 'embed_result_within_padded_region' setting.
 *
 */


#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stddef.h>

#include "convolve.h"

#ifdef _OPENMP
#include <omp.h>
#endif


void convolveNd_c(DTYPE * const result,
        const DTYPE * const f,
        const unsigned n_dim,
        const size_t * const image_shape,
        const DTYPE * const g,
        const size_t * const kernel_shape,
        const bool nan_interpolate,
        const bool embed_result_within_padded_region,
        const unsigned n_threads)
{
#ifdef NDEBUG
    if (!result || !f || !g || !image_shape || !kernel_shape)
        return;
#else
    assert(result);
    assert(f);
    assert(g);
    assert(image_shape);
    assert(kernel_shape);
#endif

    if (n_dim == 1)
        convolve1d_c(result, f,
                image_shape[0],
                g, kernel_shape[0],
                nan_interpolate,
                embed_result_within_padded_region,
                n_threads);
    else if (n_dim == 2)
        convolve2d_c(result, f,
                image_shape[0], image_shape[1],
                g, kernel_shape[0], kernel_shape[1],
                nan_interpolate,
                embed_result_within_padded_region,
                n_threads);
    else if (n_dim == 3)
        convolve3d_c(result, f,
                        image_shape[0], image_shape[1], image_shape[2],
                        g, kernel_shape[0], kernel_shape[1], kernel_shape[2],
                        nan_interpolate,
                        embed_result_within_padded_region,
                        n_threads);
    else
        assert(0); // Unimplemented: n_dim > 3
}

/*-------------------------PERFORMANCE NOTES--------------------------------
 * The function wrappers below are designed to take advantage of the following:
 * The preprocessor will inline convolve<N>d(), effectively
 * expanding the two logical branches, replacing nan_interpolate
 * for their literal equivalents. The corresponding conditionals
 * within these functions will then be optimized away, this
 * being the goal - removing the unnecessary conditionals from
 * the loops without duplicating code.
 *--------------------------------------------------------------------------
 */

void convolve1d_c(DTYPE * const result,
        const DTYPE * const f, const size_t nx,
        const DTYPE * const g, const size_t nkx,
        const bool nan_interpolate,
        const bool embed_result_within_padded_region,
        const unsigned n_threads)
{
#ifdef NDEBUG
    if (!result || !f || !g)
        return;
#else
    assert(result);
    assert(f);
    assert(g);
#endif

    if (nan_interpolate) {
      if (embed_result_within_padded_region)
        convolve1d(result, f, nx, g, nkx, true, true, n_threads);
      else
        convolve1d(result, f, nx, g, nkx, true, false, n_threads);
    } else {
      if (embed_result_within_padded_region)
        convolve1d(result, f, nx, g, nkx, false, true, n_threads);
      else
        convolve1d(result, f, nx, g, nkx, false, false, n_threads);
    }
}

void convolve2d_c(DTYPE * const result,
        const DTYPE * const f, const size_t nx, const size_t ny,
        const DTYPE * const g, const size_t nkx, const size_t nky,
        const bool nan_interpolate,
        const bool embed_result_within_padded_region,
        const unsigned n_threads)
{
#ifdef NDEBUG
    if (!result || !f || !g)
        return;
#else
    assert(result);
    assert(f);
    assert(g);
#endif

    if (nan_interpolate) {
      if (embed_result_within_padded_region)
        convolve2d(result, f, nx, ny, g, nkx, nky, true, true, n_threads);
      else
        convolve2d(result, f, nx, ny, g, nkx, nky, true, false, n_threads);
    } else {
      if (embed_result_within_padded_region)
        convolve2d(result, f, nx, ny, g, nkx, nky, false, true, n_threads);
      else
        convolve2d(result, f, nx, ny, g, nkx, nky, false, false, n_threads);
    }
}

void convolve3d_c(DTYPE * const result,
        const DTYPE * const f, const size_t nx, const size_t ny, const size_t nz,
        const DTYPE * const g, const size_t nkx, const size_t nky, const size_t nkz,
        const bool nan_interpolate,
        const bool embed_result_within_padded_region,
        const unsigned n_threads)
{
#ifdef NDEBUG
    if (!result || !f || !g)
        return;
#else
    assert(result);
    assert(f);
    assert(g);
#endif

    if (nan_interpolate) {
      if (embed_result_within_padded_region)
        convolve3d(result, f, nx, ny, nz, g, nkx, nky, nkz, true, true, n_threads);
      else
        convolve3d(result, f, nx, ny, nz, g, nkx, nky, nkz, true, false, n_threads);
    } else {
      if (embed_result_within_padded_region)
        convolve3d(result, f, nx, ny, nz, g, nkx, nky, nkz, false, true, n_threads);
      else
        convolve3d(result, f, nx, ny, nz, g, nkx, nky, nkz, false, false, n_threads);
    }
}

// 1D
FORCE_INLINE void convolve1d(DTYPE * const result,
        const DTYPE * const f, const size_t _nx,
        const DTYPE * const g, const size_t _nkx,
        const bool _nan_interpolate,
        const bool _embed_result_within_padded_region,
        const unsigned n_threads)
{
#ifdef NDEBUG
    if (!result || !f || !g)
        return;
#else
    assert(result);
    assert(f);
    assert(g);
#endif

    const size_t _wkx = _nkx / 2;

#ifdef NDEBUG
    if (!(_nx > 2*_wkx))
        return;
#else
    assert(_nx > 2*_wkx);
#endif

#ifdef _OPENMP
    omp_set_num_threads(n_threads); // Set number of threads to use
#pragma omp parallel
    { // Code within this block is threaded
#endif

    // Copy these to thread locals to allow compiler to optimize (hoist/loads licm)
    // when threaded. Without these, compile time constant conditionals may
    // not be optimized away.
    const size_t nx = _nx;
    const size_t nkx = _nkx;
    const size_t nkx_minus_1 = nkx - 1;
    const bool nan_interpolate = _nan_interpolate;
    const bool embed_result_within_padded_region = _embed_result_within_padded_region;

    // Thread locals
    const size_t wkx = _wkx;
    const omp_iter_var nx_minus_wkx = nx - wkx;
    size_t i_minus_wkx;
    size_t result_index;

    DTYPE top, bot=0., ker, val;

    {omp_iter_var i;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (i = wkx; i < nx_minus_wkx; ++i)
    {
        i_minus_wkx = i - wkx;

        top = 0.;
        if (nan_interpolate) // compile time constant
            bot = 0.;
        {omp_iter_var ii;
        for (ii = 0; ii < nkx; ++ii)
        {
            val = f[i_minus_wkx + ii];
            ker = g[nkx_minus_1 - ii];
            if (nan_interpolate) // compile time constant
            {
                if (!isnan(val))
                {
                    top += val * ker;
                    bot += ker;
                }
            }
            else
                top += val * ker;
        }}

        if (embed_result_within_padded_region) { // compile time constant
            result_index = i;
        } else {
            result_index = i_minus_wkx;
        }

        if (nan_interpolate) // compile time constant
        {
            if (bot == 0) // This should prob be np.isclose(kernel_sum, 0, atol=normalization_zero_tol)
                result[result_index]  = f[i];
            else
                result[result_index]  = top / bot;
        }
        else
            result[result_index] = top;
    }}
#ifdef _OPENMP
    }//end parallel scope
#endif
}

// 2D
FORCE_INLINE void convolve2d(DTYPE * const result,
        const DTYPE * const f, const size_t _nx, const size_t _ny,
        const DTYPE * const g, const size_t _nkx, const size_t _nky,
        const bool _nan_interpolate,
        const bool _embed_result_within_padded_region,
        const unsigned n_threads)
{
#ifdef NDEBUG
    if (!result || !f || !g)
        return;
#else
    assert(result);
    assert(f);
    assert(g);
#endif

    const size_t _wkx = _nkx / 2;
    const size_t _wky = _nky / 2;
#ifdef NDEBUG
    if (!(_nx > 2*_wkx) || !(_ny > 2*_wky))
        return;
#else
    assert(_nx > 2*_wkx);
    assert(_ny > 2*_wky);
#endif

#ifdef _OPENMP
    omp_set_num_threads(n_threads); // Set number of threads to use
#pragma omp parallel
    { // Code within this block is threaded
#endif

    // Copy these to thread locals to allow compiler to optimize (hoist/loads licm)
    // when threaded. Without these, compile time constant conditionals may
    // not be optimized away.
    const size_t nx = _nx, ny = _ny;
    const size_t nkx = _nkx, nky = _nky;
    const size_t nkx_minus_1 = nkx - 1, nky_minus_1 = nky - 1;
    const bool nan_interpolate = _nan_interpolate;
    const bool embed_result_within_padded_region = _embed_result_within_padded_region;

    // Thread locals
    const size_t wkx = _wkx;
    const size_t wky = _wky;
    const omp_iter_var nx_minus_wkx = nx - wkx;
    const omp_iter_var ny_minus_wky = ny - wky;
    const size_t ny_minus_2wky = ny - 2 * wky;
    size_t i_minus_wkx, j_minus_wky;
    size_t result_cursor;
    size_t f_cursor, g_cursor;
    size_t result_index;

    DTYPE top, bot=0., ker, val;

    {omp_iter_var i;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (i = wkx; i < nx_minus_wkx; ++i)
    {
        i_minus_wkx = i - wkx;
        result_cursor = i*ny;

        {omp_iter_var j;
        for (j = wky; j < ny_minus_wky; ++j)
        {
            j_minus_wky = j - wky;

            top = 0.;
            if (nan_interpolate) // compile time constant
                bot = 0.;
            {omp_iter_var ii;
            for (ii = 0; ii < nkx; ++ii)
            {
                f_cursor = (i_minus_wkx + ii)*ny + j_minus_wky;
                g_cursor = (nkx_minus_1 - ii)*nky + nky_minus_1;

                {omp_iter_var jj;
                for (jj = 0; jj < nky; ++jj)
                {
                    val = f[f_cursor + jj];
                    ker = g[g_cursor - jj];
                    if (nan_interpolate) // compile time constant
                    {
                        if (!isnan(val))
                        {
                            top += val * ker;
                            bot += ker;
                        }
                    }
                    else
                        top += val * ker;
                }}
            }}

            if (embed_result_within_padded_region) { // compile time constant
                result_index = result_cursor + j;
            } else {
                result_index = i_minus_wkx * ny_minus_2wky + j_minus_wky;
            }

            if (nan_interpolate) // compile time constant
            {
                if (bot == 0) // This should prob be np.isclose(kernel_sum, 0, atol=normalization_zero_tol)
                    result[result_index]  = f[result_cursor + j] ;
                else
                    result[result_index]  = top / bot;
            }
            else
                result[result_index] = top;
        }}
    }}
#ifdef _OPENMP
    }//end parallel scope
#endif
}

// 3D
FORCE_INLINE void convolve3d(DTYPE * const result,
        const DTYPE * const f, const size_t _nx, const size_t _ny, const size_t _nz,
        const DTYPE * const g, const size_t _nkx, const size_t _nky, const size_t _nkz,
        const bool _nan_interpolate,
        const bool _embed_result_within_padded_region,
        const unsigned n_threads)
{
#ifdef NDEBUG
    if (!result || !f || !g)
        return;
#else
    assert(result);
    assert(f);
    assert(g);
#endif

    const size_t _wkx = _nkx / 2;
    const size_t _wky = _nky / 2;
    const size_t _wkz = _nkz / 2;
#ifdef NDEBUG
    if (!(_nx > 2*_wkx) || !(_ny > 2*_wky) || !(_nz > 2*_wkz))
        return;
#else
    assert(_nx > 2*_wkx);
    assert(_ny > 2*_wky);
    assert(_nz > 2*_wkz);
#endif

#ifdef _OPENMP
    omp_set_num_threads(n_threads); // Set number of threads to use
#pragma omp parallel
    { // Code within this block is threaded
#endif

    // Copy these to thread locals to allow compiler to optimize (hoist/loads licm)
    // when threaded. Without these, compile time constant conditionals may
    // not be optimized away.
    const size_t nx = _nx, ny = _ny, nz = _nz;
    const size_t nkx = _nkx, nky = _nky, nkz = _nkz;
    const size_t nkx_minus_1 = nkx - 1, nky_minus_1 = nky - 1, nkz_minus_1 = nkz - 1;
    const bool nan_interpolate = _nan_interpolate;
    const bool embed_result_within_padded_region = _embed_result_within_padded_region;

    // Thread locals
    const size_t wkx = _wkx;
    const size_t wky = _wky;
    const size_t wkz = _wkz;
    const size_t nx_minus_wkx = nx - wkx;
    const omp_iter_var ny_minus_wky = ny - wky;
    const omp_iter_var nz_minus_wkz = nz - wkz;
    const size_t ny_minus_2wky = ny - 2 * wky;
    const size_t nz_minus_2wkz = nz - 2 * wkz;
    size_t i_minus_wkx, j_minus_wky, k_minus_wkz;
    size_t f_ii_cursor, g_ii_cursor;
    size_t f_cursor, g_cursor;
    size_t array_cursor, array_i_cursor;
    size_t result_index;

    DTYPE top, bot=0., ker, val;

    {omp_iter_var i;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (i = wkx; i < nx_minus_wkx; ++i)
    {
        i_minus_wkx = i - wkx;
        array_i_cursor = i*ny;

        {omp_iter_var j;
        for (j = wky; j < ny_minus_wky; ++j)
        {
            j_minus_wky = j - wky;
            array_cursor = (array_i_cursor + j)*nz;

            {omp_iter_var k;
            for (k = wkz; k < nz_minus_wkz; ++k)
            {
                k_minus_wkz = k - wkz;

                top = 0.;
                if (nan_interpolate) // compile time constant
                    bot = 0.;
                {omp_iter_var ii;
                for (ii = 0; ii < nkx; ++ii)
                {
                    f_ii_cursor = ((i_minus_wkx + ii)*ny + j_minus_wky)*nz + k_minus_wkz;
                    g_ii_cursor = ((nkx_minus_1 - ii)*nky + nky_minus_1)*nkz + nkz_minus_1;

                    {omp_iter_var jj;
                    for (jj = 0; jj < nky; ++jj)
                    {
                        f_cursor = f_ii_cursor + jj*nz;
                        g_cursor = g_ii_cursor - jj*nkz;
                        {omp_iter_var kk;
                        for (kk = 0; kk < nkz; ++kk)
                        {
                            val = f[f_cursor + kk];
                            ker = g[g_cursor - kk];
                            if (nan_interpolate) // compile time constant
                            {
                                if (!isnan(val))
                                {
                                    top += val * ker;
                                    bot += ker;
                                }
                            }
                            else
                                top += val * ker;
                        }}
                    }}
                }}

                if (embed_result_within_padded_region) { // compile time constant
                    result_index = array_cursor + k;
                } else {
                    result_index = (i_minus_wkx*ny_minus_2wky + j_minus_wky)*nz_minus_2wkz + k_minus_wkz;
                }

                if (nan_interpolate) // compile time constant
                {
                    if (bot == 0) // This should prob be np.isclose(kernel_sum, 0, atol=normalization_zero_tol)
                        result[result_index]  = f[array_cursor+ k] ;
                    else
                        result[result_index]  = top / bot;
                }
                else
                    result[result_index] = top;
            }}
        }}
    }}
#ifdef _OPENMP
    }//end parallel scope
#endif
}
