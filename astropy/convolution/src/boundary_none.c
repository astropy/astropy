// Licensed under a 3-clause BSD style license - see LICENSE.rst

/*------------------------------WARNING!------------------------------
 * The C functions below are NOT designed to be called externally to
 * the Python function astropy/astropy/convolution/convolve.py.
 * They do NOT include any of the required correct usage checking.
 *------------------------------WARNING!------------------------------
 */

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "convolve.h"

// Distutils on Windows will automatically exports ``PyInit_lib_convolve``,
// create dummy to prevent linker complaining about missing symbol.
#if defined(_MSC_VER)
void PyInit_lib_convolve(void)
{
    return;
}
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


void convolveNd_boundary_none_c(DTYPE * const result,
        const DTYPE * const f,
        const unsigned n_dim,
        const size_t * const image_shape,
        const DTYPE * const g,
        const size_t * const kernel_shape,
        const bool nan_interpolate,
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
        convolve1d_boundary_none_c(result, f,
                image_shape[0],
                g, kernel_shape[0],
                nan_interpolate, n_threads);
    else if (n_dim == 2)
        convolve2d_boundary_none_c(result, f,
                image_shape[0], image_shape[1],
                g, kernel_shape[0], kernel_shape[1],
                nan_interpolate, n_threads);
    else if (n_dim == 3)
        convolve3d_boundary_none_c(result, f,
                        image_shape[0], image_shape[1], image_shape[2],
                        g, kernel_shape[0], kernel_shape[1], kernel_shape[2],
                        nan_interpolate, n_threads);
    else
        assert(0); // Unimplemented: n_dim > 3
}

/*-------------------------PERFORMANCE NOTES--------------------------------
 * The function wrappers below are designed to take advantage of the following:
 * The preprocessor will inline convolve<N>d_boundary_none(), effectively
 * expanding the two logical branches, replacing nan_interpolate
 * for their literal equivalents. The corresponding conditionals
 * within these functions will then be optimized away, this
 * being the goal - removing the unnecessary conditionals from
 * the loops without duplicating code.
 *--------------------------------------------------------------------------
 */

void convolve1d_boundary_none_c(DTYPE * const result,
        const DTYPE * const f, const size_t nx,
        const DTYPE * const g, const size_t nkx,
        const bool nan_interpolate,
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

    if (nan_interpolate)
        convolve1d_boundary_none(result, f, nx, g, nkx, true, n_threads);
    else
        convolve1d_boundary_none(result, f, nx, g, nkx, false, n_threads);
}

void convolve2d_boundary_none_c(DTYPE * const result,
        const DTYPE * const f, const size_t nx, const size_t ny,
        const DTYPE * const g, const size_t nkx, const size_t nky,
        const bool nan_interpolate,
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

    if (nan_interpolate)
        convolve2d_boundary_none(result, f, nx, ny, g, nkx, nky, true, n_threads);
    else
        convolve2d_boundary_none(result, f, nx, ny, g, nkx, nky, false, n_threads);
}

void convolve3d_boundary_none_c(DTYPE * const result,
        const DTYPE * const f, const size_t nx, const size_t ny, const size_t nz,
        const DTYPE * const g, const size_t nkx, const size_t nky, const size_t nkz,
        const bool nan_interpolate,
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

    if (nan_interpolate)
        convolve3d_boundary_none(result, f, nx, ny, nz, g, nkx, nky, nkz, true, n_threads);
    else
        convolve3d_boundary_none(result, f, nx, ny, nz, g, nkx, nky, nkz, false, n_threads);
}

// 1D
FORCE_INLINE void convolve1d_boundary_none(DTYPE * const result,
        const DTYPE * const f, const size_t _nx,
        const DTYPE * const g, const size_t _nkx,
        const bool _nan_interpolate,
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

    const unsigned _wkx = _nkx / 2;
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
    const bool nan_interpolate = _nan_interpolate;

    // Thread locals
    const unsigned wkx = _wkx;
    const unsigned nkx_minus_1 = nkx-1;
    unsigned wkx_minus_i;
    unsigned ker_i;
    const omp_iter_var nx_minus_wkx = nx - wkx;
    unsigned i_minus_wkx;
    const unsigned wkx_plus_1 = wkx + 1;
    omp_iter_var i_plus_wkx_plus_1;
    unsigned nkx_minus_1_minus_wkx_plus_i;

    DTYPE top, bot=0., ker, val;

    {omp_iter_var i;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (i = wkx; i < nx_minus_wkx; ++i)
    {
        wkx_minus_i = wkx - i; // wkx - 1
        i_minus_wkx = i - wkx; //i - wkx
        i_plus_wkx_plus_1 = i + wkx_plus_1; // i + wkx + 1
        nkx_minus_1_minus_wkx_plus_i = nkx_minus_1 - wkx_minus_i; // nkx - 1 - (wkx - i)

        top = 0.;
        if (nan_interpolate) // compile time constant
            bot = 0.;
        {omp_iter_var ii;
        for (ii = i_minus_wkx; ii < i_plus_wkx_plus_1; ++ii)
        {
            ker_i = nkx_minus_1_minus_wkx_plus_i - ii; // nkx - 1 - (wkx + ii - i)
            val = f[ii];
            ker = g[ker_i];
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

        if (nan_interpolate) // compile time constant
        {
            if (bot == 0) // This should prob be np.isclose(kernel_sum, 0, atol=normalization_zero_tol)
                result[i]  = f[i] ;
            else
                result[i]  = top / bot;
        }
        else
            result[i] = top;
    }}
#ifdef _OPENMP
    }//end parallel scope
#endif
}

// 2D
FORCE_INLINE void convolve2d_boundary_none(DTYPE * const result,
        const DTYPE * const f, const size_t _nx, const size_t _ny,
        const DTYPE * const g, const size_t _nkx, const size_t _nky,
        const bool _nan_interpolate,
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

    const unsigned _wkx = _nkx / 2;
    const unsigned _wky = _nky / 2;
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
    const bool nan_interpolate = _nan_interpolate;

    // Thread locals
    const unsigned wkx = _wkx;
    const unsigned wky = _wky;
    const unsigned nkx_minus_1 = nkx-1, nky_minus_1 = nky-1;
    unsigned wkx_minus_i, wky_minus_j;
    unsigned ker_i, ker_j;
    const omp_iter_var nx_minus_wkx = nx - wkx;
    const omp_iter_var ny_minus_wky = ny - wky;
    unsigned i_minus_wkx, j_minus_wky;
    const unsigned wkx_plus_1 = wkx + 1;
    const unsigned wky_plus_1 = wky + 1;
    omp_iter_var i_plus_wkx_plus_1, j_plus_wky_plus_1;
    unsigned nkx_minus_1_minus_wkx_plus_i, nky_minus_1_minus_wky_plus_j;

    DTYPE top, bot=0., ker, val;

    {omp_iter_var i;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (i = wkx; i < nx_minus_wkx; ++i)
    {
        wkx_minus_i = wkx - i; // wkx - 1
        i_minus_wkx = i - wkx; //i - wkx
        i_plus_wkx_plus_1 = i + wkx_plus_1; // i + wkx + 1
        nkx_minus_1_minus_wkx_plus_i = nkx_minus_1 - wkx_minus_i; // nkx - 1 - (wkx - i)

        {omp_iter_var j;
        for (j = wky; j < ny_minus_wky; ++j)
        {
            wky_minus_j = wkx - j; // wky - j
            j_minus_wky = j - wky; // j - wky
            j_plus_wky_plus_1 = j + wky_plus_1; // j + wky + 1
            nky_minus_1_minus_wky_plus_j = nky_minus_1 - wky_minus_j; // nky - 1 - (wky - i)

            top = 0.;
            if (nan_interpolate) // compile time constant
                bot = 0.;
            {omp_iter_var ii;
            for (ii = i_minus_wkx; ii < i_plus_wkx_plus_1; ++ii)
            {
                ker_i = nkx_minus_1_minus_wkx_plus_i - ii; // nkx - 1 - (wkx + ii - i)
                {omp_iter_var jj;
                for (jj = j_minus_wky; jj < j_plus_wky_plus_1; ++jj)
                {
                    ker_j = nky_minus_1_minus_wky_plus_j - jj; // nky - 1 - (wky + jj - j)
                    val = f[ii*ny + jj]; //[ii, jj];
                    ker = g[ker_i*nky + ker_j]; // [ker_i, ker_j];
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
            if (nan_interpolate) // compile time constant
            {
                if (bot == 0) // This should prob be np.isclose(kernel_sum, 0, atol=normalization_zero_tol)
                    result[i*ny + j]  = f[i*ny + j] ;
                else
                    result[i*ny + j]  = top / bot;
            }
            else
                result[i*ny + j] = top;
        }}
    }}
#ifdef _OPENMP
    }//end parallel scope
#endif
}

// 3D
FORCE_INLINE void convolve3d_boundary_none(DTYPE * const result,
        const DTYPE * const f, const size_t _nx, const size_t _ny, const size_t _nz,
        const DTYPE * const g, const size_t _nkx, const size_t _nky, const size_t _nkz,
        const bool _nan_interpolate,
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

    const unsigned _wkx = _nkx / 2;
    const unsigned _wky = _nky / 2;
    const unsigned _wkz = _nkz / 2;
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
    const bool nan_interpolate = _nan_interpolate;

    // Thread locals
    const unsigned wkx = _wkx;
    const unsigned wky = _wky;
    const unsigned wkz = _wkz;
    const unsigned nkx_minus_1 = nkx-1, nky_minus_1 = nky-1, nkz_minus_1 = nkz-1;
    unsigned wkx_minus_i, wky_minus_j, wkz_minus_k;
    unsigned ker_i, ker_j, ker_k;
    const unsigned nx_minus_wkx = nx - wkx;
    const omp_iter_var ny_minus_wky = ny - wky;
    const omp_iter_var nz_minus_wkz = nz - wkz;
    unsigned i_minus_wkx, j_minus_wky, k_minus_wkz;
    const unsigned wkx_plus_1 = wkx + 1;
    const unsigned wky_plus_1 = wky + 1;
    const unsigned wkz_plus_1 = wkz + 1;
    omp_iter_var i_plus_wkx_plus_1, j_plus_wky_plus_1, k_plus_wkz_plus_1;
    unsigned nkx_minus_1_minus_wkx_plus_i, nky_minus_1_minus_wky_plus_j, nkz_minus_1_minus_wkz_plus_k;

    DTYPE top, bot=0., ker, val;

    {omp_iter_var i;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (i = wkx; i < nx_minus_wkx; ++i)
    {
        wkx_minus_i = wkx - i; // wkx - 1
        i_minus_wkx = i - wkx; //i - wkx
        i_plus_wkx_plus_1 = i + wkx_plus_1; // i + wkx + 1
        nkx_minus_1_minus_wkx_plus_i = nkx_minus_1 - wkx_minus_i; // nkx - 1 - (wkx - i)

        {omp_iter_var j;
        for (j = wky; j < ny_minus_wky; ++j)
        {
            wky_minus_j = wkx - j; // wky - j
            j_minus_wky = j - wky; // j - wky
            j_plus_wky_plus_1 = j + wky_plus_1; // j + wky + 1
            nky_minus_1_minus_wky_plus_j = nky_minus_1 - wky_minus_j; // nky - 1 - (wky - i)

            {omp_iter_var k;
            for (k = wkz; k < nz_minus_wkz; ++k)
            {
                wkz_minus_k = wkz - k; // wkz - k
                k_minus_wkz = k - wkz; // k - wkz
                k_plus_wkz_plus_1 = k + wkz_plus_1; // k + wkz + 1
                nkz_minus_1_minus_wkz_plus_k = nkz_minus_1 - wkz_minus_k; // nkz - 1 - (wkz - i)

                top = 0.;
                if (nan_interpolate) // compile time constant
                    bot = 0.;
                {omp_iter_var ii;
                for (ii = i_minus_wkx; ii < i_plus_wkx_plus_1; ++ii)
                {
                    ker_i = nkx_minus_1_minus_wkx_plus_i - ii; // nkx - 1 - (wkx + ii - i)
                    {omp_iter_var jj;
                    for (jj = j_minus_wky; jj < j_plus_wky_plus_1; ++jj)
                    {
                        ker_j = nky_minus_1_minus_wky_plus_j - jj; // nky - 1 - (wky + jj - j)
                        {omp_iter_var kk;
                        for (kk = k_minus_wkz; kk < k_plus_wkz_plus_1; ++kk)
                        {
                            ker_k = nkz_minus_1_minus_wkz_plus_k - kk; // nkz - 1 - (wkz + kk - k)

                            val = f[(ii*ny + jj)*nz + kk]; //[ii, jj, kk];
                            ker = g[(ker_i*nky + ker_j)*nkz + ker_k]; // [ker_i, ker_j, ker_k];
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
                if (nan_interpolate) // compile time constant
                {
                    if (bot == 0) // This should prob be np.isclose(kernel_sum, 0, atol=normalization_zero_tol)
                        result[(i*ny + j)*nz + k]  = f[(i*ny + j)*nz + k] ;
                    else
                        result[(i*ny + j)*nz + k]  = top / bot;
                }
                else
                    result[(i*ny + j)*nz + k] = top;
            }}
        }}
    }}
#ifdef _OPENMP
    }//end parallel scope
#endif
}
