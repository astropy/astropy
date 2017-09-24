# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
cimport numpy as np


DTYPE = float
ctypedef np.float_t DTYPE_t

cdef extern from "numpy/npy_math.h" nogil:
    bint npy_isnan(double x)

cimport cython


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def convolve1d_boundary_fill(np.ndarray[DTYPE_t, ndim=1] f,
                             np.ndarray[DTYPE_t, ndim=1] g,
                             float fill_value,
                             bint normalize_by_kernel
                            ):

    if g.shape[0] % 2 != 1:
        raise ValueError("Convolution kernel must have odd dimensions")

    assert f.dtype == DTYPE and g.dtype == DTYPE

    cdef int nx = f.shape[0]
    cdef int nkx = g.shape[0]
    cdef int wkx = nkx // 2
    cdef np.ndarray[DTYPE_t, ndim=1] conv = np.empty([nx], dtype=DTYPE)
    cdef unsigned int i, iii
    cdef int ii

    cdef int iimin, iimax

    cdef DTYPE_t top, bot, ker, val

    # release the GIL
    with nogil:

        # Now run the proper convolution
        for i in range(nx):
            top = 0.
            bot = 0.
            iimin = i - wkx
            iimax = i + wkx + 1
            for ii in range(iimin, iimax):
                if ii < 0 or ii > nx - 1:
                    val = fill_value
                else:
                    val = f[ii]
                ker = g[<unsigned int>(nkx - 1 - (wkx + ii - i))]
                if not npy_isnan(val):
                    top += val * ker
                    bot += ker
            if normalize_by_kernel:
                if bot == 0:
                    conv[i] = f[i]
                else:
                    conv[i] = top / bot
            else:
                conv[i] = top
    # GIL acquired again here
    return conv


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def convolve2d_boundary_fill(np.ndarray[DTYPE_t, ndim=2] f,
                             np.ndarray[DTYPE_t, ndim=2] g,
                             float fill_value,
                             bint normalize_by_kernel
                            ):

    if g.shape[0] % 2 != 1 or g.shape[1] % 2 != 1:
        raise ValueError("Convolution kernel must have odd dimensions")

    assert f.dtype == DTYPE and g.dtype == DTYPE

    cdef int nx = f.shape[0]
    cdef int ny = f.shape[1]
    cdef int nkx = g.shape[0]
    cdef int nky = g.shape[1]
    cdef int wkx = nkx // 2
    cdef int wky = nky // 2
    cdef np.ndarray[DTYPE_t, ndim=2] conv = np.empty([nx, ny], dtype=DTYPE)
    cdef unsigned int i, j, iii, jjj
    cdef int ii, jj

    cdef int iimin, iimax, jjmin, jjmax

    cdef DTYPE_t top, bot, ker, val

    # release the GIL
    with nogil:

        # now run the proper convolution
        for i in range(nx):
            for j in range(ny):
                top = 0.
                bot = 0.
                iimin = i - wkx
                iimax = i + wkx + 1
                jjmin = j - wky
                jjmax = j + wky + 1
                for ii in range(iimin, iimax):
                    for jj in range(jjmin, jjmax):
                        if ii < 0 or ii > nx - 1 or jj < 0 or jj > ny - 1:
                            val = fill_value
                        else:
                            val = f[ii, jj]
                        ker = g[<unsigned int>(nkx - 1 - (wkx + ii - i)),
                                <unsigned int>(nky - 1 - (wky + jj - j))]
                        if not npy_isnan(val):
                            top += val * ker
                            bot += ker
                if normalize_by_kernel:
                    if bot == 0:
                        conv[i, j] = f[i, j]
                    else:
                        conv[i, j] = top / bot
                else:
                    conv[i, j] = top
    # GIL acquired again here
    return conv


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def convolve3d_boundary_fill(np.ndarray[DTYPE_t, ndim=3] f,
                             np.ndarray[DTYPE_t, ndim=3] g,
                             float fill_value,
                             bint normalize_by_kernel):

    if g.shape[0] % 2 != 1 or g.shape[1] % 2 != 1 or g.shape[2] % 2 != 1:
        raise ValueError("Convolution kernel must have odd dimensions")

    assert f.dtype == DTYPE and g.dtype == DTYPE

    cdef int nx = f.shape[0]
    cdef int ny = f.shape[1]
    cdef int nz = f.shape[2]
    cdef int nkx = g.shape[0]
    cdef int nky = g.shape[1]
    cdef int nkz = g.shape[2]
    cdef int wkx = nkx // 2
    cdef int wky = nky // 2
    cdef int wkz = nkz // 2
    cdef np.ndarray[DTYPE_t, ndim=3] conv = np.empty([nx, ny, nz], dtype=DTYPE)
    cdef unsigned int i, j, k, iii, jjj, kkk
    cdef int ii, jj, kk

    cdef int iimin, iimax, jjmin, jjmax, kkmin, kkmax

    cdef DTYPE_t top, bot, ker, val

    # release the GIL
    with nogil:

        # Now run the proper convolution
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    top = 0.
                    bot = 0.
                    iimin = i - wkx
                    iimax = i + wkx + 1
                    jjmin = j - wky
                    jjmax = j + wky + 1
                    kkmin = k - wkz
                    kkmax = k + wkz + 1
                    for ii in range(iimin, iimax):
                        for jj in range(jjmin, jjmax):
                            for kk in range(kkmin, kkmax):
                                if ii < 0 or ii > nx - 1 or jj < 0 or jj > ny - 1 or kk < 0 or kk > nz - 1:
                                    val = fill_value
                                else:
                                    val = f[ii, jj, kk]
                                ker = g[<unsigned int>(nkx - 1 - (wkx + ii - i)),
                                        <unsigned int>(nky - 1 - (wky + jj - j)),
                                        <unsigned int>(nkz - 1 - (wkz + kk - k))]
                                if not npy_isnan(val):
                                    top += val * ker
                                    bot += ker
                    if normalize_by_kernel:
                        if bot == 0:
                            conv[i, j, k] = f[i, j, k]
                        else:
                            conv[i, j, k] = top / bot
                    else:
                        conv[i, j, k] = top
    # GIl acquired again here
    return conv
