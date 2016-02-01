# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
import numpy as np
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_t

cdef inline int int_max(int a, int b) nogil: return a if a >= b else b
cdef inline int int_min(int a, int b) nogil: return a if a <= b else b

cdef extern from "numpy/npy_math.h" nogil:
    bint npy_isnan(double x)

cimport cython


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def convolve1d_boundary_extend(np.ndarray[DTYPE_t, ndim=1] f,
                               np.ndarray[DTYPE_t, ndim=1] g):

    if g.shape[0] % 2 != 1:
        raise ValueError("Convolution kernel must have odd dimensions")

    assert f.dtype == DTYPE and g.dtype == DTYPE

    cdef int nx = f.shape[0]
    cdef int nkx = g.shape[0]
    cdef int wkx = nkx // 2
    cdef np.ndarray[DTYPE_t, ndim=1] fixed = np.empty([nx], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] conv = np.empty([nx], dtype=DTYPE)
    cdef unsigned int i, iii
    cdef int ii

    cdef int iimin, iimax

    cdef DTYPE_t top, bot, ker, val

    # release the GIL
    with nogil:
        # Need a first pass to replace NaN values with value convolved from
        # neighboring values
        for i in range(nx):
            if npy_isnan(f[i]):
                top = 0.
                bot = 0.
                iimin = i - wkx
                iimax = i + wkx + 1
                for ii in range(iimin, iimax):
                    iii = int_min(int_max(ii, 0), nx - 1)
                    val = f[iii]
                    if not npy_isnan(val):
                        ker = g[<unsigned int>(wkx + ii - i)]
                        top += val * ker
                        bot += ker

                if bot != 0.:
                    fixed[i] = top / bot
                else:
                    fixed[i] = f[i]
            else:
                fixed[i] = f[i]

        # Now run the proper convolution
        for i in range(nx):
            if not npy_isnan(fixed[i]):
                top = 0.
                bot = 0.
                iimin = i - wkx
                iimax = i + wkx + 1
                for ii in range(iimin, iimax):
                    iii = int_min(int_max(ii, 0), nx - 1)
                    val = fixed[iii]
                    ker = g[<unsigned int>(wkx + ii - i)]
                    if not npy_isnan(val):
                        top += val * ker
                        bot += ker
                if bot != 0:
                    conv[i] = top / bot
                else:
                    conv[i] = fixed[i]
            else:
                conv[i] = fixed[i]
    # GIL acquired again here
    return conv


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def convolve2d_boundary_extend(np.ndarray[DTYPE_t, ndim=2] f,
                               np.ndarray[DTYPE_t, ndim=2] g):

    if g.shape[0] % 2 != 1 or g.shape[1] % 2 != 1:
        raise ValueError("Convolution kernel must have odd dimensions")

    assert f.dtype == DTYPE and g.dtype == DTYPE

    cdef int nx = f.shape[0]
    cdef int ny = f.shape[1]
    cdef int nkx = g.shape[0]
    cdef int nky = g.shape[1]
    cdef int wkx = nkx // 2
    cdef int wky = nky // 2
    cdef np.ndarray[DTYPE_t, ndim=2] fixed = np.empty([nx, ny], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] conv = np.empty([nx, ny], dtype=DTYPE)
    cdef unsigned int i, j, iii, jjj
    cdef int ii, jj

    cdef int iimin, iimax, jjmin, jjmax

    cdef DTYPE_t top, bot, ker, val

    # release the GIL
    with nogil:
        # Need a first pass to replace NaN values with value convolved from
        # neighboring values
        for i in range(nx):
            for j in range(ny):
                if npy_isnan(f[i, j]):
                    top = 0.
                    bot = 0.
                    iimin = i - wkx
                    iimax = i + wkx + 1
                    jjmin = j - wky
                    jjmax = j + wky + 1
                    for ii in range(iimin, iimax):
                        for jj in range(jjmin, jjmax):
                            iii = int_min(int_max(ii, 0), nx - 1)
                            jjj = int_min(int_max(jj, 0), ny - 1)
                            val = f[iii, jjj]
                            if not npy_isnan(val):
                                ker = g[<unsigned int>(wkx + ii - i),
                                        <unsigned int>(wky + jj - j)]
                                top += val * ker
                                bot += ker

                    if bot != 0.:
                        fixed[i, j] = top / bot
                    else:
                        fixed[i, j] = f[i, j]
                else:
                    fixed[i, j] = f[i, j]

        # Now run the proper convolution
        for i in range(nx):
            for j in range(ny):
                if not npy_isnan(fixed[i, j]):
                    top = 0.
                    bot = 0.
                    iimin = i - wkx
                    iimax = i + wkx + 1
                    jjmin = j - wky
                    jjmax = j + wky + 1
                    for ii in range(iimin, iimax):
                        for jj in range(jjmin, jjmax):
                            iii = int_min(int_max(ii, 0), nx - 1)
                            jjj = int_min(int_max(jj, 0), ny - 1)
                            val = fixed[iii, jjj]
                            ker = g[<unsigned int>(wkx + ii - i),
                                    <unsigned int>(wky + jj - j)]
                            if not npy_isnan(val):
                                top += val * ker
                                bot += ker
                    if bot != 0:
                        conv[i, j] = top / bot
                    else:
                        conv[i, j] = fixed[i, j]
                else:
                    conv[i, j] = fixed[i, j]
    # GIL acquired again here
    return conv


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def convolve3d_boundary_extend(np.ndarray[DTYPE_t, ndim=3] f,
                               np.ndarray[DTYPE_t, ndim=3] g):

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
    cdef np.ndarray[DTYPE_t, ndim=3] fixed = np.empty([nx, ny, nz], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=3] conv = np.empty([nx, ny, nz], dtype=DTYPE)
    cdef unsigned int i, j, k, iii, jjj, kkk
    cdef int ii, jj, kk

    cdef int iimin, iimax, jjmin, jjmax, kkmin, kkmax

    cdef DTYPE_t top, bot, ker, val

    # release the GIL
    with nogil:
        # Need a first pass to replace NaN values with value convolved from
        # neighboring values
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if npy_isnan(f[i, j, k]):
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
                                    iii = int_min(int_max(ii, 0), nx - 1)
                                    jjj = int_min(int_max(jj, 0), ny - 1)
                                    kkk = int_min(int_max(kk, 0), nz - 1)
                                    val = f[iii, jjj, kkk]
                                    if not npy_isnan(val):
                                        ker = g[<unsigned int>(wkx + ii - i),
                                                <unsigned int>(wky + jj - j),
                                                <unsigned int>(wkz + kk - k)]
                                        top += val * ker
                                        bot += ker

                        if bot != 0.:
                            fixed[i, j, k] = top / bot
                        else:
                            fixed[i, j, k] = f[i, j, k]
                    else:
                        fixed[i, j, k] = f[i, j, k]

        # Now run the proper convolution
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if not npy_isnan(fixed[i, j, k]):
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
                                    iii = int_min(int_max(ii, 0), nx - 1)
                                    jjj = int_min(int_max(jj, 0), ny - 1)
                                    kkk = int_min(int_max(kk, 0), nz - 1)
                                    val = fixed[iii, jjj, kkk]
                                    ker = g[<unsigned int>(wkx + ii - i),
                                            <unsigned int>(wky + jj - j),
                                            <unsigned int>(wkz + kk - k)]
                                    if not npy_isnan(val):
                                        top += val * ker
                                        bot += ker
                        if bot != 0:
                            conv[i, j, k] = top / bot
                        else:
                            conv[i, j, k] = fixed[i, j, k]
                    else:
                        conv[i, j, k] = fixed[i, j, k]
    # GIL acquired again here
    return conv
