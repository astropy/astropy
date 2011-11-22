from __future__ import division
import numpy as np
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_t

cdef extern from "math.h":
    bint isnan(double x)

cimport cython


@cython.boundscheck(False)  # turn of bounds-checking for entire function
def convolve2d_edge_fill(np.ndarray[DTYPE_t, ndim=2] f,
                         np.ndarray[DTYPE_t, ndim=2] g,
                         float fill_value):

    if g.shape[0] % 2 != 1 or g.shape[1] % 2 != 1:
        raise ValueError("Only odd dimensions on filter supported")

    assert f.dtype == DTYPE and g.dtype == DTYPE

    cdef int nx = f.shape[0]
    cdef int ny = f.shape[1]
    cdef int nkx = g.shape[0]
    cdef int nky = g.shape[1]
    cdef int wkx = nkx // 2
    cdef int wky = nky // 2
    cdef np.ndarray[DTYPE_t, ndim=2] fixed = np.zeros([nx, ny], dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] conv = np.zeros([nx, ny], dtype=DTYPE)
    cdef unsigned int i, j, iii, jjj
    cdef int ii, jj

    cdef int iimin, iimax, jjmin, jjmax

    cdef DTYPE_t top, bot, ker, val

    # Need a first pass to replace NaN values with value convolved from
    # neighboring values
    for i in range(nx):
        for j in range(ny):
            if isnan(f[i, j]):
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
                        if not isnan(val):
                            ker = g[<unsigned int>(wkx + ii - i),
                                    <unsigned int>(wky + jj - j)]
                            top += val * ker
                            bot += ker
                if bot > 0.:
                    fixed[i, j] = top / bot
                else:
                    fixed[i, j] = f[i, j]
            else:
                fixed[i, j] = f[i, j]

    # Now run the proper convolution
    for i in range(nx):
        for j in range(ny):
            if not isnan(fixed[i, j]):
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
                            val = fixed[ii, jj]
                        ker = g[<unsigned int>(wkx + ii - i),
                                <unsigned int>(wky + jj - j)]
                        top += val * ker
                        bot += ker
                if bot > 0:
                    conv[i, j] = top / bot
                else:
                    conv[i, j] = fixed[i, j]
            else:
                conv[i, j] = fixed[i, j]

    return conv
