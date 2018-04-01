# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
cimport numpy as np
from .utils import KernelSizeError


DTYPE = float
ctypedef np.float_t DTYPE_t

cdef extern from "numpy/npy_math.h" nogil:
    bint npy_isnan(double x)

cimport cython

cdef _raise_if_kernel_larger_than_array(const np.npy_intp * array_shape, const np.npy_intp * kernel_shape, length):
    """
    Helper function to raise an exception due to the following:

    For boundary=None only the center space is convolved. All array indices within a
    distance kernel.shape//2 from the edge are completely ignored (zeroed).
    E.g. (1D list) only the indices len(kernel)//2 : len(array)-len(kernel)//2
    are convolved. It is therefore not possible to use this method to convolve an
    array by a kernel that is larger* than the array - as ALL pixels would be ignored
    leaving an array of only zeros.
    * For even kernels the correctness condition is array_shape > kernel_shape.
    For odd kernels it is:
    array_shape >= kernel_shape OR array_shape > kernel_shape-1 OR array_shape > 2*(kernel_shape//2).
    Since the latter is equal to the former two for even lengths, the latter condition is complete.
    """

    for i in range(length):
        if not array_shape[i] > 2*(kernel_shape[i]//2):
            raise KernelSizeError("for boundary=None all kernel axes must be smaller than array's. "
                                  "Triggered by: if not numpy.all(array.shape > 2*(kernel.shape//2)). "
                                  "Use boundary in ('fill', 'extend', 'wrap') instead.")


@cython.boundscheck(False)  # turn off bounds-checking for entire function
def convolve1d_boundary_none(np.ndarray[DTYPE_t, ndim=1] f,
                             np.ndarray[DTYPE_t, ndim=1] g,
                             bint normalize_by_kernel):

    if g.shape[0] % 2 != 1:
        raise ValueError("Convolution kernel must have odd dimensions")

    _raise_if_kernel_larger_than_array(f.shape, g.shape, 1)

    assert f.dtype == DTYPE and g.dtype == DTYPE

    cdef int nx = f.shape[0]
    cdef int nkx = g.shape[0]
    cdef int wkx = nkx // 2

    # The following need to be set to zeros rather than empty because the
    # boundary does not get reset.
    cdef np.ndarray[DTYPE_t, ndim=1] conv = np.zeros([nx], dtype=DTYPE)

    cdef unsigned int i, ii

    cdef int iimin, iimax

    cdef DTYPE_t top, bot, ker, val

    # release the GIL
    with nogil:

        # Now run the proper convolution
        for i in range(wkx, nx - wkx):
            top = 0.
            bot = 0.
            for ii in range(i - wkx, i + wkx + 1):
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
def convolve2d_boundary_none(np.ndarray[DTYPE_t, ndim=2] f,
                             np.ndarray[DTYPE_t, ndim=2] g,
                             bint normalize_by_kernel):

    if g.shape[0] % 2 != 1 or g.shape[1] % 2 != 1:
        raise ValueError("Convolution kernel must have odd dimensions")

    _raise_if_kernel_larger_than_array(f.shape, g.shape, 2)

    assert f.dtype == DTYPE and g.dtype == DTYPE

    cdef int nx = f.shape[0]
    cdef int ny = f.shape[1]
    cdef int nkx = g.shape[0]
    cdef int nky = g.shape[1]
    cdef int wkx = nkx // 2
    cdef int wky = nky // 2

    # The following need to be set to zeros rather than empty because the
    # boundary does not get reset.
    cdef np.ndarray[DTYPE_t, ndim=2] conv = np.zeros([nx, ny], dtype=DTYPE)

    cdef unsigned int i, j, ii, jj

    cdef int iimin, iimax, jjmin, jjmax

    cdef DTYPE_t top, bot, ker, val

    # release the GIL
    with nogil:

        # Now run the proper convolution
        for i in range(wkx, nx - wkx):
            for j in range(wky, ny - wky):
                top = 0.
                bot = 0.
                for ii in range(i - wkx, i + wkx + 1):
                    for jj in range(j - wky, j + wky + 1):
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
def convolve3d_boundary_none(np.ndarray[DTYPE_t, ndim=3] f,
                             np.ndarray[DTYPE_t, ndim=3] g,
                             bint normalize_by_kernel):

    if g.shape[0] % 2 != 1 or g.shape[1] % 2 != 1 or g.shape[2] % 2 != 1:
        raise ValueError("Convolution kernel must have odd dimensions")

    _raise_if_kernel_larger_than_array(f.shape, g.shape, 3)

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

    # The following need to be set to zeros rather than empty because the
    # boundary does not get reset.
    cdef np.ndarray[DTYPE_t, ndim=3] conv = np.zeros([nx, ny, nz], dtype=DTYPE)

    cdef unsigned int i, j, k, ii, jj, kk

    cdef int iimin, iimax, jjmin, jjmax, kkmin, kkmax

    cdef DTYPE_t top, bot, ker, val

    # release the GIL
    with nogil:

        # Now run the proper convolution
        for i in range(wkx, nx - wkx):
            for j in range(wky, ny - wky):
                for k in range(wkz, nz - wkz):
                    top = 0.
                    bot = 0.
                    for ii in range(i - wkx, i + wkx + 1):
                        for jj in range(j - wky, j + wky + 1):
                            for kk in range(k - wkz, k + wkz + 1):
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
    # GIL acquired again here
    return conv
