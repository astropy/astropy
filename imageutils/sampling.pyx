# Licensed under a 3-clause BSD style license - see LICENSE.rst
#cython: boundscheck=False
#cython: wraparound=False
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np

cimport cython
cimport numpy as np
DTYPE = np.float
ctypedef np.float_t DTYPE_t

__all__ = ['downsample', 'upsample']


def downsample(np.ndarray[DTYPE_t, ndim=2] array, int block_size):
    """
    downsample(np.ndarray[DTYPE_t, ndim=2] array, int block_size)

    Downsample an image by block summing image pixels.  This process
    conserves image flux.  If the dimensions of ``image`` are not a
    whole-multiple of ``block_size``, the extra rows/columns will not be
    included in the output downsampled image.

    This function is similar to the `scikit-image`_ `block_reduce`_
    function, which is more flexible, but slower than `downsample`.

    .. _scikit-image: http://scikit-image.org/
    .. _block_reduce: http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.block_reduce

    Parameters
    ----------
    image : array_like, float
        The 2D array of the image.

    block_size : int
        Downsampling integer factor along each axis.

    Returns
    -------
    output : array_like
        The 2D array of the resampled image.

    See Also
    --------
    upsample

    Examples
    --------
    >>> import numpy as np
    >>> from photutils import utils
    >>> img = np.arange(16.).reshape((4, 4))
    >>> utils.downsample(img, 2)
    array([[ 10.,  18.],
           [ 42.,  50.]])
    """

    cdef int nx = array.shape[1]
    cdef int ny = array.shape[0]
    cdef int nx_new = nx // block_size
    cdef int ny_new = ny // block_size
    cdef unsigned int i, j, ii, jj
    cdef np.ndarray[DTYPE_t, ndim=2] result = np.zeros([ny_new, nx_new],
                                                       dtype=DTYPE)
    for i in range(nx_new):
        for j in range(ny_new):
            for ii in range(block_size):
                for jj in range(block_size):
                    result[j, i] += array[j * block_size + jj,
                                          i * block_size + ii]
    return result


def upsample(np.ndarray[DTYPE_t, ndim=2] array, int block_size):
    """
    upsample(np.ndarray[DTYPE_t, ndim=2] array, int block_size)
    
    Upsample an image by block replicating image pixels.  To conserve
    image flux, the block-replicated image is then divided by
    ``block_size**2``.

    Parameters
    ----------
    image : array_like, float
        The 2D array of the image.

    block_size : int
        Upsampling integer factor along each axis.

    Returns
    -------
    output : array_like
        The 2D array of the resampled image.

    See Also
    --------
    downsample

    Examples
    --------
    >>> import numpy as np
    >>> from photutils import utils
    >>> img = np.array([[0., 1.], [2., 3.]])
    >>> utils.upsample(img, 2)
    array([[ 0.  ,  0.  ,  0.25,  0.25],
           [ 0.  ,  0.  ,  0.25,  0.25],
           [ 0.5 ,  0.5 ,  0.75,  0.75],
           [ 0.5 ,  0.5 ,  0.75,  0.75]])
    """

    cdef int nx = array.shape[1]
    cdef int ny = array.shape[0]
    cdef int nx_new = nx * block_size
    cdef int ny_new = ny * block_size
    cdef unsigned int i, j
    cdef np.ndarray[DTYPE_t, ndim=2] result = np.zeros((ny_new, nx_new),
                                                       dtype=DTYPE)
    cdef float block_size_sq = block_size * block_size
    for i in range(nx_new):
        for j in range(ny_new):
            result[j, i] += array[j // block_size, i // block_size]
    return result / block_size_sq
