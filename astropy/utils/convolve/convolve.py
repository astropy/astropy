import numpy as np

# from .core import convolve1d_core, convolve2d_core, convolve3d_core
from .edge_default import convolve2d_edge_default
from .edge_truncate import convolve2d_edge_truncate
from .edge_fill import convolve2d_edge_fill
from .edge_wrap import convolve2d_edge_wrap


def convolve(array, kernel, edge=None, fill_value=0.):
    '''
    Convolve an array with a Kernel

    This routine differs from scipy.ndimage.convolve because it includes a
    special treatment for NaN values. Rather than including NaNs in the
    convolution calculation, which causes large NaN holes in the convolved
    image, NaN values are ignored in the summations.

    Parameters
    ----------
    array, np.ndarray
        The array to convolve. This should be a 2-dimensional array.
    kernel, np.ndarray
        The convolution kernel. The number of dimensions should match those
        for the array, and the number of dimensions should be odd in all
        directions.

    Returns
    -------
    result, np.ndarray
        An array with the same dimensions as the input array, convolved with
        kernel.

    Notes
    -----
    Masked arrays are not supported at this time.
    '''

    # Check that the arguemnts are Numpy arrays
    if type(array) != np.ndarray:
        raise TypeError("array should be a Numpy array")
    if type(kernel) != np.ndarray:
        raise TypeError("kernel should be a Numpy array")

    # Check that the number of dimensions is compatible
    if array.ndim != kernel.ndim:
        raise Exception("array and kernel have differing number of dimensions")

    # The .dtype.type attribute returs the datatype without the endian. We can
    # use this to check that the arrays are 32- or 64-bit arrays
    if array.dtype.type not in [np.float32, np.float64]:
        raise TypeError("array should be a 32- or 64-bit Numpy array")
    if kernel.dtype.type not in [np.float32, np.float64]:
        raise TypeError("kernel should be a 32- or 64-bit Numpy array")

    # The cython routines are written for np.flaot, but the default endian
    # depends on platform. For that reason, we first save the original
    # array datatype, cast to np.float, then convert back
    array_dtype = array.dtype

    if array.ndim == 0:
        raise Exception("cannot convolve 0-dimensional arrays")
    elif array.ndim == 2:
        if edge == 'truncate':
            result = convolve2d_edge_truncate(array.astype(np.float),
                                              kernel.astype(np.float))
        elif edge == 'fill':
            result = convolve2d_edge_fill(array.astype(np.float),
                                          kernel.astype(np.float),
                                          float(fill_value))
        elif edge == 'wrap':
            result = convolve2d_edge_wrap(array.astype(np.float),
                                          kernel.astype(np.float))
        else:
            result = convolve2d_edge_default(array.astype(np.float),
                                             kernel.astype(np.float))
    else:
        raise NotImplemented("convolve only supports 2-dimensional arrays at this time")

    # Cast back to original dtype and return
    return result.astype(array_dtype)
