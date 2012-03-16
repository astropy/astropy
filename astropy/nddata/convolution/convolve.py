import numpy as np

from .boundary_none import convolve2d_boundary_none
from .boundary_extend import convolve2d_boundary_extend
from .boundary_fill import convolve2d_boundary_fill
from .boundary_wrap import convolve2d_boundary_wrap


def convolve(array, kernel, boundary=None, fill_value=0., normalize_kernel=False):
    '''
    Convolve an array with a kernel

    This routine differs from scipy.ndimage.convolve because it includes a
    special treatment for NaN values. Rather than including NaNs in the
    convolution calculation, which causes large NaN holes in the convolved
    image, NaN values are replaced with interpolated values using the kernel
    as an interpolation function.

    Parameters
    ----------
    array: np.ndarray
        The array to convolve. This should be a 2-dimensional array.
    kernel: np.ndarray
        The convolution kernel. The number of dimensions should match those
        for the array, and the dimensions should be odd in all directions.
    boundary: str, optional
        A flag indicating how to handle boundaries:
            * None : set the ``result`` values to zero where the kernel
                     extends beyond the edge of the array (default)
            * 'fill' : set values outside the array boundary to fill_value
            * 'wrap' : periodic boundary
            * 'extend' : set values outside the array to the nearest array
                         value
    fill_value: float, optional
        The value to use outside the array when using boundary='fill'
    normalize_kernel: bool, optional
        Whether to normalize the kernel prior to convolving

    Returns
    -------
    result, np.ndarray
        An array with the same dimensions and type as the input array,
        convolved with kernel.

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
        raise Exception('array and kernel have differing number of'
                        'dimensions')

    # The .dtype.type attribute returs the datatype without the endian. We can
    # use this to check that the arrays are 32- or 64-bit arrays
    if array.dtype.type not in [np.float32, np.float64]:
        raise TypeError("array should be a 32- or 64-bit Numpy array")
    if kernel.dtype.type not in [np.float32, np.float64]:
        raise TypeError("kernel should be a 32- or 64-bit Numpy array")

    # Because the Cython routines have to normalize the kernel on the fly, we
    # explicitly normalize the kernel here, and then scale the image at the
    # end if normalization was not requested.
    kernel_sum = np.sum(kernel)
    kernel /= kernel_sum

    # The cython routines are written for np.flaot, but the default endian
    # depends on platform. For that reason, we first save the original
    # array datatype, cast to np.float, then convert back
    array_dtype = array.dtype

    if array.ndim == 0:
        raise Exception("cannot convolve 0-dimensional arrays")
    elif array.ndim == 2:
        if boundary == 'extend':
            result = convolve2d_boundary_extend(array.astype(np.float),
                                                kernel.astype(np.float))
        elif boundary == 'fill':
            result = convolve2d_boundary_fill(array.astype(np.float),
                                              kernel.astype(np.float),
                                              float(fill_value))
        elif boundary == 'wrap':
            result = convolve2d_boundary_wrap(array.astype(np.float),
                                              kernel.astype(np.float))
        else:
            result = convolve2d_boundary_none(array.astype(np.float),
                                              kernel.astype(np.float))
    else:
        raise NotImplemented('convolve only supports 2-dimensional arrays'
                             'at this time')

    # If normalization was not requested, we need to scale the array (since
    # the kernel was normalized prior to convolution)
    if not normalize_kernel:
        result *= kernel_sum

    # Cast back to original dtype and return
    return result.astype(array_dtype)
