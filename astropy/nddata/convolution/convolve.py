import numpy as np


def convolve(array, kernel, boundary=None, fill_value=0.,
             normalize_kernel=False):
    '''
    Convolve an array with a kernel.

    This routine differs from `scipy.ndimage.filters.convolve` because
    it includes a special treatment for `NaN` values. Rather than
    including `NaN`s in the convolution calculation, which causes large
    `NaN` holes in the convolved image, `NaN` values are replaced with
    interpolated values using the kernel as an interpolation function.

    Parameters
    ----------
    array : `numpy.ndarray`
        The array to convolve. This should be a 1, 2, or 3-dimensional array
        or a list or a set of nested lists representing a 1, 2, or
        3-dimensional array.
    kernel : `numpy.ndarray`
        The convolution kernel. The number of dimensions should match those
        for the array, and the dimensions should be odd in all directions.
    boundary : str, optional
        A flag indicating how to handle boundaries:
            * `None`
                Set the `result` values to zero where the kernel
                extends beyond the edge of the array (default).
            * 'fill'
                Set values outside the array boundary to `fill_value`.
            * 'wrap'
                Periodic boundary that wrap to the other side of `array`.
            * 'extend'
                Set values outside the array to the nearest `array`
                value.
    fill_value : float, optional
        The value to use outside the array when using boundary='fill'
    normalize_kernel : bool, optional
        Whether to normalize the kernel prior to convolving

    Returns
    -------
    result : `numpy.ndarray`
        An array with the same dimensions and type as the input array,
        convolved with kernel.

    Notes
    -----
    Masked arrays are not supported at this time.
    '''
    from .boundary_none import convolve1d_boundary_none, \
                               convolve2d_boundary_none, \
                               convolve3d_boundary_none

    from .boundary_extend import convolve1d_boundary_extend, \
                                 convolve2d_boundary_extend, \
                                 convolve3d_boundary_extend

    from .boundary_fill import convolve1d_boundary_fill, \
                               convolve2d_boundary_fill, \
                               convolve3d_boundary_fill

    from .boundary_wrap import convolve1d_boundary_wrap, \
                               convolve2d_boundary_wrap, \
                               convolve3d_boundary_wrap

    # Check that the arguemnts are lists or Numpy arrays
    if type(array) == list:
        array = np.array(array, dtype=float)
    elif type(array) != np.ndarray:
        raise TypeError("array should be a list or a Numpy array")
    if type(kernel) == list:
        kernel = np.array(kernel, dtype=float)
    elif type(kernel) != np.ndarray:
        raise TypeError("kernel should be a list or a Numpy array")

    # Check that the number of dimensions is compatible
    if array.ndim != kernel.ndim:
        raise Exception('array and kernel have differing number of'
                        'dimensions')

    # The .dtype.type attribute returs the datatype without the endian. We can
    # use this to check that the arrays are 32- or 64-bit arrays
    if array.dtype.kind == 'i':
        array = array.astype(float)
    elif array.dtype.kind != 'f':
        raise TypeError('array should be an integer or a '
                        'floating-point Numpy array')
    if kernel.dtype.kind == 'i':
        kernel = kernel.astype(float)
    elif kernel.dtype.kind != 'f':
        raise TypeError('kernel should be an integer or a '
                        'floating-point Numpy array')

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
    elif array.ndim == 1:
        if boundary == 'extend':
            result = convolve1d_boundary_extend(array.astype(np.float),
                                                kernel.astype(np.float))
        elif boundary == 'fill':
            result = convolve1d_boundary_fill(array.astype(np.float),
                                              kernel.astype(np.float),
                                              float(fill_value))
        elif boundary == 'wrap':
            result = convolve1d_boundary_wrap(array.astype(np.float),
                                              kernel.astype(np.float))
        else:
            result = convolve1d_boundary_none(array.astype(np.float),
                                              kernel.astype(np.float))
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
    elif array.ndim == 3:
        if boundary == 'extend':
            result = convolve3d_boundary_extend(array.astype(np.float),
                                                kernel.astype(np.float))
        elif boundary == 'fill':
            result = convolve3d_boundary_fill(array.astype(np.float),
                                              kernel.astype(np.float),
                                              float(fill_value))
        elif boundary == 'wrap':
            result = convolve3d_boundary_wrap(array.astype(np.float),
                                              kernel.astype(np.float))
        else:
            result = convolve3d_boundary_none(array.astype(np.float),
                                              kernel.astype(np.float))
    else:
        raise NotImplemented('convolve only supports 1, 2, and 3-dimensional '
                             'arrays at this time')

    # If normalization was not requested, we need to scale the array (since
    # the kernel was normalized prior to convolution)
    if not normalize_kernel:
        result *= kernel_sum

    # Cast back to original dtype and return
    return result.astype(array_dtype)
