# Licensed under a 3-clause BSD style license - see PYFITS.rst

import warnings

import numpy as np
from ...config import ConfigurationItem


def convolve(array, kernel, boundary=None, fill_value=0.,
             normalize_kernel=False):
    '''
    Convolve an array with a kernel.

    This routine differs from `scipy.ndimage.filters.convolve` because
    it includes a special treatment for `NaN` values. Rather than
    including `NaNs` in the convolution calculation, which causes large
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
    from .boundary_none import (convolve1d_boundary_none,
                                convolve2d_boundary_none,
                                convolve3d_boundary_none)

    from .boundary_extend import (convolve1d_boundary_extend,
                                  convolve2d_boundary_extend,
                                  convolve3d_boundary_extend)

    from .boundary_fill import (convolve1d_boundary_fill,
                                convolve2d_boundary_fill,
                                convolve3d_boundary_fill)

    from .boundary_wrap import (convolve1d_boundary_wrap,
                                convolve2d_boundary_wrap,
                                convolve3d_boundary_wrap)

    # Check that the arguemnts are lists or Numpy arrays
    if isinstance(array, list):
        array = np.array(array, dtype=float)
    elif not isinstance(array, np.ndarray):
        raise TypeError("array should be a list or a Numpy array")
    if isinstance(kernel, list):
        kernel = np.array(kernel, dtype=float)
    elif not isinstance(kernel, np.ndarray):
        raise TypeError("kernel should be a list or a Numpy array")

    # Check that the number of dimensions is compatible
    if array.ndim != kernel.ndim:
        raise Exception('array and kernel have differing number of '
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

    # The cython routines are written for np.float, but the default endian
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


def convolve_fft(array, kernel, boundary='fill', fill_value=0, crop=True,
                 return_fft=False, fft_pad=True, psf_pad=False,
                 interpolate_nan=False, quiet=False, ignore_edge_zeros=False,
                 min_wt=0.0, normalize_kernel=False, fftn=np.fft.fftn,
                 ifftn=np.fft.ifftn, complex_dtype=np.complex):
    """
    Convolve an ndarray with an nd-kernel.  Returns a convolved image with
    shape = array.shape.  Assumes kernel is centered.

    `convolve_fft` differs from `scipy.signal.fftconvolve` in a few ways:

    * It can treat NaN's as zeros or interpolate over them.
    * `inf` values are treated as `NaN`
    * (optionally) It pads to the nearest 2^n size to improve FFT speed.
    * Its only valid `mode` is 'same' (i.e., the same shape array is returned)
    * It lets you use your own fft, e.g.,
      `pyFFTW <http://pypi.python.org/pypi/pyFFTW>`_ or
      `pyFFTW3 <http://pypi.python.org/pypi/PyFFTW3/0.2.1>`_ , which can lead to
      performance improvements, depending on your system configuration.  pyFFTW3
      is threaded, and therefore may yield significant performance benefits on
      multi-core machines at the cost of greater memory requirements.  Specify
      the `fftn` and `ifftn` keywords to override the default, which is
      `numpy.fft.fft` and `numpy.fft.ifft`.

    Parameters
    ----------
    array : `numpy.ndarray`
          Array to be convolved with `kernel`
    kernel : `numpy.ndarray`
          Will be normalized if `normalize_kernel` is set.  Assumed to be
          centered (i.e., shifts may result if your kernel is asymmetric)
    boundary : {'fill', 'wrap'}
        A flag indicating how to handle boundaries:

            * 'fill': set values outside the array boundary to fill_value
              (default)
            * 'wrap': periodic boundary

    interpolate_nan : bool
        The convolution will be re-weighted assuming NAN values are meant to be
        ignored, not treated as zero.  If this is off, all NaN values will be
        treated as zero.
    ignore_edge_zeros : bool
        Ignore the zero-pad-created zeros.  This will effectively decrease
        the kernel area on the edges but will not re-normalize the kernel.
        This parameter may result in 'edge-brightening' effects if you're using
        a normalized kernel
    min_wt : float
        If ignoring NANs/zeros, force all grid points with a weight less than
        this value to NAN (the weight of a grid point with *no* ignored
        neighbors is 1.0).
        If `min_wt` == 0.0, then all zero-weight points will be set to zero
        instead of NAN (which they would be otherwise, because 1/0 = nan).
        See the examples below
    normalize_kernel : function or boolean
        If specified, this is the function to divide kernel by to normalize it.
        e.g., normalize_kernel=np.sum means that kernel will be modified to be:
        kernel = kernel / np.sum(kernel).  If True, defaults to
        normalize_kernel = np.sum

    Other Parameters
    ----------------
    fft_pad : bool
        Default on.  Zero-pad image to the nearest 2^n
    psf_pad : bool
        Default off.  Zero-pad image to be at least the sum of the image sizes
        (in order to avoid edge-wrapping when smoothing)
    crop : bool
        Default on.  Return an image of the size of the largest input image.
        If the images are asymmetric in opposite directions, will return the
        largest image in both directions.
        For example, if an input image has shape [100,3] but a kernel with shape
        [6,6] is used, the output will be [100,6].
    return_fft : bool
        Return the fft(image)*fft(kernel) instead of the convolution (which is
        ifft(fft(image)*fft(kernel))).  Useful for making PSDs.
    nthreads : int
        if fftw3 is installed, can specify the number of threads to allow FFTs
        to use.  Probably only helpful for large arrays
    fftn, ifftn : functions
        The fft and inverse fft functions.  Can be overridden to use your own
        ffts, e.g. an fftw3 wrapper or scipy's fftn, e.g.
        `fftn=scipy.fftpack.fftn`
    complex_dtype : np.complex
        Which complex dtype to use.  `numpy` has a range of options, from 64 to
        256.

    See Also
    --------
    convolve : Convolve is a non-fft version of this code.

    Returns
    -------
    default : ndarray
        **array** convolved with `kernel`.
        If `return_fft` is set, returns fft(**array**) * fft(`kernel`).
        If crop is not set, returns the image, but with the fft-padded size
        instead of the input size

    Examples
    --------
    >>> convolve_fft([1, 0, 3], [1, 1, 1])
    array([ 1.,  4.,  3.])

    >>> convolve_fft([1, np.nan, 3], [1, 1, 1], quiet=True)
    array([ 1.,  4.,  3.])

    >>> convolve_fft([1, 0, 3], [0, 1, 0])
    array([ 1.,  0.,  3.])

    >>> convolve_fft([1, 2, 3], [1])
    array([ 1.,  2.,  3.])

    >>> convolve_fft([1, np.nan, 3], [0, 1, 0], interpolate_nan=True,
    ...               quiet=True)
    array([ 1.,  0.,  3.])

    >>> convolve_fft([1, np.nan, 3], [0, 1, 0], interpolate_nan=True,
    ...              min_wt=1e-8)
    array([ 1.,  nan,  3.])

    >>> convolve_fft([1, np.nan, 3], [1, 1, 1], interpolate_nan=True)
    array([ 1.,  4.,  3.])

    >>> convolve_fft([1, np.nan, 3], [1, 1, 1], interpolate_nan=True,
    ...               normalize_kernel=True, ignore_edge_zeros=True)
    array([ 1.,  2.,  3.])

    >>> import scipy.fftpack  # optional - requires scipy
    >>> convolve_fft([1, np.nan, 3], [1, 1, 1], interpolate_nan=True,
    ...               normalize_kernel=True, ignore_edge_zeros=True,
    ...               fftn=scipy.fftpack.fft, ifftn=scipy.fftpack.ifft)
    array([ 1.,  2.,  3.])

    """

    # Checking copied from convolve.py - however, since FFTs have real &
    # complex components, we change the types.  Only the real part will be
    # returned!
    # Check that the arguments are lists or Numpy arrays
    array = np.asarray(array, dtype=np.complex)
    kernel = np.asarray(kernel, dtype=np.complex)

    # Check that the number of dimensions is compatible
    if array.ndim != kernel.ndim:
        raise Exception('array and kernel have differing number of dimensions')

    # turn the arrays into 'complex' arrays
    if array.dtype.kind != 'c':
        array = array.astype(np.complex)
    if kernel.dtype.kind != 'c':
        kernel = kernel.astype(np.complex)

    # mask catching - masks must be turned into NaNs for use later
    if np.ma.is_masked(array):
        mask = array.mask
        array = np.array(array)
        array[mask] = np.nan
    if np.ma.is_masked(kernel):
        mask = kernel.mask
        kernel = np.array(kernel)
        kernel[mask] = np.nan

    # NAN and inf catching
    nanmaskarray = np.isnan(array) + np.isinf(array)
    array[nanmaskarray] = 0
    nanmaskkernel = np.isnan(kernel) + np.isinf(kernel)
    kernel[nanmaskkernel] = 0
    if ((nanmaskarray.sum() > 0 or nanmaskkernel.sum() > 0) and
            not interpolate_nan and not quiet):
        warnings.warn("NOT ignoring nan values even though they are present "
                      " (they are treated as 0)")

    if normalize_kernel is True:
        kernel = kernel / kernel.sum()
        kernel_is_normalized = True
    elif normalize_kernel:
        # try this.  If a function is not passed, the code will just crash... I
        # think type checking would be better but PEPs say otherwise...
        kernel = kernel / normalize_kernel(kernel)
        kernel_is_normalized = True
    else:
        if np.abs(kernel.sum() - 1) < 1e-8:
            kernel_is_normalized = True
        else:
            kernel_is_normalized = False
            if (interpolate_nan or ignore_edge_zeros):
                warnings.warn("Kernel is not normalized, therefore "
                              "ignore_edge_zeros and interpolate_nan will be "
                              "ignored.")

    if boundary is None:
        warnings.warn("The convolve_fft version of boundary=None is "
                      "equivalent to the convolve boundary='fill'.  There is "
                      "no FFT equivalent to convolve's "
                      "zero-if-kernel-leaves-boundary")
        psf_pad = True
    elif boundary == 'fill':
        # create a boundary region at least as large as the kernel
        psf_pad = True
    elif boundary == 'wrap':
        psf_pad = False
        fft_pad = False
        fill_value = 0  # force zero; it should not be used
    elif boundary == 'extend':
        raise NotImplementedError("The 'extend' option is not implemented "
                                  "for fft-based convolution")

    arrayshape = array.shape
    kernshape = kernel.shape
    if array.ndim != kernel.ndim:
        raise ValueError("Image and kernel must have same number of "
                         "dimensions")
    # find ideal size (power of 2) for fft.
    # Can add shapes because they are tuples
    if fft_pad:
        if psf_pad:
            # add the dimensions and then take the max (bigger)
            fsize = 2 ** np.ceil(np.log2(
                np.max(np.array(arrayshape) + np.array(kernshape))))
        else:
            # add the shape lists (max of a list of length 4) (smaller)
            # also makes the shapes square
            fsize = 2 ** np.ceil(np.log2(np.max(arrayshape + kernshape)))
        newshape = np.array([fsize for ii in range(array.ndim)])
    else:
        if psf_pad:
            # just add the biggest dimensions
            newshape = np.array(arrayshape) + np.array(kernshape)
        else:
            newshape = np.array([np.max([imsh, kernsh])
                                 for imsh, kernsh in zip(arrayshape, kernshape)])

    # separate each dimension by the padding size...  this is to determine the
    # appropriate slice size to get back to the input dimensions
    arrayslices = []
    kernslices = []
    for ii, (newdimsize, arraydimsize, kerndimsize) in enumerate(zip(newshape, arrayshape, kernshape)):
        center = newdimsize - (newdimsize + 1) // 2
        arrayslices += [slice(center - arraydimsize // 2,
                              center + (arraydimsize + 1) // 2)]
        kernslices += [slice(center - kerndimsize // 2,
                             center + (kerndimsize + 1) // 2)]

    bigarray = np.ones(newshape, dtype=complex_dtype) * fill_value
    bigkernel = np.zeros(newshape, dtype=complex_dtype)
    bigarray[arrayslices] = array
    bigkernel[kernslices] = kernel
    arrayfft = fftn(bigarray)
    # need to shift the kernel so that, e.g., [0,0,1,0] -> [1,0,0,0] = unity
    kernfft = fftn(np.fft.ifftshift(bigkernel))
    fftmult = arrayfft * kernfft
    if (interpolate_nan or ignore_edge_zeros) and kernel_is_normalized:
        if ignore_edge_zeros:
            bigimwt = np.zeros(newshape, dtype=complex_dtype)
        else:
            bigimwt = np.ones(newshape, dtype=complex_dtype)
        bigimwt[arrayslices] = 1.0 - nanmaskarray * interpolate_nan
        wtfft = fftn(bigimwt)
        # I think this one HAS to be normalized (i.e., the weights can't be
        # computed with a non-normalized kernel)
        wtfftmult = wtfft * kernfft / kernel.sum()
        wtsm = ifftn(wtfftmult)
        # need to re-zero weights outside of the image (if it is padded, we
        # still don't weight those regions)
        bigimwt[arrayslices] = wtsm.real[arrayslices]
        # curiously, at the floating-point limit, can get slightly negative numbers
        # they break the min_wt=0 "flag" and must therefore be removed
        bigimwt[bigimwt < 0] = 0
    else:
        bigimwt = 1

    if np.isnan(fftmult).any():
        # this check should be unnecessary; call it an insanity check
        raise ValueError("Encountered NaNs in convolve.  This is disallowed.")

    # restore nans in original image (they were modified inplace earlier)
    # We don't have to worry about masked arrays - if input was masked, it was
    # copied
    array[nanmaskarray] = np.nan
    kernel[nanmaskkernel] = np.nan

    if return_fft:
        return fftmult

    if interpolate_nan or ignore_edge_zeros:
        rifft = (ifftn(fftmult)) / bigimwt
        if not np.isscalar(bigimwt):
            rifft[bigimwt < min_wt] = np.nan
            if min_wt == 0.0:
                rifft[bigimwt == 0.0] = 0.0
    else:
        rifft = (ifftn(fftmult))

    if crop:
        result = rifft[arrayslices].real
        return result
    else:
        return rifft.real
