# Licensed under a 3-clause BSD style license - see LICENSE.rst


import warnings

import numpy as np
from functools import partial

from .core import Kernel, Kernel1D, Kernel2D, MAX_NORMALIZATION
from ..utils.exceptions import AstropyUserWarning
from ..utils.console import human_file_size
from ..utils.decorators import deprecated_renamed_argument
from .. import units as u
from ..nddata import support_nddata
from ..modeling.core import _make_arithmetic_operator, BINARY_OPERATORS
from ..modeling.core import _CompoundModelMeta



# Disabling all doctests in this module until a better way of handling warnings
# in doctests can be determined
__doctest_skip__ = ['*']

BOUNDARY_OPTIONS = [None, 'fill', 'wrap', 'extend']


@support_nddata(data='array')
def convolve(array, kernel, boundary='fill', fill_value=0.,
             nan_treatment='interpolate', normalize_kernel=True, mask=None,
             preserve_nan=False, normalization_zero_tol=1e-8):
    '''
    Convolve an array with a kernel.

    This routine differs from `scipy.ndimage.convolve` because
    it includes a special treatment for ``NaN`` values. Rather than
    including ``NaN`` values in the array in the convolution calculation, which
    causes large ``NaN`` holes in the convolved array, ``NaN`` values are
    replaced with interpolated values using the kernel as an interpolation
    function.

    Parameters
    ----------
    array : `numpy.ndarray` or `~astropy.nddata.NDData`
        The array to convolve. This should be a 1, 2, or 3-dimensional array
        or a list or a set of nested lists representing a 1, 2, or
        3-dimensional array.  If an `~astropy.nddata.NDData`, the ``mask`` of
        the `~astropy.nddata.NDData` will be used as the ``mask`` argument.
    kernel : `numpy.ndarray` or `~astropy.convolution.Kernel`
        The convolution kernel. The number of dimensions should match those for
        the array, and the dimensions should be odd in all directions.  If a
        masked array, the masked values will be replaced by ``fill_value``.
    boundary : str, optional
        A flag indicating how to handle boundaries:
            * `None`
                Set the ``result`` values to zero where the kernel
                extends beyond the edge of the array.
            * 'fill'
                Set values outside the array boundary to ``fill_value`` (default).
            * 'wrap'
                Periodic boundary that wrap to the other side of ``array``.
            * 'extend'
                Set values outside the array to the nearest ``array``
                value.
    fill_value : float, optional
        The value to use outside the array when using ``boundary='fill'``
    normalize_kernel : bool, optional
        Whether to normalize the kernel to have a sum of one prior to
        convolving
    nan_treatment : 'interpolate', 'fill'
        interpolate will result in renormalization of the kernel at each
        position ignoring (pixels that are NaN in the image) in both the image
        and the kernel.
        'fill' will replace the NaN pixels with a fixed numerical value (default
        zero, see ``fill_value``) prior to convolution
        Note that if the kernel has a sum equal to zero, NaN interpolation
        is not possible and will raise an exception
    preserve_nan : bool
        After performing convolution, should pixels that were originally NaN
        again become NaN?
    mask : `None` or `numpy.ndarray`
        A "mask" array.  Shape must match ``array``, and anything that is masked
        (i.e., not 0/`False`) will be set to NaN for the convolution.  If
        `None`, no masking will be performed unless ``array`` is a masked array.
        If ``mask`` is not `None` *and* ``array`` is a masked array, a pixel is
        masked of it is masked in either ``mask`` *or* ``array.mask``.
    normalization_zero_tol: float, optional
        The absolute tolerance on whether the kernel is different than zero.
        If the kernel sums to zero to within this precision, it cannot be
        normalized. Default is "1e-8".

    Returns
    -------
    result : `numpy.ndarray`
        An array with the same dimensions and as the input array,
        convolved with kernel.  The data type depends on the input
        array type.  If array is a floating point type, then the
        return array keeps the same data type, otherwise the type
        is ``numpy.float``.

    Notes
    -----
    For masked arrays, masked values are treated as NaNs.  The convolution
    is always done at ``numpy.float`` precision.
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

    if boundary not in BOUNDARY_OPTIONS:
        raise ValueError("Invalid boundary option: must be one of {0}"
                         .format(BOUNDARY_OPTIONS))

    if nan_treatment not in ('interpolate', 'fill'):
        raise ValueError("nan_treatment must be one of 'interpolate','fill'")

    # The cython routines all need float type inputs (so, a particular
    # bit size, endianness, etc.).  So we have to convert, which also
    # has the effect of making copies so we don't modify the inputs.
    # After this, the variables we work with will be array_internal, and
    # kernel_internal.  However -- we do want to keep track of what type
    # the input array was so we can cast the result to that at the end
    # if it's a floating point type.  Don't bother with this for lists --
    # just always push those as float.
    # It is always necessary to make a copy of kernel (since it is modified),
    # but, if we just so happen to be lucky enough to have the input array
    # have exactly the desired type, we just alias to array_internal

    # Check if kernel is kernel instance
    if isinstance(kernel, Kernel):
        # Check if array is also kernel instance, if so convolve and
        # return new kernel instance
        if isinstance(array, Kernel):
            if isinstance(array, Kernel1D) and isinstance(kernel, Kernel1D):
                new_array = convolve1d_boundary_fill(array.array, kernel.array,
                                                     0, True)
                new_kernel = Kernel1D(array=new_array)
            elif isinstance(array, Kernel2D) and isinstance(kernel, Kernel2D):
                new_array = convolve2d_boundary_fill(array.array, kernel.array,
                                                     0, True)
                new_kernel = Kernel2D(array=new_array)
            else:
                raise Exception("Can't convolve 1D and 2D kernel.")
            new_kernel._separable = kernel._separable and array._separable
            new_kernel._is_bool = False
            return new_kernel
        kernel = kernel.array

    # Check that the arguments are lists or Numpy arrays

    if isinstance(array, list):
        array_internal = np.array(array, dtype=float)
        array_dtype = array_internal.dtype
    elif isinstance(array, np.ndarray):
        # Note this won't copy if it doesn't have to -- which is okay
        # because none of what follows modifies array_internal.
        array_dtype = array.dtype
        array_internal = array.astype(float, copy=False)
    else:
        raise TypeError("array should be a list or a Numpy array")

    if isinstance(kernel, list):
        kernel_internal = np.array(kernel, dtype=float)
    elif isinstance(kernel, np.ndarray):
        # Note this always makes a copy, since we will be modifying it
        kernel_internal = kernel.astype(float)
    else:
        raise TypeError("kernel should be a list or a Numpy array")

    # Check that the number of dimensions is compatible
    if array_internal.ndim != kernel_internal.ndim:
        raise Exception('array and kernel have differing number of '
                        'dimensions.')

    # anything that's masked must be turned into NaNs for the interpolation.
    # This requires copying the array_internal
    array_internal_copied = False
    if np.ma.is_masked(array):
        array_internal = array_internal.filled(np.nan)
        array_internal_copied = True
    if mask is not None:
        if not array_internal_copied:
            array_internal = array_internal.copy()
            array_internal_copied = True
        # mask != 0 yields a bool mask for all ints/floats/bool
        array_internal[mask != 0] = np.nan
    if np.ma.is_masked(kernel):
        # *kernel* doesn't support NaN interpolation, so instead we just fill it
        kernel_internal = kernel_internal.filled(fill_value)

    # Mark the NaN values so we can replace them later if interpolate_nan is
    # not set
    if preserve_nan:
        badvals = np.isnan(array_internal)

    if nan_treatment == 'fill':
        initially_nan = np.isnan(array_internal)
        array_internal[initially_nan] = fill_value

    # Because the Cython routines have to normalize the kernel on the fly, we
    # explicitly normalize the kernel here, and then scale the image at the
    # end if normalization was not requested.
    kernel_sum = kernel_internal.sum()
    kernel_sums_to_zero = np.isclose(kernel_sum, 0, atol=normalization_zero_tol)

    if (kernel_sum < 1. / MAX_NORMALIZATION or kernel_sums_to_zero) and normalize_kernel:
        raise Exception("The kernel can't be normalized, because its sum is "
                        "close to zero. The sum of the given kernel is < {0}"
                        .format(1. / MAX_NORMALIZATION))

    if not kernel_sums_to_zero:
        kernel_internal /= kernel_sum

    renormalize_by_kernel = not kernel_sums_to_zero

    if array_internal.ndim == 0:
        raise Exception("cannot convolve 0-dimensional arrays")
    elif array_internal.ndim == 1:
        if boundary == 'extend':
            result = convolve1d_boundary_extend(array_internal,
                                                kernel_internal,
                                                renormalize_by_kernel)
        elif boundary == 'fill':
            result = convolve1d_boundary_fill(array_internal,
                                              kernel_internal,
                                              float(fill_value),
                                              renormalize_by_kernel)
        elif boundary == 'wrap':
            result = convolve1d_boundary_wrap(array_internal,
                                              kernel_internal,
                                              renormalize_by_kernel)
        elif boundary is None:
            result = convolve1d_boundary_none(array_internal,
                                              kernel_internal,
                                              renormalize_by_kernel)
    elif array_internal.ndim == 2:
        if boundary == 'extend':
            result = convolve2d_boundary_extend(array_internal,
                                                kernel_internal,
                                                renormalize_by_kernel,
                                               )
        elif boundary == 'fill':
            result = convolve2d_boundary_fill(array_internal,
                                              kernel_internal,
                                              float(fill_value),
                                              renormalize_by_kernel,
                                             )
        elif boundary == 'wrap':
            result = convolve2d_boundary_wrap(array_internal,
                                              kernel_internal,
                                              renormalize_by_kernel,
                                             )
        elif boundary is None:
            result = convolve2d_boundary_none(array_internal,
                                              kernel_internal,
                                              renormalize_by_kernel,
                                             )
    elif array_internal.ndim == 3:
        if boundary == 'extend':
            result = convolve3d_boundary_extend(array_internal,
                                                kernel_internal,
                                                renormalize_by_kernel)
        elif boundary == 'fill':
            result = convolve3d_boundary_fill(array_internal,
                                              kernel_internal,
                                              float(fill_value),
                                              renormalize_by_kernel)
        elif boundary == 'wrap':
            result = convolve3d_boundary_wrap(array_internal,
                                              kernel_internal,
                                              renormalize_by_kernel)
        elif boundary is None:
            result = convolve3d_boundary_none(array_internal,
                                              kernel_internal,
                                              renormalize_by_kernel)
    else:
        raise NotImplementedError('convolve only supports 1, 2, and 3-dimensional '
                                  'arrays at this time')

    # If normalization was not requested, we need to scale the array (since
    # the kernel is effectively normalized within the cython functions)
    if not normalize_kernel and not kernel_sums_to_zero:
        result *= kernel_sum

    if preserve_nan:
        result[badvals] = np.nan

    if nan_treatment == 'fill':
        array_internal[initially_nan] = np.nan

    # Try to preserve the input type if it's a floating point type
    if array_dtype.kind == 'f':
        # Avoid making another copy if possible
        try:
            return result.astype(array_dtype, copy=False)
        except TypeError:
            return result.astype(array_dtype)
    else:
        return result


@deprecated_renamed_argument('interpolate_nan', 'nan_treatment', 'v2.0.0')
@support_nddata(data='array')
def convolve_fft(array, kernel, boundary='fill', fill_value=0.,
                 nan_treatment='interpolate', normalize_kernel=True,
                 normalization_zero_tol=1e-8,
                 preserve_nan=False, mask=None, crop=True, return_fft=False,
                 fft_pad=None, psf_pad=None, quiet=False,
                 min_wt=0.0, allow_huge=False,
                 fftn=np.fft.fftn, ifftn=np.fft.ifftn,
                 complex_dtype=complex):
    """
    Convolve an ndarray with an nd-kernel.  Returns a convolved image with
    ``shape = array.shape``.  Assumes kernel is centered.

    `convolve_fft` is very similar to `convolve` in that it replaces ``NaN``
    values in the original image with interpolated values using the kernel as
    an interpolation function.  However, it also includes many additional
    options specific to the implementation.

    `convolve_fft` differs from `scipy.signal.fftconvolve` in a few ways:

    * It can treat ``NaN`` values as zeros or interpolate over them.
    * ``inf`` values are treated as ``NaN``
    * (optionally) It pads to the nearest 2^n size to improve FFT speed.
    * Its only valid ``mode`` is 'same' (i.e., the same shape array is returned)
    * It lets you use your own fft, e.g.,
      `pyFFTW <https://pypi.python.org/pypi/pyFFTW>`_ or
      `pyFFTW3 <https://pypi.python.org/pypi/PyFFTW3/0.2.1>`_ , which can lead to
      performance improvements, depending on your system configuration.  pyFFTW3
      is threaded, and therefore may yield significant performance benefits on
      multi-core machines at the cost of greater memory requirements.  Specify
      the ``fftn`` and ``ifftn`` keywords to override the default, which is
      `numpy.fft.fft` and `numpy.fft.ifft`.

    Parameters
    ----------
    array : `numpy.ndarray`
        Array to be convolved with ``kernel``.  It can be of any
        dimensionality, though only 1, 2, and 3d arrays have been tested.
    kernel : `numpy.ndarray` or `astropy.convolution.Kernel`
        The convolution kernel. The number of dimensions should match those
        for the array.  The dimensions *do not* have to be odd in all directions,
        unlike in the non-fft `convolve` function.  The kernel will be
        normalized if ``normalize_kernel`` is set.  It is assumed to be centered
        (i.e., shifts may result if your kernel is asymmetric)
    boundary : {'fill', 'wrap'}, optional
        A flag indicating how to handle boundaries:

            * 'fill': set values outside the array boundary to fill_value
              (default)
            * 'wrap': periodic boundary

        The `None` and 'extend' parameters are not supported for FFT-based
        convolution
    fill_value : float, optional
        The value to use outside the array when using boundary='fill'
    nan_treatment : 'interpolate', 'fill'
        ``interpolate`` will result in renormalization of the kernel at each
        position ignoring (pixels that are NaN in the image) in both the image
        and the kernel.  ``fill`` will replace the NaN pixels with a fixed
        numerical value (default zero, see ``fill_value``) prior to
        convolution.  Note that if the kernel has a sum equal to zero, NaN
        interpolation is not possible and will raise an exception.
    normalize_kernel : function or boolean, optional
        If specified, this is the function to divide kernel by to normalize it.
        e.g., ``normalize_kernel=np.sum`` means that kernel will be modified to be:
        ``kernel = kernel / np.sum(kernel)``.  If True, defaults to
        ``normalize_kernel = np.sum``.
    normalization_zero_tol: float, optional
        The absolute tolerance on whether the kernel is different than zero.
        If the kernel sums to zero to within this precision, it cannot be
        normalized. Default is "1e-8".
    preserve_nan : bool
        After performing convolution, should pixels that were originally NaN
        again become NaN?
    mask : `None` or `numpy.ndarray`
        A "mask" array.  Shape must match ``array``, and anything that is masked
        (i.e., not 0/`False`) will be set to NaN for the convolution.  If
        `None`, no masking will be performed unless ``array`` is a masked array.
        If ``mask`` is not `None` *and* ``array`` is a masked array, a pixel is
        masked of it is masked in either ``mask`` *or* ``array.mask``.


    Other Parameters
    ----------------
    min_wt : float, optional
        If ignoring ``NaN`` / zeros, force all grid points with a weight less than
        this value to ``NaN`` (the weight of a grid point with *no* ignored
        neighbors is 1.0).
        If ``min_wt`` is zero, then all zero-weight points will be set to zero
        instead of ``NaN`` (which they would be otherwise, because 1/0 = nan).
        See the examples below
    fft_pad : bool, optional
        Default on.  Zero-pad image to the nearest 2^n.  With
        ``boundary='wrap'``, this will be disabled.
    psf_pad : bool, optional
        Zero-pad image to be at least the sum of the image sizes to avoid
        edge-wrapping when smoothing.  This is enabled by default with
        ``boundary='fill'``, but it can be overridden with a boolean option.
        ``boundary='wrap'`` and ``psf_pad=True`` are not compatible.
    crop : bool, optional
        Default on.  Return an image of the size of the larger of the input
        image and the kernel.
        If the image and kernel are asymmetric in opposite directions, will
        return the largest image in both directions.
        For example, if an input image has shape [100,3] but a kernel with shape
        [6,6] is used, the output will be [100,6].
    return_fft : bool, optional
        Return the ``fft(image)*fft(kernel)`` instead of the convolution (which is
        ``ifft(fft(image)*fft(kernel))``).  Useful for making PSDs.
    fftn, ifftn : functions, optional
        The fft and inverse fft functions.  Can be overridden to use your own
        ffts, e.g. an fftw3 wrapper or scipy's fftn,
        ``fft=scipy.fftpack.fftn``
    complex_dtype : numpy.complex, optional
        Which complex dtype to use.  `numpy` has a range of options, from 64 to
        256.
    quiet : bool, optional
        Silence warning message about NaN interpolation
    allow_huge : bool, optional
        Allow huge arrays in the FFT?  If False, will raise an exception if the
        array or kernel size is >1 GB

    Raises
    ------
    ValueError:
        If the array is bigger than 1 GB after padding, will raise this exception
        unless ``allow_huge`` is True

    See Also
    --------
    convolve:
        Convolve is a non-fft version of this code.  It is more memory
        efficient and for small kernels can be faster.

    Returns
    -------
    default : ndarray
        ``array`` convolved with ``kernel``.  If ``return_fft`` is set, returns
        ``fft(array) * fft(kernel)``.  If crop is not set, returns the
        image, but with the fft-padded size instead of the input size

    Notes
    -----
        With ``psf_pad=True`` and a large PSF, the resulting data can become
        very large and consume a lot of memory.  See Issue
        https://github.com/astropy/astropy/pull/4366 for further detail.

    Examples
    --------
    >>> convolve_fft([1, 0, 3], [1, 1, 1])
    array([ 1.,  4.,  3.])

    >>> convolve_fft([1, np.nan, 3], [1, 1, 1])
    array([ 1.,  4.,  3.])

    >>> convolve_fft([1, 0, 3], [0, 1, 0])
    array([ 1.,  0.,  3.])

    >>> convolve_fft([1, 2, 3], [1])
    array([ 1.,  2.,  3.])

    >>> convolve_fft([1, np.nan, 3], [0, 1, 0], nan_treatment='interpolate')
    ...
    array([ 1.,  0.,  3.])

    >>> convolve_fft([1, np.nan, 3], [0, 1, 0], nan_treatment='interpolate',
    ...              min_wt=1e-8)
    array([ 1.,  nan,  3.])

    >>> convolve_fft([1, np.nan, 3], [1, 1, 1], nan_treatment='interpolate')
    array([ 1.,  4.,  3.])

    >>> convolve_fft([1, np.nan, 3], [1, 1, 1], nan_treatment='interpolate',
    ...               normalize_kernel=True)
    array([ 1.,  2.,  3.])

    >>> import scipy.fftpack  # optional - requires scipy
    >>> convolve_fft([1, np.nan, 3], [1, 1, 1], nan_treatment='interpolate',
    ...               normalize_kernel=True,
    ...               fftn=scipy.fftpack.fft, ifftn=scipy.fftpack.ifft)
    array([ 1.,  2.,  3.])

    """
    # Checking copied from convolve.py - however, since FFTs have real &
    # complex components, we change the types.  Only the real part will be
    # returned! Note that this always makes a copy.

    # Check kernel is kernel instance
    if isinstance(kernel, Kernel):
        kernel = kernel.array
        if isinstance(array, Kernel):
            raise TypeError("Can't convolve two kernels with convolve_fft.  "
                            "Use convolve instead.")

    if nan_treatment not in ('interpolate', 'fill'):
        raise ValueError("nan_treatment must be one of 'interpolate','fill'")

    # Convert array dtype to complex
    # and ensure that list inputs become arrays
    array = np.asarray(array, dtype=complex)
    kernel = np.asarray(kernel, dtype=complex)

    # Check that the number of dimensions is compatible
    if array.ndim != kernel.ndim:
        raise ValueError("Image and kernel must have same number of "
                         "dimensions")

    arrayshape = array.shape
    kernshape = kernel.shape

    array_size_B = (np.product(arrayshape, dtype=np.int64) *
                    np.dtype(complex_dtype).itemsize)*u.byte
    if array_size_B > 1*u.GB and not allow_huge:
        raise ValueError("Size Error: Arrays will be {}.  Use "
                         "allow_huge=True to override this exception."
                         .format(human_file_size(array_size_B.to_value(u.byte))))

    # mask catching - masks must be turned into NaNs for use later in the image
    if np.ma.is_masked(array):
        mamask = array.mask
        array = np.array(array)
        array[mamask] = np.nan
    elif mask is not None:
        # copying here because we have to mask it below.  But no need to copy
        # if mask is None because we won't modify it.
        array = np.array(array)
    if mask is not None:
        # mask != 0 yields a bool mask for all ints/floats/bool
        array[mask != 0] = np.nan
    # the *kernel* doesn't support NaN interpolation, so instead we just fill it
    if np.ma.is_masked(kernel):
        kernel = kernel.filled(0)

    # NaN and inf catching
    nanmaskarray = np.isnan(array) | np.isinf(array)
    array[nanmaskarray] = 0
    nanmaskkernel = np.isnan(kernel) | np.isinf(kernel)
    kernel[nanmaskkernel] = 0

    if normalize_kernel is True:
        if kernel.sum() < 1. / MAX_NORMALIZATION:
            raise Exception("The kernel can't be normalized, because its sum is "
                            "close to zero. The sum of the given kernel is < {0}"
                            .format(1. / MAX_NORMALIZATION))
        kernel_scale = kernel.sum()
        normalized_kernel = kernel / kernel_scale
        kernel_scale = 1  # if we want to normalize it, leave it normed!
    elif normalize_kernel:
        # try this.  If a function is not passed, the code will just crash... I
        # think type checking would be better but PEPs say otherwise...
        kernel_scale = normalize_kernel(kernel)
        normalized_kernel = kernel / kernel_scale
    else:
        kernel_scale = kernel.sum()
        if np.abs(kernel_scale) < normalization_zero_tol:
            if nan_treatment == 'interpolate':
                raise ValueError('Cannot interpolate NaNs with an unnormalizable kernel')
            else:
                # the kernel's sum is near-zero, so it can't be scaled
                kernel_scale = 1
                normalized_kernel = kernel
        else:
            # the kernel is normalizable; we'll temporarily normalize it
            # now and undo the normalization later.
            normalized_kernel = kernel / kernel_scale

    if boundary is None:
        warnings.warn("The convolve_fft version of boundary=None is "
                      "equivalent to the convolve boundary='fill'.  There is "
                      "no FFT equivalent to convolve's "
                      "zero-if-kernel-leaves-boundary", AstropyUserWarning)
        if psf_pad is None:
            psf_pad = True
        if fft_pad is None:
            fft_pad = True
    elif boundary == 'fill':
        # create a boundary region at least as large as the kernel
        if psf_pad is False:
            warnings.warn("psf_pad was set to {0}, which overrides the "
                          "boundary='fill' setting.".format(psf_pad),
                          AstropyUserWarning)
        else:
            psf_pad = True
        if fft_pad is None:
            # default is 'True' according to the docstring
            fft_pad = True
    elif boundary == 'wrap':
        if psf_pad:
            raise ValueError("With boundary='wrap', psf_pad cannot be enabled.")
        psf_pad = False
        if fft_pad:
            raise ValueError("With boundary='wrap', fft_pad cannot be enabled.")
        fft_pad = False
        fill_value = 0  # force zero; it should not be used
    elif boundary == 'extend':
        raise NotImplementedError("The 'extend' option is not implemented "
                                  "for fft-based convolution")

    # find ideal size (power of 2) for fft.
    # Can add shapes because they are tuples
    if fft_pad:  # default=True
        if psf_pad:  # default=False
            # add the dimensions and then take the max (bigger)
            fsize = 2 ** np.ceil(np.log2(
                np.max(np.array(arrayshape) + np.array(kernshape))))
        else:
            # add the shape lists (max of a list of length 4) (smaller)
            # also makes the shapes square
            fsize = 2 ** np.ceil(np.log2(np.max(arrayshape + kernshape)))
        newshape = np.array([fsize for ii in range(array.ndim)], dtype=int)
    else:
        if psf_pad:
            # just add the biggest dimensions
            newshape = np.array(arrayshape) + np.array(kernshape)
        else:
            newshape = np.array([np.max([imsh, kernsh])
                                 for imsh, kernsh in zip(arrayshape, kernshape)])

    # perform a second check after padding
    array_size_C = (np.product(newshape, dtype=np.int64) *
                    np.dtype(complex_dtype).itemsize)*u.byte
    if array_size_C > 1*u.GB and not allow_huge:
        raise ValueError("Size Error: Arrays will be {}.  Use "
                         "allow_huge=True to override this exception."
                         .format(human_file_size(array_size_C)))

    # For future reference, this can be used to predict "almost exactly"
    # how much *additional* memory will be used.
    # size * (array + kernel + kernelfft + arrayfft +
    #         (kernel*array)fft +
    #         optional(weight image + weight_fft + weight_ifft) +
    #         optional(returned_fft))
    # total_memory_used_GB = (np.product(newshape)*np.dtype(complex_dtype).itemsize
    #                        * (5 + 3*((interpolate_nan or ) and kernel_is_normalized))
    #                        + (1 + (not return_fft)) *
    #                          np.product(arrayshape)*np.dtype(complex_dtype).itemsize
    #                        + np.product(arrayshape)*np.dtype(bool).itemsize
    #                        + np.product(kernshape)*np.dtype(bool).itemsize)
    #                        ) / 1024.**3

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
    arrayslices = tuple(arrayslices)
    kernslices = tuple(kernslices)

    if not np.all(newshape == arrayshape):
        if np.isfinite(fill_value):
            bigarray = np.ones(newshape, dtype=complex_dtype) * fill_value
        else:
            bigarray = np.zeros(newshape, dtype=complex_dtype)
        bigarray[arrayslices] = array
    else:
        bigarray = array

    if not np.all(newshape == kernshape):
        bigkernel = np.zeros(newshape, dtype=complex_dtype)
        bigkernel[kernslices] = normalized_kernel
    else:
        bigkernel = normalized_kernel

    arrayfft = fftn(bigarray)
    # need to shift the kernel so that, e.g., [0,0,1,0] -> [1,0,0,0] = unity
    kernfft = fftn(np.fft.ifftshift(bigkernel))
    fftmult = arrayfft * kernfft

    interpolate_nan = (nan_treatment == 'interpolate')
    if interpolate_nan:
        if not np.isfinite(fill_value):
            bigimwt = np.zeros(newshape, dtype=complex_dtype)
        else:
            bigimwt = np.ones(newshape, dtype=complex_dtype)

        bigimwt[arrayslices] = 1.0 - nanmaskarray * interpolate_nan
        wtfft = fftn(bigimwt)

        # You can only get to this point if kernel_is_normalized
        wtfftmult = wtfft * kernfft
        wtsm = ifftn(wtfftmult)
        # need to re-zero weights outside of the image (if it is padded, we
        # still don't weight those regions)
        bigimwt[arrayslices] = wtsm.real[arrayslices]
    else:
        bigimwt = 1

    if np.isnan(fftmult).any():
        # this check should be unnecessary; call it an insanity check
        raise ValueError("Encountered NaNs in convolve.  This is disallowed.")

    # restore NaNs in original image (they were modified inplace earlier)
    # We don't have to worry about masked arrays - if input was masked, it was
    # copied
    array[nanmaskarray] = np.nan
    kernel[nanmaskkernel] = np.nan

    fftmult *= kernel_scale

    if return_fft:
        return fftmult

    if interpolate_nan:
        rifft = (ifftn(fftmult)) / bigimwt
        if not np.isscalar(bigimwt):
            if min_wt > 0.:
                rifft[bigimwt < min_wt] = np.nan
            else:
                # Set anything with no weight to zero (taking into account
                # slight offsets due to floating-point errors).
                rifft[bigimwt < 10 * np.finfo(bigimwt.dtype).eps] = 0.0
    else:
        rifft = ifftn(fftmult)

    if preserve_nan:
        rifft[arrayslices][nanmaskarray] = np.nan

    if crop:
        result = rifft[arrayslices].real
        return result
    else:
        return rifft.real


def interpolate_replace_nans(array, kernel, convolve=convolve, **kwargs):
    """
    Given a data set containing NaNs, replace the NaNs by interpolating from
    neighboring data points with a given kernel.

    Parameters
    ----------
    array : `numpy.ndarray`
        Array to be convolved with ``kernel``.  It can be of any
        dimensionality, though only 1, 2, and 3d arrays have been tested.
    kernel : `numpy.ndarray` or `astropy.convolution.Kernel`
        The convolution kernel. The number of dimensions should match those
        for the array.  The dimensions *do not* have to be odd in all directions,
        unlike in the non-fft `convolve` function.  The kernel will be
        normalized if ``normalize_kernel`` is set.  It is assumed to be centered
        (i.e., shifts may result if your kernel is asymmetric).  The kernel
        *must be normalizable* (i.e., its sum cannot be zero).
    convolve : `convolve` or `convolve_fft`
        One of the two convolution functions defined in this package.

    Returns
    -------
    newarray : `numpy.ndarray`
        A copy of the original array with NaN pixels replaced with their
        interpolated counterparts
    """

    if not np.any(np.isnan(array)):
        return array.copy()

    newarray = array.copy()

    convolved = convolve(array, kernel, nan_treatment='interpolate',
                         normalize_kernel=True, **kwargs)

    isnan = np.isnan(array)
    newarray[isnan] = convolved[isnan]

    return newarray


def convolve_models(model, kernel, mode='convolve_fft', **kwargs):
    """
    Convolve two models using `~astropy.convolution.convolve_fft`.

    Parameters
    ----------
    model : `~astropy.modeling.core.Model`
        Functional model
    kernel : `~astropy.modeling.core.Model`
        Convolution kernel
    mode : str
        Keyword representing which function to use for convolution.
            * 'convolve_fft' : use `~astropy.convolution.convolve_fft` function.
            * 'convolve' : use `~astropy.convolution.convolve`.
    kwargs : dict
        Keyword arguments to me passed either to `~astropy.convolution.convolve`
        or `~astropy.convolution.convolve_fft` depending on ``mode``.

    Returns
    -------
    default : CompoundModel
        Convolved model
    """

    if mode == 'convolve_fft':
        BINARY_OPERATORS['convolve_fft'] = _make_arithmetic_operator(partial(convolve_fft, **kwargs))
    elif mode == 'convolve':
        BINARY_OPERATORS['convolve'] = _make_arithmetic_operator(partial(convolve, **kwargs))
    else:
        raise ValueError('Mode {} is not supported.'.format(mode))

    return _CompoundModelMeta._from_operator(mode, model, kernel)
