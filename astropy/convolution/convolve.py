# Licensed under a 3-clause BSD style license - see LICENSE.rst

import ctypes
import os
import warnings
from functools import partial

import numpy as np
from numpy.ctypeslib import load_library, ndpointer

from astropy import units as u
from astropy.modeling.convolution import Convolution
from astropy.modeling.core import SPECIAL_OPERATORS, CompoundModel
from astropy.nddata import support_nddata
from astropy.utils.console import human_file_size
from astropy.utils.exceptions import AstropyUserWarning

from .core import MAX_NORMALIZATION, Kernel, Kernel1D, Kernel2D
from .utils import KernelSizeError, has_even_axis, raise_even_kernel_exception

LIBRARY_PATH = os.path.dirname(__file__)

try:
    with warnings.catch_warnings():
        # numpy.distutils is deprecated since numpy 1.23
        # see https://github.com/astropy/astropy/issues/12865
        warnings.simplefilter('ignore', DeprecationWarning)
        _convolve = load_library("_convolve", LIBRARY_PATH)
except Exception:
    raise ImportError("Convolution C extension is missing. Try re-building astropy.")

# The GIL is automatically released by default when calling functions imported
# from libraries loaded by ctypes.cdll.LoadLibrary(<path>)

# Declare prototypes
# Boundary None
_convolveNd_c = _convolve.convolveNd_c
_convolveNd_c.restype = None
_convolveNd_c.argtypes = [ndpointer(ctypes.c_double, flags={"C_CONTIGUOUS", "WRITEABLE"}),  # return array
                          ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # input array
                          ctypes.c_uint,  # N dim
                          # size array for input and result unless
                          # embed_result_within_padded_region is False,
                          # in which case the result array is assumed to be
                          # input.shape - 2*(kernel.shape//2). Note: integer division.
                          ndpointer(ctypes.c_size_t, flags="C_CONTIGUOUS"),
                          ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # kernel array
                          ndpointer(ctypes.c_size_t, flags="C_CONTIGUOUS"),  # size array for kernel
                          ctypes.c_bool,  # nan_interpolate
                          ctypes.c_bool,  # embed_result_within_padded_region
                          ctypes.c_uint]  # n_threads

# np.unique([scipy.fft.next_fast_len(i, real=True) for i in range(10000)])
_good_sizes = np.array([   0,    1,    2,    3,    4,    5,    6,    8,    9,   10,   12,  # noqa E201
                          15,   16,   18,   20,   24,   25,   27,   30,   32,   36,   40,
                          45,   48,   50,   54,   60,   64,   72,   75,   80,   81,   90,
                          96,  100,  108,  120,  125,  128,  135,  144,  150,  160,  162,
                         180,  192,  200,  216,  225,  240,  243,  250,  256,  270,  288,
                         300,  320,  324,  360,  375,  384,  400,  405,  432,  450,  480,
                         486,  500,  512,  540,  576,  600,  625,  640,  648,  675,  720,
                         729,  750,  768,  800,  810,  864,  900,  960,  972, 1000, 1024,
                        1080, 1125, 1152, 1200, 1215, 1250, 1280, 1296, 1350, 1440, 1458,
                        1500, 1536, 1600, 1620, 1728, 1800, 1875, 1920, 1944, 2000, 2025,
                        2048, 2160, 2187, 2250, 2304, 2400, 2430, 2500, 2560, 2592, 2700,
                        2880, 2916, 3000, 3072, 3125, 3200, 3240, 3375, 3456, 3600, 3645,
                        3750, 3840, 3888, 4000, 4050, 4096, 4320, 4374, 4500, 4608, 4800,
                        4860, 5000, 5120, 5184, 5400, 5625, 5760, 5832, 6000, 6075, 6144,
                        6250, 6400, 6480, 6561, 6750, 6912, 7200, 7290, 7500, 7680, 7776,
                        8000, 8100, 8192, 8640, 8748, 9000, 9216, 9375, 9600, 9720, 10000])
_good_range = int(np.log10(_good_sizes[-1]))

# Disabling doctests when scipy isn't present.
__doctest_requires__ = {('convolve_fft',): ['scipy.fft']}

BOUNDARY_OPTIONS = [None, 'fill', 'wrap', 'extend']


def _next_fast_lengths(shape):
    """
    Find optimal or good sizes to pad an array of ``shape`` to for better
    performance with `numpy.fft.*fft` and `scipy.fft.*fft`.
    Calculated directly with `scipy.fft.next_fast_len`, if available; otherwise
    looked up from list and scaled by powers of 10, if necessary.
    """

    try:
        import scipy.fft
        return np.array([scipy.fft.next_fast_len(j) for j in shape])
    except ImportError:
        pass

    newshape = np.empty(len(np.atleast_1d(shape)), dtype=int)
    for i, j in enumerate(shape):
        scale = 10 ** max(int(np.ceil(np.log10(j))) - _good_range, 0)
        for n in _good_sizes:
            if n * scale >= j:
                newshape[i] = n * scale
                break
        else:
            raise ValueError(f'No next fast length for {j} found in list of _good_sizes '
                             f'<= {_good_sizes[-1] * scale}.')
    return newshape


def _copy_input_if_needed(input, dtype=float, order='C', nan_treatment=None,
                          mask=None, fill_value=None):
    # Alias input
    input = input.array if isinstance(input, Kernel) else input
    # strip quantity attributes
    if hasattr(input, 'unit'):
        input = input.value
    output = input
    # Copy input
    try:
        # Anything that's masked must be turned into NaNs for the interpolation.
        # This requires copying. A copy is also needed for nan_treatment == 'fill'
        # A copy prevents possible function side-effects of the input array.
        if nan_treatment == 'fill' or np.ma.is_masked(input) or mask is not None:
            if np.ma.is_masked(input):
                # ``np.ma.maskedarray.filled()`` returns a copy, however there
                # is no way to specify the return type or order etc. In addition
                # ``np.nan`` is a ``float`` and there is no conversion to an
                # ``int`` type. Therefore, a pre-fill copy is needed for non
                # ``float`` masked arrays. ``subok=True`` is needed to retain
                # ``np.ma.maskedarray.filled()``. ``copy=False`` allows the fill
                # to act as the copy if type and order are already correct.
                output = np.array(input, dtype=dtype, copy=False, order=order, subok=True)
                output = output.filled(fill_value)
            else:
                # Since we're making a copy, we might as well use `subok=False` to save,
                # what is probably, a negligible amount of memory.
                output = np.array(input, dtype=dtype, copy=True, order=order, subok=False)

            if mask is not None:
                # mask != 0 yields a bool mask for all ints/floats/bool
                output[mask != 0] = fill_value
        else:
            # The call below is synonymous with np.asanyarray(array, ftype=float, order='C')
            # The advantage of `subok=True` is that it won't copy when array is an ndarray subclass. If it
            # is and `subok=False` (default), then it will copy even if `copy=False`. This uses less memory
            # when ndarray subclasses are passed in.
            output = np.array(input, dtype=dtype, copy=False, order=order, subok=True)
    except (TypeError, ValueError) as e:
        raise TypeError('input should be a Numpy array or something '
                        'convertible into a float array', e)
    return output


@support_nddata(data='array')
def convolve(array, kernel, boundary='fill', fill_value=0.,
             nan_treatment='interpolate', normalize_kernel=True, mask=None,
             preserve_nan=False, normalization_zero_tol=1e-8):
    """
    Convolve an array with a kernel.

    This routine differs from `scipy.ndimage.convolve` because
    it includes a special treatment for ``NaN`` values. Rather than
    including ``NaN`` values in the array in the convolution calculation, which
    causes large ``NaN`` holes in the convolved array, ``NaN`` values are
    replaced with interpolated values using the kernel as an interpolation
    function.

    Parameters
    ----------
    array : `~astropy.nddata.NDData` or array-like
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
        The value to use outside the array when using ``boundary='fill'``.
    normalize_kernel : bool, optional
        Whether to normalize the kernel to have a sum of one.
    nan_treatment : {'interpolate', 'fill'}, optional
        The method used to handle NaNs in the input ``array``:
            * ``'interpolate'``: ``NaN`` values are replaced with
              interpolated values using the kernel as an interpolation
              function. Note that if the kernel has a sum equal to
              zero, NaN interpolation is not possible and will raise an
              exception.
            * ``'fill'``: ``NaN`` values are replaced by ``fill_value``
              prior to convolution.
    preserve_nan : bool, optional
        After performing convolution, should pixels that were originally NaN
        again become NaN?
    mask : None or ndarray, optional
        A "mask" array.  Shape must match ``array``, and anything that is masked
        (i.e., not 0/`False`) will be set to NaN for the convolution.  If
        `None`, no masking will be performed unless ``array`` is a masked array.
        If ``mask`` is not `None` *and* ``array`` is a masked array, a pixel is
        masked of it is masked in either ``mask`` *or* ``array.mask``.
    normalization_zero_tol : float, optional
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
    """

    if boundary not in BOUNDARY_OPTIONS:
        raise ValueError(f"Invalid boundary option: must be one of {BOUNDARY_OPTIONS}")

    if nan_treatment not in ('interpolate', 'fill'):
        raise ValueError("nan_treatment must be one of 'interpolate','fill'")

    # OpenMP support is disabled at the C src code level, changing this will have
    # no effect.
    n_threads = 1

    # Keep refs to originals
    passed_kernel = kernel
    passed_array = array

    # The C routines all need float type inputs (so, a particular
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
    # Convert kernel to ndarray if not already

    # Copy or alias array to array_internal
    array_internal = _copy_input_if_needed(passed_array, dtype=float, order='C',
                                           nan_treatment=nan_treatment, mask=mask,
                                           fill_value=np.nan)
    array_dtype = getattr(passed_array, 'dtype', array_internal.dtype)
    # Copy or alias kernel to kernel_internal
    kernel_internal = _copy_input_if_needed(passed_kernel, dtype=float, order='C',
                                            nan_treatment=None, mask=None,
                                            fill_value=fill_value)

    # Make sure kernel has all odd axes
    if has_even_axis(kernel_internal):
        raise_even_kernel_exception()

    # If both image array and kernel are Kernel instances
    # constrain convolution method
    # This must occur before the main alias/copy of ``passed_kernel`` to
    # ``kernel_internal`` as it is used for filling masked kernels.
    if isinstance(passed_array, Kernel) and isinstance(passed_kernel, Kernel):
        warnings.warn("Both array and kernel are Kernel instances, hardwiring "
                      "the following parameters: boundary='fill', fill_value=0,"
                      " normalize_Kernel=True, nan_treatment='interpolate'",
                      AstropyUserWarning)
        boundary = 'fill'
        fill_value = 0
        normalize_kernel = True
        nan_treatment = 'interpolate'

    # -----------------------------------------------------------------------
    # From this point onwards refer only to ``array_internal`` and
    # ``kernel_internal``.
    # Assume both are base np.ndarrays and NOT subclasses e.g. NOT
    # ``Kernel`` nor ``np.ma.maskedarray`` classes.
    # -----------------------------------------------------------------------

    # Check dimensionality
    if array_internal.ndim == 0:
        raise Exception("cannot convolve 0-dimensional arrays")
    elif array_internal.ndim > 3:
        raise NotImplementedError('convolve only supports 1, 2, and 3-dimensional '
                                  'arrays at this time')
    elif array_internal.ndim != kernel_internal.ndim:
        raise Exception('array and kernel have differing number of '
                        'dimensions.')

    array_shape = np.array(array_internal.shape)
    kernel_shape = np.array(kernel_internal.shape)
    pad_width = kernel_shape//2

    # For boundary=None only the center space is convolved. All array indices within a
    # distance kernel.shape//2 from the edge are completely ignored (zeroed).
    # E.g. (1D list) only the indices len(kernel)//2 : len(array)-len(kernel)//2
    # are convolved. It is therefore not possible to use this method to convolve an
    # array by a kernel that is larger (see note below) than the array - as ALL pixels would be ignored
    # leaving an array of only zeros.
    # Note: For even kernels the correctness condition is array_shape > kernel_shape.
    # For odd kernels it is:
    # array_shape >= kernel_shape OR array_shape > kernel_shape-1 OR array_shape > 2*(kernel_shape//2).
    # Since the latter is equal to the former two for even lengths, the latter condition is complete.
    if boundary is None and not np.all(array_shape > 2*pad_width):
        raise KernelSizeError("for boundary=None all kernel axes must be smaller than array's - "
                              "use boundary in ['fill', 'extend', 'wrap'] instead.")

    # NaN interpolation significantly slows down the C convolution
    # computation. Since nan_treatment = 'interpolate', is the default
    # check whether it is even needed, if not, don't interpolate.
    # NB: np.isnan(array_internal.sum()) is faster than np.isnan(array_internal).any()
    nan_interpolate = (nan_treatment == 'interpolate') and np.isnan(array_internal.sum())

    # Check if kernel is normalizable
    if normalize_kernel or nan_interpolate:
        kernel_sum = kernel_internal.sum()
        kernel_sums_to_zero = np.isclose(kernel_sum, 0,
                                         atol=normalization_zero_tol)

        if kernel_sum < 1. / MAX_NORMALIZATION or kernel_sums_to_zero:
            if nan_interpolate:
                raise ValueError("Setting nan_treatment='interpolate' "
                                 "requires the kernel to be normalized, "
                                 "but the input kernel has a sum close "
                                 "to zero. For a zero-sum kernel and "
                                 "data with NaNs, set nan_treatment='fill'.")
            else:
                raise ValueError("The kernel can't be normalized, because "
                                 "its sum is close to zero. The sum of the "
                                 "given kernel is < {}"
                                 .format(1. / MAX_NORMALIZATION))

    # Mark the NaN values so we can replace them later if interpolate_nan is
    # not set
    if preserve_nan or nan_treatment == 'fill':
        initially_nan = np.isnan(array_internal)
        if nan_treatment == 'fill':
            array_internal[initially_nan] = fill_value

    # Avoid any memory allocation within the C code. Allocate output array
    # here and pass through instead.
    result = np.zeros(array_internal.shape, dtype=float, order='C')

    embed_result_within_padded_region = True
    array_to_convolve = array_internal
    if boundary in ('fill', 'extend', 'wrap'):
        embed_result_within_padded_region = False
        if boundary == 'fill':
            # This method is faster than using numpy.pad(..., mode='constant')
            array_to_convolve = np.full(array_shape + 2*pad_width, fill_value=fill_value, dtype=float, order='C')
            # Use bounds [pad_width[0]:array_shape[0]+pad_width[0]] instead of [pad_width[0]:-pad_width[0]]
            # to account for when the kernel has size of 1 making pad_width = 0.
            if array_internal.ndim == 1:
                array_to_convolve[pad_width[0]:array_shape[0]+pad_width[0]] = array_internal
            elif array_internal.ndim == 2:
                array_to_convolve[pad_width[0]:array_shape[0]+pad_width[0],
                                  pad_width[1]:array_shape[1]+pad_width[1]] = array_internal
            else:
                array_to_convolve[pad_width[0]:array_shape[0]+pad_width[0],
                                  pad_width[1]:array_shape[1]+pad_width[1],
                                  pad_width[2]:array_shape[2]+pad_width[2]] = array_internal
        else:
            np_pad_mode_dict = {'fill': 'constant', 'extend': 'edge', 'wrap': 'wrap'}
            np_pad_mode = np_pad_mode_dict[boundary]
            pad_width = kernel_shape // 2

            if array_internal.ndim == 1:
                np_pad_width = (pad_width[0],)
            elif array_internal.ndim == 2:
                np_pad_width = ((pad_width[0],), (pad_width[1],))
            else:
                np_pad_width = ((pad_width[0],), (pad_width[1],), (pad_width[2],))

            array_to_convolve = np.pad(array_internal, pad_width=np_pad_width,
                                       mode=np_pad_mode)

    _convolveNd_c(result, array_to_convolve,
                  array_to_convolve.ndim,
                  np.array(array_to_convolve.shape, dtype=ctypes.c_size_t, order='C'),
                  kernel_internal,
                  np.array(kernel_shape, dtype=ctypes.c_size_t, order='C'),
                  nan_interpolate, embed_result_within_padded_region,
                  n_threads)

    # So far, normalization has only occurred for nan_treatment == 'interpolate'
    # because this had to happen within the C extension so as to ignore
    # any NaNs
    if normalize_kernel:
        if not nan_interpolate:
            result /= kernel_sum
    elif nan_interpolate:
        result *= kernel_sum

    if nan_interpolate and not preserve_nan and np.isnan(result.sum()):
        warnings.warn("nan_treatment='interpolate', however, NaN values detected "
                      "post convolution. A contiguous region of NaN values, larger "
                      "than the kernel size, are present in the input array. "
                      "Increase the kernel size to avoid this.", AstropyUserWarning)

    if preserve_nan:
        result[initially_nan] = np.nan

    # Convert result to original data type
    array_unit = getattr(passed_array, "unit", None)
    if array_unit is not None:
        result <<= array_unit

    if isinstance(passed_array, Kernel):
        if isinstance(passed_array, Kernel1D):
            new_result = Kernel1D(array=result)
        elif isinstance(passed_array, Kernel2D):
            new_result = Kernel2D(array=result)
        else:
            raise TypeError("Only 1D and 2D Kernels are supported.")
        new_result._is_bool = False
        new_result._separable = passed_array._separable
        if isinstance(passed_kernel, Kernel):
            new_result._separable = new_result._separable and passed_kernel._separable
        return new_result
    elif array_dtype.kind == 'f':
        # Try to preserve the input type if it's a floating point type
        # Avoid making another copy if possible
        try:
            return result.astype(array_dtype, copy=False)
        except TypeError:
            return result.astype(array_dtype)
    else:
        return result


@support_nddata(data='array')
def convolve_fft(array, kernel, boundary='fill', fill_value=0.,
                 nan_treatment='interpolate', normalize_kernel=True,
                 normalization_zero_tol=1e-8,
                 preserve_nan=False, mask=None, crop=True, return_fft=False,
                 fft_pad=None, psf_pad=None, min_wt=0.0, allow_huge=False,
                 fftn=np.fft.fftn, ifftn=np.fft.ifftn,
                 complex_dtype=complex, dealias=False):
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
    * It optionally pads to the nearest faster sizes to improve FFT speed.
      These sizes are optimized for the numpy and scipy implementations, and
      ``fftconvolve`` uses them by default as well; when using other external
      functions (see below), results may vary.
    * Its only valid ``mode`` is 'same' (i.e., the same shape array is returned)
    * It lets you use your own fft, e.g.,
      `pyFFTW <https://pypi.org/project/pyFFTW/>`_ or
      `pyFFTW3 <https://pypi.org/project/PyFFTW3/0.2.1/>`_ , which can lead to
      performance improvements, depending on your system configuration.  pyFFTW3
      is threaded, and therefore may yield significant performance benefits on
      multi-core machines at the cost of greater memory requirements.  Specify
      the ``fftn`` and ``ifftn`` keywords to override the default, which is
      `numpy.fft.fftn` and `numpy.fft.ifftn`.  The `scipy.fft` functions also
      offer somewhat better performance and a multi-threaded option.

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
        convolution.
    fill_value : float, optional
        The value to use outside the array when using boundary='fill'.
    nan_treatment : {'interpolate', 'fill'}, optional
        The method used to handle NaNs in the input ``array``:
            * ``'interpolate'``: ``NaN`` values are replaced with
              interpolated values using the kernel as an interpolation
              function. Note that if the kernel has a sum equal to
              zero, NaN interpolation is not possible and will raise an
              exception.
            * ``'fill'``: ``NaN`` values are replaced by ``fill_value``
              prior to convolution.
    normalize_kernel : callable or boolean, optional
        If specified, this is the function to divide kernel by to normalize it.
        e.g., ``normalize_kernel=np.sum`` means that kernel will be modified to be:
        ``kernel = kernel / np.sum(kernel)``.  If True, defaults to
        ``normalize_kernel = np.sum``.
    normalization_zero_tol : float, optional
        The absolute tolerance on whether the kernel is different than zero.
        If the kernel sums to zero to within this precision, it cannot be
        normalized. Default is "1e-8".
    preserve_nan : bool, optional
        After performing convolution, should pixels that were originally NaN
        again become NaN?
    mask : None or ndarray, optional
        A "mask" array.  Shape must match ``array``, and anything that is masked
        (i.e., not 0/`False`) will be set to NaN for the convolution.  If
        `None`, no masking will be performed unless ``array`` is a masked array.
        If ``mask`` is not `None` *and* ``array`` is a masked array, a pixel is
        masked of it is masked in either ``mask`` *or* ``array.mask``.
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
    fft_pad : bool, optional
        Default on.  Zero-pad image to the nearest size supporting more efficient
        execution of the FFT, generally values factorizable into the first 3-5
        prime numbers.  With ``boundary='wrap'``, this will be disabled.
    psf_pad : bool, optional
        Zero-pad image to be at least the sum of the image sizes to avoid
        edge-wrapping when smoothing.  This is enabled by default with
        ``boundary='fill'``, but it can be overridden with a boolean option.
        ``boundary='wrap'`` and ``psf_pad=True`` are not compatible.
    min_wt : float, optional
        If ignoring ``NaN`` / zeros, force all grid points with a weight less than
        this value to ``NaN`` (the weight of a grid point with *no* ignored
        neighbors is 1.0).
        If ``min_wt`` is zero, then all zero-weight points will be set to zero
        instead of ``NaN`` (which they would be otherwise, because 1/0 = nan).
        See the examples below.
    allow_huge : bool, optional
        Allow huge arrays in the FFT?  If False, will raise an exception if the
        array or kernel size is >1 GB.
    fftn : callable, optional
        The fft function.  Can be overridden to use your own ffts,
        e.g. an fftw3 wrapper or scipy's fftn, ``fft=scipy.fftpack.fftn``.
    ifftn : callable, optional
        The inverse fft function. Can be overridden the same way ``fttn``.
    complex_dtype : complex type, optional
        Which complex dtype to use.  `numpy` has a range of options, from 64 to
        256.
    dealias: bool, optional
        Default off. Zero-pad image to enable explicit dealiasing
        of convolution. With ``boundary='wrap'``, this will be disabled.
        Note that for an input of nd dimensions this will increase
        the size of the temporary arrays by at least ``1.5**nd``.
        This may result in significantly more memory usage.

    Returns
    -------
    default : ndarray
        ``array`` convolved with ``kernel``.  If ``return_fft`` is set, returns
        ``fft(array) * fft(kernel)``.  If crop is not set, returns the
        image, but with the fft-padded size instead of the input size.

    Raises
    ------
    `ValueError`
        If the array is bigger than 1 GB after padding, will raise this
        exception unless ``allow_huge`` is True.

    See Also
    --------
    convolve:
        Convolve is a non-fft version of this code.  It is more memory
        efficient and for small kernels can be faster.

    Notes
    -----
    With ``psf_pad=True`` and a large PSF, the resulting data
    can become large and consume a lot of memory. See Issue
    https://github.com/astropy/astropy/pull/4366 and the update in
    https://github.com/astropy/astropy/pull/11533 for further details.

    Dealiasing of pseudospectral convolutions is necessary for
    numerical stability of the underlying algorithms. A common
    method for handling this is to zero pad the image by at least
    1/2 to eliminate the wavenumbers which have been aliased
    by convolution. This is so that the aliased 1/3 of the
    results of the convolution computation can be thrown out. See
    https://doi.org/10.1175/1520-0469(1971)028%3C1074:OTEOAI%3E2.0.CO;2
    https://iopscience.iop.org/article/10.1088/1742-6596/318/7/072037

    Note that if dealiasing is necessary to your application, but your
    process is memory constrained, you may want to consider using
    FFTW++: https://github.com/dealias/fftwpp. It includes python
    wrappers for a pseudospectral convolution which will implicitly
    dealias your convolution without the need for additional padding.
    Note that one cannot use FFTW++'s convlution directly in this
    method as in handles the entire convolution process internally.
    Additionally, FFTW++ includes other useful pseudospectral methods to
    consider.

    Examples
    --------
    >>> convolve_fft([1, 0, 3], [1, 1, 1])
    array([0.33333333, 1.33333333, 1.        ])

    >>> convolve_fft([1, np.nan, 3], [1, 1, 1])
    array([0.5, 2. , 1.5])

    >>> convolve_fft([1, 0, 3], [0, 1, 0])  # doctest: +FLOAT_CMP
    array([ 1.00000000e+00, -3.70074342e-17,  3.00000000e+00])

    >>> convolve_fft([1, 2, 3], [1])
    array([1., 2., 3.])

    >>> convolve_fft([1, np.nan, 3], [0, 1, 0], nan_treatment='interpolate')
    array([1., 0., 3.])

    >>> convolve_fft([1, np.nan, 3], [0, 1, 0], nan_treatment='interpolate',
    ...              min_wt=1e-8)
    array([ 1., nan,  3.])

    >>> convolve_fft([1, np.nan, 3], [1, 1, 1], nan_treatment='interpolate')
    array([0.5, 2. , 1.5])

    >>> convolve_fft([1, np.nan, 3], [1, 1, 1], nan_treatment='interpolate',
    ...               normalize_kernel=True)
    array([0.5, 2. , 1.5])

    >>> import scipy.fft  # optional - requires scipy
    >>> convolve_fft([1, np.nan, 3], [1, 1, 1], nan_treatment='interpolate',
    ...               normalize_kernel=True,
    ...               fftn=scipy.fft.fftn, ifftn=scipy.fft.ifftn)
    array([0.5, 2. , 1.5])

    >>> fft_mp = lambda a: scipy.fft.fftn(a, workers=-1)  # use all available cores
    >>> ifft_mp = lambda a: scipy.fft.ifftn(a, workers=-1)
    >>> convolve_fft([1, np.nan, 3], [1, 1, 1], nan_treatment='interpolate',
    ...               normalize_kernel=True, fftn=fft_mp, ifftn=ifft_mp)
    array([0.5, 2. , 1.5])
    """
    # Checking copied from convolve.py - however, since FFTs have real &
    # complex components, we change the types.  Only the real part will be
    # returned! Note that this always makes a copy.

    # Check kernel is kernel instance
    if isinstance(kernel, Kernel):
        kernel = kernel.array
        if isinstance(array, Kernel):
            raise TypeError("Can't convolve two kernels with convolve_fft.  Use convolve instead.")

    if nan_treatment not in ('interpolate', 'fill'):
        raise ValueError("nan_treatment must be one of 'interpolate','fill'")

    # Get array quantity if it exists
    array_unit = getattr(array, "unit", None)

    # Convert array dtype to complex
    # and ensure that list inputs become arrays
    array = _copy_input_if_needed(array, dtype=complex, order='C',
                                  nan_treatment=nan_treatment, mask=mask,
                                  fill_value=np.nan)
    kernel = _copy_input_if_needed(kernel, dtype=complex, order='C',
                                   nan_treatment=None, mask=None,
                                   fill_value=0)

    # Check that the number of dimensions is compatible
    if array.ndim != kernel.ndim:
        raise ValueError("Image and kernel must have same number of dimensions")

    arrayshape = array.shape
    kernshape = kernel.shape

    array_size_B = (np.product(arrayshape, dtype=np.int64) *
                    np.dtype(complex_dtype).itemsize) * u.byte
    if array_size_B > 1 * u.GB and not allow_huge:
        raise ValueError(f"Size Error: Arrays will be {human_file_size(array_size_B)}.  "
                         f"Use allow_huge=True to override this exception.")

    # NaN and inf catching
    nanmaskarray = np.isnan(array) | np.isinf(array)
    if nan_treatment == 'fill':
        array[nanmaskarray] = fill_value
    else:
        array[nanmaskarray] = 0
    nanmaskkernel = np.isnan(kernel) | np.isinf(kernel)
    kernel[nanmaskkernel] = 0

    if normalize_kernel is True:
        if kernel.sum() < 1. / MAX_NORMALIZATION:
            raise Exception("The kernel can't be normalized, because its sum is "
                            "close to zero. The sum of the given kernel is < {}"
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
            warnings.warn(f"psf_pad was set to {psf_pad}, which overrides the "
                          f"boundary='fill' setting.", AstropyUserWarning)
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
        if dealias:
            raise ValueError("With boundary='wrap', dealias cannot be enabled.")
        fill_value = 0  # force zero; it should not be used
    elif boundary == 'extend':
        raise NotImplementedError("The 'extend' option is not implemented "
                                  "for fft-based convolution")

    # Add shapes elementwise for psf_pad.
    if psf_pad:  # default=False
        # add the sizes along each dimension (bigger)
        newshape = np.array(arrayshape) + np.array(kernshape)
    else:
        # take the larger shape in each dimension (smaller)
        newshape = np.maximum(arrayshape, kernshape)

    if dealias:
        # Extend shape by 1/2 for dealiasing
        newshape += np.ceil(newshape / 2).astype(int)

    # Find ideal size for fft (was power of 2, now any powers of prime factors 2, 3, 5).
    if fft_pad:  # default=True
        # Get optimized sizes from scipy.
        newshape = _next_fast_lengths(newshape)

    # perform a second check after padding
    array_size_C = (np.product(newshape, dtype=np.int64) *
                    np.dtype(complex_dtype).itemsize) * u.byte
    if array_size_C > 1 * u.GB and not allow_huge:
        raise ValueError(f"Size Error: Arrays will be {human_file_size(array_size_C)}.  "
                         f"Use allow_huge=True to override this exception.")

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

    fftmult *= kernel_scale

    if array_unit is not None:
        fftmult <<= array_unit

    if return_fft:
        return fftmult

    if interpolate_nan:
        with np.errstate(divide='ignore', invalid='ignore'):
            # divide by zeros are expected here; if the weight is zero, we want
            # the output to be nan or inf
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
                         normalize_kernel=True, preserve_nan=False, **kwargs)

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
    **kwargs : dict
        Keyword arguments to me passed either to `~astropy.convolution.convolve`
        or `~astropy.convolution.convolve_fft` depending on ``mode``.

    Returns
    -------
    default : `~astropy.modeling.core.CompoundModel`
        Convolved model
    """

    if mode == 'convolve_fft':
        operator = SPECIAL_OPERATORS.add('convolve_fft', partial(convolve_fft, **kwargs))
    elif mode == 'convolve':
        operator = SPECIAL_OPERATORS.add('convolve', partial(convolve, **kwargs))
    else:
        raise ValueError(f'Mode {mode} is not supported.')

    return CompoundModel(operator, model, kernel)


def convolve_models_fft(model, kernel, bounding_box, resolution, cache=True, **kwargs):
    """
    Convolve two models using `~astropy.convolution.convolve_fft`.

    Parameters
    ----------
    model : `~astropy.modeling.core.Model`
        Functional model
    kernel : `~astropy.modeling.core.Model`
        Convolution kernel
    bounding_box : tuple
        The bounding box which encompasses enough of the support of both
        the ``model`` and ``kernel`` so that an accurate convolution can be
        computed.
    resolution : float
        The resolution that one wishes to approximate the convolution
        integral at.
    cache : optional, bool
        Default value True. Allow for the storage of the convolution
        computation for later reuse.
    **kwargs : dict
        Keyword arguments to be passed either to `~astropy.convolution.convolve`
        or `~astropy.convolution.convolve_fft` depending on ``mode``.

    Returns
    -------
    default : `~astropy.modeling.core.CompoundModel`
        Convolved model
    """

    operator = SPECIAL_OPERATORS.add('convolve_fft', partial(convolve_fft, **kwargs))

    return Convolution(operator, model, kernel, bounding_box, resolution, cache)
