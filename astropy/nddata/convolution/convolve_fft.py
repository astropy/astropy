import numpy as np
import warnings

try:
    import fftw3
    has_fftw = True

    def fftwn(array, nthreads=1):
        array = array.astype('complex').copy()
        outarray = array.copy()
        fft_forward = fftw3.Plan(array, outarray, direction='forward',
                flags=['estimate'], nthreads=nthreads)
        fft_forward.execute()
        return outarray

    def ifftwn(array, nthreads=1):
        array = array.astype('complex').copy()
        outarray = array.copy()
        fft_backward = fftw3.Plan(array, outarray, direction='backward',
                flags=['estimate'], nthreads=nthreads)
        fft_backward.execute()
        return outarray / np.size(array)
except ImportError:
    fftn = np.fft.fftn
    ifftn = np.fft.ifftn
    has_fftw = False
# I performed some fft speed tests and found that scipy is slower than numpy
# http://code.google.com/p/agpy/source/browse/trunk/tests/test_ffts.py However,
# the speed varied on machines - YMMV.  If someone finds that scipy's fft is
# faster, we should add that as an option here... not sure how exactly


def convolve_fft(array, kernel, boundary='fill', fill_value=0,
        crop=True, return_fft=False, fft_pad=True,
        psf_pad=False, interpolate_nan=False, quiet=False,
        ignore_edge_zeros=False, min_wt=0.0, normalize_kernel=False,
        use_numpy_fft=not has_fftw, nthreads=1):
    """
    Convolve an ndarray with an nd-kernel.  Returns a convolved image with
    shape = array.shape.  Assumes kernel is centered.

    convolve_fft differs from `scipy.signal.fftconvolve` in a few ways:
    * can treat NaN's as zeros or interpolate over them 
    * defaults to using the faster FFTW algorithm if installed
    * (optionally) pads to the nearest 2^n size to improve FFT speed
    * only operates in mode='same' (i.e., the same shape array is returned) mode

    Parameters
    ----------
    array : `numpy.ndarray`
          Array to be convolved with *kernel*
    kernel : `numpy.ndarray`
          Will be normalized if *normalize_kernel* is set.  Assumed to be
          centered (i.e., shifts may result if your kernel is asymmetric)
    boundary : {'fill', 'wrap'}
        A flag indicating how to handle boundaries:
            * 'fill' : set values outside the array boundary to fill_value
                       (default)
            * 'wrap' : periodic boundary
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
    use_numpy_fft : bool
        Force the code to use the numpy FFTs instead of FFTW even if FFTW is
        installed

    Returns
    -------
    default : ndarray
        `array` convolved with `kernel`.
        If `return_fft` is set, returns fft(`array`) * fft(`kernel`).
        If crop is not set, returns the image, but with the fft-padded size
        instead of the input size

    Examples
    --------
    >>> convolve_fft([1,0,3],[1,1,1])
    array([ 1.,  4.,  3.])

    >>> convolve_fft([1,np.nan,3],[1,1,1])
    array([ 1.,  4.,  3.])

    >>> convolve_fft([1,0,3],[0,1,0])
    array([ 1.,  0.,  3.])

    >>> convolve_fft([1,2,3],[1])
    array([ 1.,  2.,  3.])

    >>> convolve_fft([1,np.nan,3],[0,1,0], interpolate_nan=True)
    array([ 1.,  0.,  3.])

    >>> convolve_fft([1,np.nan,3],[0,1,0], interpolate_nan=True, min_wt=1e-8)
    array([ 1.,  nan,  3.])

    >>> convolve_fft([1,np.nan,3],[1,1,1], interpolate_nan=True)
    array([ 1.,  4.,  3.])

    >>> convolve_fft([1,np.nan,3],[1,1,1], interpolate_nan=True, normalize_kernel=True)
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
        raise Exception('array and kernel have differing number of'
                        'dimensions')

    # store the dtype for conversion back later
    array_dtype = array.dtype
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

    # replace fftn if has_fftw so that nthreads can be passed
    global fftn, ifftn
    if has_fftw and not use_numpy_fft:
        def fftn(*args, **kwargs):
            return fftwn(*args, nthreads=nthreads, **kwargs)

        def ifftn(*args, **kwargs):
            return ifftwn(*args, nthreads=nthreads, **kwargs)
    elif use_numpy_fft:
        fftn = np.fft.fftn
        ifftn = np.fft.ifftn


    # NAN catching
    nanmaskarray = (array != array)
    array[nanmaskarray] = 0
    nanmaskkernel = (kernel != kernel)
    kernel[nanmaskkernel] = 0
    if ((nanmaskarray.sum() > 0 or nanmaskkernel.sum() > 0) and not interpolate_nan
            and not quiet):
        warnings.warn("NOT ignoring nan values even though they are present" +
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


    if boundary is None:
        WARNING = ("The convolve_fft version of boundary=None is equivalent" +
                " to the convolve boundary='fill'.  There is no FFT " +
                " equivalent to convolve's zero-if-kernel-leaves-boundary" )
        warnings.warn(WARNING)
        psf_pad = True
    elif boundary == 'fill':
        # create a boundary region at least as large as the kernel
        psf_pad = True
    elif boundary == 'wrap':
        psf_pad = False
        fft_pad = False
        fill_value = 0 # force zero; it should not be used
    elif boundary == 'extend':
        raise NotImplementedError("The 'extend' option is not implemented " +
                "for fft-based convolution")

    arrayshape = array.shape
    kernshape = kernel.shape
    ndim = len(array.shape)
    if ndim != len(kernshape):
        raise ValueError("Image and kernel must " +
            "have same number of dimensions")
    # find ideal size (power of 2) for fft.
    # Can add shapes because they are tuples
    if fft_pad:
        if psf_pad:
            # add the dimensions and then take the max (bigger)
            fsize = 2**np.ceil(np.log2(
                np.max(np.array(arrayshape) + np.array(kernshape))))
        else:
            # add the shape lists (max of a list of length 4) (smaller)
            # also makes the shapes square
            fsize = 2**np.ceil(np.log2(np.max(arrayshape+kernshape)))
        newshape = np.array([fsize for ii in range(ndim)])
    else:
        if psf_pad:
            # just add the biggest dimensions
            newshape = np.array(arrayshape)+np.array(kernshape)
        else:
            newshape = np.array([np.max([imsh, kernsh])
                for imsh, kernsh in zip(arrayshape, kernshape)])


    # separate each dimension by the padding size...  this is to determine the
    # appropriate slice size to get back to the input dimensions
    arrayslices = []
    kernslices = []
    for ii, (newdimsize, arraydimsize, kerndimsize) in enumerate(zip(newshape, arrayshape, kernshape)):
        center = newdimsize - (newdimsize+1)//2
        arrayslices += [slice(center - arraydimsize//2,
            center + (arraydimsize+1)//2)]
        kernslices += [slice(center - kerndimsize//2,
            center + (kerndimsize+1)//2)]

    bigarray = np.ones(newshape, dtype=np.complex128) * fill_value
    bigkernel = np.zeros(newshape, dtype=np.complex128)
    bigarray[arrayslices] = array
    bigkernel[kernslices] = kernel
    arrayfft = fftn(bigarray)
    # need to shift the kernel so that, e.g., [0,0,1,0] -> [1,0,0,0] = unity
    kernfft = fftn(np.fft.ifftshift(bigkernel))
    fftmult = arrayfft*kernfft
    if (interpolate_nan or ignore_edge_zeros) and kernel_is_normalized:
        if ignore_edge_zeros:
            bigimwt = np.zeros(newshape, dtype=np.complex128)
        else:
            bigimwt = np.ones(newshape, dtype=np.complex128)
        bigimwt[arrayslices] = 1.0-nanmaskarray*interpolate_nan
        wtfft = fftn(bigimwt)
        # I think this one HAS to be normalized (i.e., the weights can't be
        # computed with a non-normalized kernel)
        wtfftmult = wtfft*kernfft/kernel.sum()
        wtsm = ifftn(wtfftmult)
        # need to re-zero weights outside of the image (if it is padded, we
        # still don't weight those regions)
        bigimwt[arrayslices] = wtsm.real[arrayslices]
        # curiously, at the floating-point limit, can get slightly negative numbers
        # they break the min_wt=0 "flag" and must therefore be removed
        bigimwt[bigimwt<0] = 0
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
