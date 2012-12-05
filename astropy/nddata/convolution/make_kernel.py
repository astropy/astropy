# Licensed under a 3-clause BSD style license - see PYFITS.rst

import numpy as np

def make_kernel(kernelshape, kernelwidth=3, kerneltype='gaussian',
        trapslope=None, normalize_kernel=np.sum, force_odd=False):
    """
    Create a smoothing kernel for use with `convolve` or `convolve_fft`.

    Parameters
    ----------
    kernelshape : n-tuple
        A tuple (or list or array) defining the shape of the kernel.  The
        length of kernelshape determines the dimensionality of the resulting
        kernel
    kernelwidth : float
        Width of kernel in pixels  (see definitions under `kerneltype`)
    kerneltype : {'gaussian', 'boxcar', 'tophat', 'brickwall', 'airy', 'trapezoid'}
        Defines the type of kernel to be generated. The following types are
        available:

        * 'gaussian'
            Uses a gaussian kernel with sigma = `kernelwidth` (in pixels),
            i.e. kernel = exp(-r**2 / (2*sigma**2)) where r is the radius.
        * 'boxcar'
            A `kernelwidth` x `kernelwidth` square kernel, i.e.,
            kernel = (x < `kernelwidth`) * (y < `kernelwidth`)
        * 'tophat'
            A flat circle  with radius = `kernelwidth`,
            i.e., kernel = (r < `kernelwidth`)
        * 'brickwall' or 'airy'
            A kernel using the airy function from optics. It requires
            `scipy.special` for the bessel function. See e.g.,
            http://en.wikipedia.org/wiki/Airy_disk.
        * 'trapezoid'
            A kernel like 'tophat' but with sloped edges. It is
            effectively a cone chopped off at the `kernelwidth` radius.

    trapslope : float
        Slope of the trapezoid kernel.  Only used if `kerneltype` == 'trapezoid'
    normalize_kernel : function
        Function to use for kernel normalization
    force_odd : boolean
        If set, forces the kernel to have odd dimensions (needed for convolve
        w/o ffts)

    Returns
    -------
    kernel : ndarray
        An N-dimensional float array

    Examples
    --------

    >>> make_kernel([3,3],1,'boxcar')
    array([[ 0.  0.  0.]
           [ 0.  1.  0.]
           [ 0.  0.  0.]])

    >>> make_kernel([9],1) # Gaussian by default
    array([  1.33830625e-04   4.43186162e-03   5.39911274e-02   2.41971446e-01
             3.98943469e-01   2.41971446e-01   5.39911274e-02   4.43186162e-03
             1.33830625e-04])

    >>> make_kernel([3,3],3,'boxcar')
    array([[ 0.11111111,  0.11111111,  0.11111111],
           [ 0.11111111,  0.11111111,  0.11111111],
           [ 0.11111111,  0.11111111,  0.11111111]])

    >>> make_kernel([3,3],1.4,'tophat')
    array([[ 0. ,  0.2,  0. ],
           [ 0.2,  0.2,  0.2],
           [ 0. ,  0.2,  0. ]])


    """

    if force_odd:
        kernelshape = [n-1 if (n%2==0) else n for n in kernelshape]

    if normalize_kernel is True:
        normalize_kernel = np.sum

    if kerneltype == 'gaussian':
        rr = np.sum([(x-(x.max()+1)//2)**2 for x in np.indices(kernelshape)],axis=0)**0.5
        kernel = np.exp(-(rr**2)/(2.*kernelwidth**2))
        kernel /= normalize_kernel(kernel) #/ (kernelwidth**2 * (2*np.pi))

    elif kerneltype == 'boxcar':
        kernel = np.zeros(kernelshape,dtype='float64')
        kernelslices = []
        for dimsize in kernelshape:
            center = dimsize - (dimsize+1)//2
            kernelslices += [slice(center - (kernelwidth)//2, center + (kernelwidth+1)//2)]
        kernel[kernelslices] = 1.0
        kernel /= normalize_kernel(kernel)
    elif kerneltype == 'tophat':
        rr = np.sum([(x-(x.max())/2.)**2 for x in np.indices(kernelshape)],axis=0)**0.5
        kernel = np.zeros(kernelshape,dtype='float64')
        kernel[rr<kernelwidth] = 1.0
        # normalize
        kernel /= normalize_kernel(kernel)
    elif kerneltype in ('brickwall','airy'):
        try:
            import scipy.special
        except ImportError:
            raise ImportError("Could not import scipy.special; cannot create an "+
                    "airy kernel without this (need the bessel function)")
        rr = np.sum([(x-(x.max())/2.)**2 for x in np.indices(kernelshape)],axis=0)**0.5
        # airy function is first bessel(x) / x  [like the sinc]
        kernel = j1(rr/kernelwidth) / (rr/kernelwidth)
        # fix NAN @ center
        kernel[rr==0] = 0.5
        kernel /= normalize_kernel(kernel)
    elif kerneltype == 'trapezoid':
        rr = np.sum([(x-(x.max())/2.)**2 for x in np.indices(kernelshape)],axis=0)**0.5
        if trapslope:
            zz = rr.max()-(rr*trapslope)
            zz[zz<0] = 0
            zz[rr<kernelwidth] = 1.0
            kernel = zz/zz.sum()
        else:
            raise ValueError("Must specify a slope for kerneltype='trapezoid'")

    return kernel
