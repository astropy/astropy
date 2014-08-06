# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from astropy.stats import sigma_clip

__all__ = ['find_cutlevels', 'normalize_image', 'scale_image',
           'sigmaclip_stats']


_MINMAX_PARAMS = \
    """
    min_cut : float, optional
        The pixel value of the minimum cut level.  Data values less than
        ``min_cut`` will set to ``min_cut`` before scaling the image.

    max_cut : float, optional
        The pixel value of the maximum cut level.  Data values greater
        than ``min_cut`` will set to ``min_cut`` before scaling the
        image.
    """.strip()


_PERCENTILE_PARAMS = \
    """
    min_percent : float, optional
        The percentile value used to determine the pixel value of
        minimum cut level.  The default is 0.0.  ``min_percent``
        overrides ``percent``.

    max_percent : float, optional
        The percentile value used to determine the pixel value of
        maximum cut level.  The default is 100.0.  ``max_percent``
        overrides ``percent``.

    percent : float, optional
        The percentage of the image values used to determine the pixel
        values of the minimum and maximum cut levels.  The lower cut
        level will set at the ``(100 - percent) / 2`` percentile, while
        the upper cut level will be set at the ``(100 + percent) / 2``
        percentile.  The default is 100.0.  ``percent`` is ignored if
        either ``min_percent`` or ``max_percent`` is input.
    """.strip()


def _insert_cutlevel_params(func):
    """Insert the cutlevel parameters into the function documentation."""
    func.__doc__ = func.__doc__.format(minmax_params=_MINMAX_PARAMS,
                                       percentile_params=_PERCENTILE_PARAMS)
    return func


def _insert_percentile_params(func):
    """Insert the cutlevel parameters into the function documentation."""
    func.__doc__ = func.__doc__.format(percentile_params=_PERCENTILE_PARAMS)
    return func


@_insert_percentile_params
def find_cutlevels(image, min_percent=None, max_percent=None, percent=None):
    """
    Find pixel values of the minimum and maximum image cut levels from
    percentiles of the image values.

    Parameters
    ----------
    image : array_like
        The 2D array of the image.

    {percentile_params}

    Returns
    -------
    cutlevels : 2-tuple of floats
        The pixel values of the minimum and maximum cut levels.
    """

    if percent is not None:
        assert (percent >= 0) and (percent <= 100.0), \
            'percent must be >= 0 and <= 100.0'
        if min_percent is None and max_percent is None:
            min_percent = (100.0 - float(percent)) / 2.0
            max_percent = 100.0 - min_percent
    if min_percent is None:
        min_percent = 0.0
    assert min_percent >= 0, 'min_percent must be >= 0'
    min_cut = np.percentile(image, min_percent)
    if max_percent is None:
        max_percent = 100.0
    assert max_percent <= 100.0, 'max_percent must be <= 100.0'
    max_cut = np.percentile(image, max_percent)
    assert min_percent <= max_percent, 'min_percent must be <= max_percent'
    return min_cut, max_cut


@_insert_cutlevel_params
def normalize_image(image, min_cut=None, max_cut=None, **kwargs):
    """
    Rescale image values between minimum and maximum cut levels to
    values between 0 and 1, inclusive.

    Parameters
    ----------
    image : array_like
        The 2D array of the image.

    {minmax_params}

    {percentile_params}

    Returns
    -------
    image : ndarray
        The normalized image with a minimum of 0.0 and a maximum of 1.0.

    cutlevels : 2-tuple of floats
        The pixel values of the minimum and maximum cut levels.
    """

    from skimage import exposure
    image = np.asarray(image)
    if min_cut is not None and max_cut is not None:
        cutlevels = (min_cut, max_cut)
    else:
        cutlevels = find_cutlevels(image, **kwargs)
    # now override percentiles if either min_cut or max_cut is set
    if min_cut is not None:
        cutlevels[0] = min_cut
    if max_cut is not None:
        cutlevels[1] = max_cut
    assert cutlevels[0] <= cutlevels[1], 'minimum cut must be <= maximum cut'
    image_norm = exposure.rescale_intensity(image, in_range=cutlevels,
                                            out_range=(0, 1))
    return image_norm, cutlevels


@_insert_cutlevel_params
def scale_image(image, scale='linear', power=1.0, noise_level=None,
                min_cut=None, max_cut=None, **kwargs):
    """
    Perform scaling/stretching of an image between minimum and maximum
    cut levels.

    Parameters
    ----------
    image : array_like
        The 2D array of the image.

    scale : {{'linear', 'sqrt', 'power', log', 'asinh'}}
        The scaling/stretch function to apply to the image.  The default
        is 'linear'.

        ``scaling='power'`` requires input of the ``power`` keyword.

        ``scaling='asinh``` can use the ``noise_level`` keyword.

    power : float, optional
        The power index for the image scaling.  The default is 1.0.

    noise_level: float, optional
        The noise level of the image.  Pixel values less than
        ``noise_level`` will approximately be linearly scaled, while
        pixel values greater than ``noise_level`` will approximately be
        logarithmically scaled.  If ``noise_level`` is not input, an
        estimate of the 2-sigma level above the background will be used.

    {minmax_params}

    {percentile_params}

    Returns
    -------
    image : ndarray
        The 2D array of the scaled/stretched image with a minimum of 0.0
        and a maximum of 1.0.
    """

    image_norm, cutlevels = normalize_image(image, min_cut=min_cut,
                                            max_cut=max_cut, **kwargs)
    if scale == 'linear':
        return image_norm
    elif scale == 'sqrt':
        return np.sqrt(image_norm)
    elif scale == 'power':
        return image_norm ** power
    elif scale == 'log':
        return np.log10(image_norm + 1.0) / np.log10(2.0)
    elif scale == 'asinh':
        if noise_level is None:
            mean, median, stddev = sigmaclip_stats(image_norm, sigma=3.0)
            noise_level = mean + (2.0 * stddev)   # 2 sigma above background
        min_cut, max_cut = cutlevels
        z = (noise_level - min_cut) / (max_cut - min_cut)
        if z == 0.0:    # avoid zero division
            z = 1.e-2
        return np.arcsinh(image_norm / z) / np.arcsinh(1.0 / z)
    else:
        raise ValueError('scale type is unknown')


def sigmaclip_stats(image, image_mask=None, mask_val=None, sigma=3.0,
                    iters=None):
    """
    Calculate sigma-clipped statistics of an image.  For example,
    sigma-clipped statistics can be used to estimate the background and
    background noise in an image.

    Parameters
    ----------
    image : array_like
        The 2D array of the image.

    image_mask : array_like, bool, optional
        A boolean mask with the same shape as ``image``, where a `True`
        value indicates the corresponding element of ``image`` is
        invalid.  Masked pixels are ignored when computing the image
        statistics.

    mask_val : float, optional
        An image data value (e.g., ``0.0``) that is ignored when
        computing the image statistics.  ``mask_val`` will be ignored if
        ``image_mask`` is input.

    sigma : float, optional
        The number of standard deviations to use as the clipping limit.

    iters : float, optional
       The number of clipping iterations to perform, or `None` to clip
       until convergence is achieved (i.e. continue until the last
       iteration clips nothing).

    Returns
    -------
    mean, median, stddev : floats
        The mean, median, and standard deviation of the sigma-clipped
        image.
    """

    if image_mask is not None:
        assert image_mask.dtype == np.bool, \
            'image_mask must be a boolean ndarray'
        image = image[~image_mask]
    if mask_val is not None and image_mask is None:
        idx = (image != mask_val).nonzero()
        image = image[idx]
    image_clip = sigma_clip(image, sig=sigma, iters=iters)
    goodvals = image_clip.data[~image_clip.mask]
    return np.mean(goodvals), np.median(goodvals), np.std(goodvals)
