from .interval import *
from .stretch import *


def scale_image(image, scale='linear', power=1.0, noise_level=None,
                min_cut=None, max_cut=None, min_percent=None,
                max_percent=None, percent=None):
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

    power : float, optional
        The power index for the image scaling.  The default is 1.0.

    min_cut : float, optional
        The pixel value of the minimum cut level.  Data values less than
        ``min_cut`` will set to ``min_cut`` before scaling the image.
        The default is the image minimum.  ``min_cut`` overrides
        ``min_percent``.

    max_cut : float, optional
        The pixel value of the maximum cut level.  Data values greater
        than ``min_cut`` will set to ``min_cut`` before scaling the
        image.  The default is the image maximum.  ``max_cut`` overrides
        ``max_percent``.

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

    Returns
    -------
    image : ndarray
        The 2D array of the scaled/stretched image with a minimum of 0.0
        and a maximum of 1.0.
    """

    if percent is not None:
        interval = PercentileInterval(percent)
    elif min_percent is not None or max_percent is not None:
        interval = AsymmetricPercentileInterval(min_percent or 0.,
                                                max_percent or 100.)
    elif min_cut is not None or max_cut is not None:
        interval = ManualInterval(min_cut or np.min(image),
                                  max_cut or np.max(image))
    else:
        interval = MinMaxInterval()

    if scale == 'linear':
        stretch = LinearStretch()
    elif scale == 'sqrt':
        stretch = SqrtStretch()
    elif scale == 'power':
        stretch = PowerStretch(power)
    elif scale == 'log':
        stretch = LogStretch()
    else:
        raise ValueError('Unknown scale: {0}'.format(scale))

    return (stretch + interval)(image)
