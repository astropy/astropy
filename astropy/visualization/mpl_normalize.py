"""
Normalization class for Matplotlib that can be used to produce
colorbars.
"""

import inspect

import numpy as np
from numpy import ma

from .interval import (PercentileInterval, AsymmetricPercentileInterval,
                       ManualInterval, MinMaxInterval, BaseInterval)
from .stretch import (LinearStretch, SqrtStretch, PowerStretch, LogStretch,
                      AsinhStretch, BaseStretch)

try:
    import matplotlib  # pylint: disable=W0611
    from matplotlib.colors import Normalize
    from matplotlib import pyplot as plt
except ImportError:
    class Normalize:
        def __init__(self, *args, **kwargs):
            raise ImportError('matplotlib is required in order to use this '
                              'class.')


__all__ = ['ImageNormalize', 'simple_norm', 'imshow_norm']

__doctest_requires__ = {'*': ['matplotlib']}


class ImageNormalize(Normalize):
    """
    Normalization class to be used with Matplotlib.

    Parameters
    ----------
    data : `~numpy.ndarray`, optional
        The image array.  This input is used only if ``interval`` is
        also input.  ``data`` and ``interval`` are used to compute the
        vmin and/or vmax values only if ``vmin`` or ``vmax`` are not
        input.
    interval : `~astropy.visualization.BaseInterval` subclass instance, optional
        The interval object to apply to the input ``data`` to determine
        the ``vmin`` and ``vmax`` values.  This input is used only if
        ``data`` is also input.  ``data`` and ``interval`` are used to
        compute the vmin and/or vmax values only if ``vmin`` or ``vmax``
        are not input.
    vmin, vmax : float, optional
        The minimum and maximum levels to show for the data.  The
        ``vmin`` and ``vmax`` inputs override any calculated values from
        the ``interval`` and ``data`` inputs.
    stretch : `~astropy.visualization.BaseStretch` subclass instance
        The stretch object to apply to the data.  The default is
        `~astropy.visualization.LinearStretch`.
    clip : bool, optional
        If `True`, data values outside the [0:1] range are clipped to
        the [0:1] range.
    """

    def __init__(self, data=None, interval=None, vmin=None, vmax=None,
                 stretch=LinearStretch(), clip=False):
        # this super call checks for matplotlib
        super().__init__(vmin=vmin, vmax=vmax, clip=clip)

        self.vmin = vmin
        self.vmax = vmax

        if stretch is None:
            raise ValueError('stretch must be input')
        if not isinstance(stretch, BaseStretch):
            raise TypeError('stretch must be an instance of a BaseStretch '
                            'subclass')
        self.stretch = stretch

        if interval is not None and not isinstance(interval, BaseInterval):
            raise TypeError('interval must be an instance of a BaseInterval '
                            'subclass')
        self.interval = interval

        self.inverse_stretch = stretch.inverse
        self.clip = clip

        # Define vmin and vmax if not None and data was input
        if data is not None:
            self._set_limits(data)

    def _set_limits(self, data):
        if self.vmin is not None and self.vmax is not None:
            return

        # Define vmin and vmax from the interval class if not None
        if self.interval is None:
            if self.vmin is None:
                self.vmin = np.min(data[np.isfinite(data)])
            if self.vmax is None:
                self.vmax = np.max(data[np.isfinite(data)])
        else:
            _vmin, _vmax = self.interval.get_limits(data)
            if self.vmin is None:
                self.vmin = _vmin
            if self.vmax is None:
                self.vmax = _vmax

    def __call__(self, values, clip=None):
        if clip is None:
            clip = self.clip

        if isinstance(values, ma.MaskedArray):
            if clip:
                mask = False
            else:
                mask = values.mask
            values = values.filled(self.vmax)
        else:
            mask = False

        # Make sure scalars get broadcast to 1-d
        if np.isscalar(values):
            values = np.array([values], dtype=float)
        else:
            # copy because of in-place operations after
            values = np.array(values, copy=True, dtype=float)

        # Define vmin and vmax if not None
        self._set_limits(values)

        # Normalize based on vmin and vmax
        np.subtract(values, self.vmin, out=values)
        np.true_divide(values, self.vmax - self.vmin, out=values)

        # Clip to the 0 to 1 range
        if clip:
            values = np.clip(values, 0., 1., out=values)

        # Stretch values
        values = self.stretch(values, out=values, clip=False)

        # Convert to masked array for matplotlib
        return ma.array(values, mask=mask)

    def inverse(self, values):
        # Find unstretched values in range 0 to 1
        values_norm = self.inverse_stretch(values, clip=False)

        # Scale to original range
        return values_norm * (self.vmax - self.vmin) + self.vmin


def simple_norm(data, stretch='linear', power=1.0, asinh_a=0.1, min_cut=None,
                max_cut=None, min_percent=None, max_percent=None,
                percent=None, clip=False, log_a=1000):
    """
    Return a Normalization class that can be used for displaying images
    with Matplotlib.

    This function enables only a subset of image stretching functions
    available in `~astropy.visualization.mpl_normalize.ImageNormalize`.

    This function is used by the
    ``astropy.visualization.scripts.fits2bitmap`` script.

    Parameters
    ----------
    data : `~numpy.ndarray`
        The image array.

    stretch : {'linear', 'sqrt', 'power', log', 'asinh'}, optional
        The stretch function to apply to the image.  The default is
        'linear'.

    power : float, optional
        The power index for ``stretch='power'``.  The default is 1.0.

    asinh_a : float, optional
        For ``stretch='asinh'``, the value where the asinh curve
        transitions from linear to logarithmic behavior, expressed as a
        fraction of the normalized image.  Must be in the range between
        0 and 1.  The default is 0.1.

    min_cut : float, optional
        The pixel value of the minimum cut level.  Data values less than
        ``min_cut`` will set to ``min_cut`` before stretching the image.
        The default is the image minimum.  ``min_cut`` overrides
        ``min_percent``.

    max_cut : float, optional
        The pixel value of the maximum cut level.  Data values greater
        than ``min_cut`` will set to ``min_cut`` before stretching the
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

    clip : bool, optional
        If `True`, data values outside the [0:1] range are clipped to
        the [0:1] range.

    log_a : float, optional
        The log index for ``stretch='log'``. The default is 1000.

    Returns
    -------
    result : `ImageNormalize` instance
        An `ImageNormalize` instance that can be used for displaying
        images with Matplotlib.
    """

    if percent is not None:
        interval = PercentileInterval(percent)
    elif min_percent is not None or max_percent is not None:
        interval = AsymmetricPercentileInterval(min_percent or 0.,
                                                max_percent or 100.)
    elif min_cut is not None or max_cut is not None:
        interval = ManualInterval(min_cut, max_cut)
    else:
        interval = MinMaxInterval()

    if stretch == 'linear':
        stretch = LinearStretch()
    elif stretch == 'sqrt':
        stretch = SqrtStretch()
    elif stretch == 'power':
        stretch = PowerStretch(power)
    elif stretch == 'log':
        stretch = LogStretch(log_a)
    elif stretch == 'asinh':
        stretch = AsinhStretch(asinh_a)
    else:
        raise ValueError(f'Unknown stretch: {stretch}.')

    vmin, vmax = interval.get_limits(data)

    return ImageNormalize(vmin=vmin, vmax=vmax, stretch=stretch, clip=clip)


# used in imshow_norm
_norm_sig = inspect.signature(ImageNormalize)


def imshow_norm(data, ax=None, imshow_only_kwargs={}, **kwargs):
    """ A convenience function to call matplotlib's `matplotlib.pyplot.imshow`
    function, using an `ImageNormalize` object as the normalization.

    Parameters
    ----------
    data : 2D or 3D array_like - see `~matplotlib.pyplot.imshow`
        The data to show. Can be whatever `~matplotlib.pyplot.imshow` and
        `ImageNormalize` both accept.
    ax : None or `~matplotlib.axes.Axes`, optional
        If None, use pyplot's imshow.  Otherwise, calls ``imshow`` method of the
        supplied axes.
    imshow_only_kwargs : dict, optional
        Arguments to be passed directly to `~matplotlib.pyplot.imshow` without
        first trying `ImageNormalize`.  This is only for keywords that have the
        same name in both `ImageNormalize` and `~matplotlib.pyplot.imshow` - if
        you want to set the `~matplotlib.pyplot.imshow` keywords only, supply
        them in this dictionary.
    kwargs : dict, optional
        All other keyword arguments are parsed first by the
        `ImageNormalize` initializer, then to
        `~matplotlib.pyplot.imshow`.

    Returns
    -------
    result : tuple
        A tuple containing the `~matplotlib.image.AxesImage` generated
        by `~matplotlib.pyplot.imshow` as well as the `ImageNormalize`
        instance.

    Notes
    -----
    The ``norm`` matplotlib keyword is not supported.

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        from astropy.visualization import (imshow_norm, MinMaxInterval,
                                           SqrtStretch)

        # Generate and display a test image
        image = np.arange(65536).reshape((256, 256))
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        im, norm = imshow_norm(image, ax, origin='lower',
                               interval=MinMaxInterval(),
                               stretch=SqrtStretch())
        fig.colorbar(im)
    """

    if 'X' in kwargs:
        raise ValueError('Cannot give both ``X`` and ``data``')

    if 'norm' in kwargs:
        raise ValueError('There is no point in using imshow_norm if you give '
                         'the ``norm`` keyword - use imshow directly if you '
                         'want that.')

    imshow_kwargs = dict(kwargs)

    norm_kwargs = {'data': data}
    for pname in _norm_sig.parameters:
        if pname in kwargs:
            norm_kwargs[pname] = imshow_kwargs.pop(pname)

    for k, v in imshow_only_kwargs.items():
        if k not in _norm_sig.parameters:
            # the below is not strictly "has to be true", but is here so that
            # users don't start using both imshow_only_kwargs *and* keyword
            # arguments to this function, as that makes for more confusing
            # user code
            raise ValueError('You provided a keyword to imshow_only_kwargs '
                             '({}) that is not a keyword for ImageNormalize. '
                             'This is not supported. Instead you should '
                             'pass the keyword directly into imshow_norm'
                             .format(k))
        imshow_kwargs[k] = v

    imshow_kwargs['norm'] = ImageNormalize(**norm_kwargs)

    if ax is None:
        imshow_result = plt.imshow(data, **imshow_kwargs)
    else:
        imshow_result = ax.imshow(data, **imshow_kwargs)

    return imshow_result, imshow_kwargs['norm']
