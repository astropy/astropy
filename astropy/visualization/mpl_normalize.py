"""
Normalization class for Matplotlib that can be used to produce
colorbars.
"""

from __future__ import division, print_function
import numpy as np
from numpy import ma
from .stretch import LinearStretch

try:
    import matplotlib    # pylint: disable=W0611
    from matplotlib.colors import Normalize

    # On older versions of matplotlib Normalize is an old-style class
    if not isinstance(Normalize, type):
        class Normalize(Normalize, object):
            pass
except ImportError:
    class Normalize(object):
        def __init__(self, *args, **kwargs):
            raise ImportError('matplotlib is required in order to use this '
                              'class')


__all__ = ['ImageNormalize']


class ImageNormalize(Normalize):
    """
    Normalization class to be used with Matplotlib.

    Parameters
    ----------
    data : `~numpy.ndarray`, optional
        The image array.  This input is used only if ``interval`` is
        also input.  ``data`` and ``interval`` are used to compute the
        vmin and/or vmax values only if ``vmin`` or ``vmax`` are not input.
    interval : `~astropy.visualization.BaseInterval` subclass instance, optional
        The interval object to apply to the input ``data`` to determine
        the ``vmin`` and ``vmax`` values.  This input is used only if
        ``data`` is also input.  ``data`` and ``interval`` are used to
        compute the vmin and/or vmax values only if ``vmin`` or ``vmax``
        are not input.
    vmin, vmax : float
        The minimum and maximum levels to show for the data.  The
        ``vmin`` and ``vmax`` inputs override any calculated values from
        the ``interval`` and ``data`` inputs.
    stretch : `~astropy.visualization.BaseStretch` subclass instance, optional
        The stretch object to apply to the data.  The default is
        `~astropy.visualization.LinearStretch`.
    clip : bool, optional
        If `True` (default), data values outside the [0:1] range are
        clipped to the [0:1] range.
    """

    def __init__(self, data=None, interval=None, vmin=None, vmax=None,
                 stretch=LinearStretch(), clip=False):
        # this super call checks for matplotlib
        super(ImageNormalize, self).__init__(vmin=vmin, vmax=vmax, clip=clip)

        self.vmin = vmin
        self.vmax = vmax
        if data is not None and interval is not None:
            _vmin, _vmax = interval.get_limits(data)
            if self.vmin is None:
                self.vmin = _vmin
            if self.vmax is None:
                self.vmax = _vmax

        self.stretch = stretch
        self.interval = interval
        self.inverse_stretch = stretch.inverse
        self.clip = clip

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

        # Set default values for vmin and vmax if not specified
        self.autoscale_None(values)

        # Normalize based on vmin and vmax
        np.subtract(values, self.vmin, out=values)
        np.true_divide(values, self.vmax - self.vmin, out=values)

        # Clip to the 0 to 1 range
        if self.clip:
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
