"""
Normalization class for Matplotlib that can be used to produce colorbars.
"""

from __future__ import division, print_function

from distutils.version import LooseVersion

import numpy as np
from numpy import ma

try:
    import matplotlib
    MATPLOTLIB_INSTALLED = True
    MATPLOTLIB_LT_12 = LooseVersion(matplotlib.__version__) < LooseVersion("1.2.0")
    from matplotlib.colors import Normalize
except ImportError:
    MATPLOTLIB_INSTALLED = False
    class Normalize(object):
        def __init__(self, *args, **kwargs):
            raise ImportError("matplotlib is required in order to use this class")

__all__ = ['ImageNormalize']


class ImageNormalize(Normalize):
    """
    Normalization class to be used with Matplotlib.

    Parameters
    ----------
    vmin, vmax : float
        The minimum and maximum levels to show for the data
    stretch : :class:`~astropy.visualization.BaseStretch` instance
        The stretch to use for the normalization
    clip : bool, optional
        Whether to clip the output values to the [0:1] range
    """

    def __init__(self, vmin=None, vmax=None, stretch=None, clip=False):

        if MATPLOTLIB_INSTALLED and MATPLOTLIB_LT_12:
            # Normalize is an old-style class
            Normalize.__init__(self, vmin=vmin, vmax=vmax, clip=clip)
        else:  # Normalize is a new-style class
            super(ImageNormalize, self).__init__(vmin=vmin, vmax=vmax, clip=clip)

        self.vmin = vmin
        self.vmax = vmax
        self.stretch = stretch
        self.inverse_stretch = stretch.inverse

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
