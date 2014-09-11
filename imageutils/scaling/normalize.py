"""
Normalization class for Matplotlib that can be used to produce colorbars.
"""

import numpy as np
from numpy import ma

from matplotlib.colors import Normalize

__all__ = ['ImageNormalize']


class ImageNormalize(Normalize):

    def __init__(self, vmin=None, vmax=None, stretch=None, clip=False):

        super(ImageNormalize, self).__init__(vmin=vmin, vmax=vmax, clip=clip)

        self.vmin = vmin
        self.vmax = vmax
        self.stretch = stretch
        self.inverse_stretch = stretch.inverted()

    def __call__(self, values, clip=None):

        if clip is None:
            clip = self.clip

        # Convert to masked array and make sure scalars get broadcast to 1-d
        values = ma.atleast_1d(values)

        # Normalize to 0 -> 1 range
        values_norm = (values - self.vmin) / (self.vmax - self.vmin)

        # Clip if required
        if clip:
            values_norm = np.clip(values_norm, 0., 1.)

        # Stretch values
        new_values = self.stretch(values_norm)
        
        # Don't assume stretch returned a masked array
        return ma.asarray(new_values)

    def inverse(self, values):

        # Find unstretched values in range 0 to 1
        values_norm = self.inverse_stretch(values)

        # Scale to original range
        return values_norm * (self.vmax - self.vmin) + self.vmin
