"""
Normalization class for Matplotlib that can be used to produce colorbars.
"""

import numpy as np
from numpy import ma

from matplotlib.colors import Normalize


__all__ = ['ImageNormalize']


class ImageNormalize(Normalize):

    def __init__(self, vmin=None, vmax=None, stretch=None, clip=True):

        super(ImageNormalize, self).__init__(vmin=vmin, vmax=vmax, clip=clip)

        self.vmin = vmin
        self.vmax = vmax
        self.stretch = stretch
        self.inverse_stretch = stretch.inverted()

    def __call__(self, values, clip=None):

        if clip is None:
            clip = self.clip

        # Make sure scalars get broadcast to 1-d
        if np.isscalar(values):
            values = np.atleast_1d(values)
        else:
            values = values.copy()  # copy because of in-place operations after

        # Normalize based on vmin and vmax
        np.subtract(values, self.vmin, out=values)
        np.divide(values, self.vmax - self.vmin, out=values)

        # Clip to the 0 to 1 range
        if self.clip:
            values = np.clip(values, 0., 1., out=values)

        # Stretch values
        values = self.stretch(values, out=values)

        # Convert to masked array for matplotlib
        return ma.array(values)

    def inverse(self, values):

        # Find unstretched values in range 0 to 1
        values_norm = self.inverse_stretch(values)

        # Scale to original range
        return values_norm * (self.vmax - self.vmin) + self.vmin
