# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module implements functions and classes for spatial statistics.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import math

class RipleysKEstimate(object):
    """
    This class implements an estimator for Ripley's K function in a two
    dimensional space.

    Parameters
    ----------
    data : 2D array
        Set of observed points in as a n by 2 array which will be used to
        estimate Ripley's K function.
    area : float
        Area of study from which the points where observed.
    x_max, y_max : float, float, optional
        Maximum rectangular dimensions of the area of study.
        Required if ``mode == 'translation'``.
    lratio : float, optional
        Greatest ratio between lengths of the rectangular window of study,
        i.e., ``lratio = max(width/height, height/width)``.
        Required if ``mode == 'ohser'``.

    Examples
    --------
    >>> import numpy as np
    >>> from matplotlib import pyplot as plt
    >>> from astropy.stats import RipleysKEstimate
    >>> z = np.random.uniform(low=5, high=10, size=(100, 2))
    >>> area = 25
    >>> Kest = RipleysKEstimate(data=z, area=area, max_height=10, max_width=10)
    >>> r = np.linspace(0, 2.5, 100)
    >>> plt.plot(r, Kest.poisson(r)) # doctest: +SKIP
    >>> plt.plot(r, Kest(r, mode='none')) # doctest: +SKIP
    >>> plt.plot(r, Kest(r, mode='translation')) # doctest: +SKIP

    References
    ----------
    .. [1] Spatial descriptive statistics.
       <https://en.wikipedia.org/wiki/Spatial_descriptive_statistics>
    .. [2] Package spatstat.
       <https://cran.r-project.org/web/packages/spatstat/spatstat.pdf>
    .. [3] Cressie, N.A.C. (1991). Statistics for Spatial Data,
       Wiley, New York.
    .. [4] Stoyan, D., Stoyan, H. (1992). Fractals, Random Shapes and
       Point Fields, Akademie Verlag GmbH, Chichester.
    """

    def __init__(self, data, area, x_max=None, y_max=None, lratio=None):
        self.data = data
        self.area = area
        self.x_max = x_max
        self.y_max = y_max
        self.lratio = lratio

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        value = np.asarray(value)

        if value.shape[1] == 2:
            self._data = value
        else:
            raise ValueError('data must be an n by 2 array, where n is the '
                             'number of observed points.')

    @property
    def area(self):
        return self._area

    @area.setter
    def area(self, value):
        if isinstance(value, (float, int)) and value > 0:
            self._area = value
        else:
            raise ValueError('area is expected to be a positive number. '
                             'Got {}.'.format(value))

    @property
    def y_max(self):
        return self._y_max

    @y_max.setter
    def y_max(self, value):
        if value is None or (isinstance(value, (float, int)) and value > 0):
            self._y_max = value
        else:
            raise ValueError('y_max is expected to be a positive number '
                             'or None. Got {}.'.format(value))

    @property
    def x_max(self):
        return self._x_max

    @x_max.setter
    def x_max(self, value):
        if value is None or (isinstance(value, (float, int)) and value > 0):
            self._x_max = value
        else:
            raise ValueError('x_max is expected to be a positive number '
                             'or None. Got {}.'.format(value))

    @property
    def lratio(self):
        return self._lratio

    @lratio.setter
    def lratio(self, value):
        if value is None or (isinstance(value, (float, int)) and value >= 1):
            self._lratio = value
        else:
            raise ValueError('lratio is expected to be a positive number'
                             ' (>= 1) or None. Got {}.'.format(value))

    def __call__(self, radii, mode='none'):
        return self.evaluate(radii=radii, mode=mode)

    def poisson(self, radii):
        """
        Evaluates the Ripley K function for the homongeneous Poisson process,
        also known as Complete State of Randomness (CSR).

        Parameters
        ----------
        radii : 1D array
            Set of distances in which Ripley's K function will be evaluated.

        Returns
        -------
        output : 1D array
            Ripley's K function evaluated at ``radii``.
        """

        return np.pi * radii * radii

    def evaluate(self, radii, mode='none'):
        """
        Evaluates the Ripley K estimator for a given set of values ``radii``.

        Parameters
        ----------
        radii : 1D array
            Set of distances in which Ripley's K estimator will be evaluated.
            Usually, it's common to consider max(radii) < (area/2)**0.5.
        mode : str
            Keyword which indicates the method for edge effects correction.
            Available methods are {'none', 'translation', 'ohser'}.

            * 'none' : this method does not take into account any edge effects
                whatsoever.
            * 'translation' : computes the intersection of rectangular areas
                centered at the given points provided the upper bounds of the
                dimensions of the rectangular area of study.
        Returns
        -------
        ripley : 1D array
            Ripley's K function estimator evaluated at ``radii``.
        """

        self.data = np.asarray(self.data)
        npts = len(self.data)
        ripley = np.zeros(len(radii))
        distances = np.zeros((npts * (npts - 1)) // 2, dtype=np.double)
        diff = np.zeros(shape=(npts * (npts - 1) // 2, 2), dtype=np.double)

        k = 0
        for i in range(npts - 1):
            for j in range(i + 1, npts):
                diff[k] = abs(self.data[i] - self.data[j])
                distances[k] = math.sqrt((diff[k] * diff[k]).sum())
                k = k + 1

        if mode == 'none':
            for idx in range(len(radii)):
                ripley[idx] = (distances < radii[idx]).sum()

            ripley = self.area * 2. * ripley / (npts * (npts - 1))
        # eq. 15.11 Stoyan book page 283
        elif mode == 'translation':
            x_diff, y_diff = diff[:][:, 0], diff[:][:, 1]
            intersec_area = ((self.x_max - x_diff) * (self.y_max - y_diff))

            for idx in range(len(radii)):
                dist_indicator = distances < radii[idx]
                ripley[idx] = ((1 / intersec_area) * dist_indicator).sum()

            ripley = (self.area**2 / (npts * (npts - 1))) * 2 * ripley
        elif mode == 'ohser':
            # Stoyan book page 123
            # Stoyan book eq 15.13
            a, b = self.area, self.lratio
            x = distances / math.sqrt(a / b)
            u = np.sqrt((x * x - 1) * (x > 1))
            v = np.sqrt((x * x - b ** 2) * (x < math.sqrt(b**2 + 1)) * (x > b))
            c1 = np.pi - 2 * x * (1 + 1 / b) + x * x / b
            c2 = 2 * np.arcsin((1 / x) * (x > 1)) - 1 / b - 2 * (x - u)
            c3 = (2 * np.arcsin(((b - u * v) / (x * x))
                                * (x > b) * (x < math.sqrt(b**2 + 1)))
                  + 2 * u + 2 * v / b - b - (1 + x * x) / b)

            cov_func = ((a / np.pi) * (c1 * (x >= 0) * (x <= 1)
                                       + c2 * (x > 1) * (x <= b)
                                       + c3 * (b < x) * (x < math.sqrt(b ** 2 + 1))))

            for idx in range(len(radii)):
                dist_indicator = distances < radii[idx]
                ripley[idx] = ((1 / cov_func) * dist_indicator).sum()

            ripley = (self.area**2 / (npts * (npts - 1))) * 2 * ripley
        elif mode == 'ripley':
            # Stoyan book eq 15.14
            pass

        else:
            raise ValueError('mode {} is not implemented'.format(mode))

        return ripley
