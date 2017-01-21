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
    area : float
        Area of study from which the points where observed.
    x_max, y_max : float, float, optional
        Maximum rectangular coordinates of the area of study.
        Required if ``mode == 'translation'``.
    x_min, y_min : float, float, optional
        Minimum rectangular coordinates of the area of study.
        Required if ``mode == 'variable-width'``.
    lratio : float, optional
        Greatest ratio between lengths of the rectangular window of study,
        i.e., ``lratio = max(width/height, height/width)``.
        Required if ``mode == 'ohser'``.

    Examples
    --------
    >>> import numpy as np
    >>> from matplotlib import pyplot as plt # doctest: +SKIP
    >>> from astropy.stats import RipleysKEstimate
    >>> z = np.random.uniform(low=5, high=10, size=(100, 2))
    >>> Kest = RipleysKEstimate(data=z, area=25, x_max=10, y_max=10,
    ... x_min=5, y_min=5, lratio=1)
    >>> r = np.linspace(0, 2.5, 100)
    >>> plt.plot(r, Kest.poisson(r)) # doctest: +SKIP
    >>> plt.plot(r, Kest(r, mode='none')) # doctest: +SKIP
    >>> plt.plot(r, Kest(r, mode='translation')) # doctest: +SKIP
    >>> plt.plot(r, Kest(r, mode='ohser')) # doctest: +SKIP
    >>> plt.plot(r, Kest(r, mode='ohser')) # doctest: +SKIP

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

    def __init__(self, area, x_max=None, y_max=None, x_min=0, y_min=0,
                 lratio=None):
        self.area = area
        self.x_max = x_max
        self.y_max = y_max
        self.x_min = x_min
        self.y_min = y_min
        self.lratio = lratio



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
    def y_min(self):
        return self._y_min

    @y_min.setter
    def y_min(self, value):
        if isinstance(value, (float, int)) and value >= 0:
            self._y_min = value
        else:
            raise ValueError('y_min is expected to be a nonnegative number. '
                             'Got {}.'.format(value))

    @property
    def x_min(self):
        return self._x_min

    @x_min.setter
    def x_min(self, value):
        if isinstance(value, (float, int)) and value >= 0:
            self._x_min = value
        else:
            raise ValueError('x_min is expected to be a nonnegative number. '
                             'Got {}.'.format(value))

    @property
    def lratio(self):
        return self._lratio

    @lratio.setter
    def lratio(self, value):
        if value is None or (isinstance(value, (float, int)) and value >= 1):
            self._lratio = value
        else:
            raise ValueError('lratio is expected to be a real number'
                             ' (>= 1) or None. Got {}.'.format(value))

    def __call__(self, data, radii, mode='none'):
        return self.evaluate(data=data, radii=radii, mode=mode)

    def _pairwise_diffs(self, data):
        npts = len(data)
        diff = np.zeros(shape=(npts * (npts - 1) // 2, 2), dtype=np.double)
        k = 0
        for i in range(npts - 1):
            for j in range(i + 1, npts):
                diff[k] = abs(data[i] - data[j])
                k = k + 1

        return diff

    def _pairwise_distances(self, data):
        npts = len(data)
        distances = np.zeros((npts * (npts - 1)) // 2, dtype=np.double)
        diff = self._pairwise_diffs(data)
        for k in range(len(distances)):
            distances[k] = math.sqrt((diff[k] * diff[k]).sum())

        return distances


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

    def evaluate(self, data, radii, mode='none'):
        """
        Evaluates the Ripley K estimator for a given set of values ``radii``.

        Parameters
        ----------
        data : 2D array
            Set of observed points in as a n by 2 array which will be used to
            estimate Ripley's K function.
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
                dimensions of the rectangular area of study. It assumes that
                all the points lie in a bounded rectangular region satisfying
                x_min < x_i < x_max; y_min < y_i < y_max.
        Returns
        -------
        ripley : 1D array
            Ripley's K function estimator evaluated at ``radii``.
        """

        data = np.asarray(data)

        if not data.shape[1] == 2:
            raise ValueError('data must be an n by 2 array, where n is the '
                             'number of observed points.')

        npts = len(data)
        ripley = np.zeros(len(radii))

        if mode == 'none':
            distances = self._pairwise_distances(data)
            for r in range(len(radii)):
                ripley[r] = (distances < radii[r]).sum()

            ripley = self.area * 2. * ripley / (npts * (npts - 1))
        # eq. 15.11 Stoyan book page 283
        elif mode == 'translation':
            diff = self._pairwise_diffs(data)
            distances = self._pairwise_distances(data)
            intersec_area = (((self.x_max - self.x_min) - diff[:][:, 0]) *
                             ((self.y_max - self.y_min) - diff[:][:, 1]))

            for r in range(len(radii)):
                dist_indicator = distances < radii[r]
                ripley[r] = ((1 / intersec_area) * dist_indicator).sum()

            ripley = (self.area**2 / (npts * (npts - 1))) * 2 * ripley
        # Stoyan book page 123 and eq 15.13
        elif mode == 'ohser':
            distances = self._pairwise_distances(data)
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

            for r in range(len(radii)):
                dist_indicator = distances < radii[r]
                ripley[r] = ((1 / cov_func) * dist_indicator).sum()

            ripley = (self.area**2 / (npts * (npts - 1))) * 2 * ripley
        # Cressie book eq 8.2.20 page 616
        elif mode == 'var-width':
            lt_dist = np.zeros(npts, dtype=np.double)
            for k in range(npts):
                lt_dist[k] = min(self.x_max - data[k][0],
                                 self.y_max - data[k][1],
                                 data[k][0] - self.x_min,
                                 data[k][1] - self.y_min)

            for r in range(len(radii)):
                for i in range(npts):
                    for j in range(npts):
                        if i != j:
                            diff = abs(data[i] - data[j])
                            dist = math.sqrt((diff * diff).sum())
                            if dist < radii[r] and lt_dist[i] > radii[r]:
                                ripley[r] = ripley[r] + 1
                if not (lt_dist > radii[r]).sum() == 0:
                    ripley[r] = ripley[r] / (lt_dist > radii[r]).sum()

            ripley = self.area * ripley / npts
        else:
            raise ValueError('mode {} is not implemented.'.format(mode))

        return ripley
