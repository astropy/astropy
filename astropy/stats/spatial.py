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

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from astropy.stats import RipleysKEstimate
    >>> x = np.random.uniform(low=5, high=10, size=100)
    >>> y = np.random.uniform(low=5, high=10, size=100)
    >>> z = np.array([x,y]).T
    >>> area = 25
    >>> Kest = RipleysKEstimate(data=z, area=area)
    >>> r = np.linspace(0, 2.5, 100)
    >>> plt.plot(r, Kest(r))

    References
    ----------
    .. [1] Spatial descriptive statistics.
       <https://en.wikipedia.org/wiki/Spatial_descriptive_statistics>
    .. [2] Package spatstat.
       <https://cran.r-project.org/web/packages/spatstat/spatstat.pdf>
    .. [3] Cressie, N.A.C. (1991). Statistics for Spatial Data,
       Wiley, New York.
    """

    def __init__(self, data, area):
        self.data = data
        self.area = area

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
        if not isinstance(value, (float, int)) or value < 0:
            raise ValueError('density is expected to be positive either '
                             'float or int. Got {}.'.format(type(value)))
        elif value < 0:
            raise ValueError('crit_separation is expected to be a positive '
                             'real number. Got {}'.format(value))
        else:
            self._area = value

    def __call__(self, radii):
        return self.evaluate(radii=radii)

    def evaluate(self, radii):
        """
        Evaluates the Ripley K estimator for a given set of values ``radii``.

        Parameters
        ----------
        radii : 1D array
            Set of distances in which Ripley's K estimator will be evaluated.
            Usually, it's common to consider max(radii) < (area/2)**0.5.

        Returns
        -------
        output : 1D array
            Ripley's K function estimator evaluated at ``radii``.
        """

        self.data = np.asarray(self.data)
        npts = len(self.data)
        distances = np.zeros((npts * (npts - 1)) // 2, dtype=np.double)

        k = 0
        for i in range(0, npts - 1):
            for j in range(i + 1, npts):
                diff = abs(self.data[i] - self.data[j])
                distances[k] = math.sqrt((diff * diff).sum())
                k = k + 1

        ripley = np.zeros(len(radii))

        for idx in range(len(radii)):
            ripley[idx] = (distances < radii[idx]).sum()


        return self.area * 2. * ripley / (npts * (npts - 1))
