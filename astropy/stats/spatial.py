# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module implements functions and classes for spatial statistics.
"""


import numpy as np
import math


class RipleysKEstimator:
    """
    Estimators for Ripley's K function for two-dimensional spatial data.
    See [1]_, [2]_, [3]_, [4]_, [5]_ for detailed mathematical and
    practical aspects of those estimators.

    Parameters
    ----------
    area : float
        Area of study from which the points where observed.
    x_max, y_max : float, float, optional
        Maximum rectangular coordinates of the area of study.
        Required if ``mode == 'translation'`` or ``mode == ohser``.
    x_min, y_min : float, float, optional
        Minimum rectangular coordinates of the area of study.
        Required if ``mode == 'variable-width'`` or ``mode == ohser``.

    Examples
    --------
    >>> import numpy as np
    >>> from matplotlib import pyplot as plt # doctest: +SKIP
    >>> from astropy.stats import RipleysKEstimator
    >>> z = np.random.uniform(low=5, high=10, size=(100, 2))
    >>> Kest = RipleysKEstimator(area=25, x_max=10, y_max=10,
    ... x_min=5, y_min=5)
    >>> r = np.linspace(0, 2.5, 100)
    >>> plt.plot(r, Kest.poisson(r)) # doctest: +SKIP
    >>> plt.plot(r, Kest(data=z, radii=r, mode='none')) # doctest: +SKIP
    >>> plt.plot(r, Kest(data=z, radii=r, mode='translation')) # doctest: +SKIP
    >>> plt.plot(r, Kest(data=z, radii=r, mode='ohser')) # doctest: +SKIP
    >>> plt.plot(r, Kest(data=z, radii=r, mode='var-width')) # doctest: +SKIP
    >>> plt.plot(r, Kest(data=z, radii=r, mode='ripley')) # doctest: +SKIP

    References
    ----------
    .. [1] Peebles, P.J.E. *The large scale structure of the universe*.
       <http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=1980lssu.book.....P&db_key=AST>
    .. [2] Spatial descriptive statistics.
       <https://en.wikipedia.org/wiki/Spatial_descriptive_statistics>
    .. [3] Package spatstat.
       <https://cran.r-project.org/web/packages/spatstat/spatstat.pdf>
    .. [4] Cressie, N.A.C. (1991). Statistics for Spatial Data,
       Wiley, New York.
    .. [5] Stoyan, D., Stoyan, H. (1992). Fractals, Random Shapes and
       Point Fields, Akademie Verlag GmbH, Chichester.
    """

    def __init__(self, area, x_max=None, y_max=None, x_min=None, y_min=None):
        self.area = area
        self.x_max = x_max
        self.y_max = y_max
        self.x_min = x_min
        self.y_min = y_min

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
        if value is None or isinstance(value, (float, int)):
            self._y_max = value
        else:
            raise ValueError('y_max is expected to be a real number '
                             'or None. Got {}.'.format(value))

    @property
    def x_max(self):
        return self._x_max

    @x_max.setter
    def x_max(self, value):
        if value is None or isinstance(value, (float, int)):
            self._x_max = value
        else:
            raise ValueError('x_max is expected to be a real number '
                             'or None. Got {}.'.format(value))

    @property
    def y_min(self):
        return self._y_min

    @y_min.setter
    def y_min(self, value):
        if value is None or isinstance(value, (float, int)):
            self._y_min = value
        else:
            raise ValueError('y_min is expected to be a real number. '
                             'Got {}.'.format(value))

    @property
    def x_min(self):
        return self._x_min

    @x_min.setter
    def x_min(self, value):
        if value is None or isinstance(value, (float, int)):
            self._x_min = value
        else:
            raise ValueError('x_min is expected to be a real number. '
                             'Got {}.'.format(value))

    def __call__(self, data, radii, mode='none'):
        return self.evaluate(data=data, radii=radii, mode=mode)

    def _pairwise_diffs(self, data):
        npts = len(data)
        diff = np.zeros(shape=(npts * (npts - 1) // 2, 2), dtype=np.double)
        k = 0
        for i in range(npts - 1):
            size = npts - i - 1
            diff[k:k + size] = abs(data[i] - data[i+1:])
            k += size

        return diff

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

    def Lfunction(self, data, radii, mode='none'):
        """
        Evaluates the L function at ``radii``. For parameter description
        see ``evaluate`` method.
        """

        return np.sqrt(self.evaluate(data, radii, mode=mode) / np.pi)

    def Hfunction(self, data, radii, mode='none'):
        """
        Evaluates the H function at ``radii``. For parameter description
        see ``evaluate`` method.
        """

        return self.Lfunction(data, radii, mode=mode) - radii

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
            Available methods are 'none', 'translation', 'ohser', 'var-width',
            and 'ripley'.

            * 'none'
                this method does not take into account any edge effects
                whatsoever.
            * 'translation'
                computes the intersection of rectangular areas centered at
                the given points provided the upper bounds of the
                dimensions of the rectangular area of study. It assumes that
                all the points lie in a bounded rectangular region satisfying
                x_min < x_i < x_max; y_min < y_i < y_max. A detailed
                description of this method can be found on ref [4].
            * 'ohser'
                this method uses the isotropized set covariance function of
                the window of study as a weigth to correct for
                edge-effects. A detailed description of this method can be
                found on ref [4].
            * 'var-width'
                this method considers the distance of each observed point to
                the nearest boundary of the study window as a factor to
                account for edge-effects. See [3] for a brief description of
                this method.
            * 'ripley'
                this method is known as Ripley's edge-corrected estimator.
                The weight for edge-correction is a function of the
                proportions of circumferences centered at each data point
                which crosses another data point of interest. See [3] for
                a detailed description of this method.

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
            diff = self._pairwise_diffs(data)
            distances = np.hypot(diff[:, 0], diff[:, 1])
            for r in range(len(radii)):
                ripley[r] = (distances < radii[r]).sum()

            ripley = self.area * 2. * ripley / (npts * (npts - 1))
        # eq. 15.11 Stoyan book page 283
        elif mode == 'translation':
            diff = self._pairwise_diffs(data)
            distances = np.hypot(diff[:, 0], diff[:, 1])
            intersec_area = (((self.x_max - self.x_min) - diff[:, 0]) *
                             ((self.y_max - self.y_min) - diff[:, 1]))

            for r in range(len(radii)):
                dist_indicator = distances < radii[r]
                ripley[r] = ((1 / intersec_area) * dist_indicator).sum()

            ripley = (self.area**2 / (npts * (npts - 1))) * 2 * ripley
        # Stoyan book page 123 and eq 15.13
        elif mode == 'ohser':
            diff = self._pairwise_diffs(data)
            distances = np.hypot(diff[:, 0], diff[:, 1])
            a = self.area
            b = max((self.y_max - self.y_min) / (self.x_max - self.x_min),
                    (self.x_max - self.x_min) / (self.y_max - self.y_min))
            x = distances / math.sqrt(a / b)
            u = np.sqrt((x * x - 1) * (x > 1))
            v = np.sqrt((x * x - b ** 2) * (x < math.sqrt(b ** 2 + 1)) * (x > b))
            c1 = np.pi - 2 * x * (1 + 1 / b) + x * x / b
            c2 = 2 * np.arcsin((1 / x) * (x > 1)) - 1 / b - 2 * (x - u)
            c3 = (2 * np.arcsin(((b - u * v) / (x * x))
                                * (x > b) * (x < math.sqrt(b ** 2 + 1)))
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
            lt_dist = np.minimum(np.minimum(self.x_max - data[:, 0], self.y_max - data[:, 1]),
                                 np.minimum(data[:, 0] - self.x_min, data[:, 1] - self.y_min))

            for r in range(len(radii)):
                for i in range(npts):
                    for j in range(npts):
                        if i != j:
                            diff = abs(data[i] - data[j])
                            dist = math.sqrt((diff * diff).sum())
                            if dist < radii[r] < lt_dist[i]:
                                ripley[r] = ripley[r] + 1
                lt_dist_sum = (lt_dist > radii[r]).sum()
                if not lt_dist_sum == 0:
                    ripley[r] = ripley[r] / lt_dist_sum

            ripley = self.area * ripley / npts
        # Cressie book eq 8.4.22 page 640
        elif mode == 'ripley':
            hor_dist = np.zeros(shape=(npts * (npts - 1)) // 2,
                                dtype=np.double)
            ver_dist = np.zeros(shape=(npts * (npts - 1)) // 2,
                                dtype=np.double)

            for k in range(npts - 1):
                min_hor_dist = min(self.x_max - data[k][0],
                                   data[k][0] - self.x_min)
                min_ver_dist = min(self.y_max - data[k][1],
                                   data[k][1] - self.y_min)
                start = (k * (2 * (npts - 1) - (k - 1))) // 2
                end = ((k + 1) * (2 * (npts - 1) - k)) // 2
                hor_dist[start: end] = min_hor_dist * np.ones(npts - 1 - k)
                ver_dist[start: end] = min_ver_dist * np.ones(npts - 1 - k)

            diff = self._pairwise_diffs(data)
            dist = np.hypot(diff[:, 0], diff[:, 1])
            dist_ind = dist <= np.hypot(hor_dist, ver_dist)

            w1 = (1 - (np.arccos(np.minimum(ver_dist, dist) / dist) +
                       np.arccos(np.minimum(hor_dist, dist) / dist)) / np.pi)
            w2 = (3 / 4 - 0.5 * (np.arccos(ver_dist / dist * ~dist_ind) +
                            np.arccos(hor_dist / dist * ~dist_ind)) / np.pi)

            weight = dist_ind * w1 + ~dist_ind * w2

            for r in range(len(radii)):
                ripley[r] = ((dist < radii[r]) / weight).sum()

            ripley = self.area * 2. * ripley / (npts * (npts - 1))
        else:
            raise ValueError('mode {} is not implemented.'.format(mode))

        return ripley
