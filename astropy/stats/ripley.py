class RipleysKEstimate(object):
    """
    This class implements an estimator for Ripley's K function in a two
    dimensional space.

    Parameters
    ----------
    data : 2D array
        Set of observed points in a 2D space which will be used to estimate
        Ripley's K function.
    area : float
        Area of study from which the points where observed.

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> a = np.array([[1, 4], [2, 5], [3, 6]])
    >>> area = 100.
    >>> Kest = RipleysKEstimate(data=a, area=area)
    >>> r = np.linspace(0, 2.5, 100)
    >>> plt.plot(r, Kest(r))
    >>> plt.plot(r, np.pi * r**2)

    References
    ----------
        [1] Spatial descriptive statistics.
        https://en.wikipedia.org/wiki/Spatial_descriptive_statistics
        [2] Package spatstat.
        https://cran.r-project.org/web/packages/spatstat/spatstat.pdf
        [3] Cressie, N.A.C. (1991). Statistics for Spatial Data,
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

    def __call__(self, radii, **kwargs):
        return self.evaluate(radii=radii, **kwargs)

    def evaluate(self, radii, metric='euclidean', **kwargs):
        """
        Evaluates the Ripley K estimator for a given set of distances
        ``radii``.

        Parameters
        ----------
        radii : 1D array
            Set of distances in which Ripley's K estimator will be evaluated.
        metric : str
            Type of distance.
        kwargs : dict
            Keyword argument for scipy's pdist function.

        Returns
        -------
        output : 1D array
            Ripley's K function estimator evaluated at ``radii``.
        """

        from scipy.spatial.distance import pdist

        distances = pdist(self.data, metric=metric, **kwargs)
        ripley = np.zeros(len(radii))

        for idx in range(len(radii)):
            ripley[idx] = (distances < radii[idx]).sum()

        npts = len(self.data)

        return self.area * 2. * ripley / (npts * (npts - 1))
