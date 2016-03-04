# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains simple functions for dealing with circular statistics, for
instance, mean, variance, standard deviation, correlation coefficient, and so
on. This module also cover tests of uniformity, e.g., the Rayleigh and V tests.
The Maximum Likelihood Estimator for the Von Mises distribution along with the
Cramer-Rao Lower Bounds are also implemented. Almost all of the implementations
are based on reference [1], which is also the basis for the R package
'CircStats' [2].
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from astropy.units import Quantity
from astropy import units as u

__all__ = ['circmean', 'circvar', 'circmoment', 'circcorrcoef', 'rayleightest',
           'vtest', 'vonmisesmle']
__doctest_requires__ = {'vtest': ['scipy.stats']}


def _components(data, w=None, p=1, phi=0.0, axis=None):
    # Utility function for computing the generalized rectangular components
    # of the circular data.
    if w is None:
        w = np.ones(data.shape)
    elif (w.shape != data.shape):
        raise ValueError('shape of w %s is not equals shape of data %s.'
                         % (w.shape, data.shape,))

    C = np.sum(w * np.cos(p * (data - phi)), axis)/np.sum(w, axis)
    S = np.sum(w * np.sin(p * (data - phi)), axis)/np.sum(w, axis)

    return C, S


def _angle(data, w=None, p=1, phi=0.0, axis=None):
    # Utility function for computing the generalized sample mean angle
    C, S = _components(data, w, p, phi, axis)
    
    # theta will be an angle in the interval [-np.pi, np.pi)
    # [-180, 180)*u.deg in case data is a Quantity
    theta = np.arctan2(S, C)
    
    if isinstance(data, Quantity):
        if data.unit == u'deg':
            theta = np.rad2deg(theta)

    return theta


def _length(data, w=None, p=1, phi=0.0, axis=None):
    # Utility function for computing the generalized sample length
    C, S = _components(data, w, p, phi, axis)
    return np.hypot(S, C)


def circmean(data, w=None, axis=None):
    """ Computes the circular mean angle of an array of circular data.

    Parameters
    ----------
    data : numpy.ndarray
        Array of circular (directional) data.
    w : numpy.ndarray, optional
        In case of grouped data, the i-th element of ``w`` represents a
        weighting factor for each group such that sum(w,axis) equals the number
        of observations. See [1], remark 1.4, page 22, for detailed
        explanation.
    axis : int, optional
        Axis along which circular means are computed. The default is to compute
        the mean of the flattened array.

    Returns
    -------
    circmean : Quantity
        Circular mean.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.stats import circmean
    >>> from astropy import units as u
    >>> data = np.array([51, 67, 40, 109, 31, 358])*u.deg
    >>> circmean(data) # doctest: +FLOAT_CMP
    <Quantity 48.62718088722989 deg>

    References
    ----------
    .. [1] S. R. Jammalamadaka, A. SenGupta. "Topics in Circular Statistics".
       Series on Multivariate Analysis, Vol. 5, 2001.
    .. [2] C. Agostinelli, U. Lund. "Circular Statistics from 'Topics in
       Circular Statistics (2001)'". 2015.
       <https://cran.r-project.org/web/packages/CircStats/CircStats.pdf>
    """
    return _angle(data, w, 1, 0.0, axis)


def circvar(data, w=None, axis=None):
    """ Computes the circular variance of an array of circular data.

    There are some concepts for defining measures of dispersion for circular
    data. The variance implemented here is based on the definition given by [1]
    , which is also the same used by the R package 'CircStats' [2].

    Parameters
    ----------
    data : numpy.ndarray
        Array of circular (directional) data.
    w : numpy.ndarray, optional
        In case of grouped data, the i-th element of ``w`` represents a
        weighting factor for each group such that sum(w,axis) equals the number
        of observations. See [1], remark 1.4, page 22, for detailed
        explanation.
    axis : int, optional
        Axis along which circular variances are computed. The default is to
        compute the variance of the flattened array.

    Returns
    -------
    circvar : Quantity
        Circular variance.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.stats import circvar
    >>> from astropy import units as u
    >>> data = np.array([51, 67, 40, 109, 31, 358])*u.deg
    >>> circvar(data) # doctest: +FLOAT_CMP
    <Quantity 0.16356352748437508>

    References
    ----------
    .. [1] S. R. Jammalamadaka, A. SenGupta. "Topics in Circular Statistics".
       Series on Multivariate Analysis, Vol. 5, 2001.
    .. [2] C. Agostinelli, U. Lund. "Circular Statistics from 'Topics in
       Circular Statistics (2001)'". 2015.
       <https://cran.r-project.org/web/packages/CircStats/CircStats.pdf>
    
    Note
    ----
    The definition used here differs from the one in scipy.stats.circvar.
    Precisely, Scipy circvar uses an approximation based on the limit of small
    angles which approaches the linear variance.
    """

    return 1.0 - _length(data, w, 1, 0.0, axis)


def circmoment(data, w=None, p=1.0, centered=False, axis=None):
    """ Computes the ``p``-th trigonometric circular moment for an array
    of circular data.

    Parameters
    ----------
    data : numpy.ndarray
        Array of circular (directional) data.
    w : numpy.ndarray, optional
        In case of grouped data, the i-th element of ``w`` represents a
        weighting factor for each group such that sum(w,axis) equals the number
        of observations. See [1], remark 1.4, page 22, for detailed
        explanation.
    p : float, optional
        Order of the circular moment.
    centered : Boolean, optinal
        If ``True``, central circular moments are computed. Default value is
        ``False``.
    axis : int, optional
        Axis along which circular moments are computed. The default is to
        compute the circular moment of the flattened array.

    Returns
    -------
    circmoment: Quantity
        The first and second elements correspond to the direction and length of
        the ``p``-th circular moment, respectively.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.stats import circmoment
    >>> from astropy import units as u
    >>> data = np.array([51, 67, 40, 109, 31, 358])*u.deg
    >>> circmoment(data, p=2) # doctest: +FLOAT_CMP
    (<Quantity 90.99263082432564 deg>, <Quantity 0.48004283892950717>)

    References
    ----------
    .. [1] S. R. Jammalamadaka, A. SenGupta. "Topics in Circular Statistics".
       Series on Multivariate Analysis, Vol. 5, 2001.
    .. [2] C. Agostinelli, U. Lund. "Circular Statistics from 'Topics in
       Circular Statistics (2001)'". 2015.
       <https://cran.r-project.org/web/packages/CircStats/CircStats.pdf>
    """
    phi = 0.0
    if centered:
        phi = circmean(data, w, axis)

    return _angle(data, w, p, phi, axis), _length(data, w, p, phi, axis)


def circcorrcoef(alpha, beta, w_alpha=None, w_beta=None, ax_alpha=None,
                 ax_beta=None):
    """ Computes the circular correlation coefficient between two array of
    circular data.

    Parameters
    ----------
    alpha : numpy.ndarray
        Array of circular (directional) data.
    beta : numpy.ndarray
        Array of circular (directional) data.
    w_alpha : numpy.ndarray, optional
        In case of grouped data, the i-th element of ``w_alpha`` represents a
        weighting factor for each group such that sum(w,axis) equals the number
        of observations. See [1], remark 1.4, page 22, for detailed
        explanation.
    w_beta : numpy.ndarray, optional
        See description of ``w_alpha``.
    ax_alpha : int, optional
        Axis along which circular correlation coefficients are computed.
        The default is the compute the circular correlation coefficient of the
        flattened array.
    ax_ beta : int, optional
        See description of ``ax_alpha``.

    Returns
    -------
    rho : Quantity
        Circular correlation coefficient.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.stats import circcorrcoef
    >>> from astropy import units as u
    >>> alpha = np.array([356, 97, 211, 232, 343, 292, 157, 302, 335, 302,
    ...                   324, 85, 324, 340, 157, 238, 254, 146, 232, 122,
    ...                   329])*u.deg
    >>> beta = np.array([119, 162, 221, 259, 270, 29, 97, 292, 40, 313, 94,
    ...                  45, 47, 108, 221, 270, 119, 248, 270, 45, 23])*u.deg
    >>> circcorrcoef(alpha, beta) # doctest: +FLOAT_CMP
    <Quantity 0.2704648826748831>

    References
    ----------
    .. [1] S. R. Jammalamadaka, A. SenGupta. "Topics in Circular Statistics".
       Series on Multivariate Analysis, Vol. 5, 2001.
    .. [2] C. Agostinelli, U. Lund. "Circular Statistics from 'Topics in
       Circular Statistics (2001)'". 2015.
       <https://cran.r-project.org/web/packages/CircStats/CircStats.pdf>
    """
    if(np.size(alpha, axis=ax_alpha) != np.size(beta, axis=ax_beta)):
        raise ValueError("alpha and beta must be arrays of the same size")

    mu_a = circmean(alpha, w_alpha, ax_alpha)
    mu_b = circmean(beta, w_beta, ax_beta)

    sin_a = np.sin(alpha - mu_a)
    sin_b = np.sin(beta - mu_b)
    rho = np.sum(sin_a*sin_b)/np.sqrt(np.sum(sin_a*sin_a)*np.sum(sin_b*sin_b))

    return rho


def rayleightest(data, w=None, axis=None):
    """ Performs the Rayleigh test of uniformity.

    This test is  used to indentify a non-uniform distribution, i.e. it is
    designed for detecting an unimodal deviation from uniformity. More
    precisely, it assumes the following hypotheses:
    - H0 (null hypothesis): The population is distributed uniformly around the
    circle.
    - H1 (alternative hypothesis): The population is not distributed uniformly
    around the circle.
    Small p-values suggest to reject the null hypothesis.

    Parameters
    ----------
    data : numpy.ndarray
        Array of circular (directional) data.
    w : numpy.ndarray, optional
        In case of grouped data, the i-th element of ``w`` represents a
        weighting factor for each group such that sum(w,axis) equals the number
        of observations. See [1], remark 1.4, page 22, for detailed
        explanation.
    axis : int, optional
        Axis along which the Rayleigh test will be performed.

    Returns
    -------
    p-value : float
        p-value.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.stats import rayleightest
    >>> from astropy import units as u
    >>> data = np.array([130, 90, 0, 145])*u.deg
    >>> rayleightest(data) # doctest: +FLOAT_CMP
    <Quantity 0.2563487733797317>

    References
    ----------
    .. [1] S. R. Jammalamadaka, A. SenGupta. "Topics in Circular Statistics".
       Series on Multivariate Analysis, Vol. 5, 2001.
    .. [2] C. Agostinelli, U. Lund. "Circular Statistics from 'Topics in
       Circular Statistics (2001)'". 2015.
       <https://cran.r-project.org/web/packages/CircStats/CircStats.pdf>
    .. [3] M. Chirstman., C. Miller. "Testing a Sample of Directions for
       Uniformity." Lecture Notes, STA 6934/5805. University of Florida, 2007.
    .. [4] D. Wilkie. "Rayleigh Test for Randomness of Circular Data". Applied
       Statistics. 1983.
       <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.211.4762>
    """
    n = np.size(data, axis=axis)
    Rbar = _length(data, w, 1, 0.0, axis)
    z = n*Rbar*Rbar

    # see [3] and [4] for the formulae below
    tmp = 1.0
    if(n < 50):
        tmp = 1.0 + (2.0*z - z*z)/(4.0*n) - (24.0*z - 132.0*z**2.0 +
                                             76.0*z**3.0 - 9.0*z**4.0)/(288.0 *
                                                                        n * n)

    p_value = np.exp(-z)*tmp
    return p_value


def vtest(data, w=None, mu=0.0, axis=None):
    """ Performs the Rayleigh test of uniformity where the alternative
    hypothesis H1 is assumed to have a known mean angle ``mu``.

    Parameters
    ----------
    data : numpy.ndarray
        Array of circular (directional) data.
    w : numpy.ndarray, optional
        In case of grouped data, the i-th element of ``w`` represents a
        weighting factor for each group such that sum(w,axis) equals the number
        of observations. See [1], remark 1.4, page 22, for detailed
        explanation.
    mu : float, optional
        Mean angle. Assumed to be known.
    axis : int, optional
        Axis along which the V test will be performed.

    Returns
    -------
    p-value : float
        p-value.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.stats import vtest
    >>> from astropy import units as u
    >>> data = np.array([130, 90, 0, 145])*u.deg
    >>> vtest(data) # doctest: +FLOAT_CMP
    <Quantity 0.6223678199713766>

    References
    ----------
    .. [1] S. R. Jammalamadaka, A. SenGupta. "Topics in Circular Statistics".
       Series on Multivariate Analysis, Vol. 5, 2001.
    .. [2] C. Agostinelli, U. Lund. "Circular Statistics from 'Topics in
       Circular Statistics (2001)'". 2015.
       <https://cran.r-project.org/web/packages/CircStats/CircStats.pdf>
    .. [3] M. Chirstman., C. Miller. "Testing a Sample of Directions for
       Uniformity." Lecture Notes, STA 6934/5805. University of Florida, 2007.
    """
    from scipy.stats import norm

    if w is None:
        w = np.ones(data.shape)
    elif (w.shape != data.shape):
        raise ValueError('shape of w %s is not equals shape of data %s.'
                         % (w.shape, data.shape,))

    n = np.size(data, axis=axis)

    R0bar = np.sum(w*np.cos(data-mu), axis)/np.sum(w, axis)
    z = np.sqrt(2.0*n)*R0bar

    pz = norm.cdf(z)
    fz = norm.pdf(z)

    # see reference [3]
    p_value = 1 - pz + fz*((3*z - z**3)/(16.0*n) +
                           (15*z + 305*z**3 - 125*z**5 + 9*z**7)/(4608.0*n*n))
    return p_value


def _A1inv(x):
    # Approximation for _A1inv(x) according R Package 'CircStats'
    # See http://www.scienceasia.org/2012.38.n1/scias38_118.pdf, equation (4)
    if(x >= 0 and x < 0.53):
        return 2.0*x + x*x*x + (5.0*x**5)/6.0
    elif x < 0.85:
        return -0.4 + 1.39*x + 0.43/(1.0 - x)
    else:
        return 1.0/(x*x*x - 4.0*x*x + 3.0*x)


def vonmisesmle(data, axis=None):
    """ Computes the Maximum Likelihood Estimator (MLE) for the parameters of
    the von Mises distribution.

    Parameters
    ----------
    data : numpy.ndarray
        Array of circular (directional) data.
    axis : int, optional
        Axis along which the mle will be computed.

    Returns
    -------
    mu : float
        the mean (aka location parameter).
    kappa : float
        the concentration parameter.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.stats import vonmisesmle
    >>> from astropy import units as u
    >>> data = np.array([130, 90, 0, 145])*u.deg
    >>> vonmisesmle(data) # doctest: +FLOAT_CMP
    (<Quantity 101.16894320013179 deg>, <Quantity 1.49358958737054>)

    References
    ----------
    .. [1] S. R. Jammalamadaka, A. SenGupta. "Topics in Circular Statistics".
       Series on Multivariate Analysis, Vol. 5, 2001.
    .. [2] C. Agostinelli, U. Lund. "Circular Statistics from 'Topics in
       Circular Statistics (2001)'". 2015.
       <https://cran.r-project.org/web/packages/CircStats/CircStats.pdf>
    """
    mu = circmean(data, axis=None)

    kappa = _A1inv(np.mean(np.cos(data - mu), axis))
    return mu, kappa
