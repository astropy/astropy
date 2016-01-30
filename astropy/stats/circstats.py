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


def _components(data, deg=False, w=None, p=1, phi=0.0, axis=None):
    # Utility function for computing the generalized rectangular components
    # of the circular data.
    if w is None:
        w = np.ones(data.shape)
    elif (w.shape != data.shape):
        raise ValueError('shape of w %s is not equals shape of data %s.'
                         % (w.shape, data.shape,))

    if deg:
        data = Quantity(data, u.deg)
    else:
        data = Quantity(data, u.rad)

    C = np.sum(w * np.cos(p * (data - phi)), axis)/np.sum(w, axis)
    S = np.sum(w * np.sin(p * (data - phi)), axis)/np.sum(w, axis)

    return C, S


def _angle(data, deg=False, w=None, p=1, phi=0.0, axis=None):
    # Utility function for computing the generalized sample mean angle
    C, S = _components(data, deg, w, p, phi, axis)
    theta = np.arctan2(S, C)

    mask = (S == .0)*(C == .0)
    if np.any(mask):
        theta[mask] = np.nan

    if deg:
        theta = np.rad2deg(theta)
    return theta


def _length(data, deg=False, w=None, p=1, phi=0.0, axis=None):
    # Utility function for computing the generalized sample length
    C, S = _components(data, deg, w, p, phi, axis)
    return np.hypot(S, C)


def circmean(data, deg=False, w=None, axis=None):
    """ Computes the circular mean angle of an array of circular data.

    Parameters
    ----------
    data : numpy.ndarray
        Array of circular (directional) data.
    deg : boolean, optional
        If ``data`` is expressed in degrees, this should be ``True``.
        If ``deg`` is ``False``, the data is assumed to be measured in radians.
        The returned circular mean value will be expressed in radians or
        degrees depending on ``deg``.
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
    circmean : float
        Circular mean.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.stats import circmean
    >>> data_deg = np.array([51, 67, 40, 109, 31, 358])
    >>> data_rad = np.deg2rad(data_deg)
    >>> circmean(data_rad) # doctest: +FLOAT_CMP
    <Quantity 0.848704412445019 rad>
    >>> circmean(data_deg, deg=True) # doctest: +FLOAT_CMP
    <Quantity 48.62718088722989 deg>
    >>> from astropy import units as u
    >>> circmean(data_deg*u.deg, deg=True) # doctest: +FLOAT_CMP
    <Quantity 48.62718088722989 deg>
    >>> circmean(data_rad*u.rad) # doctest: +FLOAT_CMP
    <Quantity 0.848704412445019 rad>

    References
    ----------
    .. [1] S. R. Jammalamadaka, A. SenGupta. "Topics in Circular Statistics".
       Series on Multivariate Analysis, Vol. 5, 2001.
    .. [2] C. Agostinelli, U. Lund. "Circular Statistics from 'Topics in
       Circular Statistics (2001)'". 2015.
       <https://cran.r-project.org/web/packages/CircStats/CircStats.pdf>
    """
    return _angle(data, deg, w, 1, 0.0, axis)


def circvar(data, deg=False, w=None, axis=None):
    """ Computes the circular variance of an array of circular data.

    There are some concepts for defining measures of dispersion for circular
    data. The variance implemented here is based on the definition given by [1]
    , which is also the same used by the R package 'CircStats' [2].

    Parameters
    ----------
    data : numpy.ndarray
        Array of circular (directional) data.
    deg : boolean, optional
        If ``data`` is expressed in degrees, this should be ``True``.
        If ``deg`` is ``False``, the data is assumed to be measured in radians.
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
    circvar : float
        Circular variance.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.stats import circvar
    >>> data_deg = np.array([51, 67, 40, 109, 31, 358])
    >>> data_rad = np.deg2rad(data_deg)
    >>> circvar(data_rad) # doctest: +FLOAT_CMP
    <Quantity 0.16356352748437508>
    >>> circvar(data_deg, deg=True) # doctest: +FLOAT_CMP
    <Quantity 0.16356352748437508>
    >>> from astropy import units as u
    >>> circvar(data_rad*u.rad) # doctest: +FLOAT_CMP
    <Quantity 0.16356352748437508>
    >>> circvar(data_deg*u.deg, deg=True) # doctest: +FLOAT_CMP
    <Quantity 0.16356352748437508>

    References
    ----------
    .. [1] S. R. Jammalamadaka, A. SenGupta. "Topics in Circular Statistics".
       Series on Multivariate Analysis, Vol. 5, 2001.
    .. [2] C. Agostinelli, U. Lund. "Circular Statistics from 'Topics in
       Circular Statistics (2001)'". 2015.
       <https://cran.r-project.org/web/packages/CircStats/CircStats.pdf>
    """
    return 1.0 - _length(data, deg, w, 1, 0.0, axis)


def circmoment(data, deg=False, w=None, p=1.0, centered=False, axis=None):
    """ Computes the ``p``-th trigonometric circular moment for an array
    of circular data.

    Parameters
    ----------
    data : numpy.ndarray
        Array of circular (directional) data.
    deg : boolean, optional
        If ``data`` is expressed in degrees, this should be ``True``.
        If ``deg`` is ``False``, the data is assumed to be measured in radians.
        The returned circular moment direction will be expressed in radians or
        degrees depending on ``deg``.
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
    circmoment: np.ndarray
        The first are second elements correspond to the direction and length of
        the ``p``-th circular moment.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.stats import circmoment
    >>> data_deg = np.array([51, 67, 40, 109, 31, 358])
    >>> data_rad = np.deg2rad(data_deg)
    >>> circmoment(data_rad, p=2) # doctest: +FLOAT_CMP
    (<Quantity 1.5881210029361645 rad>, <Quantity 0.48004283892950717>)
    >>> circmoment(data_deg,deg=True,p=2) # doctest: +FLOAT_CMP
    (<Quantity 90.99263082432564 deg>, <Quantity 0.48004283892950717>)
    >>> from astropy import units as u
    >>> circmoment(data_rad*u.rad,p=2) # doctest: +FLOAT_CMP
    (<Quantity 1.5881210029361645 rad>, <Quantity 0.48004283892950717>)
    >>> circmoment(data_deg*u.deg,deg=True,p=2) # doctest: +FLOAT_CMP
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
        phi = circmean(data, deg, w, axis)

    return _angle(data, deg, w, p, phi, axis), _length(data, deg, w, p, phi,
                                                       axis)


def circcorrcoef(alpha, beta, deg=False, w_alpha=None, w_beta=None,
                 ax_alpha=None, ax_beta=None):
    """ Computes the circular correlation coefficient between two array of
    circular data.

    Parameters
    ----------
    alpha : numpy.ndarray
        Array of circular (directional) data.
    beta : numpy.ndarray
        Array of circular (directional) data.
    deg : boolean, optional
        If ``data`` is expressed in degrees, this should be ``True``.
        If ``deg`` is ``False``, the data is assumed to be measured in radians.
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
    rho : float
        Circular correlation coefficient.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.stats import circcorrcoef
    >>> alpha_deg = np.array([356, 97, 211, 232, 343, 292, 157, 302, 335, 302,
    ...                       324, 85, 324, 340, 157, 238, 254, 146, 232, 122,
    ...                       329])
    >>> beta_deg = np.array([119, 162, 221, 259, 270, 29, 97, 292, 40, 313, 94,
    ...                      45, 47, 108, 221, 270, 119, 248, 270, 45, 23])
    >>> alpha_rad = np.deg2rad(alpha_deg)
    >>> beta_rad = np.deg2rad(beta_deg)
    >>> circcorrcoef(alpha_rad, beta_rad) # doctest: +FLOAT_CMP
    <Quantity 0.2704648826748831>
    >>> circcorrcoef(alpha_deg, beta_deg, deg=True) # doctest: +FLOAT_CMP
    <Quantity 0.2704648826748831>
    >>> from astropy import units as u
    >>> circcorrcoef(alpha_rad*u.rad, beta_rad*u.rad) # doctest: +FLOAT_CMP
    <Quantity 0.2704648826748831>
    >>> circcorrcoef(alpha_deg*u.deg, beta_deg*u.deg,
    ...              deg=True) # doctest: +FLOAT_CMP
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

    if deg:
        alpha = Quantity(alpha, u.deg)
        beta = Quantity(beta, u.deg)
    else:
        alpha = Quantity(alpha, u.rad)
        beta = Quantity(beta, u.rad)

    mu_a = circmean(alpha, deg, w_alpha, ax_alpha)
    mu_b = circmean(beta, deg, w_beta, ax_beta)

    sin_a = np.sin(alpha - mu_a)
    sin_b = np.sin(beta - mu_b)
    rho = np.sum(sin_a*sin_b)/np.sqrt(np.sum(sin_a*sin_a)*np.sum(sin_b*sin_b))

    return rho


def rayleightest(data, deg=False, w=None, axis=None):
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
    deg : boolean, optional
        If ``data`` is expressed in degrees, this should be ``True``.
        If ``deg`` is ``False``, the data is assumed to be measured in radians.
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
    >>> data_rad = np.array([3.8955614, 3.1700932, 3.1804325, 2.7887885,
    ...                      2.9513829, 3.1140079, 4.5052974, 3.7984484,
    ...                      3.7510773, 2.9764646, 2.2100108, 2.9529469,
    ...                      2.9722527, 3.5099678, 2.6844710, 4.5020762,
    ...                      2.9285673, 2.2871037, 3.5805262, 3.4247286,
    ...                      2.4355773, 3.7549614, 4.1346585, 2.8426607,
    ...                      3.3291801, 2.5847791, 2.7217270, 3.4323084,
    ...                      3.3058256, 3.2147845, 2.4918005, 3.4849414,
    ...                      0.8249985, 3.4881397, 3.2070389, 3.0859854,
    ...                      4.0566486, 2.3463984, 3.7370984, 5.6853310,
    ...                      5.8108009, 3.3921987, 3.2166975, 3.6827617,
    ...                      4.5055291, 3.5295258, 3.0162183, 3.2317242,
    ...                      3.2778354, 3.1713455])
    >>> data_deg = np.rad2deg(data_rad)
    >>> rayleightest(data_rad) # doctest: +FLOAT_CMP
    <Quantity 2.726928722464598e-13>
    >>> rayleightest(data_deg, deg=True) # doctest: +FLOAT_CMP
    <Quantity 2.726928722464598e-13>
    >>> from astropy import units as u
    >>> rayleightest(data_rad*u.rad) # doctest: +FLOAT_CMP
    <Quantity 2.726928722464598e-13>
    >>> rayleightest(data_deg*u.deg,deg=True) # doctest: +FLOAT_CMP
    <Quantity 2.726928722464598e-13>

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
    Rbar = _length(data, deg, w, 1, 0.0, axis)
    z = n*Rbar*Rbar

    # see [3] and [4] for the formulae below
    tmp = 1.0
    if(n < 50):
        tmp = 1.0 + (2.0*z - z*z)/(4.0*n) - (24.0*z - 132.0*z**2.0 +
                                             76.0*z**3.0 - 9.0*z**4.0)/(288.0 *
                                                                        n * n)

    p_value = np.exp(-z)*tmp
    return p_value


def vtest(data, deg=False, w=None, mu=0.0, axis=None):
    """ Performs the Rayleigh test of uniformity where the alternative
    hypothesis H1 is assumed to have a known mean angle ``mu``.

    Parameters
    ----------
    data : numpy.ndarray
        Array of circular (directional) data.
    deg : boolean, optional
        If ``data`` is expressed in degrees, this should be ``True``.
        If ``deg`` is ``False``, the data is assumed to be measured in radians.
    w : numpy.ndarray, optional
        In case of grouped data, the i-th element of ``w`` represents a
        weighting factor for each group such that sum(w,axis) equals the number
        of observations. See [1], remark 1.4, page 22, for detailed
        explanation.
    mu : float, optional
        Assumed known mean angle.
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
    >>> data_rad = np.array([1.316075, 4.439193, 3.096231, 4.807068, 2.986021,
    ...                      1.756324, 3.046718, 3.299150, 3.360557,
    ...                      4.842499])*u.rad
    >>> data_deg = np.rad2deg(data_rad)
    >>> vtest(data_rad) # doctest: +FLOAT_CMP
    <Quantity 0.9871420365205558>
    >>> vtest(data_deg, deg=True) # doctest +FLOAT_CMP
    <Quantity 0.9871420365205558>

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


def vonmisesmle(data, deg=False, axis=None):
    """ Computes the Maximum Likelihood Estimator (MLE) for the parameters of
    the von Mises distribution.

    Parameters
    ----------
    data : numpy.ndarray
        Array of circular (directional) data.
    deg : boolean, optional
        If ``data`` is expressed in degrees, this should be ``True``.
        If ``deg`` is ``False``, the data is assumed to be measured in radians.
        The estimated circular mean (``mu``) will be expressed in radians or
        degrees depending on ``deg``.
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
    >>> data_rad = np.array([3.3699057, 4.0411630, 0.5014477, 2.6223103,
    ...                      3.7336524, 1.8136389, 4.1566039, 2.7806317,
    ...                      2.4672173, 2.8493644])
    >>> data_deg = np.rad2deg(data_rad)
    >>> vonmisesmle(data_rad) # doctest: +FLOAT_CMP
    (<Quantity 3.0065143178219063 rad>, <Quantity 1.474132390391715>)
    >>> vonmisesmle(data_deg, deg=True) # doctest: +FLOAT_CMP
    (<Quantity 172.26058145684905 deg>, <Quantity 1.474132390391715>)
    >>> from astropy import units as u
    >>> vonmisesmle(data_rad*u.rad) # doctest: +FLOAT_CMP
    (<Quantity 3.0065143178219063 rad>, <Quantity 1.474132390391715>)
    >>> vonmisesmle(data_deg*u.deg, deg=True) # doctest: +FLOAT_CMP
    (<Quantity 172.26058145684905 deg>, <Quantity 1.474132390391715>)

    References
    ----------
    .. [1] S. R. Jammalamadaka, A. SenGupta. "Topics in Circular Statistics".
       Series on Multivariate Analysis, Vol. 5, 2001.
    .. [2] C. Agostinelli, U. Lund. "Circular Statistics from 'Topics in
       Circular Statistics (2001)'". 2015.
       <https://cran.r-project.org/web/packages/CircStats/CircStats.pdf>
    """
    mu = circmean(data, deg, axis=None)

    if deg:
        #data = np.deg2rad(data)
        #mu_rad = np.deg2rad(mu)
        data = Quantity(data, u.deg)
        mu = Quantity(mu, u.deg)
    else:
        data = Quantity(data, u.rad)
        mu = Quantity(mu, u.rad)

    kappa = _A1inv(np.mean(np.cos(data - mu), axis))
    return mu, kappa
