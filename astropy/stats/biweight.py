# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains functions for computing robust statistics using
Tukey's biweight function.
"""


import numpy as np

from .funcs import median_absolute_deviation
from ..utils.decorators import deprecated_renamed_argument


__all__ = ['biweight_location', 'biweight_scale', 'biweight_midvariance',
           'biweight_midcovariance', 'biweight_midcorrelation']


@deprecated_renamed_argument('a', 'data', '2.0')
def biweight_location(data, c=6.0, M=None, axis=None):
    r"""
    Compute the biweight location.

    The biweight location is a robust statistic for determining the
    central location of a distribution.  It is given by:

    .. math::

        \zeta_{biloc}= M + \frac{\Sigma_{|u_i|<1} \ (x_i - M) (1 - u_i^2)^2}
            {\Sigma_{|u_i|<1} \ (1 - u_i^2)^2}

    where :math:`x` is the input data, :math:`M` is the sample median
    (or the input initial location guess) and :math:`u_i` is given by:

    .. math::

        u_{i} = \frac{(x_i - M)}{c * MAD}

    where :math:`c` is the tuning constant and :math:`MAD` is the
    `median absolute deviation
    <https://en.wikipedia.org/wiki/Median_absolute_deviation>`_.  The
    biweight location tuning constant ``c`` is typically 6.0 (the
    default).

    Parameters
    ----------
    data : array-like
        Input array or object that can be converted to an array.
    c : float, optional
        Tuning constant for the biweight estimator (default = 6.0).
    M : float or array-like, optional
        Initial guess for the location.  If ``M`` is a scalar value,
        then its value will be used for the entire array (or along each
        ``axis``, if specified).  If ``M`` is an array, then its must be
        an array containing the initial location estimate along each
        ``axis`` of the input array.  If `None` (default), then the
        median of the input array will be used (or along each ``axis``,
        if specified).
    axis : int, optional
        The axis along which the biweight locations are computed.  If
        `None` (default), then the biweight location of the flattened
        input array will be computed.

    Returns
    -------
    biweight_location : float or `~numpy.ndarray`
        The biweight location of the input data.  If ``axis`` is `None`
        then a scalar will be returned, otherwise a `~numpy.ndarray`
        will be returned.

    See Also
    --------
    biweight_scale, biweight_midvariance, biweight_midcovariance

    References
    ----------
    .. [1] Beers, Flynn, and Gebhardt (1990; AJ 100, 32) (http://adsabs.harvard.edu/abs/1990AJ....100...32B)

    .. [2] http://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwloc.htm

    Examples
    --------
    Generate random variates from a Gaussian distribution and return the
    biweight location of the distribution:

    >>> import numpy as np
    >>> from astropy.stats import biweight_location
    >>> rand = np.random.RandomState(12345)
    >>> biloc = biweight_location(rand.randn(1000))
    >>> print(biloc)    # doctest: +FLOAT_CMP
    -0.0175741540445
    """

    data = np.asanyarray(data).astype(np.float64)

    if M is None:
        M = np.median(data, axis=axis)
    if axis is not None:
        M = np.expand_dims(M, axis=axis)

    # set up the differences
    d = data - M

    # set up the weighting
    mad = median_absolute_deviation(data, axis=axis)
    if axis is not None:
        mad = np.expand_dims(mad, axis=axis)
    u = d / (c * mad)

    # now remove the outlier points
    mask = (np.abs(u) >= 1)
    u = (1 - u ** 2) ** 2
    u[mask] = 0

    return M.squeeze() + (d * u).sum(axis=axis) / u.sum(axis=axis)


def biweight_scale(data, c=9.0, M=None, axis=None, modify_sample_size=False):
    r"""
    Compute the biweight scale.

    The biweight scale is a robust statistic for determining the
    standard deviation of a distribution.  It is the square root of the
    `biweight midvariance
    <https://en.wikipedia.org/wiki/Robust_measures_of_scale#The_biweight_midvariance>`_.
    It is given by:

    .. math::

        \zeta_{biscl} = \sqrt{n} \ \frac{\sqrt{\Sigma_{|u_i| < 1} \
            (x_i - M)^2 (1 - u_i^2)^4}} {|(\Sigma_{|u_i| < 1} \
            (1 - u_i^2) (1 - 5u_i^2))|}

    where :math:`x` is the input data, :math:`M` is the sample median
    (or the input location) and :math:`u_i` is given by:

    .. math::

        u_{i} = \frac{(x_i - M)}{c * MAD}

    where :math:`c` is the tuning constant and :math:`MAD` is the
    `median absolute deviation
    <https://en.wikipedia.org/wiki/Median_absolute_deviation>`_.  The
    biweight midvariance tuning constant ``c`` is typically 9.0 (the
    default).

    For the standard definition of biweight scale, :math:`n` is the
    total number of points in the array (or along the input ``axis``, if
    specified).  That definition is used if ``modify_sample_size`` is
    `False`, which is the default.

    However, if ``modify_sample_size = True``, then :math:`n` is the
    number of points for which :math:`|u_i| < 1` (i.e. the total number
    of non-rejected values), i.e.

    .. math::

        n = \Sigma_{|u_i| < 1} \ 1

    which results in a value closer to the true standard deviation for
    small sample sizes or for a large number of rejected values.

    Parameters
    ----------
    data : array-like
        Input array or object that can be converted to an array.
    c : float, optional
        Tuning constant for the biweight estimator (default = 9.0).
    M : float or array-like, optional
        The location estimate.  If ``M`` is a scalar value, then its
        value will be used for the entire array (or along each ``axis``,
        if specified).  If ``M`` is an array, then its must be an array
        containing the location estimate along each ``axis`` of the
        input array.  If `None` (default), then the median of the input
        array will be used (or along each ``axis``, if specified).
    axis : int, optional
        The axis along which the biweight scales are computed.  If
        `None` (default), then the biweight scale of the flattened input
        array will be computed.
    modify_sample_size : bool, optional
        If `False` (default), then the sample size used is the total
        number of elements in the array (or along the input ``axis``, if
        specified), which follows the standard definition of biweight
        scale.  If `True`, then the sample size is reduced to correct
        for any rejected values (i.e. the sample size used includes only
        the non-rejected values), which results in a value closer to the
        true standard deviation for small sample sizes or for a large
        number of rejected values.

    Returns
    -------
    biweight_scale : float or `~numpy.ndarray`
        The biweight scale of the input data.  If ``axis`` is `None`
        then a scalar will be returned, otherwise a `~numpy.ndarray`
        will be returned.

    See Also
    --------
    biweight_midvariance, biweight_midcovariance, biweight_location, astropy.stats.mad_std, astropy.stats.median_absolute_deviation

    References
    ----------
    .. [1] Beers, Flynn, and Gebhardt (1990; AJ 100, 32) (http://adsabs.harvard.edu/abs/1990AJ....100...32B)

    .. [2] http://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwscale.htm

    Examples
    --------
    Generate random variates from a Gaussian distribution and return the
    biweight scale of the distribution:

    >>> import numpy as np
    >>> from astropy.stats import biweight_scale
    >>> rand = np.random.RandomState(12345)
    >>> biscl = biweight_scale(rand.randn(1000))
    >>> print(biscl)    # doctest: +FLOAT_CMP
    0.986726249291
    """

    return np.sqrt(
        biweight_midvariance(data, c=c, M=M, axis=axis,
                             modify_sample_size=modify_sample_size))


@deprecated_renamed_argument('a', 'data', '2.0')
def biweight_midvariance(data, c=9.0, M=None, axis=None,
                         modify_sample_size=False):
    r"""
    Compute the biweight midvariance.

    The biweight midvariance is a robust statistic for determining the
    variance of a distribution.  Its square root is a robust estimator
    of scale (i.e. standard deviation).  It is given by:

    .. math::

        \zeta_{bivar} = n \ \frac{\Sigma_{|u_i| < 1} \
            (x_i - M)^2 (1 - u_i^2)^4} {(\Sigma_{|u_i| < 1} \
            (1 - u_i^2) (1 - 5u_i^2))^2}

    where :math:`x` is the input data, :math:`M` is the sample median
    (or the input location) and :math:`u_i` is given by:

    .. math::

        u_{i} = \frac{(x_i - M)}{c * MAD}

    where :math:`c` is the tuning constant and :math:`MAD` is the
    `median absolute deviation
    <https://en.wikipedia.org/wiki/Median_absolute_deviation>`_.  The
    biweight midvariance tuning constant ``c`` is typically 9.0 (the
    default).

    For the standard definition of `biweight midvariance
    <https://en.wikipedia.org/wiki/Robust_measures_of_scale#The_biweight_midvariance>`_,
    :math:`n` is the total number of points in the array (or along the
    input ``axis``, if specified).  That definition is used if
    ``modify_sample_size`` is `False`, which is the default.

    However, if ``modify_sample_size = True``, then :math:`n` is the
    number of points for which :math:`|u_i| < 1` (i.e. the total number
    of non-rejected values), i.e.

    .. math::

        n = \Sigma_{|u_i| < 1} \ 1

    which results in a value closer to the true variance for small
    sample sizes or for a large number of rejected values.

    Parameters
    ----------
    data : array-like
        Input array or object that can be converted to an array.
    c : float, optional
        Tuning constant for the biweight estimator (default = 9.0).
    M : float or array-like, optional
        The location estimate.  If ``M`` is a scalar value, then its
        value will be used for the entire array (or along each ``axis``,
        if specified).  If ``M`` is an array, then its must be an array
        containing the location estimate along each ``axis`` of the
        input array.  If `None` (default), then the median of the input
        array will be used (or along each ``axis``, if specified).
    axis : int, optional
        The axis along which the biweight midvariances are computed.  If
        `None` (default), then the biweight midvariance of the flattened
        input array will be computed.
    modify_sample_size : bool, optional
        If `False` (default), then the sample size used is the total
        number of elements in the array (or along the input ``axis``, if
        specified), which follows the standard definition of biweight
        midvariance.  If `True`, then the sample size is reduced to
        correct for any rejected values (i.e. the sample size used
        includes only the non-rejected values), which results in a value
        closer to the true variance for small sample sizes or for a
        large number of rejected values.

    Returns
    -------
    biweight_midvariance : float or `~numpy.ndarray`
        The biweight midvariance of the input data.  If ``axis`` is
        `None` then a scalar will be returned, otherwise a
        `~numpy.ndarray` will be returned.

    See Also
    --------
    biweight_midcovariance, biweight_midcorrelation, astropy.stats.mad_std, astropy.stats.median_absolute_deviation

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Robust_measures_of_scale#The_biweight_midvariance

    .. [2] Beers, Flynn, and Gebhardt (1990; AJ 100, 32) (http://adsabs.harvard.edu/abs/1990AJ....100...32B)

    Examples
    --------
    Generate random variates from a Gaussian distribution and return the
    biweight midvariance of the distribution:

    >>> import numpy as np
    >>> from astropy.stats import biweight_midvariance
    >>> rand = np.random.RandomState(12345)
    >>> bivar = biweight_midvariance(rand.randn(1000))
    >>> print(bivar)    # doctest: +FLOAT_CMP
    0.97362869104
    """

    data = np.asanyarray(data).astype(np.float64)

    if M is None:
        M = np.median(data, axis=axis)
    if axis is not None:
        M = np.expand_dims(M, axis=axis)

    # set up the differences
    d = data - M

    # set up the weighting
    mad = median_absolute_deviation(data, axis=axis)
    if axis is not None:
        mad = np.expand_dims(mad, axis=axis)
    u = d / (c * mad)

    # now remove the outlier points
    mask = np.abs(u) < 1
    u = u ** 2

    if modify_sample_size:
        n = mask.sum(axis=axis)
    else:
        if axis is None:
            n = data.size
        else:
            n = data.shape[axis]

    f1 = d * d * (1. - u)**4
    f1[~mask] = 0.
    f1 = f1.sum(axis=axis)
    f2 = (1. - u) * (1. - 5.*u)
    f2[~mask] = 0.
    f2 = np.abs(f2.sum(axis=axis))**2

    return n * f1 / f2


@deprecated_renamed_argument('a', 'data', '2.0')
def biweight_midcovariance(data, c=9.0, M=None, modify_sample_size=False):
    r"""
    Compute the biweight midcovariance between pairs of multiple
    variables.

    The biweight midcovariance is a robust and resistant estimator of
    the covariance between two variables.

    This function computes the biweight midcovariance between all pairs
    of the input variables (rows) in the input data.  The output array
    will have a shape of (N_variables, N_variables).  The diagonal
    elements will be the biweight midvariances of each input variable
    (see :func:`biweight_midvariance`).  The off-diagonal elements will
    be the biweight midcovariances between each pair of input variables.

    For example, if the input array ``data`` contains three variables
    (rows) ``x``, ``y``, and ``z``, the output `~numpy.ndarray`
    midcovariance matrix will be:

    .. math::

         \begin{pmatrix}
         \zeta_{xx}  & \zeta_{xy}  & \zeta_{xz} \\
         \zeta_{yx}  & \zeta_{yy}  & \zeta_{yz} \\
         \zeta_{zx}  & \zeta_{zy}  & \zeta_{zz}
         \end{pmatrix}

    where :math:`\zeta_{xx}`, :math:`\zeta_{yy}`, and :math:`\zeta_{zz}`
    are the biweight midvariances of each variable.  The biweight
    midcovariance between :math:`x` and :math:`y` is :math:`\zeta_{xy}`
    (:math:`= \zeta_{yx}`).  The biweight midcovariance between
    :math:`x` and :math:`z` is :math:`\zeta_{xz}` (:math:`=
    \zeta_{zx}`).  The biweight midcovariance between :math:`y` and
    :math:`z` is :math:`\zeta_{yz}` (:math:`= \zeta_{zy}`).

    The biweight midcovariance between two variables :math:`x` and
    :math:`y` is given by:

    .. math::

        \zeta_{xy} = n \ \frac{\Sigma_{|u_i| < 1, \ |v_i| < 1} \
            (x_i - M_x) (1 - u_i^2)^2 (y_i - M_y) (1 - v_i^2)^2}
            {(\Sigma_{|u_i| < 1} \ (1 - u_i^2) (1 - 5u_i^2))
            (\Sigma_{|v_i| < 1} \ (1 - v_i^2) (1 - 5v_i^2))}

    where :math:`M_x` and :math:`M_y` are the medians (or the input
    locations) of the two variables and :math:`u_i` and :math:`v_i` are
    given by:

    .. math::

        u_{i} = \frac{(x_i - M_x)}{c * MAD_x}

        v_{i} = \frac{(y_i - M_y)}{c * MAD_y}

    where :math:`c` is the biweight tuning constant and :math:`MAD_x`
    and :math:`MAD_y` are the `median absolute deviation
    <https://en.wikipedia.org/wiki/Median_absolute_deviation>`_ of the
    :math:`x` and :math:`y` variables.  The biweight midvariance tuning
    constant ``c`` is typically 9.0 (the default).

    For the standard definition of biweight midcovariance :math:`n` is
    the total number of observations of each variable.  That definition
    is used if ``modify_sample_size`` is `False`, which is the default.

    However, if ``modify_sample_size = True``, then :math:`n` is the
    number of observations for which :math:`|u_i| < 1` and :math:`|v_i|
    < 1`, i.e.

    .. math::

        n = \Sigma_{|u_i| < 1, \ |v_i| < 1} \ 1

    which results in a value closer to the true variance for small
    sample sizes or for a large number of rejected values.

    Parameters
    ----------
    data : 2D or 1D array-like
        Input data either as a 2D or 1D array.  For a 2D array, it
        should have a shape (N_variables, N_observations).  A 1D array
        may be input for observations of a single variable, in which
        case the biweight midvariance will be calculated (no
        covariance).  Each row of ``data`` represents a variable, and
        each column a single observation of all those variables (same as
        the `numpy.cov` convention).

    c : float, optional
        Tuning constant for the biweight estimator (default = 9.0).

    M : float or 1D array-like, optional
        The location estimate of each variable, either as a scalar or
        array.  If ``M`` is an array, then its must be a 1D array
        containing the location estimate of each row (i.e. ``a.ndim``
        elements).  If ``M`` is a scalar value, then its value will be
        used for each variable (row).  If `None` (default), then the
        median of each variable (row) will be used.

    modify_sample_size : bool, optional
        If `False` (default), then the sample size used is the total
        number of observations of each variable, which follows the
        standard definition of biweight midcovariance.  If `True`, then
        the sample size is reduced to correct for any rejected values
        (see formula above), which results in a value closer to the true
        covariance for small sample sizes or for a large number of
        rejected values.

    Returns
    -------
    biweight_midcovariance : `~numpy.ndarray`
        A 2D array representing the biweight midcovariances between each
        pair of the variables (rows) in the input array.  The output
        array will have a shape of (N_variables, N_variables).  The
        diagonal elements will be the biweight midvariances of each
        input variable.  The off-diagonal elements will be the biweight
        midcovariances between each pair of input variables.

    See Also
    --------
    biweight_midvariance, biweight_midcorrelation, biweight_scale, biweight_location

    References
    ----------
    .. [1] http://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidc.htm

    Examples
    --------
    Compute the biweight midcovariance between two random variables:

    >>> import numpy as np
    >>> from astropy.stats import biweight_midcovariance
    >>> # Generate two random variables x and y
    >>> rng = np.random.RandomState(1)
    >>> x = rng.normal(0, 1, 200)
    >>> y = rng.normal(0, 3, 200)
    >>> # Introduce an obvious outlier
    >>> x[0] = 30.0
    >>> # Calculate the biweight midcovariances between x and y
    >>> bicov = biweight_midcovariance([x, y])
    >>> print(bicov)    # doctest: +FLOAT_CMP
    [[ 0.82483155 -0.18961219]
     [-0.18961219 9.80265764]]
    >>> # Print standard deviation estimates
    >>> print(np.sqrt(bicov.diagonal()))    # doctest: +FLOAT_CMP
    [ 0.90820237  3.13091961]
    """

    data = np.asanyarray(data).astype(np.float64)

    # ensure data is 2D
    if data.ndim == 1:
        data = data[np.newaxis, :]
    if data.ndim != 2:
        raise ValueError('The input array must be 2D or 1D.')

    # estimate location if not given
    if M is None:
        M = np.median(data, axis=1)
    M = np.asanyarray(M)
    if M.ndim > 1:
        raise ValueError('M must be a scalar or 1D array.')

    # set up the differences
    d = (data.T - M).T

    # set up the weighting
    mad = median_absolute_deviation(data, axis=1)
    u = (d.T / (c * mad)).T

    # now remove the outlier points
    mask = np.abs(u) < 1
    u = u ** 2

    if modify_sample_size:
        maskf = mask.astype(float)
        n = np.inner(maskf, maskf)
    else:
        n = data[0].size

    usub1 = (1. - u)
    usub5 = (1. - 5. * u)
    usub1[~mask] = 0.

    numerator = d * usub1 ** 2
    denominator = (usub1 * usub5).sum(axis=1)[:, np.newaxis]
    numerator_matrix = np.dot(numerator, numerator.T)
    denominator_matrix = np.dot(denominator, denominator.T)

    return n * (numerator_matrix / denominator_matrix)


def biweight_midcorrelation(x, y, c=9.0, M=None, modify_sample_size=False):
    r"""
    Compute the biweight midcorrelation between two variables.

    The `biweight midcorrelation
    <https://en.wikipedia.org/wiki/Biweight_midcorrelation>`_ is a
    measure of similarity between samples.  It is given by:

    .. math::

        r_{bicorr} = \frac{\zeta_{xy}}{\sqrt{\zeta_{xx} \ \zeta_{yy}}}

    where :math:`\zeta_{xx}` is the biweight midvariance of :math:`x`,
    :math:`\zeta_{yy}` is the biweight midvariance of :math:`y`, and
    :math:`\zeta_{xy}` is the biweight midcovariance of :math:`x` and
    :math:`y`.

    Parameters
    ----------
    x, y : 1D array-like
        Input arrays for the two variables.  ``x`` and ``y`` must be 1D
        arrays and have the same number of elements.
    c : float, optional
        Tuning constant for the biweight estimator (default = 9.0).  See
        `biweight_midcovariance` for more details.
    M : float or array-like, optional
        The location estimate.  If ``M`` is a scalar value, then its
        value will be used for the entire array (or along each ``axis``,
        if specified).  If ``M`` is an array, then its must be an array
        containing the location estimate along each ``axis`` of the
        input array.  If `None` (default), then the median of the input
        array will be used (or along each ``axis``, if specified).  See
        `biweight_midcovariance` for more details.
    modify_sample_size : bool, optional
        If `False` (default), then the sample size used is the total
        number of elements in the array (or along the input ``axis``, if
        specified), which follows the standard definition of biweight
        midcovariance.  If `True`, then the sample size is reduced to
        correct for any rejected values (i.e. the sample size used
        includes only the non-rejected values), which results in a value
        closer to the true midcovariance for small sample sizes or for a
        large number of rejected values.  See `biweight_midcovariance`
        for more details.

    Returns
    -------
    biweight_midcorrelation : float
        The biweight midcorrelation between ``x`` and ``y``.

    See Also
    --------
    biweight_scale, biweight_midvariance, biweight_midcovariance, biweight_location

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Biweight_midcorrelation

    Examples
    --------
    Calculate the biweight midcorrelation between two variables:

    >>> import numpy as np
    >>> from astropy.stats import biweight_midcorrelation
    >>> rng = np.random.RandomState(12345)
    >>> x = rng.normal(0, 1, 200)
    >>> y = rng.normal(0, 3, 200)
    >>> # Introduce an obvious outlier
    >>> x[0] = 30.0
    >>> bicorr = biweight_midcorrelation(x, y)
    >>> print(bicorr)    # doctest: +FLOAT_CMP
    -0.0495780713907
    """

    x = np.asanyarray(x)
    y = np.asanyarray(y)
    if x.ndim != 1:
        raise ValueError('x must be a 1D array.')
    if y.ndim != 1:
        raise ValueError('y must be a 1D array.')
    if x.shape != y.shape:
        raise ValueError('x and y must have the same shape.')

    bicorr = biweight_midcovariance([x, y], c=c, M=M,
                                    modify_sample_size=modify_sample_size)

    return bicorr[0, 1] / (np.sqrt(bicorr[0, 0] * bicorr[1, 1]))
