# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains simple functions for model selection.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np

__all__ = ['bayesian_info_criterion', 'bayesian_info_criterion_lsq',
           'akaike_info_criterion', 'akaike_info_criterion_lsq']

__doctest_requires__ = {'bayesian_info_criterion_lsq': ['scipy'],
                        'akaike_info_criterion_lsq': ['scipy']}


def bayesian_info_criterion(log_likelihood, n_params, n_samples):
    r""" Computes the Bayesian Information Criterion (BIC) given the log of the
    likelihood function evaluated at the estimated (or analytically derived)
    parameters, the number of parameters, and the number of samples.

    The BIC is usually applied to decide whether increasing the number of free
    parameters (hence, increasing the model complexity) yields significantly
    better fittings. The decision is in favor of the model with the lowest
    BIC.

    BIC is given as

    .. math::

        \mathrm{BIC} = k \ln(n) - 2L,

    in which :math:`n` is the sample size, :math:`k` is the number of free
    parameters, and :math:`L` is the log likelihood function of the model
    evaluated at the maximum likelihood estimate (i. e., the parameters for
    which L is maximized).

    When comparing two models define
    :math:`\Delta \mathrm{BIC} = \mathrm{BIC}_h - \mathrm{BIC}_l`, in which
    :math:`\mathrm{BIC}_h` is the higher BIC, and :math:`\mathrm{BIC}_l` is
    the lower BIC. The higher is :math:`\Delta \mathrm{BIC}` the stronger is
    the evidence against the model with higher BIC.

    The general rule of thumb is:

    :math:`0 < \Delta\mathrm{BIC} \leq 2`: weak evidence that model low is
    better

    :math:`2 < \Delta\mathrm{BIC} \leq 6`: moderate evidence that model low is
    better

    :math:`6 < \Delta\mathrm{BIC} \leq 10`: strong evidence that model low is
    better

    :math:`\Delta\mathrm{BIC} > 10`: very strong evidence that model low is
    better

    For a detailed explanation, see [1]_ - [5]_.

    Parameters
    ----------
    log_likelihood : float
        Logarithm of the likelihood function of the model evaluated at the
        point of maxima (with respect to the parameter space).
    n_params : int
        Number of free parameters of the model, i.e., dimension of the
        parameter space.
    n_samples : int
        Number of observations.

    Returns
    -------
    bic : float
        Bayesian Information Criterion.

    Examples
    --------
    The following example was originally presented in [1]_. Consider a
    Gaussian model (mu, sigma) and a t-Student model (mu, sigma, delta).
    In addition, assume that the t model has presented a higher likelihood.
    The question that the BIC is proposed to answer is: "Is the increase in
    likelihood due to larger number of parameters?"

    >>> from astropy.stats.info_theory import bayesian_info_criterion
    >>> lnL_g = -176.4
    >>> lnL_t = -173.0
    >>> n_params_g = 2
    >>> n_params_t = 3
    >>> n_samples = 100
    >>> bic_g = bayesian_info_criterion(lnL_g, n_params_g, n_samples)
    >>> bic_t = bayesian_info_criterion(lnL_t, n_params_t, n_samples)
    >>> bic_g - bic_t # doctest: +FLOAT_CMP
    2.1948298140119391

    Therefore, there exist a moderate evidence that the increasing in
    likelihood for t-Student model is due to the larger number of parameters.

    References
    ----------
    .. [1] Richards, D. Maximum Likelihood Estimation and the Bayesian
       Information Criterion.
       <https://hea-www.harvard.edu/astrostat/Stat310_0910/dr_20100323_mle.pdf>
    .. [2] Wikipedia. Bayesian Information Criterion.
       <https://en.wikipedia.org/wiki/Bayesian_information_criterion>
    .. [3] Origin Lab. Comparing Two Fitting Functions.
       <http://www.originlab.com/doc/Origin-Help/PostFit-CompareFitFunc>
    .. [4] Liddle, A. R. Information Criteria for Astrophysical Model
       Selection. 2008. <https://arxiv.org/pdf/astro-ph/0701113v2.pdf>
    .. [5] Liddle, A. R. How many cosmological parameters? 2008.
       <https://arxiv.org/pdf/astro-ph/0401198v3.pdf>
    """

    return n_params*np.log(n_samples) - 2.0*log_likelihood


def bayesian_info_criterion_lsq(ssr, n_params, n_samples):
    r"""
    Computes the Bayesian Information Criterion (BIC) assuming that the
    observations come from a Gaussian distribution.

    In this case, BIC is given as

    .. math::

        \mathrm{BIC} = n\ln\left(\dfrac{\mathrm{SSR}}{n}\right) + k\ln(n)

    in which :math:`n` is the sample size, :math:`k` is the number of free
    parameters and :math:`\mathrm{SSR}` stands for the sum of squared residuals
    between model and data.

    This is applicable, for instance, when the parameters of a model are
    estimated using the least squares statistic. See [1]_ and [2]_.

    Parameters
    ----------
    ssr : float
        Sum of squared residuals (SSR) between model and data.
    n_params : int
        Number of free parameters of the model, i.e., dimension of the
        parameter space.
    n_samples : int
        Number of observations.

    Returns
    -------
    bic : float

    Examples
    --------
    Consider the simple 1-D fitting example presented in the Astropy
    modeling webpage [3]_. There, two models (Box and Gaussian) were fitted to
    a source flux using the least squares statistic. However, the fittings
    themselves do not tell much about which model better represents this
    hypothetical source. Therefore, we are going to apply to BIC in order to
    decide in favor of a model.

    >>> import numpy as np
    >>> from astropy.modeling import models, fitting
    >>> from astropy.stats.info_theory import bayesian_info_criterion_lsq
    >>> # Generate fake data
    >>> np.random.seed(0)
    >>> x = np.linspace(-5., 5., 200)
    >>> y = 3 * np.exp(-0.5 * (x - 1.3)**2 / 0.8**2)
    >>> y += np.random.normal(0., 0.2, x.shape)
    >>> # Fit the data using a Box model.
    >>> # Bounds are not really needed but included here to demonstrate usage.
    >>> t_init = models.Trapezoid1D(amplitude=1., x_0=0., width=1., slope=0.5,
    ...                             bounds={"x_0": (-5., 5.)})
    >>> fit_t = fitting.LevMarLSQFitter()
    >>> t = fit_t(t_init, x, y)
    >>> # Fit the data using a Gaussian
    >>> g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
    >>> fit_g = fitting.LevMarLSQFitter()
    >>> g = fit_g(g_init, x, y)
    >>> # Compute the mean squared errors
    >>> ssr_t = np.sum((t(x) - y)*(t(x) - y))
    >>> ssr_g = np.sum((g(x) - y)*(g(x) - y))
    >>> # Compute the bics
    >>> bic_t = bayesian_info_criterion_lsq(ssr_t, 4, x.shape[0])
    >>> bic_g = bayesian_info_criterion_lsq(ssr_g, 3, x.shape[0])
    >>> bic_t - bic_g # doctest: +FLOAT_CMP
    30.644474706065466

    Hence, there is a very strong evidence that the Gaussian model has a
    significantly better representation of the data than the Box model. This
    is, obviously, expected since the true model is Gaussian.

    References
    ----------
    .. [1] Wikipedia. Bayesian Information Criterion.
       <https://en.wikipedia.org/wiki/Bayesian_information_criterion>
    .. [2] Origin Lab. Comparing Two Fitting Functions.
       <http://www.originlab.com/doc/Origin-Help/PostFit-CompareFitFunc>
    .. [3] Astropy Models and Fitting
        <http://docs.astropy.org/en/stable/modeling>
    """

    return bayesian_info_criterion(-0.5 * n_samples * np.log(ssr / n_samples),
                                   n_params, n_samples)


def akaike_info_criterion(log_likelihood, n_params, n_samples):
    r"""
    Computes the Akaike Information Criterion (AIC).

    Like the Bayesian Information Criterion, the AIC is a measure of
    relative fitting quality which is used for fitting evaluation and model
    selection. The decision is in favor of the model with the lowest AIC.

    AIC is given as

    .. math::

        \mathrm{AIC} = 2(k - L)

    in which :math:`n` is the sample size, :math:`k` is the number of free
    parameters, and :math:`L` is the log likelihood function of the model
    evaluated at the maximum likelihood estimate (i. e., the parameters for
    which L is maximized).

    In case that the sample size is not "large enough" a correction is
    applied, i.e.

    .. math::

        \mathrm{AIC} = 2(k - L) + \dfrac{2k(k+1)}{n - k - 1}

    Rule of thumb [1]_:

    :math:`\Delta\mathrm{AIC}_i = \mathrm{AIC}_i - \mathrm{AIC}_{min}`

    :math:`\Delta\mathrm{AIC}_i < 2`: substantial support for model i

    :math:`3 < \Delta\mathrm{AIC}_i < 7`: considerably less support for model i

    :math:`\Delta\mathrm{AIC}_i > 10`: essentially none support for model i

    in which :math:`\mathrm{AIC}_{min}` stands for the lower AIC among the
    models which are being compared.

    For detailed explanations see [1]_-[6]_.

    Parameters
    ----------
    log_likelihood : float
        Logarithm of the likelihood function of the model evaluated at the
        point of maxima (with respect to the parameter space).
    n_params : int
        Number of free parameters of the model, i.e., dimension of the
        parameter space.
    n_samples : int
        Number of observations.

    Returns
    -------
    aic : float
        Akaike Information Criterion.

    Examples
    --------
    The following example was originally presented in [2]_. Basically, two
    models are being compared. One with six parameters (model 1) and another
    with five parameters (model 2). Despite of the fact that model 2 has a
    lower AIC, we could decide in favor of model 1 since the difference (in
    AIC)  between them is only about 1.0.

    >>> n_samples = 121
    >>> lnL1 = -3.54
    >>> n1_params = 6
    >>> lnL2 = -4.17
    >>> n2_params = 5
    >>> aic1 = akaike_info_criterion(lnL1, n1_params, n_samples)
    >>> aic2 = akaike_info_criterion(lnL2, n2_params, n_samples)
    >>> aic1 - aic2 # doctest: +FLOAT_CMP
    0.9551029748283746

    Therefore, we can strongly support the model 1 with the advantage that
    it has more free parameters.

    References
    ----------
    .. [1] Cavanaugh, J. E.  Model Selection Lecture II: The Akaike
       Information Criterion.
       <http://machinelearning102.pbworks.com/w/file/fetch/47699383/ms_lec_2_ho.pdf>
    .. [2] Mazerolle, M. J. Making sense out of Akaike's Information
       Criterion (AIC): its use and interpretation in model selection and
       inference from ecological data.
       <http://theses.ulaval.ca/archimede/fichiers/21842/apa.html>
    .. [3] Wikipedia. Akaike Information Criterion.
       <https://en.wikipedia.org/wiki/Akaike_information_criterion>
    .. [4] Origin Lab. Comparing Two Fitting Functions.
       <http://www.originlab.com/doc/Origin-Help/PostFit-CompareFitFunc>
    .. [5] Liddle, A. R. Information Criteria for Astrophysical Model
       Selection. 2008. <https://arxiv.org/pdf/astro-ph/0701113v2.pdf>
    .. [6] Liddle, A. R. How many cosmological parameters? 2008.
       <https://arxiv.org/pdf/astro-ph/0401198v3.pdf>
    """
    # Correction in case of small number of observations
    if n_samples/float(n_params) >= 40.0:
        aic = 2.0 * (n_params - log_likelihood)
    else:
        aic = (2.0 * (n_params - log_likelihood) +
               2.0 * n_params * (n_params + 1.0) /
               (n_samples - n_params - 1.0))
    return aic


def akaike_info_criterion_lsq(ssr, n_params, n_samples):
    r"""
    Computes the Akaike Information Criterion assuming that the observations
    are Gaussian distributed.

    In this case, AIC is given as

    .. math::

        \mathrm{AIC} = n\ln\left(\dfrac{\mathrm{SSR}}{n}\right) + 2k

    In case that the sample size is not "large enough", a correction is
    applied, i.e.

    .. math::

        \mathrm{AIC} = n\ln\left(\dfrac{\mathrm{SSR}}{n}\right) + 2k +
                       \dfrac{2k(k+1)}{n-k-1}


    in which :math:`n` is the sample size, :math:`k` is the number of free
    parameters and :math:`\mathrm{SSR}` stands for the sum of squared residuals
    between model and data.

    This is applicable, for instance, when the parameters of a model are
    estimated using the least squares statistic.

    Parameters
    ----------
    ssr : float
        Sum of squared residuals (SSR) between model and data.
    n_params : int
        Number of free parameters of the model, i.e.,  the dimension of the
        parameter space.
    n_samples : int
        Number of observations.

    Returns
    -------
    aic : float
        Akaike Information Criterion.

    Examples
    --------
    This example is based on Astropy Modeling webpage, Compound models
    section.

    >>> import numpy as np
    >>> from astropy.modeling import models, fitting
    >>> from astropy.stats.info_theory import akaike_info_criterion_lsq
    >>> np.random.seed(42)
    >>> # Generate fake data
    >>> g1 = models.Gaussian1D(.1, 0, 0.2) # changed this to noise level
    >>> g2 = models.Gaussian1D(.1, 0.3, 0.2) # and added another Gaussian
    >>> g3 = models.Gaussian1D(2.5, 0.5, 0.1)
    >>> x = np.linspace(-1, 1, 200)
    >>> y = g1(x) + g2(x) + g3(x) + np.random.normal(0., 0.2, x.shape)
    >>> # Fit with three Gaussians
    >>> g3_init = (models.Gaussian1D(.1, 0, 0.1)
    ...            + models.Gaussian1D(.1, 0.2, 0.15)
    ...            + models.Gaussian1D(2., .4, 0.1))
    >>> fitter = fitting.LevMarLSQFitter()
    >>> g3_fit = fitter(g3_init, x, y)
    >>> # Fit with two Gaussians
    >>> g2_init = (models.Gaussian1D(.1, 0, 0.1) +
    ...            models.Gaussian1D(2, 0.5, 0.1))
    >>> g2_fit = fitter(g2_init, x, y)
    >>> # Fit with only one Gaussian
    >>> g1_init = models.Gaussian1D(amplitude=2., mean=0.3, stddev=.5)
    >>> g1_fit = fitter(g1_init, x, y)
    >>> # Compute the mean squared errors
    >>> ssr_g3 = np.sum((g3_fit(x) - y)**2.0)
    >>> ssr_g2 = np.sum((g2_fit(x) - y)**2.0)
    >>> ssr_g1 = np.sum((g1_fit(x) - y)**2.0)
    >>> akaike_info_criterion_lsq(ssr_g3, 9, x.shape[0]) # doctest: +FLOAT_CMP
    -660.41075962620482
    >>> akaike_info_criterion_lsq(ssr_g2, 6, x.shape[0]) # doctest: +FLOAT_CMP
    -662.83834510232043
    >>> akaike_info_criterion_lsq(ssr_g1, 3, x.shape[0]) # doctest: +FLOAT_CMP
    -647.47312032659499

    Hence, from the AIC values, we would prefer to choose the model g2_fit.
    However, we can considerably support the model g3_fit, since the
    difference in AIC is about 2.4. We should reject the model g1_fit.

    References
    ----------
    .. [1] Akaike Information Criteria
       <http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf>
    .. [2] Hu, S. Akaike Information Criterion.
       <http://www4.ncsu.edu/~shu3/Presentation/AIC.pdf>
    .. [3] Origin Lab. Comparing Two Fitting Functions.
       <http://www.originlab.com/doc/Origin-Help/PostFit-CompareFitFunc>
    """

    return akaike_info_criterion(-0.5 * n_samples * np.log(ssr / n_samples),
                                 n_params, n_samples)
