# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains simple functions for model selection.
"""

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
       <https://www.originlab.com/doc/Origin-Help/PostFit-CompareFitFunc>
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
    >>> # Set up fake data
    >>> x = np.array(
    ...     [-5.        , -4.94974874, -4.89949749, -4.84924623, -4.79899497,
    ...      -4.74874372, -4.69849246, -4.64824121, -4.59798995, -4.54773869,
    ...      -4.49748744, -4.44723618, -4.39698492, -4.34673367, -4.29648241,
    ...      -4.24623116, -4.1959799 , -4.14572864, -4.09547739, -4.04522613,
    ...      -3.99497487, -3.94472362, -3.89447236, -3.84422111, -3.79396985,
    ...      -3.74371859, -3.69346734, -3.64321608, -3.59296482, -3.54271357,
    ...      -3.49246231, -3.44221106, -3.3919598 , -3.34170854, -3.29145729,
    ...      -3.24120603, -3.19095477, -3.14070352, -3.09045226, -3.04020101,
    ...      -2.98994975, -2.93969849, -2.88944724, -2.83919598, -2.78894472,
    ...      -2.73869347, -2.68844221, -2.63819095, -2.5879397 , -2.53768844,
    ...      -2.48743719, -2.43718593, -2.38693467, -2.33668342, -2.28643216,
    ...      -2.2361809 , -2.18592965, -2.13567839, -2.08542714, -2.03517588,
    ...      -1.98492462, -1.93467337, -1.88442211, -1.83417085, -1.7839196 ,
    ...      -1.73366834, -1.68341709, -1.63316583, -1.58291457, -1.53266332,
    ...      -1.48241206, -1.4321608 , -1.38190955, -1.33165829, -1.28140704,
    ...      -1.23115578, -1.18090452, -1.13065327, -1.08040201, -1.03015075,
    ...      -0.9798995 , -0.92964824, -0.87939698, -0.82914573, -0.77889447,
    ...      -0.72864322, -0.67839196, -0.6281407 , -0.57788945, -0.52763819,
    ...      -0.47738693, -0.42713568, -0.37688442, -0.32663317, -0.27638191,
    ...      -0.22613065, -0.1758794 , -0.12562814, -0.07537688, -0.02512563,
    ...       0.02512563,  0.07537688,  0.12562814,  0.1758794 ,  0.22613065,
    ...       0.27638191,  0.32663317,  0.37688442,  0.42713568,  0.47738693,
    ...       0.52763819,  0.57788945,  0.6281407 ,  0.67839196,  0.72864322,
    ...       0.77889447,  0.82914573,  0.87939698,  0.92964824,  0.9798995 ,
    ...       1.03015075,  1.08040201,  1.13065327,  1.18090452,  1.23115578,
    ...       1.28140704,  1.33165829,  1.38190955,  1.4321608 ,  1.48241206,
    ...       1.53266332,  1.58291457,  1.63316583,  1.68341709,  1.73366834,
    ...       1.7839196 ,  1.83417085,  1.88442211,  1.93467337,  1.98492462,
    ...       2.03517588,  2.08542714,  2.13567839,  2.18592965,  2.2361809 ,
    ...       2.28643216,  2.33668342,  2.38693467,  2.43718593,  2.48743719,
    ...       2.53768844,  2.5879397 ,  2.63819095,  2.68844221,  2.73869347,
    ...       2.78894472,  2.83919598,  2.88944724,  2.93969849,  2.98994975,
    ...       3.04020101,  3.09045226,  3.14070352,  3.19095477,  3.24120603,
    ...       3.29145729,  3.34170854,  3.3919598 ,  3.44221106,  3.49246231,
    ...       3.54271357,  3.59296482,  3.64321608,  3.69346734,  3.74371859,
    ...       3.79396985,  3.84422111,  3.89447236,  3.94472362,  3.99497487,
    ...       4.04522613,  4.09547739,  4.14572864,  4.1959799 ,  4.24623116,
    ...       4.29648241,  4.34673367,  4.39698492,  4.44723618,  4.49748744,
    ...       4.54773869,  4.59798995,  4.64824121,  4.69849246,  4.74874372,
    ...       4.79899497,  4.84924623,  4.89949749,  4.94974874,  5.        ])
    >>> y = np.array(
    ...     [ 0.35281047,  0.08003144,  0.1957476 ,  0.44817864,  0.3735116 ,
    ...      -0.19545558,  0.19001768, -0.03027144, -0.02064377,  0.0821197 ,
    ...       0.02880871,  0.2908547 ,  0.15220755,  0.024335  ,  0.08877265,
    ...       0.06673487,  0.29881581, -0.04103165,  0.06261354, -0.17081915,
    ...      -0.51059796,  0.13072372,  0.17288724, -0.148433  ,  0.45395093,
    ...      -0.29087313,  0.00915171, -0.03743675,  0.30655587,  0.29387179,
    ...       0.03098953,  0.07563257, -0.17755705, -0.39615915, -0.06958222,
    ...       0.0312701 ,  0.24605857,  0.24047658, -0.0774645 , -0.06045933,
    ...      -0.20970888, -0.2840012 , -0.34125071,  0.39015969, -0.10192406,
    ...      -0.08760609, -0.25054706,  0.15551447, -0.32275727, -0.04251785,
    ...      -0.17905257,  0.07743525, -0.10208775, -0.23602875, -0.00550674,
    ...       0.08583793,  0.01352946,  0.06079097, -0.12647676, -0.07204353,
    ...      -0.13383759, -0.07106523, -0.16154155, -0.34386259,  0.03726457,
    ...      -0.07809382, -0.32317448,  0.0961708 , -0.17691823,  0.01607295,
    ...       0.15290377,  0.03459497,  0.23876223, -0.23355889,  0.09691925,
    ...      -0.11685449, -0.1496793 , -0.08608369, -0.02645275,  0.05437479,
    ...      -0.18132905,  0.24187956,  0.16650936, -0.22034883,  0.40016036,
    ...       0.49962565,  0.37672321,  0.12834784, -0.02333372,  0.43158683,
    ...       0.17361443,  0.53624   ,  0.37512015,  0.57497049,  0.50178606,
    ...       0.62759138,  0.54920088,  0.97028319,  0.70975807,  0.84131575,
    ...       1.21931942,  0.66002139,  0.76728021,  1.31170785,  0.98394576,
    ...       1.71188479,  1.34834625,  1.39219417,  2.03889851,  2.06427535,
    ...       2.25594134,  2.1773871 ,  1.93621948,  2.60031042,  2.2710606 ,
    ...       2.58702332,  2.71234334,  2.58174921,  2.81797085,  2.95365119,
    ...       2.90938034,  2.66919909,  2.99318021,  3.23221756,  2.8499987 ,
    ...       2.96926297,  2.91062119,  3.35416926,  3.09379993,  3.00451097,
    ...       2.72178974,  2.92599951,  2.61593945,  2.68086431,  2.4628984 ,
    ...       2.63371179,  2.51585306,  2.25573778,  2.26923466,  1.86085328,
    ...       1.66845337,  1.94061247,  1.77182972,  1.75188237,  1.98933131,
    ...       1.59162293,  1.11306047,  1.4153912 ,  0.82914102,  0.90473171,
    ...       0.89285259,  1.16359943,  0.59155982,  0.50005075,  0.57575265,
    ...       0.39809643,  0.69662455,  0.20083703,  0.13770019,  0.23463588,
    ...       0.1819994 ,  0.63106354,  0.40246973,  0.2011258 , -0.08711821,
    ...       0.30424189, -0.084497  , -0.21071687,  0.32079853,  0.13356295,
    ...       0.24313206,  0.11308869,  0.2124982 , -0.09605276, -0.1786032 ,
    ...       0.15958697, -0.14158957, -0.12230561, -0.07840316,  0.01379677,
    ...      -0.06247881, -0.26829675, -0.12335989, -0.44039932,  0.12845012,
    ...      -0.317716  , -0.21875044,  0.01210356, -0.14660524,  0.30962202,
    ...      -0.2577801 ,  0.05402214, -0.00738514, -0.23325698,  0.10493179,
    ...      -0.0340988 ,  0.1545177 ,  0.16482137,  0.43273787,  0.26737354])
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
       <https://www.originlab.com/doc/Origin-Help/PostFit-CompareFitFunc>
    .. [3] Astropy Models and Fitting
        <https://docs.astropy.org/en/stable/modeling>
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
       <https://corpus.ulaval.ca/jspui/handle/20.500.11794/17461>
    .. [3] Wikipedia. Akaike Information Criterion.
       <https://en.wikipedia.org/wiki/Akaike_information_criterion>
    .. [4] Origin Lab. Comparing Two Fitting Functions.
       <https://www.originlab.com/doc/Origin-Help/PostFit-CompareFitFunc>
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
    ...            + models.Gaussian1D(2.4, .4, 0.1))
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
    -634.5257517810961
    >>> akaike_info_criterion_lsq(ssr_g2, 6, x.shape[0]) # doctest: +FLOAT_CMP
    -662.83834510232043
    >>> akaike_info_criterion_lsq(ssr_g1, 3, x.shape[0]) # doctest: +FLOAT_CMP
    -647.47312032659499

    Hence, from the AIC values, we would prefer to choose the model g2_fit.
    However, we can considerably support the model g3_fit, since the
    difference in AIC is about 2.4. We should reject the model g1_fit.

    References
    ----------
    .. [1] Akaike Information Criterion.
       <https://en.wikipedia.org/wiki/Akaike_information_criterion>
    .. [2] Origin Lab. Comparing Two Fitting Functions.
       <https://www.originlab.com/doc/Origin-Help/PostFit-CompareFitFunc>
    """

    return akaike_info_criterion(-0.5 * n_samples * np.log(ssr / n_samples),
                                 n_params, n_samples)
