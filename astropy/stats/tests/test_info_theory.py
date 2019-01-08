
from numpy.testing import assert_allclose

from astropy.stats.info_theory import bayesian_info_criterion, bayesian_info_criterion_lsq
from astropy.stats.info_theory import akaike_info_criterion, akaike_info_criterion_lsq


def test_bayesian_info_criterion():
    # This test is from an example presented in Ref [1]
    lnL = (-176.4, -173.0)
    n_params = (2, 3)
    n_samples = 100
    answer = 2.195
    bic_g = bayesian_info_criterion(lnL[0], n_params[0], n_samples)
    bic_t = bayesian_info_criterion(lnL[1], n_params[1], n_samples)
    assert_allclose(answer, bic_g - bic_t, atol=1e-1)


def test_akaike_info_criterion():
    # This test is from an example presented in Ref [2]
    n_samples = 121
    lnL = (-3.54, -4.17)
    n_params = (6, 5)
    answer = 0.95
    aic_1 = akaike_info_criterion(lnL[0], n_params[0], n_samples)
    aic_2 = akaike_info_criterion(lnL[1], n_params[1], n_samples)
    assert_allclose(answer, aic_1 - aic_2, atol=1e-2)


def test_akaike_info_criterion_lsq():
    # This test is from an example presented in Ref [1]
    n_samples = 100
    n_params = (4, 3, 3)
    ssr = (25.0, 26.0, 27.0)
    answer = (-130.21, -128.46, -124.68)

    assert_allclose(answer[0],
                    akaike_info_criterion_lsq(ssr[0], n_params[0], n_samples),
                    atol=1e-2)
    assert_allclose(answer[1],
                    akaike_info_criterion_lsq(ssr[1], n_params[1], n_samples),
                    atol=1e-2)
    assert_allclose(answer[2],
                    akaike_info_criterion_lsq(ssr[2], n_params[2], n_samples),
                    atol=1e-2)


def test_bayesian_info_criterion_lsq():
    """This test is from:
    http://www.statoek.wiso.uni-goettingen.de/veranstaltungen/non_semi_models/
    AkaikeLsg.pdf
    Note that in there, they compute a "normalized BIC". Therefore, the
    answers presented here are recalculated versions based on their values.
    """

    n_samples = 25
    n_params = (1, 2, 1)
    ssr = (48959, 32512, 37980)
    answer = (192.706, 185.706, 186.360)

    assert_allclose(answer[0], bayesian_info_criterion_lsq(ssr[0],
                                                           n_params[0],
                                                           n_samples),
                    atol=1e-2)
    assert_allclose(answer[1], bayesian_info_criterion_lsq(ssr[1],
                                                           n_params[1],
                                                           n_samples),
                    atol=1e-2)
    assert_allclose(answer[2], bayesian_info_criterion_lsq(ssr[2],
                                                           n_params[2],
                                                           n_samples),
                    atol=1e-2)
