import pytest
import numpy as np
from numpy.testing import assert_allclose

from .... import units
from ....tests.helper import assert_quantity_allclose
from .. import LombScargle


ALL_METHODS = LombScargle.available_methods
ALL_METHODS_NO_AUTO = [method for method in ALL_METHODS if method != 'auto']
FAST_METHODS = [method for method in ALL_METHODS if 'fast' in method]
NTERMS_METHODS = [method for method in ALL_METHODS if 'chi2' in method]
NORMALIZATIONS = ['standard', 'psd', 'log', 'model']


@pytest.fixture
def data(N=100, period=1, theta=[10, 2, 3], dy=1, rseed=0):
    """Generate some data for testing"""
    rng = np.random.RandomState(rseed)
    t = 20 * period * rng.rand(N)
    omega = 2 * np.pi / period
    y = theta[0] + theta[1] * np.sin(omega * t) + theta[2] * np.cos(omega * t)
    dy = dy * (0.5 + rng.rand(N))
    y += dy * rng.randn(N)

    return t, y, dy


@pytest.mark.parametrize('minimum_frequency', [None, 1.0])
@pytest.mark.parametrize('maximum_frequency', [None, 5.0])
@pytest.mark.parametrize('nyquist_factor', [1, 10])
@pytest.mark.parametrize('samples_per_peak', [1, 5])
def test_autofrequency(data, minimum_frequency, maximum_frequency,
                       nyquist_factor, samples_per_peak):
    t, y, dy = data
    baseline = t.max() - t.min()

    freq = LombScargle(t, y, dy).autofrequency(samples_per_peak,
                                               nyquist_factor,
                                               minimum_frequency,
                                               maximum_frequency)
    df = freq[1] - freq[0]

    # Check sample spacing
    assert_allclose(df, 1. / baseline / samples_per_peak)

    # Check minimum frequency
    if minimum_frequency is None:
        assert_allclose(freq[0], 0.5 * df)
    else:
        assert_allclose(freq[0], minimum_frequency)

    if maximum_frequency is None:
        avg_nyquist = 0.5 * len(t) / baseline
        assert_allclose(freq[-1], avg_nyquist * nyquist_factor, atol=0.5*df)
    else:
        assert_allclose(freq[-1], maximum_frequency, atol=0.5*df)


@pytest.mark.parametrize('method', ALL_METHODS_NO_AUTO)
@pytest.mark.parametrize('center_data', [True, False])
@pytest.mark.parametrize('fit_mean', [True, False])
@pytest.mark.parametrize('with_errors', [True, False])
@pytest.mark.parametrize('with_units', [True, False])
@pytest.mark.parametrize('normalization', NORMALIZATIONS)
def test_all_methods(data, method, center_data, fit_mean,
                     with_errors, with_units, normalization):
    if method == 'scipy' and (fit_mean or with_errors):
        return

    t, y, dy = data
    frequency = 0.8 + 0.01 * np.arange(40)
    if with_units:
        t = t * units.day
        y = y * units.mag
        dy = dy * units.mag
        frequency = frequency / t.unit
    if not with_errors:
        dy = None

    kwds = dict(normalization=normalization)
    ls = LombScargle(t, y, dy, center_data=center_data, fit_mean=fit_mean)
    P_expected = ls.power(frequency, **kwds)

    # don't use the fft approximation here; we'll test this elsewhere
    if method in FAST_METHODS:
        kwds['method_kwds'] = dict(use_fft=False)
    P_method = ls.power(frequency, method=method, **kwds)

    if with_units:
        if normalization == 'psd' and not with_errors:
            assert P_method.unit == y.unit ** 2
        else:
            assert P_method.unit == units.dimensionless_unscaled
    else:
        assert not hasattr(P_method, 'unit')

    assert_quantity_allclose(P_expected, P_method)


@pytest.mark.parametrize('method', ALL_METHODS_NO_AUTO)
@pytest.mark.parametrize('center_data', [True, False])
@pytest.mark.parametrize('fit_mean', [True, False])
@pytest.mark.parametrize('with_errors', [True, False])
@pytest.mark.parametrize('normalization', NORMALIZATIONS)
def test_integer_inputs(data, method, center_data, fit_mean, with_errors,
                        normalization):
    if method == 'scipy' and (fit_mean or with_errors):
        return

    t, y, dy = data

    t = np.floor(100 * t)
    t_int = t.astype(int)

    y = np.floor(100 * y)
    y_int = y.astype(int)

    dy = np.floor(100 * dy)
    dy_int = dy.astype('int32')

    frequency = 1E-2 * (0.8 + 0.01 * np.arange(40))

    if not with_errors:
        dy = None
        dy_int = None

    kwds = dict(center_data=center_data,
                fit_mean=fit_mean)
    P_float = LombScargle(t, y, dy, **kwds).power(frequency,
                                                  method=method,
                                                  normalization=normalization)
    P_int = LombScargle(t_int, y_int, dy_int,
                        **kwds).power(frequency,
                                      method=method,
                                      normalization=normalization)
    assert_allclose(P_float, P_int)


@pytest.mark.parametrize('method', NTERMS_METHODS)
@pytest.mark.parametrize('center_data', [True, False])
@pytest.mark.parametrize('fit_mean', [True, False])
@pytest.mark.parametrize('with_errors', [True, False])
@pytest.mark.parametrize('nterms', [0, 2, 4])
@pytest.mark.parametrize('normalization', NORMALIZATIONS)
def test_nterms_methods(method, center_data, fit_mean, with_errors,
                        nterms, normalization, data):
    t, y, dy = data
    frequency = 0.8 + 0.01 * np.arange(40)
    if not with_errors:
        dy = None

    ls = LombScargle(t, y, dy, center_data=center_data,
                     fit_mean=fit_mean, nterms=nterms)

    kwds = dict(normalization=normalization)

    if nterms == 0 and not fit_mean:
        with pytest.raises(ValueError) as err:
            ls.power(frequency, method=method, **kwds)
        assert 'nterms' in str(err.value) and 'bias' in str(err.value)
    else:
        P_expected = ls.power(frequency, **kwds)

        # don't use fast fft approximations here
        if 'fast' in method:
            kwds['method_kwds'] = dict(use_fft=False)
        P_method = ls.power(frequency, method=method, **kwds)

        assert_allclose(P_expected, P_method, rtol=1E-7, atol=1E-25)


@pytest.mark.parametrize('method', FAST_METHODS)
@pytest.mark.parametrize('center_data', [True, False])
@pytest.mark.parametrize('fit_mean', [True, False])
@pytest.mark.parametrize('with_errors', [True, False])
@pytest.mark.parametrize('nterms', [0, 1, 2])
def test_fast_approximations(method, center_data, fit_mean,
                             with_errors, nterms, data):
    t, y, dy = data
    frequency = 0.8 + 0.01 * np.arange(40)
    if not with_errors:
        dy = None

    ls = LombScargle(t, y, dy, center_data=center_data,
                     fit_mean=fit_mean, nterms=nterms)

    # use only standard normalization because we compare via absolute tolerance
    kwds = dict(method=method, normalization='standard')

    if method == 'fast' and nterms != 1:
        with pytest.raises(ValueError) as err:
            ls.power(frequency, **kwds)
        assert 'nterms' in str(err.value)

    elif nterms == 0 and not fit_mean:
        with pytest.raises(ValueError) as err:
            ls.power(frequency, **kwds)
        assert 'nterms' in str(err.value) and 'bias' in str(err.value)

    else:
        P_fast = ls.power(frequency, **kwds)
        kwds['method_kwds'] = dict(use_fft=False)
        P_slow = ls.power(frequency, **kwds)

        assert_allclose(P_fast, P_slow, atol=0.008)


@pytest.mark.parametrize('method', LombScargle.available_methods)
@pytest.mark.parametrize('shape', [(), (1,), (2,), (3,), (2, 3)])
def test_output_shapes(method, shape, data):
    t, y, dy = data
    freq = np.asarray(np.zeros(shape))
    freq.flat = np.arange(1, freq.size + 1)
    PLS = LombScargle(t, y, fit_mean=False).power(freq, method=method)
    assert PLS.shape == shape


@pytest.mark.parametrize('method', LombScargle.available_methods)
def test_errors_on_unit_mismatch(method, data):
    t, y, dy = data

    t = t * units.second
    y = y * units.mag
    frequency = np.linspace(0.5, 1.5, 10)

    # this should fail because frequency and 1/t units do not match
    with pytest.raises(ValueError) as err:
        LombScargle(t, y, fit_mean=False).power(frequency, method=method)
    assert str(err.value).startswith('Units of frequency not equivalent')

    # this should fail because dy and y units do not match
    with pytest.raises(ValueError) as err:
        LombScargle(t, y, dy, fit_mean=False).power(frequency / t.unit)
    assert str(err.value).startswith('Units of dy not equivalent')


# we don't test all normalizations here because they are tested above
# only test method='auto' because unit handling does not depend on method
@pytest.mark.parametrize('fit_mean', [True, False])
@pytest.mark.parametrize('center_data', [True, False])
@pytest.mark.parametrize('normalization', ['standard', 'psd'])
@pytest.mark.parametrize('with_error', [True, False])
def test_unit_conversions(data, fit_mean, center_data,
                          normalization, with_error):
    t, y, dy = data

    t_day = t * units.day
    t_hour = units.Quantity(t_day, 'hour')

    y_meter = y * units.meter
    y_millimeter = units.Quantity(y_meter, 'millimeter')

    # sanity check on inputs
    assert_quantity_allclose(t_day, t_hour)
    assert_quantity_allclose(y_meter, y_millimeter)

    if with_error:
        dy = dy * units.meter
    else:
        dy = None

    freq_day, P1 = LombScargle(t_day, y_meter, dy).autopower()
    freq_hour, P2 = LombScargle(t_hour, y_millimeter, dy).autopower()

    # Check units of frequency
    assert freq_day.unit == 1. / units.day
    assert freq_hour.unit == 1. / units.hour

    # Check that results match
    assert_quantity_allclose(freq_day, freq_hour)
    assert_quantity_allclose(P1, P2)

    # Check that switching frequency units doesn't change things
    P3 = LombScargle(t_day, y_meter, dy).power(freq_hour)
    P4 = LombScargle(t_hour, y_meter, dy).power(freq_day)
    assert_quantity_allclose(P3, P4)


@pytest.mark.parametrize('fit_mean', [True, False])
@pytest.mark.parametrize('with_units', [True, False])
@pytest.mark.parametrize('freq', [1.0, 2.0])
def test_model(fit_mean, with_units, freq):
    rand = np.random.RandomState(0)
    t = 10 * rand.rand(40)
    params = 10 * rand.rand(3)

    y = np.zeros_like(t)
    if fit_mean:
        y += params[0]
    y += params[1] * np.sin(2 * np.pi * freq * (t - params[2]))

    if with_units:
        t = t * units.day
        y = y * units.mag
        freq = freq / units.day

    ls = LombScargle(t, y, center_data=False, fit_mean=fit_mean)
    y_fit = ls.model(t, freq)
    assert_quantity_allclose(y_fit, y)


@pytest.mark.parametrize('t_unit', [units.second, units.day])
@pytest.mark.parametrize('frequency_unit', [units.Hz, 1. / units.second])
@pytest.mark.parametrize('y_unit', [units.mag, units.jansky])
def test_model_units_match(data, t_unit, frequency_unit, y_unit):
    t, y, dy = data
    t_fit = t[:5]
    frequency = 1.0

    t = t * t_unit
    t_fit = t_fit * t_unit
    y = y * y_unit
    dy = dy * y_unit
    frequency = frequency * frequency_unit

    ls = LombScargle(t, y, dy)
    y_fit = ls.model(t_fit, frequency)
    assert y_fit.unit == y_unit


def test_model_units_mismatch(data):
    t, y, dy = data
    frequency = 1.0
    t_fit = t[:5]

    t = t * units.second
    t_fit = t_fit * units.second
    y = y * units.mag
    frequency = 1.0 / t.unit

    # this should fail because frequency and 1/t units do not match
    with pytest.raises(ValueError) as err:
        LombScargle(t, y).model(t_fit, frequency=1.0)
    assert str(err.value).startswith('Units of frequency not equivalent')

    # this should fail because t and t_fit units do not match
    with pytest.raises(ValueError) as err:
        LombScargle(t, y).model([1, 2], frequency)
    assert str(err.value).startswith('Units of t not equivalent')

    # this should fail because dy and y units do not match
    with pytest.raises(ValueError) as err:
        LombScargle(t, y, dy).model(t_fit, frequency)
    assert str(err.value).startswith('Units of dy not equivalent')


def test_autopower(data):
    t, y, dy = data
    ls = LombScargle(t, y, dy)
    kwargs = dict(samples_per_peak=6, nyquist_factor=2,
                  minimum_frequency=2, maximum_frequency=None)
    freq1 = ls.autofrequency(**kwargs)
    power1 = ls.power(freq1)
    freq2, power2 = ls.autopower(**kwargs)

    assert_allclose(freq1, freq2)
    assert_allclose(power1, power2)


@pytest.fixture
def null_data(N=1000, dy=1, rseed=0):
    """Generate null hypothesis data"""
    rng = np.random.RandomState(rseed)
    t = 100 * rng.rand(N)
    dy = 0.5 * dy * (1 + rng.rand(N))
    y = dy * rng.randn(N)
    return t, y, dy


@pytest.mark.parametrize('normalization', NORMALIZATIONS)
def test_distribution(null_data, normalization):
    t, y, dy = null_data
    N = len(t)
    ls = LombScargle(t, y, dy)
    freq, power = ls.autopower(normalization=normalization,
                               maximum_frequency=40)
    z = np.linspace(0, power.max(), 1000)

    # Test that pdf and cdf are consistent
    dz = z[1] - z[0]
    z_mid = z[:-1] + 0.5 * dz
    pdf = _lombscargle_pdf(z_mid, N, normalization=normalization)
    cdf = _lombscargle_cdf(z, N, normalization=normalization)
    assert_allclose(pdf, np.diff(cdf) / dz, rtol=1E-5, atol=1E-8)

    # Test that observed power is distributed according to the theoretical pdf
    hist, bins = np.histogram(power, 30, normed=True)
    midpoints = 0.5 * (bins[1:] + bins[:-1])
    pdf = _lombscargle_pdf(midpoints, N, normalization=normalization)
    assert_allclose(hist, pdf, rtol=0.05, atol=0.05 * pdf[0])


# The following are convenience functions used to compute statistics of the
# periodogram under various normalizations; they are used in the preceding
# test.


def _lombscargle_pdf(z, N, normalization, dH=1, dK=3):
    """Probability density function for Lomb-Scargle periodogram

    Compute the expected probability density function of the periodogram
    for the null hypothesis - i.e. data consisting of Gaussian noise.

    Parameters
    ----------
    z : array-like
        the periodogram value
    N : int
        the number of data points from which the periodogram was computed
    normalization : string
        The periodogram normalization. Must be one of
        ['standard', 'model', 'log', 'psd']
    dH, dK : integers (optional)
        The number of parameters in the null hypothesis and the model

    Returns
    -------
    pdf : np.ndarray
        The expected probability density function

    Notes
    -----
    For normalization='psd', the distribution can only be computed for
    periodograms constructed with errors specified.
    All expressions used here are adapted from Table 1 of Baluev 2008 [1]_.

    References
    ----------
    .. [1] Baluev, R.V. MNRAS 385, 1279 (2008)
    """
    if dK - dH != 2:
        raise NotImplementedError("Degrees of freedom != 2")
    Nk = N - dK

    if normalization == 'psd':
        return np.exp(-z)
    elif normalization == 'standard':
        return 0.5 * Nk * (1 + z) ** (-0.5 * Nk - 1)
    elif normalization == 'model':
        return 0.5 * Nk * (1 - z) ** (0.5 * Nk - 1)
    elif normalization == 'log':
        return 0.5 * Nk * np.exp(-0.5 * Nk * z)
    else:
        raise ValueError("normalization='{0}' is not recognized"
                         "".format(normalization))


def _lombscargle_cdf(z, N, normalization, dH=1, dK=3):
    """Cumulative distribution for the Lomb-Scargle periodogram

    Compute the expected cumulative distribution of the periodogram
    for the null hypothesis - i.e. data consisting of Gaussian noise.

    Parameters
    ----------
    z : array-like
        the periodogram value
    N : int
        the number of data points from which the periodogram was computed
    normalization : string
        The periodogram normalization. Must be one of
        ['standard', 'model', 'log', 'psd']
    dH, dK : integers (optional)
        The number of parameters in the null hypothesis and the model

    Returns
    -------
    cdf : np.ndarray
        The expected cumulative distribution function

    Notes
    -----
    For normalization='psd', the distribution can only be computed for
    periodograms constructed with errors specified.
    All expressions used here are adapted from Table 1 of Baluev 2008 [1]_.

    References
    ----------
    .. [1] Baluev, R.V. MNRAS 385, 1279 (2008)
    """
    if dK - dH != 2:
        raise NotImplementedError("Degrees of freedom != 2")
    Nk = N - dK

    if normalization == 'psd':
        return 1 - np.exp(-z)
    elif normalization == 'standard':
        return 1 - (1 + z) ** (-0.5 * Nk)
    elif normalization == 'model':
        return 1 - (1 - z) ** (0.5 * Nk)
    elif normalization == 'log':
        return 1 - np.exp(-0.5 * Nk * z)
    else:
        raise ValueError("normalization='{0}' is not recognized"
                         "".format(normalization))
