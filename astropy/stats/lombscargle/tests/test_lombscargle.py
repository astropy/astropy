import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal

from astropy import units

from .. import LombScargle
from ..implementations import lombscargle_slow, lombscargle


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


@pytest.mark.parametrize('method', LombScargle.available_methods)
@pytest.mark.parametrize('shape', [(), (1,), (2,), (3,), (2, 3)])
def test_output_shapes(method, shape, data):
    t, y, dy = data
    freq = np.asarray(np.random.rand(*shape))
    freq.flat = np.arange(1, freq.size + 1)
    PLS = lombscargle(t, y, frequency=freq,
                      fit_bias=False, method=method)
    assert_equal(PLS.shape, shape)


@pytest.mark.parametrize('method', LombScargle.available_methods)
@pytest.mark.parametrize('t_unit', [units.second, units.day])
@pytest.mark.parametrize('frequency_unit', [units.Hz, 1. / units.second])
@pytest.mark.parametrize('y_unit', [units.mag, units.jansky])
def test_units_match(method, t_unit, frequency_unit, y_unit, data):
    t, y, dy = data
    dy = dy.mean()  # scipy only supports constant errors

    t = t * t_unit
    y = y * y_unit
    dy = dy * y_unit
    frequency = np.linspace(0.5, 1.5, 10) * frequency_unit
    PLS = lombscargle(t, y, frequency=frequency,
                      fit_bias=False, method=method)
    assert_equal(PLS.unit, units.dimensionless_unscaled)

    PLS = lombscargle(t, y, dy,
                      frequency=frequency,
                      fit_bias=False, method=method)
    assert_equal(PLS.unit, units.dimensionless_unscaled)


@pytest.mark.parametrize('method', LombScargle.available_methods)
def test_units_mismatch(method, data):
    t, y, dy = data

    t = t * units.second
    y = y * units.mag
    frequency = np.linspace(0.5, 1.5, 10)

    # this should fail because frequency and 1/t units do not match
    with pytest.raises(ValueError) as err:
        lombscargle(t, y, frequency=frequency,
                    method=method, fit_bias=False)
    assert str(err.value).startswith('Units of frequency not equivalent')

    # this should fail because dy and y units do not match
    with pytest.raises(ValueError) as err:
        lombscargle(t, y, dy, frequency / t.unit,
                    method=method, fit_bias=False)
    assert str(err.value).startswith('Units of y not equivalent')


@pytest.mark.parametrize('method', LombScargle.available_methods)
@pytest.mark.parametrize('center_data', [True, False])
@pytest.mark.parametrize('freq', [0.8 + 0.01 * np.arange(40)])
def test_common_interface(method, center_data, freq, data):
    t, y, dy = data

    expected_PLS = lombscargle_slow(t, y, dy=None, frequency=freq,
                                    fit_bias=False, center_data=center_data)
    PLS = lombscargle(t, y, frequency=freq, method=method,
                      fit_bias=False, center_data=center_data)

    if method in ['fastchi2', 'fast', 'auto']:
        atol = 0.005
    else:
        atol = 0
    assert_allclose(PLS, expected_PLS, atol=atol)


@pytest.mark.parametrize('method', LombScargle.available_methods)
@pytest.mark.parametrize('center_data', [True, False])
@pytest.mark.parametrize('fit_bias', [True, False])
@pytest.mark.parametrize('freq', [0.8 + 0.01 * np.arange(40)])
def test_object_interface_power(data, method, center_data, fit_bias, freq):
    t, y, dy = data
    if method == 'scipy' and fit_bias:
        return
    if method == 'scipy':
        dy = None
    expected_PLS = lombscargle(t, y, dy,
                               frequency=freq,
                               method=method,
                               fit_bias=fit_bias,
                               center_data=center_data)
    ls = LombScargle(t, y, dy, fit_bias=fit_bias, center_data=center_data)
    PLS = ls.power(freq, method=method)
    assert_allclose(PLS, expected_PLS)


@pytest.mark.parametrize('method', LombScargle.available_methods)
@pytest.mark.parametrize('center_data', [True, False])
@pytest.mark.parametrize('fit_bias', [True, False])
def test_object_interface_autopower(data, method, center_data, fit_bias):
    t, y, dy = data
    if method == 'scipy' and fit_bias:
        return
    if method == 'scipy':
        dy = None

    ls = LombScargle(t, y, dy, fit_bias=fit_bias, center_data=center_data)
    freq, PLS = ls.autopower(method=method)

    expected_PLS = lombscargle(t, y, dy, freq,
                               method=method,
                               fit_bias=fit_bias,
                               center_data=center_data)

    assert_allclose(PLS, expected_PLS)


@pytest.mark.parametrize('fit_bias', [True, False])
@pytest.mark.parametrize('freq', [1.0, 2.0])
def test_object_interface_model(fit_bias, freq):
    rand = np.random.RandomState(0)
    t = 10 * rand.rand(40)
    params = 10 * rand.rand(3)

    y = np.zeros_like(t)
    if fit_bias:
        y += params[0]
    y += params[1] * np.sin(2 * np.pi * freq * (t - params[2]))

    ls = LombScargle(t, y, center_data=False, fit_bias=fit_bias)
    y_fit = ls.model(t, freq)
    assert_allclose(y_fit, y)


def test_object_interface_bad_input(data):
    t, y, dy = data
    ls = LombScargle(t, y, dy)
    # this should fail because frequency and 1/t unitsdo not match
    with pytest.raises(ValueError) as err:
        ls.power(frequency=None)
    assert str(err.value).startswith('Must supply a valid frequency')


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
    assert_equal(y_fit.unit, y_unit)


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
    assert str(err.value).startswith('Units of t_fit are not equivalent')

    # this should fail because dy and y units do not match
    with pytest.raises(ValueError) as err:
        LombScargle(t, y, dy).model(t_fit, frequency)
    assert str(err.value).startswith('Units of y not equivalent')
