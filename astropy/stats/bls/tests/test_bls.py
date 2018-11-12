# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
from numpy.testing import assert_allclose

from .... import units
from ....tests.helper import assert_quantity_allclose
from .. import BoxLeastSquares
from ...lombscargle.core import has_units


def assert_allclose_blsresults(blsresult, other, **kwargs):
    """Assert that another BoxLeastSquaresResults object is consistent

    This method loops over all attributes and compares the values using
    :func:`~astropy.tests.helper.assert_quantity_allclose` function.

    Parameters
    ----------
    other : BoxLeastSquaresResults
        The other results object to compare.

    """
    for k, v in blsresult.items():
        if k not in other:
            raise AssertionError("missing key '{0}'".format(k))
        if k == "objective":
            assert v == other[k], (
                "Mismatched objectives. Expected '{0}', got '{1}'"
                .format(v, other[k])
            )
            continue
        assert_quantity_allclose(v, other[k], **kwargs)


@pytest.fixture
def data():
    rand = np.random.RandomState(123)
    t = rand.uniform(0, 10, 500)
    y = np.ones_like(t)
    dy = rand.uniform(0.005, 0.01, len(t))
    period = 2.0
    transit_time = 0.5
    duration = 0.16
    depth = 0.2
    m = np.abs((t-transit_time+0.5*period) % period-0.5*period) < 0.5*duration
    y[m] = 1.0 - depth
    y += dy * rand.randn(len(t))
    return t, y, dy, dict(period=period, transit_time=transit_time,
                          duration=duration, depth=depth)


def test_32bit_bug():
    rand = np.random.RandomState(42)
    t = rand.uniform(0, 10, 500)
    y = np.ones_like(t)
    y[np.abs((t + 1.0) % 2.0-1) < 0.08] = 1.0 - 0.1
    y += 0.01 * rand.randn(len(t))

    model = BoxLeastSquares(t, y)
    results = model.autopower(0.16)
    assert np.allclose(results.period[np.argmax(results.power)],
                       1.9923406038842544)
    periods = np.linspace(1.9, 2.1, 5)
    results = model.power(periods, 0.16)
    assert np.allclose(
        results.power,
        np.array([0.01421067, 0.02842475, 0.10867671, 0.05117755, 0.01783253])
    )


@pytest.mark.parametrize("objective", ["likelihood", "snr"])
def test_correct_model(data, objective):
    t, y, dy, params = data
    model = BoxLeastSquares(t, y, dy)
    periods = np.exp(np.linspace(np.log(params["period"]) - 0.1,
                                 np.log(params["period"]) + 0.1, 1000))
    results = model.power(periods, params["duration"], objective=objective)
    ind = np.argmax(results.power)
    for k, v in params.items():
        assert_allclose(results[k][ind], v, atol=0.01)
    chi = (results.depth[ind]-params["depth"]) / results.depth_err[ind]
    assert np.abs(chi) < 1


@pytest.mark.parametrize("objective", ["likelihood", "snr"])
@pytest.mark.parametrize("offset", [False, True])
def test_fast_method(data, objective, offset):
    t, y, dy, params = data
    if offset:
        t = t - params["transit_time"] + params["period"]
    model = BoxLeastSquares(t, y, dy)
    periods = np.exp(np.linspace(np.log(params["period"]) - 1,
                                 np.log(params["period"]) + 1, 10))
    durations = params["duration"]
    results = model.power(periods, durations, objective=objective)

    assert_allclose_blsresults(results, model.power(periods, durations,
                                                    method="slow",
                                                    objective=objective))


def test_input_units(data):
    t, y, dy, params = data

    t_unit = units.day
    y_unit = units.mag

    with pytest.raises(units.UnitConversionError):
        BoxLeastSquares(t * t_unit, y * y_unit, dy * units.one)
    with pytest.raises(units.UnitConversionError):
        BoxLeastSquares(t * t_unit, y * units.one, dy * y_unit)
    with pytest.raises(units.UnitConversionError):
        BoxLeastSquares(t * t_unit, y, dy * y_unit)
    model = BoxLeastSquares(t*t_unit, y * units.one, dy)
    assert model.dy.unit == model.y.unit
    model = BoxLeastSquares(t*t_unit, y * y_unit, dy)
    assert model.dy.unit == model.y.unit
    model = BoxLeastSquares(t*t_unit, y*y_unit)
    assert model.dy is None


def test_period_units(data):
    t, y, dy, params = data
    t_unit = units.day
    y_unit = units.mag
    model = BoxLeastSquares(t * t_unit, y * y_unit, dy)

    p = model.autoperiod(params["duration"])
    assert p.unit == t_unit
    p = model.autoperiod(params["duration"] * 24 * units.hour)
    assert p.unit == t_unit
    with pytest.raises(units.UnitConversionError):
        model.autoperiod(params["duration"] * units.mag)

    p = model.autoperiod(params["duration"], minimum_period=0.5)
    assert p.unit == t_unit
    with pytest.raises(units.UnitConversionError):
        p = model.autoperiod(params["duration"], minimum_period=0.5*units.mag)

    p = model.autoperiod(params["duration"], maximum_period=0.5)
    assert p.unit == t_unit
    with pytest.raises(units.UnitConversionError):
        p = model.autoperiod(params["duration"], maximum_period=0.5*units.mag)

    p = model.autoperiod(params["duration"], minimum_period=0.5,
                         maximum_period=1.5)
    p2 = model.autoperiod(params["duration"], maximum_period=0.5,
                          minimum_period=1.5)
    assert_quantity_allclose(p, p2)


@pytest.mark.parametrize("method", ["fast", "slow"])
@pytest.mark.parametrize("with_err", [True, False])
@pytest.mark.parametrize("t_unit", [None, units.day])
@pytest.mark.parametrize("y_unit", [None, units.mag])
@pytest.mark.parametrize("objective", ["likelihood", "snr"])
def test_results_units(data, method, with_err, t_unit, y_unit, objective):
    t, y, dy, params = data

    periods = np.linspace(params["period"]-1.0, params["period"]+1.0, 3)

    if t_unit is not None:
        t = t * t_unit
    if y_unit is not None:
        y = y * y_unit
        dy = dy * y_unit
    if not with_err:
        dy = None

    model = BoxLeastSquares(t, y, dy)
    results = model.power(periods, params["duration"], method=method,
                          objective=objective)

    if t_unit is None:
        assert not has_units(results.period)
        assert not has_units(results.duration)
        assert not has_units(results.transit_time)
    else:
        assert results.period.unit == t_unit
        assert results.duration.unit == t_unit
        assert results.transit_time.unit == t_unit

    if y_unit is None:
        assert not has_units(results.power)
        assert not has_units(results.depth)
        assert not has_units(results.depth_err)
        assert not has_units(results.depth_snr)
        assert not has_units(results.log_likelihood)
    else:
        assert results.depth.unit == y_unit
        assert results.depth_err.unit == y_unit
        assert results.depth_snr.unit == units.one

        if dy is None:
            assert results.log_likelihood.unit == y_unit * y_unit
            if objective == "snr":
                assert results.power.unit == units.one
            else:
                assert results.power.unit == y_unit * y_unit
        else:
            assert results.log_likelihood.unit == units.one
            assert results.power.unit == units.one


def test_autopower(data):
    t, y, dy, params = data
    duration = params["duration"] + np.linspace(-0.1, 0.1, 3)

    model = BoxLeastSquares(t, y, dy)
    period = model.autoperiod(duration)
    results1 = model.power(period, duration)
    results2 = model.autopower(duration)

    assert_allclose_blsresults(results1, results2)


@pytest.mark.parametrize("with_units", [True, False])
def test_model(data, with_units):
    t, y, dy, params = data

    # Compute the model using linear regression
    A = np.zeros((len(t), 2))
    p = params["period"]
    dt = np.abs((t-params["transit_time"]+0.5*p) % p-0.5*p)
    m_in = dt < 0.5*params["duration"]
    A[~m_in, 0] = 1.0
    A[m_in, 1] = 1.0
    w = np.linalg.solve(np.dot(A.T, A / dy[:, None]**2),
                        np.dot(A.T, y / dy**2))
    model_true = np.dot(A, w)

    if with_units:
        t = t * units.day
        y = y * units.mag
        dy = dy * units.mag
        model_true = model_true * units.mag

    # Compute the model using the periodogram
    pgram = BoxLeastSquares(t, y, dy)
    model = pgram.model(t, p, params["duration"], params["transit_time"])

    # Make sure that the transit mask is consistent with the model
    transit_mask = pgram.transit_mask(t, p, params["duration"],
                                      params["transit_time"])
    transit_mask0 = (model - model.max()) < 0.0
    assert_allclose(transit_mask, transit_mask0)

    assert_quantity_allclose(model, model_true)


@pytest.mark.parametrize("shape", [(1,), (2,), (3,), (2, 3)])
def test_shapes(data, shape):
    t, y, dy, params = data
    duration = params["duration"]
    model = BoxLeastSquares(t, y, dy)

    period = np.empty(shape)
    period.flat = np.linspace(params["period"]-1, params["period"]+1,
                              period.size)
    if len(period.shape) > 1:
        with pytest.raises(ValueError):
            results = model.power(period, duration)
    else:
        results = model.power(period, duration)
        for k, v in results.items():
            if k == "objective":
                continue
            assert v.shape == shape


@pytest.mark.parametrize("with_units", [True, False])
@pytest.mark.parametrize("with_err", [True, False])
def test_compute_stats(data, with_units, with_err):
    t, y, dy, params = data

    y_unit = 1
    if with_units:
        y_unit = units.mag
        t = t * units.day
        y = y * units.mag
        dy = dy * units.mag
        params["period"] = params["period"] * units.day
        params["duration"] = params["duration"] * units.day
        params["transit_time"] = params["transit_time"] * units.day
        params["depth"] = params["depth"] * units.mag
    if not with_err:
        dy = None

    model = BoxLeastSquares(t, y, dy)
    results = model.power(params["period"], params["duration"],
                          oversample=1000)
    stats = model.compute_stats(params["period"], params["duration"],
                                params["transit_time"])

    # Test the calculated transit times
    tt = params["period"] * np.arange(int(t.max() / params["period"]) + 1)
    tt += params["transit_time"]
    assert_quantity_allclose(tt, stats["transit_times"])

    # Test that the other parameters are consistent with the periodogram
    assert_allclose(stats["per_transit_count"], np.array([9, 7, 7, 7, 8]))
    assert_quantity_allclose(np.sum(stats["per_transit_log_likelihood"]),
                             results["log_likelihood"])
    assert_quantity_allclose(stats["depth"][0], results["depth"])

    # Check the half period result
    results_half = model.power(0.5*params["period"], params["duration"],
                               oversample=1000)
    assert_quantity_allclose(stats["depth_half"][0], results_half["depth"])

    # Skip the uncertainty tests when the input errors are None
    if not with_err:
        assert_quantity_allclose(stats["harmonic_amplitude"],
                                 0.029945029964964204 * y_unit)
        assert_quantity_allclose(stats["harmonic_delta_log_likelihood"],
                                 -0.5875918155223113 * y_unit * y_unit)
        return

    assert_quantity_allclose(stats["harmonic_amplitude"],
                             0.033027988742275853 * y_unit)
    assert_quantity_allclose(stats["harmonic_delta_log_likelihood"],
                             -12407.505922833765)

    assert_quantity_allclose(stats["depth"][1], results["depth_err"])
    assert_quantity_allclose(stats["depth_half"][1], results_half["depth_err"])
    for f, k in zip((1.0, 1.0, 1.0, 0.0),
                    ("depth", "depth_even", "depth_odd", "depth_phased")):
        assert np.abs((stats[k][0]-f*params["depth"]) / stats[k][1]) < 1.0


def test_negative_times(data):
    t, y, dy, params = data
    mu = np.mean(t)
    duration = params["duration"] + np.linspace(-0.1, 0.1, 3)

    model1 = BoxLeastSquares(t, y, dy)
    results1 = model1.autopower(duration)

    # Compute the periodogram with offset (negative) times
    model2 = BoxLeastSquares(t - mu, y, dy)
    results2 = model2.autopower(duration)

    # Shift the transit times back into the unshifted coordinates
    results2.transit_time = (results2.transit_time + mu) % results2.period

    assert_allclose_blsresults(results1, results2)
