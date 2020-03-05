# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal

from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.tests.helper import assert_quantity_allclose
from astropy.timeseries.periodograms.bls import BoxLeastSquares
from astropy.timeseries.periodograms.lombscargle.core import has_units


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
            raise AssertionError(f"missing key '{k}'")
        if k == "objective":
            assert v == other[k], (
                "Mismatched objectives. Expected '{}', got '{}'"
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

    t_unit = u.day
    y_unit = u.mag

    with pytest.raises(u.UnitConversionError):
        BoxLeastSquares(t * t_unit, y * y_unit, dy * u.one)
    with pytest.raises(u.UnitConversionError):
        BoxLeastSquares(t * t_unit, y * u.one, dy * y_unit)
    with pytest.raises(u.UnitConversionError):
        BoxLeastSquares(t * t_unit, y, dy * y_unit)
    model = BoxLeastSquares(t*t_unit, y * u.one, dy)
    assert model.dy.unit == model.y.unit
    model = BoxLeastSquares(t*t_unit, y * y_unit, dy)
    assert model.dy.unit == model.y.unit
    model = BoxLeastSquares(t*t_unit, y*y_unit)
    assert model.dy is None


def test_period_units(data):
    t, y, dy, params = data
    t_unit = u.day
    y_unit = u.mag
    model = BoxLeastSquares(t * t_unit, y * y_unit, dy)

    p = model.autoperiod(params["duration"])
    assert p.unit == t_unit
    p = model.autoperiod(params["duration"] * 24 * u.hour)
    assert p.unit == t_unit
    with pytest.raises(u.UnitConversionError):
        model.autoperiod(params["duration"] * u.mag)

    p = model.autoperiod(params["duration"], minimum_period=0.5)
    assert p.unit == t_unit
    with pytest.raises(u.UnitConversionError):
        p = model.autoperiod(params["duration"], minimum_period=0.5*u.mag)

    p = model.autoperiod(params["duration"], maximum_period=0.5)
    assert p.unit == t_unit
    with pytest.raises(u.UnitConversionError):
        p = model.autoperiod(params["duration"], maximum_period=0.5*u.mag)

    p = model.autoperiod(params["duration"], minimum_period=0.5,
                         maximum_period=1.5)
    p2 = model.autoperiod(params["duration"], maximum_period=0.5,
                          minimum_period=1.5)
    assert_quantity_allclose(p, p2)


@pytest.mark.parametrize("method", ["fast", "slow"])
@pytest.mark.parametrize("with_err", [True, False])
@pytest.mark.parametrize("t_unit", [None, u.day])
@pytest.mark.parametrize("y_unit", [None, u.mag])
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
        assert results.depth_snr.unit == u.one

        if dy is None:
            assert results.log_likelihood.unit == y_unit * y_unit
            if objective == "snr":
                assert results.power.unit == u.one
            else:
                assert results.power.unit == y_unit * y_unit
        else:
            assert results.log_likelihood.unit == u.one
            assert results.power.unit == u.one


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
        t = t * u.day
        y = y * u.mag
        dy = dy * u.mag
        model_true = model_true * u.mag

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
        y_unit = u.mag
        t = t * u.day
        y = y * u.mag
        dy = dy * u.mag
        params["period"] = params["period"] * u.day
        params["duration"] = params["duration"] * u.day
        params["transit_time"] = params["transit_time"] * u.day
        params["depth"] = params["depth"] * u.mag
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


@pytest.mark.parametrize('timedelta', [False, True])
def test_absolute_times(data, timedelta):

    # Make sure that we handle absolute times correctly. We also check that
    # TimeDelta works properly when timedelta is True.

    # The example data uses relative times
    t, y, dy, params = data

    # Add units
    t = t * u.day
    y = y * u.mag
    dy = dy * u.mag

    # We now construct a set of absolute times but keeping the rest the same.
    start = Time('2019-05-04T12:34:56')
    trel = TimeDelta(t) if timedelta else t
    t = trel + start

    # and we set up two instances of BoxLeastSquares, one with absolute and one
    # with relative times.
    bls1 = BoxLeastSquares(t, y, dy)
    bls2 = BoxLeastSquares(trel, y, dy)

    results1 = bls1.autopower(0.16 * u.day)
    results2 = bls2.autopower(0.16 * u.day)

    # All the results should match except transit time which should be
    # absolute instead of relative in the first case.

    for key in results1:
        if key == 'transit_time':
            assert_quantity_allclose((results1[key] - start).to(u.day), results2[key])
        elif key == 'objective':
            assert results1[key] == results2[key]
        else:
            assert_allclose(results1[key], results2[key])

    # Check that model evaluation works fine

    model1 = bls1.model(t, 0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    model2 = bls2.model(trel, 0.2 * u.day, 0.05 * u.day, TimeDelta(1 * u.day))
    assert_quantity_allclose(model1, model2)

    # Check model validation

    with pytest.raises(TypeError) as exc:
        bls1.model(t, 0.2 * u.day, 0.05 * u.day, 1 * u.day)
    assert exc.value.args[0] == ('transit_time was provided as a relative time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with absolute times.')

    with pytest.raises(TypeError) as exc:
        bls1.model(trel, 0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    assert exc.value.args[0] == ('t_model was provided as a relative time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with absolute times.')

    with pytest.raises(TypeError) as exc:
        bls2.model(trel, 0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    assert exc.value.args[0] == ('transit_time was provided as an absolute time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with relative times.')

    with pytest.raises(TypeError) as exc:
        bls2.model(t, 0.2 * u.day, 0.05 * u.day, 1 * u.day)
    assert exc.value.args[0] == ('t_model was provided as an absolute time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with relative times.')

    # Check compute_stats

    stats1 = bls1.compute_stats(0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    stats2 = bls2.compute_stats(0.2 * u.day, 0.05 * u.day, 1 * u.day)

    for key in stats1:
        if key == 'transit_times':
            assert_quantity_allclose((stats1[key] - start).to(u.day), stats2[key], atol=1e-10 * u.day)
        elif key.startswith('depth'):
            for value1, value2 in zip(stats1[key], stats2[key]):
                assert_quantity_allclose(value1, value2)
        else:
            assert_allclose(stats1[key], stats2[key])

    # Check compute_stats validation

    with pytest.raises(TypeError) as exc:
        bls1.compute_stats(0.2 * u.day, 0.05 * u.day, 1 * u.day)
    assert exc.value.args[0] == ('transit_time was provided as a relative time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with absolute times.')

    with pytest.raises(TypeError) as exc:
        bls2.compute_stats(0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    assert exc.value.args[0] == ('transit_time was provided as an absolute time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with relative times.')

    # Check transit_mask

    mask1 = bls1.transit_mask(t, 0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    mask2 = bls2.transit_mask(trel, 0.2 * u.day, 0.05 * u.day, 1 * u.day)

    assert_equal(mask1, mask2)

    # Check transit_mask validation

    with pytest.raises(TypeError) as exc:
        bls1.transit_mask(t, 0.2 * u.day, 0.05 * u.day, 1 * u.day)
    assert exc.value.args[0] == ('transit_time was provided as a relative time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with absolute times.')

    with pytest.raises(TypeError) as exc:
        bls1.transit_mask(trel, 0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    assert exc.value.args[0] == ('t was provided as a relative time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with absolute times.')

    with pytest.raises(TypeError) as exc:
        bls2.transit_mask(trel, 0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    assert exc.value.args[0] == ('transit_time was provided as an absolute time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with relative times.')

    with pytest.raises(TypeError) as exc:
        bls2.transit_mask(t, 0.2 * u.day, 0.05 * u.day, 1 * u.day)
    assert exc.value.args[0] == ('t was provided as an absolute time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with relative times.')


def test_transit_time_in_range(data):
    t, y, dy, params = data

    t_ref = 10230.0
    t2 = t + t_ref
    bls1 = BoxLeastSquares(t, y, dy)
    bls2 = BoxLeastSquares(t2, y, dy)

    results1 = bls1.autopower(0.16)
    results2 = bls2.autopower(0.16)

    assert np.allclose(results1.transit_time, results2.transit_time - t_ref)
    assert np.all(results1.transit_time >= t.min())
    assert np.all(results1.transit_time <= t.max())
    assert np.all(results2.transit_time >= t2.min())
    assert np.all(results2.transit_time <= t2.max())
