import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time, TimeDelta
from astropy.timeseries.periodograms.lombscargle import LombScargle

ALL_METHODS = LombScargle.available_methods
ALL_METHODS_NO_AUTO = [method for method in ALL_METHODS if method != "auto"]
FAST_METHODS = [method for method in ALL_METHODS if "fast" in method]
NTERMS_METHODS = [method for method in ALL_METHODS if "chi2" in method]
NORMALIZATIONS = ["standard", "psd", "log", "model"]


@pytest.fixture
def data(N=100, period=1, theta=[10, 2, 3], dy=1, rseed=0):
    """Generate some data for testing"""
    rng = np.random.default_rng(rseed)
    t = 20 * period * rng.random(N)
    omega = 2 * np.pi / period
    y = theta[0] + theta[1] * np.sin(omega * t) + theta[2] * np.cos(omega * t)
    dy = dy * (0.5 + rng.random(N))
    y += dy * rng.standard_normal(N)

    return t, y, dy


@pytest.mark.parametrize("minimum_frequency", [None, 1.0])
@pytest.mark.parametrize("maximum_frequency", [None, 5.0])
@pytest.mark.parametrize("nyquist_factor", [1, 10])
@pytest.mark.parametrize("samples_per_peak", [1, 5])
def test_autofrequency(
    data, minimum_frequency, maximum_frequency, nyquist_factor, samples_per_peak
):
    t, y, dy = data
    baseline = t.max() - t.min()

    freq = LombScargle(t, y, dy).autofrequency(
        samples_per_peak, nyquist_factor, minimum_frequency, maximum_frequency
    )
    df = freq[1] - freq[0]

    # Check sample spacing
    assert_allclose(df, 1.0 / baseline / samples_per_peak)

    # Check minimum frequency
    if minimum_frequency is None:
        assert_allclose(freq[0], 0.5 * df)
    else:
        assert_allclose(freq[0], minimum_frequency)

    if maximum_frequency is None:
        avg_nyquist = 0.5 * len(t) / baseline
        assert_allclose(freq[-1], avg_nyquist * nyquist_factor, atol=0.5 * df)
    else:
        assert_allclose(freq[-1], maximum_frequency, atol=0.5 * df)


@pytest.mark.parametrize("method", ALL_METHODS_NO_AUTO)
@pytest.mark.parametrize("center_data", [True, False])
@pytest.mark.parametrize("fit_mean", [True, False])
@pytest.mark.parametrize("errors", ["none", "partial", "full"])
@pytest.mark.parametrize("with_units", [True, False])
@pytest.mark.parametrize("normalization", NORMALIZATIONS)
def test_all_methods(
    data, method, center_data, fit_mean, errors, with_units, normalization
):
    if method == "scipy" and (fit_mean or errors != "none"):
        return

    t, y, dy = data
    frequency = 0.8 + 0.01 * np.arange(40)
    if with_units:
        t = t * u.day
        y = y * u.mag
        dy = dy * u.mag
        frequency = frequency / t.unit

    if errors == "none":
        dy = None
    elif errors == "partial":
        dy = dy[0]
    elif errors == "full":
        pass
    else:
        raise ValueError(f"Unrecognized error type: '{errors}'")

    kwds = {}

    ls = LombScargle(
        t,
        y,
        dy,
        center_data=center_data,
        fit_mean=fit_mean,
        normalization=normalization,
    )
    P_expected = ls.power(frequency)

    # don't use the fft approximation here; we'll test this elsewhere
    if method in FAST_METHODS:
        kwds["method_kwds"] = dict(use_fft=False)
    P_method = ls.power(frequency, method=method, **kwds)

    if with_units:
        if normalization == "psd" and errors == "none":
            assert P_method.unit == y.unit**2
        else:
            assert P_method.unit == u.dimensionless_unscaled
    else:
        assert not hasattr(P_method, "unit")

    assert_quantity_allclose(P_expected, P_method)


@pytest.mark.parametrize("method", ALL_METHODS_NO_AUTO)
@pytest.mark.parametrize("center_data", [True, False])
@pytest.mark.parametrize("fit_mean", [True, False])
@pytest.mark.parametrize("with_errors", [True, False])
@pytest.mark.parametrize("normalization", NORMALIZATIONS)
def test_integer_inputs(
    data, method, center_data, fit_mean, with_errors, normalization
):
    if method == "scipy" and (fit_mean or with_errors):
        return

    t, y, dy = data

    t = np.floor(100 * t)
    t_int = t.astype(int)

    y = np.floor(100 * y)
    y_int = y.astype(int)

    dy = np.floor(100 * dy)
    dy_int = dy.astype("int32")

    frequency = 1e-2 * (0.8 + 0.01 * np.arange(40))

    if not with_errors:
        dy = None
        dy_int = None

    kwds = dict(center_data=center_data, fit_mean=fit_mean, normalization=normalization)
    P_float = LombScargle(t, y, dy, **kwds).power(frequency, method=method)
    P_int = LombScargle(t_int, y_int, dy_int, **kwds).power(frequency, method=method)
    assert_allclose(P_float, P_int)


@pytest.mark.parametrize("method", NTERMS_METHODS)
@pytest.mark.parametrize("center_data", [True, False])
@pytest.mark.parametrize("fit_mean", [True, False])
@pytest.mark.parametrize("errors", ["none", "partial", "full"])
@pytest.mark.parametrize("nterms", [0, 2, 4])
@pytest.mark.parametrize("normalization", NORMALIZATIONS)
def test_nterms_methods(
    method, center_data, fit_mean, errors, nterms, normalization, data
):
    t, y, dy = data
    frequency = 0.8 + 0.01 * np.arange(40)

    if errors == "none":
        dy = None
    elif errors == "partial":
        dy = dy[0]
    elif errors == "full":
        pass
    else:
        raise ValueError(f"Unrecognized error type: '{errors}'")

    ls = LombScargle(
        t,
        y,
        dy,
        center_data=center_data,
        fit_mean=fit_mean,
        nterms=nterms,
        normalization=normalization,
    )

    if nterms == 0 and not fit_mean:
        with pytest.raises(ValueError, match=r"[nterms, blas]"):
            ls.power(frequency, method=method)
    else:
        P_expected = ls.power(frequency)

        # don't use fast fft approximations here
        kwds = {}
        if "fast" in method:
            kwds["method_kwds"] = dict(use_fft=False)
        P_method = ls.power(frequency, method=method, **kwds)

        assert_allclose(P_expected, P_method, rtol=1e-7, atol=1e-25)


@pytest.mark.parametrize("method", FAST_METHODS)
@pytest.mark.parametrize("center_data", [True, False])
@pytest.mark.parametrize("fit_mean", [True, False])
@pytest.mark.parametrize("errors", ["none", "partial", "full"])
@pytest.mark.parametrize("nterms", [0, 1, 2])
def test_fast_approximations(method, center_data, fit_mean, errors, nterms, data):
    t, y, dy = data
    frequency = 0.8 + 0.01 * np.arange(40)

    if errors == "none":
        dy = None
    elif errors == "partial":
        dy = dy[0]
    elif errors == "full":
        pass
    else:
        raise ValueError(f"Unrecognized error type: '{errors}'")

    ls = LombScargle(
        t,
        y,
        dy,
        center_data=center_data,
        fit_mean=fit_mean,
        nterms=nterms,
        normalization="standard",
    )

    # use only standard normalization because we compare via absolute tolerance
    kwds = dict(method=method)

    if method == "fast" and nterms != 1:
        with pytest.raises(ValueError, match=r"nterms"):
            ls.power(frequency, **kwds)

    elif nterms == 0 and not fit_mean:
        with pytest.raises(ValueError, match=r"[nterms, blas]"):
            ls.power(frequency, **kwds)

    else:
        P_fast = ls.power(frequency, **kwds)
        kwds["method_kwds"] = dict(use_fft=False)
        P_slow = ls.power(frequency, **kwds)

        assert_allclose(P_fast, P_slow, atol=0.008)


@pytest.mark.parametrize("method", LombScargle.available_methods)
@pytest.mark.parametrize("shape", [(), (1,), (2,), (3,), (2, 3)])
def test_output_shapes(method, shape, data):
    t, y, dy = data
    freq = np.asarray(np.zeros(shape))
    freq.flat = np.arange(1, freq.size + 1)
    PLS = LombScargle(t, y, fit_mean=False).power(freq, method=method)
    assert PLS.shape == shape


@pytest.mark.parametrize("method", LombScargle.available_methods)
def test_errors_on_unit_mismatch(method, data):
    t, y, dy = data

    t = t * u.second
    y = y * u.mag
    frequency = np.linspace(0.5, 1.5, 10)

    # this should fail because frequency and 1/t units do not match
    MESSAGE = r"Units of {} not equivalent"

    with pytest.raises(ValueError, match=MESSAGE.format("frequency")):
        LombScargle(t, y, fit_mean=False).power(frequency, method=method)

    # this should fail because dy and y units do not match
    with pytest.raises(ValueError, match=MESSAGE.format("dy")):
        LombScargle(t, y, dy, fit_mean=False).power(frequency / t.unit)


# we don't test all normalizations here because they are tested above
# only test method='auto' because unit handling does not depend on method
@pytest.mark.parametrize("with_error", [True, False])
def test_unit_conversions(data, with_error):
    t, y, dy = data

    t_day = t * u.day
    t_hour = u.Quantity(t_day, "hour")

    y_meter = y * u.meter
    y_millimeter = u.Quantity(y_meter, "millimeter")

    # sanity check on inputs
    assert_quantity_allclose(t_day, t_hour)
    assert_quantity_allclose(y_meter, y_millimeter)

    if with_error:
        dy = dy * u.meter
    else:
        dy = None

    freq_day, P1 = LombScargle(t_day, y_meter, dy).autopower()
    freq_hour, P2 = LombScargle(t_hour, y_millimeter, dy).autopower()

    # Check units of frequency
    assert freq_day.unit == 1.0 / u.day
    assert freq_hour.unit == 1.0 / u.hour

    # Check that results match
    assert_quantity_allclose(freq_day, freq_hour)
    assert_quantity_allclose(P1, P2)

    # Check that switching frequency units doesn't change things
    P3 = LombScargle(t_day, y_meter, dy).power(freq_hour)
    P4 = LombScargle(t_hour, y_meter, dy).power(freq_day)
    assert_quantity_allclose(P3, P4)


@pytest.mark.parametrize("fit_mean", [True, False])
@pytest.mark.parametrize("with_units", [True, False])
@pytest.mark.parametrize("freq", [1.0, 2.0])
def test_model(fit_mean, with_units, freq):
    rand = np.random.default_rng(0)
    t = 10 * rand.random(40)
    params = 10 * rand.random(3)

    y = np.zeros_like(t)
    if fit_mean:
        y += params[0]
    y += params[1] * np.sin(2 * np.pi * freq * (t - params[2]))

    if with_units:
        t = t * u.day
        y = y * u.mag
        freq = freq / u.day

    ls = LombScargle(t, y, center_data=False, fit_mean=fit_mean)
    y_fit = ls.model(t, freq)
    assert_quantity_allclose(y_fit, y)


@pytest.mark.parametrize("t_unit", [u.second, u.day])
@pytest.mark.parametrize("frequency_unit", [u.Hz, 1.0 / u.second])
@pytest.mark.parametrize("y_unit", [u.mag, u.jansky])
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

    t = t * u.second
    t_fit = t_fit * u.second
    y = y * u.mag
    frequency = 1.0 / t.unit

    # this should fail because frequency and 1/t units do not match
    MESSAGE = r"Units of {} not equivalent"

    with pytest.raises(ValueError, match=MESSAGE.format("frequency")):
        LombScargle(t, y).model(t_fit, frequency=1.0)

    # this should fail because t and t_fit units do not match
    with pytest.raises(ValueError, match=MESSAGE.format("t")):
        LombScargle(t, y).model([1, 2], frequency)

    # this should fail because dy and y units do not match
    with pytest.raises(ValueError, match=MESSAGE.format("dy")):
        LombScargle(t, y, dy).model(t_fit, frequency)


def test_autopower(data):
    t, y, dy = data
    ls = LombScargle(t, y, dy)
    kwargs = dict(
        samples_per_peak=6,
        nyquist_factor=2,
        minimum_frequency=2,
        maximum_frequency=None,
    )
    freq1 = ls.autofrequency(**kwargs)
    power1 = ls.power(freq1)
    freq2, power2 = ls.autopower(**kwargs)

    assert_allclose(freq1, freq2)
    assert_allclose(power1, power2)


@pytest.mark.parametrize("with_units", [True, False])
@pytest.mark.parametrize("errors", ["none", "partial", "full"])
@pytest.mark.parametrize("center_data", [True, False])
@pytest.mark.parametrize("fit_mean", [True, False])
@pytest.mark.parametrize("nterms", [0, 1, 2])
def test_model_parameters(data, nterms, fit_mean, center_data, errors, with_units):
    if nterms == 0 and not fit_mean:
        return

    t, y, dy = data
    frequency = 1.5
    if with_units:
        t = t * u.day
        y = y * u.mag
        dy = dy * u.mag
        frequency = frequency / t.unit

    if errors == "none":
        dy = None
    elif errors == "partial":
        dy = dy[0]
    elif errors == "full":
        pass
    else:
        raise ValueError(f"Unrecognized error type: '{errors}'")

    ls = LombScargle(
        t, y, dy, nterms=nterms, fit_mean=fit_mean, center_data=center_data
    )
    tfit = np.linspace(0, 20, 10)
    if with_units:
        tfit = tfit * u.day

    model = ls.model(tfit, frequency)
    params = ls.model_parameters(frequency)
    design = ls.design_matrix(frequency, t=tfit)
    offset = ls.offset()

    assert len(params) == int(fit_mean) + 2 * nterms

    assert_quantity_allclose(offset + design.dot(params), model)


@pytest.mark.parametrize("timedelta", [False, True])
def test_absolute_times(data, timedelta):
    # Make sure that we handle absolute times correctly. We also check that
    # TimeDelta works properly when timedelta is True.

    # The example data uses relative times
    t, y, dy = data

    # FIXME: There seems to be a numerical stability issue in that if we run
    # the algorithm with the same values but offset in time, the transit_time
    # is not offset by a fixed amount. To avoid this issue in this test, we
    # make sure the first time is also the smallest so that internally the
    # values of the relative time should be the same.
    t[0] = 0.0

    # Add units
    t = t * u.day
    y = y * u.mag
    dy = dy * u.mag

    # We now construct a set of absolute times but keeping the rest the same
    start = Time("2019-05-04T12:34:56")
    trel = TimeDelta(t) if timedelta else t
    t = trel + start

    # and we set up two instances of LombScargle, one with absolute and one
    # with relative times.
    ls1 = LombScargle(t, y, dy)
    ls2 = LombScargle(trel, y, dy)

    kwargs = dict(
        samples_per_peak=6,
        nyquist_factor=2,
        minimum_frequency=2 / u.day,
        maximum_frequency=None,
    )

    freq1 = ls1.autofrequency(**kwargs)
    freq2 = ls2.autofrequency(**kwargs)
    assert_quantity_allclose(freq1, freq2)

    power1 = ls1.power(freq1)
    power2 = ls2.power(freq2)
    assert_quantity_allclose(power1, power2)

    freq1, power1 = ls1.autopower(**kwargs)
    freq2, power2 = ls2.autopower(**kwargs)
    assert_quantity_allclose(freq1, freq2)
    assert_quantity_allclose(power1, power2)

    model1 = ls1.model(t, 2 / u.day)
    model2 = ls2.model(trel, 2 / u.day)
    assert_quantity_allclose(model1, model2)

    # Check model validation
    MESSAGE = (
        r"t was provided as {} time but the LombScargle class was initialized with {}"
        r" times."
    )

    with pytest.raises(TypeError, match=MESSAGE.format("a relative", "absolute")):
        ls1.model(trel, 2 / u.day)

    with pytest.raises(TypeError, match=MESSAGE.format("an absolute", "relative")):
        ls2.model(t, 2 / u.day)

    # Check design matrix

    design1 = ls1.design_matrix(2 / u.day, t=t)
    design2 = ls2.design_matrix(2 / u.day, t=trel)
    assert_quantity_allclose(design1, design2)

    # Check design matrix validation

    with pytest.raises(TypeError, match=MESSAGE.format("a relative", "absolute")):
        ls1.design_matrix(2 / u.day, t=trel)

    with pytest.raises(TypeError, match=MESSAGE.format("an absolute", "relative")):
        ls2.design_matrix(2 / u.day, t=t)
