import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy import units as u
from astropy.table import MaskedColumn
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time, TimeDelta
from astropy.timeseries import TimeSeries
from astropy.timeseries.periodograms.lombscargle import LombScargle
from astropy.timeseries.periodograms.lombscargle_multiband import LombScargleMultiband

ALL_METHODS = LombScargleMultiband.available_methods
ALL_SB_METHODS = LombScargle.available_methods
NORMALIZATIONS = ["standard", "psd", "log", "model"]


@pytest.fixture
def data(N=100, period=1, theta=[10, 2, 3], nbands=3, dy=1, rseed=0):
    """Generate some data for testing"""
    t_arr = []
    y_arr = []
    band_arr = []
    dy_arr = []

    for band in range(nbands):
        rng = np.random.default_rng(rseed + band)
        t_band = 20 * period * rng.random(N)
        omega = 2 * np.pi / period
        y_band = (
            theta[0]
            + theta[1] * np.sin(omega * t_band)
            + theta[2] * np.cos(omega * t_band)
        )
        dy_band = dy * (0.5 + rng.random(N))
        y_band += dy_band * rng.standard_normal(N)

        t_arr += list(t_band)
        y_arr += list(y_band)
        dy_arr += list(dy_band)
        band_arr += ["a" * (band + 1)] * N  # labels bands as "a","aa","aaa",....

    t_arr = np.array(t_arr)
    y_arr = np.array(y_arr)
    band_arr = np.array(band_arr)
    dy_arr = np.array(dy_arr)

    return t_arr, y_arr, band_arr, dy_arr


@pytest.fixture
def timeseries_data():
    """Generate an astropy.timeseries.TimeSeries table"""
    rng = np.random.default_rng(1)
    deltas = 240 * rng.random(180)
    ts1 = TimeSeries(time_start="2011-01-01T00:00:00", time_delta=deltas * u.minute)

    # g band fluxes
    g_flux = [0] * 180 * u.mJy
    g_err = [0] * 180 * u.mJy
    y_g = np.round(3 + 2 * np.sin(10 * np.pi * ts1["time"].mjd[0:60]), 3)
    dy_g = np.round(0.01 * (0.5 + rng.random(60)), 3)  # uncertainties
    g_flux.value[0:60] = y_g
    g_err.value[0:60] = dy_g
    ts1["g_flux"] = MaskedColumn(g_flux, mask=[False] * 60 + [True] * 120)
    ts1["g_err"] = MaskedColumn(g_err, mask=[False] * 60 + [True] * 120)
    # r band fluxes
    r_flux = [0] * 180 * u.mJy
    r_err = [0] * 180 * u.mJy
    y_r = np.round(3 + 2 * np.sin(10 * np.pi * ts1["time"].mjd[60:120]), 3)
    dy_r = np.round(0.01 * (0.5 + rng.random(60)), 3)  # uncertainties
    r_flux.value[60:120] = y_r
    r_err.value[60:120] = dy_r
    ts1["r_flux"] = MaskedColumn(r_flux, mask=[True] * 60 + [False] * 60 + [True] * 60)
    ts1["r_err"] = MaskedColumn(r_err, mask=[True] * 60 + [False] * 60 + [True] * 60)
    # i band fluxes
    i_flux = [0] * 180 * u.mJy
    i_err = [0] * 180 * u.mJy
    y_i = np.round(3 + 2 * np.sin(10 * np.pi * ts1["time"].mjd[120:]), 3)
    dy_i = np.round(0.01 * (0.5 + rng.random(60)), 3)  # uncertainties
    i_flux.value[120:] = y_i
    i_err.value[120:] = dy_i
    ts1["i_flux"] = MaskedColumn(i_flux, mask=[True] * 120 + [False] * 60)
    ts1["i_err"] = MaskedColumn(i_err, mask=[True] * 120 + [False] * 60)

    return ts1


@pytest.mark.parametrize("minimum_frequency", [None, 1.0])
@pytest.mark.parametrize("maximum_frequency", [None, 5.0])
@pytest.mark.parametrize("nyquist_factor", [1, 10])
@pytest.mark.parametrize("samples_per_peak", [1, 5])
def test_autofrequency(
    data, minimum_frequency, maximum_frequency, nyquist_factor, samples_per_peak
):
    t, y, band, dy = data
    baseline = t.max() - t.min()

    freq = LombScargleMultiband(t, y, band, dy).autofrequency(
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


@pytest.mark.parametrize("method", ALL_METHODS)
@pytest.mark.parametrize("nterms_base", [1, 3])
@pytest.mark.parametrize("nterms_band", [0, 1])
@pytest.mark.parametrize("center_data", [True, False])
@pytest.mark.parametrize("errors", ["none", "partial", "full"])
@pytest.mark.parametrize("with_units", [True, False])
@pytest.mark.parametrize("normalization", NORMALIZATIONS)
def test_all_methods(
    data,
    method,
    nterms_base,
    nterms_band,
    center_data,
    errors,
    with_units,
    normalization,
):
    t, y, band, dy = data
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

    ls = LombScargleMultiband(
        t,
        y,
        band,
        dy,
        nterms_base=nterms_base,
        nterms_band=nterms_band,
        center_data=center_data,
        normalization=normalization,
    )

    P_expected = ls.power(frequency, method=method)

    P_method = ls.power(frequency, method=method, **kwds)
    freq_maxpower = frequency[np.argmax(P_method)]
    if with_units:
        assert P_method.unit == u.dimensionless_unscaled
        assert np.isclose(
            freq_maxpower.value, 1.0, rtol=1e-2
        )  # period=1 check peak frequency
    else:
        assert not hasattr(P_method, "unit")
        assert np.isclose(
            freq_maxpower, 1.0, rtol=1e-2
        )  # period=1, check peak frequency

    assert_quantity_allclose(P_expected, P_method)


@pytest.mark.parametrize("method", ALL_METHODS)
@pytest.mark.parametrize("nterms_base", [1, 3])
@pytest.mark.parametrize("nterms_band", [0, 1])
@pytest.mark.parametrize("center_data", [True, False])
@pytest.mark.parametrize("with_errors", [True, False])
@pytest.mark.parametrize("normalization", NORMALIZATIONS)
def test_integer_inputs(
    data, method, nterms_base, nterms_band, center_data, with_errors, normalization
):
    if method == "scipy" and with_errors:
        return

    t, y, band, dy = data

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

    kwds = dict(center_data=center_data, normalization=normalization)
    P_float = LombScargleMultiband(t, y, band, dy, **kwds).power(
        frequency, method=method
    )
    P_int = LombScargleMultiband(t_int, y_int, band, dy_int, **kwds).power(
        frequency, method=method
    )
    assert_allclose(P_float, P_int)


@pytest.mark.parametrize("method", ["flexible"])
@pytest.mark.parametrize("nterms_base", [0, 1, 2, 3])
@pytest.mark.parametrize("nterms_band", [0, 1, 2])
@pytest.mark.parametrize("center_data", [True, False])
@pytest.mark.parametrize("errors", ["none", "partial", "full"])
@pytest.mark.parametrize("normalization", NORMALIZATIONS)
def test_nterms_methods(
    method, nterms_base, nterms_band, center_data, errors, normalization, data
):
    t, y, band, dy = data
    frequency = 0.8 + 0.01 * np.arange(40)

    if errors == "none":
        dy = None
    elif errors == "partial":
        dy = dy[0]
    elif errors == "full":
        pass
    else:
        raise ValueError(f"Unrecognized error type: '{errors}'")

    ls = LombScargleMultiband(
        t,
        y,
        band,
        dy,
        center_data=center_data,
        nterms_base=nterms_base,
        nterms_band=nterms_band,
        normalization=normalization,
    )

    if (nterms_base == 0) and (nterms_band == 0):
        with pytest.raises(ValueError) as err:
            ls.power(frequency, method=method)
        assert "nterms_base" in str(err.value)
    else:
        P_expected = ls.power(frequency)

        # don't use fast fft approximations here
        kwds = {}
        if "fast" in method:
            kwds["method_kwds"] = dict(use_fft=False)
        P_method = ls.power(frequency, method=method, **kwds)

        assert_allclose(P_expected, P_method, rtol=1e-7, atol=1e-25)


@pytest.mark.parametrize("method", ALL_METHODS)
@pytest.mark.parametrize("shape", [(), (1,), (2,), (3,), (2, 3)])
def test_output_shapes(method, shape, data):
    t, y, band, dy = data
    freq = np.asarray(np.zeros(shape))
    freq.flat = np.arange(1, freq.size + 1)
    PLS = LombScargleMultiband(t, y, band).power(freq, method=method)
    assert PLS.shape == shape


@pytest.mark.parametrize("method", ALL_METHODS)
def test_errors_on_unit_mismatch(method, data):
    t, y, band, dy = data

    t = t * u.second
    y = y * u.mag
    frequency = np.linspace(0.5, 1.5, 10)

    # this should fail because frequency and 1/t units do not match
    with pytest.raises(ValueError) as err:
        LombScargleMultiband(t, y, band).power(frequency, method=method)
    assert str(err.value).startswith("Units of frequency not equivalent")

    # this should fail because dy and y units do not match
    with pytest.raises(ValueError) as err:
        LombScargleMultiband(t, y, band, dy).power(frequency / t.unit)
    assert str(err.value).startswith("Units of dy not equivalent")


@pytest.mark.parametrize("method", ALL_METHODS)
@pytest.mark.parametrize("with_error", [True, False])
def test_unit_conversions(data, method, with_error):
    t, y, band, dy = data

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

    freq_day, P1 = LombScargleMultiband(t_day, y_meter, band, dy).autopower(
        method=method
    )
    freq_hour, P2 = LombScargleMultiband(t_hour, y_millimeter, band, dy).autopower(
        method=method
    )

    # Check units of frequency
    assert freq_day.unit == 1.0 / u.day
    assert freq_hour.unit == 1.0 / u.hour

    # Check that results match
    assert_quantity_allclose(freq_day, freq_hour)
    assert_quantity_allclose(P1, P2)

    # Check that switching frequency units doesn't change things
    P3 = LombScargleMultiband(t_day, y_meter, band, dy).power(freq_hour, method=method)
    P4 = LombScargleMultiband(t_hour, y_meter, band, dy).power(freq_day, method=method)
    assert_quantity_allclose(P3, P4)


@pytest.mark.parametrize("with_units", [True, False])
@pytest.mark.parametrize("freq", [1.0, 2.0])
def test_model(with_units, freq):
    rand = np.random.default_rng(0)
    t = 10 * rand.random(120)
    band = 40 * ["a"] + 40 * ["b"] + 40 * ["c"]
    params = 10 * rand.random(3)

    y = np.zeros_like(t)
    y += params[1] * np.sin(2 * np.pi * freq * (t - params[2]))

    if with_units:
        t = t * u.day
        y = y * u.mag
        freq = freq / u.day

    ls = LombScargleMultiband(t, y, band, center_data=False)
    y_fit = ls.model(t, freq, bands_fit=None)
    assert_quantity_allclose(y_fit[0][0:40], y[0:40])
    assert_quantity_allclose(y_fit[1][40:80], y[40:80])
    assert_quantity_allclose(y_fit[2][80:], y[80:])


@pytest.mark.parametrize("t_unit", [u.second, u.day])
@pytest.mark.parametrize("frequency_unit", [u.Hz, 1.0 / u.second])
@pytest.mark.parametrize("y_unit", [u.mag, u.jansky])
def test_model_units_match(data, t_unit, frequency_unit, y_unit):
    t, y, band, dy = data
    t_fit = t[:5]
    frequency = 1.0

    t = t * t_unit
    t_fit = t_fit * t_unit
    y = y * y_unit
    dy = dy * y_unit
    frequency = frequency * frequency_unit

    ls = LombScargleMultiband(t, y, band, dy)
    y_fit = ls.model(t_fit, frequency)
    assert y_fit.unit == y_unit


def test_model_units_mismatch(data):
    t, y, band, dy = data
    frequency = 1.0
    t_fit = t[:5]

    t = t * u.second
    t_fit = t_fit * u.second
    y = y * u.mag
    frequency = 1.0 / t.unit

    # this should fail because frequency and 1/t units do not match
    with pytest.raises(ValueError) as err:
        LombScargleMultiband(t, y, band).model(t_fit, frequency=1.0)
    assert str(err.value).startswith("Units of frequency not equivalent")

    # this should fail because t and t_fit units do not match
    with pytest.raises(ValueError) as err:
        LombScargleMultiband(t, y, band).model([1, 2], frequency)
    assert str(err.value).startswith("Units of t not equivalent")

    # this should fail because dy and y units do not match
    with pytest.raises(ValueError) as err:
        LombScargleMultiband(t, y, band, dy).model(t_fit, frequency)
    assert str(err.value).startswith("Units of dy not equivalent")


def test_autopower(data):
    t, y, band, dy = data
    ls = LombScargleMultiband(t, y, band, dy)
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
@pytest.mark.parametrize("nterms_base", [0, 1, 2])
@pytest.mark.parametrize("nterms_band", [0, 1])
def test_model_parameters(
    data, nterms_base, nterms_band, center_data, errors, with_units
):
    if (nterms_base == 0) and (nterms_band == 0):
        return

    t, y, band, dy = data
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

    ls = LombScargleMultiband(
        t, y, band, dy, nterms_base=nterms_base, center_data=center_data
    )
    tfit = np.linspace(0, 20, 10)
    if with_units:
        tfit = tfit * u.day

    model = ls.model(tfit, frequency)
    params = ls.model_parameters(frequency)
    design = ls.design_matrix(frequency, t_fit=tfit, bands_fit=None)
    offset = ls.offset(t_fit=tfit)

    if nterms_band == 0:
        nterms_band = 1
    assert len(params) == 1 + 2 * nterms_base + len(np.unique(band)) * (
        2 * nterms_band + 1
    )

    from_funcs = offset + design.dot(params)
    from_funcs = from_funcs.reshape((len(np.unique(band)), len(tfit)))
    assert_quantity_allclose(from_funcs, model)


@pytest.mark.parametrize("timedelta", [False, True])
def test_absolute_times(data, timedelta):
    # Make sure that we handle absolute times correctly. We also check that
    # TimeDelta works properly when timedelta is True.

    # The example data uses relative times
    t, y, band, dy = data

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
    ls1 = LombScargleMultiband(t, y, band, dy)
    ls2 = LombScargleMultiband(trel, y, band, dy)

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

    with pytest.raises(TypeError) as exc:
        ls1.model(trel, 2 / u.day)
    assert exc.value.args[0] == (
        "t was provided as a relative time but the "
        "LombScargle class was initialized with "
        "absolute times."
    )

    with pytest.raises(TypeError) as exc:
        ls2.model(t, 2 / u.day)
    assert exc.value.args[0] == (
        "t was provided as an absolute time but the "
        "LombScargle class was initialized with "
        "relative times."
    )

    # Check design matrix

    design1 = ls1.design_matrix(2 / u.day, t_fit=t)
    design2 = ls2.design_matrix(2 / u.day, t_fit=trel)
    assert_quantity_allclose(design1, design2)

    # Check design matrix validation

    with pytest.raises(TypeError) as exc:
        ls1.design_matrix(2 / u.day, t_fit=trel)
    assert exc.value.args[0] == (
        "t was provided as a relative time but the "
        "LombScargle class was initialized with "
        "absolute times."
    )

    with pytest.raises(TypeError) as exc:
        ls2.design_matrix(2 / u.day, t_fit=t)
    assert exc.value.args[0] == (
        "t was provided as an absolute time but the "
        "LombScargle class was initialized with "
        "relative times."
    )


@pytest.mark.parametrize("uncertainty_column", [None, ["g_err", "r_err", "i_err"]])
@pytest.mark.parametrize("band_labels", [None, ["g", "r", "i"]])
def test_from_timeseries(timeseries_data, uncertainty_column, band_labels):
    ts = timeseries_data

    ls = LombScargleMultiband.from_timeseries(
        ts,
        signal_column=["g_flux", "r_flux", "i_flux"],
        uncertainty_column=["g_err", "r_err", "i_err"],
        band_labels=["g", "r", "i"],
    )

    frequency, power = ls.autopower()
    freq_maxpower = frequency[np.argmax(power)]

    assert_allclose(freq_maxpower.value, 5, rtol=0.01)


@pytest.mark.parametrize("errors", ["none", "partial", "full"])
@pytest.mark.parametrize("with_units", [True, False])
@pytest.mark.parametrize("sb_method", ALL_SB_METHODS)
def test_single_band_equivalence(data, with_units, errors, sb_method):
    fit_mean = True
    if sb_method == "scipy":
        fit_mean = False

    t, y, band, dy = data
    frequency = 0.8 + 0.01 * np.arange(40)

    # select just one band
    a_mask = band == "a"
    t = t[a_mask]
    y = y[a_mask]
    band = band[a_mask]
    dy = dy[a_mask]

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
        if sb_method == "scipy":
            return
        pass
    else:
        raise ValueError(f"Unrecognized error type: '{errors}'")

    lsmb = LombScargleMultiband(t, y, band, dy, fit_mean=fit_mean)
    P_lsmb = lsmb.power(frequency, method="fast", sb_method=sb_method)

    ls = LombScargle(t, y, dy, fit_mean=fit_mean)
    P_ls = ls.power(frequency, method=sb_method)

    # test to see if lombscargle multiband and lombscargle are equivalent
    assert_quantity_allclose(P_lsmb, P_ls)
