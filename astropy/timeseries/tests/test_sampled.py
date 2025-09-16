# Licensed under a 3-clause BSD style license - see LICENSE.rst

from datetime import datetime

import pytest
from numpy.testing import assert_allclose, assert_equal

from astropy import units as u
from astropy.table import Column, Table
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time, TimeDelta
from astropy.timeseries.periodograms import BoxLeastSquares, LombScargle
from astropy.timeseries.sampled import TimeSeries
from astropy.units import Quantity, UnitsWarning
from astropy.utils.data import get_pkg_data_filename

INPUT_TIME = Time(["2016-03-22T12:30:31", "2015-01-21T12:30:32", "2016-03-22T12:30:40"])
PLAIN_TABLE = Table([[1, 2, 11], [3, 4, 1], [1, 1, 1]], names=["a", "b", "c"])

CSV_FILE = get_pkg_data_filename("data/sampled.csv")


def test_empty_initialization():
    ts = TimeSeries()
    ts["time"] = Time([50001, 50002, 50003], format="mjd")


def test_empty_initialization_invalid():
    # Make sure things crash when the first column added is not a time column
    ts = TimeSeries()
    with pytest.raises(
        ValueError,
        match=(
            r"TimeSeries object is invalid - expected 'time' as the first column but"
            r" found 'flux'"
        ),
    ):
        ts["flux"] = [1, 2, 3]


def test_initialize_only_time():
    ts = TimeSeries(time=INPUT_TIME)
    assert ts["time"] is ts.time
    # NOTE: the object in the table is a copy
    assert_equal(ts.time.isot, INPUT_TIME.isot)


def test_initialization_with_data():
    ts = TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=["a", "b"])
    assert_equal(ts.time.isot, INPUT_TIME.isot)
    assert_equal(ts["a"], [10, 2, 3])
    assert_equal(ts["b"], [4, 5, 6])


def test_initialize_only_data():
    with pytest.raises(
        TypeError, match=r"Either 'time' or 'time_start' should be specified"
    ):
        TimeSeries(data=[[10, 2, 3], [4, 5, 6]], names=["a", "b"])


def test_initialization_with_table():
    ts = TimeSeries(time=INPUT_TIME, data=PLAIN_TABLE)
    assert ts.colnames == ["time", "a", "b", "c"]


def test_initialization_with_time_delta():
    ts = TimeSeries(
        time_start=datetime(2018, 7, 1, 10, 10, 10),
        time_delta=TimeDelta(3, format="sec"),
        data=[[10, 2, 3], [4, 5, 6]],
        names=["a", "b"],
    )
    assert_equal(
        ts.time.isot,
        [
            "2018-07-01T10:10:10.000",
            "2018-07-01T10:10:13.000",
            "2018-07-01T10:10:16.000",
        ],
    )


def test_initialization_missing_time_delta():
    with pytest.raises(
        TypeError, match=r"'time' is scalar, so 'time_delta' is required"
    ):
        TimeSeries(
            time_start=datetime(2018, 7, 1, 10, 10, 10),
            data=[[10, 2, 3], [4, 5, 6]],
            names=["a", "b"],
        )


def test_initialization_invalid_time_and_time_start():
    with pytest.raises(TypeError, match=r"Cannot specify both 'time' and 'time_start'"):
        TimeSeries(
            time=INPUT_TIME,
            time_start=datetime(2018, 7, 1, 10, 10, 10),
            data=[[10, 2, 3], [4, 5, 6]],
            names=["a", "b"],
        )


def test_initialization_invalid_time_delta():
    with pytest.raises(
        TypeError, match=r"'time_delta' should be a Quantity or a TimeDelta"
    ):
        TimeSeries(
            time_start=datetime(2018, 7, 1, 10, 10, 10),
            time_delta=[1, 4, 3],
            data=[[10, 2, 3], [4, 5, 6]],
            names=["a", "b"],
        )


def test_initialization_with_time_in_data():
    data = PLAIN_TABLE.copy()
    data["time"] = INPUT_TIME

    ts1 = TimeSeries(data=data)

    assert set(ts1.colnames) == {"time", "a", "b", "c"}
    assert all(ts1.time == INPUT_TIME)

    ts2 = TimeSeries(data=[[10, 2, 3], INPUT_TIME], names=["a", "time"])
    assert set(ts2.colnames) == {"time", "a"}
    assert all(ts2.time == INPUT_TIME)

    MESSAGE = r"'time' has been given both in the table and as a keyword argument"

    with pytest.raises(TypeError, match=MESSAGE):
        # Don't allow ambiguous cases of passing multiple 'time' columns
        TimeSeries(data=data, time=INPUT_TIME)

    with pytest.raises(TypeError, match=MESSAGE):
        # 'time' is a protected name, don't allow ambiguous cases
        TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], INPUT_TIME], names=["a", "time"])


def test_initialization_n_samples():
    # Make sure things crash with incorrect n_samples

    with pytest.raises(
        TypeError,
        match=(
            r"'n_samples' has been given both and it is not the same length as the"
            r" input data."
        ),
    ):
        TimeSeries(time=INPUT_TIME, data=PLAIN_TABLE, n_samples=1000)


def test_initialization_length_mismatch():
    with pytest.raises(
        ValueError, match=r"Length of 'time' \(3\) should match data length \(2\)"
    ):
        TimeSeries(time=INPUT_TIME, data=[[10, 2], [4, 5]], names=["a", "b"])


def test_initialization_invalid_both_time_and_time_delta():
    with pytest.raises(
        TypeError,
        match=r"'time_delta' should not be specified since 'time' is an array",
    ):
        TimeSeries(time=INPUT_TIME, time_delta=TimeDelta(3, format="sec"))


def test_fold():
    times = Time([1, 2, 3, 8, 9, 12], format="unix")

    ts = TimeSeries(time=times)
    ts["flux"] = [1, 4, 4, 3, 2, 3]

    # Try without epoch time, as it should default to the first time and
    # wrapping at half the period.
    tsf = ts.fold(period=3.2 * u.s)
    assert isinstance(tsf.time, TimeDelta)
    assert_allclose(tsf.time.sec, [0, 1, -1.2, 0.6, -1.6, 1.4], rtol=1e-6)

    # Try with epoch time
    tsf = ts.fold(period=3.2 * u.s, epoch_time=Time(1.6, format="unix"))
    assert isinstance(tsf.time, TimeDelta)
    assert_allclose(tsf.time.sec, [-0.6, 0.4, 1.4, 0.0, 1.0, 0.8], rtol=1e-6, atol=1e-6)

    # Now with wrap_phase set to the full period
    tsf = ts.fold(period=3.2 * u.s, wrap_phase=3.2 * u.s)
    assert isinstance(tsf.time, TimeDelta)
    assert_allclose(tsf.time.sec, [0, 1, 2, 0.6, 1.6, 1.4], rtol=1e-6)

    # Now set epoch_phase to be 1/4 of the way through the phase
    tsf = ts.fold(period=3.2 * u.s, epoch_phase=0.8 * u.s)
    assert isinstance(tsf.time, TimeDelta)
    assert_allclose(tsf.time.sec, [0.8, -1.4, -0.4, 1.4, -0.8, -1.0], rtol=1e-6)

    # And combining epoch_phase and wrap_phase
    tsf = ts.fold(period=3.2 * u.s, epoch_phase=0.8 * u.s, wrap_phase=3.2 * u.s)
    assert isinstance(tsf.time, TimeDelta)
    assert_allclose(tsf.time.sec, [0.8, 1.8, 2.8, 1.4, 2.4, 2.2], rtol=1e-6)

    # Now repeat the above tests but with normalization applied

    # Try without epoch time, as it should default to the first time and
    # wrapping at half the period.
    tsf = ts.fold(period=3.2 * u.s, normalize_phase=True)
    assert isinstance(tsf.time, Quantity)
    assert_allclose(
        tsf.time.to_value(u.one),
        [0, 1 / 3.2, -1.2 / 3.2, 0.6 / 3.2, -1.6 / 3.2, 1.4 / 3.2],
        rtol=1e-6,
    )

    # Try with epoch time
    tsf = ts.fold(
        period=3.2 * u.s, epoch_time=Time(1.6, format="unix"), normalize_phase=True
    )
    assert isinstance(tsf.time, Quantity)
    assert_allclose(
        tsf.time.to_value(u.one),
        [-0.6 / 3.2, 0.4 / 3.2, 1.4 / 3.2, 0.0 / 3.2, 1.0 / 3.2, 0.8 / 3.2],
        rtol=1e-6,
        atol=1e-6,
    )

    # Now with wrap_phase set to the full period
    tsf = ts.fold(period=3.2 * u.s, wrap_phase=1, normalize_phase=True)
    assert isinstance(tsf.time, Quantity)
    assert_allclose(
        tsf.time.to_value(u.one),
        [0, 1 / 3.2, 2 / 3.2, 0.6 / 3.2, 1.6 / 3.2, 1.4 / 3.2],
        rtol=1e-6,
    )

    # Now set epoch_phase to be 1/4 of the way through the phase
    tsf = ts.fold(period=3.2 * u.s, epoch_phase=0.25, normalize_phase=True)
    assert isinstance(tsf.time, Quantity)
    assert_allclose(
        tsf.time.to_value(u.one),
        [0.8 / 3.2, -1.4 / 3.2, -0.4 / 3.2, 1.4 / 3.2, -0.8 / 3.2, -1.0 / 3.2],
        rtol=1e-6,
    )

    # And combining epoch_phase and wrap_phase
    tsf = ts.fold(
        period=3.2 * u.s, epoch_phase=0.25, wrap_phase=1, normalize_phase=True
    )
    assert isinstance(tsf.time, Quantity)
    assert_allclose(
        tsf.time.to_value(u.one),
        [0.8 / 3.2, 1.8 / 3.2, 2.8 / 3.2, 1.4 / 3.2, 2.4 / 3.2, 2.2 / 3.2],
        rtol=1e-6,
    )


def test_fold_invalid_options():
    times = Time([1, 2, 3, 8, 9, 12], format="unix")

    ts = TimeSeries(time=times)
    ts["flux"] = [1, 4, 4, 3, 2, 3]

    with pytest.raises(
        u.UnitsError, match="period should be a Quantity in units of time"
    ):
        ts.fold(period=3.2)

    with pytest.raises(
        u.UnitsError, match="period should be a Quantity in units of time"
    ):
        ts.fold(period=3.2 * u.m)

    with pytest.raises(
        u.UnitsError,
        match=(
            "epoch_phase should be a Quantity in units of "
            "time when normalize_phase=False"
        ),
    ):
        ts.fold(period=3.2 * u.s, epoch_phase=0.2)

    with pytest.raises(
        u.UnitsError,
        match=(
            "epoch_phase should be a dimensionless Quantity "
            "or a float when normalize_phase=True"
        ),
    ):
        ts.fold(period=3.2 * u.s, epoch_phase=0.2 * u.s, normalize_phase=True)

    with pytest.raises(
        u.UnitsError,
        match=(
            "wrap_phase should be a Quantity in units of "
            "time when normalize_phase=False"
        ),
    ):
        ts.fold(period=3.2 * u.s, wrap_phase=0.2)

    with pytest.raises(
        u.UnitsError,
        match="wrap_phase should be dimensionless when normalize_phase=True",
    ):
        ts.fold(period=3.2 * u.s, wrap_phase=0.2 * u.s, normalize_phase=True)

    with pytest.raises(
        ValueError, match="wrap_phase should be between 0 and the period"
    ):
        ts.fold(period=3.2 * u.s, wrap_phase=-0.1 * u.s)

    with pytest.raises(
        ValueError, match="wrap_phase should be between 0 and the period"
    ):
        ts.fold(period=3.2 * u.s, wrap_phase=-4.2 * u.s)

    with pytest.raises(ValueError, match="wrap_phase should be between 0 and 1"):
        ts.fold(period=3.2 * u.s, wrap_phase=-0.1, normalize_phase=True)

    with pytest.raises(ValueError, match="wrap_phase should be between 0 and 1"):
        ts.fold(period=3.2 * u.s, wrap_phase=2.2, normalize_phase=True)


def test_pandas():
    pandas = pytest.importorskip("pandas")

    df1 = pandas.DataFrame()
    df1["a"] = [1, 2, 3]
    df1.set_index(pandas.DatetimeIndex(INPUT_TIME.datetime64), inplace=True)

    ts = TimeSeries.from_pandas(df1)
    assert_equal(ts.time.isot, INPUT_TIME.isot)
    assert ts.colnames == ["time", "a"]
    assert len(ts.indices) == 1
    assert (ts.indices["time"].columns[0] == INPUT_TIME).all()

    ts_tcb = TimeSeries.from_pandas(df1, time_scale="tcb")
    assert ts_tcb.time.scale == "tcb"

    df2 = ts.to_pandas()
    assert (df2.index.values == pandas.Index(INPUT_TIME.datetime64).values).all()
    assert df2.columns == pandas.Index(["a"])
    assert (df1["a"] == df2["a"]).all()

    with pytest.raises(TypeError, match=r"Input should be a pandas DataFrame"):
        TimeSeries.from_pandas(None)

    df4 = pandas.DataFrame()
    df4["a"] = [1, 2, 3]

    with pytest.raises(TypeError, match=r"DataFrame does not have a DatetimeIndex"):
        TimeSeries.from_pandas(df4)


def test_read_time_missing():
    with pytest.raises(
        ValueError,
        match=(
            r"``time_column`` should be provided since the default Table readers are"
            r" being used\."
        ),
    ):
        TimeSeries.read(CSV_FILE, format="csv")


def test_read_time_wrong():
    with pytest.raises(
        ValueError, match=r"Time column 'abc' not found in the input data\."
    ):
        TimeSeries.read(CSV_FILE, time_column="abc", format="csv")


def test_read():
    timeseries = TimeSeries.read(CSV_FILE, time_column="Date", format="csv")
    assert timeseries.colnames == ["time", "A", "B", "C", "D", "E", "F", "G"]
    assert len(timeseries) == 11
    assert timeseries["time"].format == "iso"
    assert timeseries["A"].sum() == 266.5


@pytest.mark.remote_data(source="astropy")
def test_kepler_astropy():
    from astropy.units import UnitsWarning

    filename = get_pkg_data_filename("timeseries/kplr010666592-2009131110544_slc.fits")

    with pytest.warns(UnitsWarning):
        timeseries = TimeSeries.read(filename, format="kepler.fits")

    assert timeseries["time"].format == "isot"
    assert timeseries["time"].scale == "tdb"
    assert timeseries["sap_flux"].unit.to_string() == "electron / s"
    assert len(timeseries) == 14280
    assert len(timeseries.columns) == 20


@pytest.mark.remote_data(source="astropy")
def test_tess_astropy():
    filename = get_pkg_data_filename(
        "timeseries/hlsp_tess-data-alerts_tess_phot_00025155310-s01_tess_v1_lc.fits"
    )
    with pytest.warns((UserWarning, UnitsWarning)) as record:
        timeseries = TimeSeries.read(filename, format="tess.fits")
    # we might hit some warnings more than once, but the exact sequence probably
    # does not matter too much, so we'll just try to match the *set* of unique warnings
    unique_warnings = {(wm.category, wm.message.args[0]) for wm in record}
    assert len(record) != len(unique_warnings)
    expected = {
        (
            UnitsWarning,
            "'BJD - 2457000, days' did not parse as fits unit: "
            "At col 0, Unit 'BJD' not supported by the FITS standard.  "
            "If this is meant to be a custom unit, define it with 'u.def_unit'. "
            "To have it recognized inside a file reader or "
            "other code, enable it with 'u.add_enabled_units'. For details, see "
            "https://docs.astropy.org/en/latest/units/combining_and_defining.html",
        ),
        (
            UnitsWarning,
            "'e-/s' did not parse as fits unit: "
            "At col 0, Unit 'e' not supported by the FITS standard.  "
            "If this is meant to be a custom unit, define it with 'u.def_unit'. "
            "To have it recognized inside a file reader or other code, "
            "enable it with 'u.add_enabled_units'. For details, see "
            "https://docs.astropy.org/en/latest/units/combining_and_defining.html",
        ),
        (
            UnitsWarning,
            "'pixels' did not parse as fits unit: "
            "At col 0, Unit 'pixels' not supported by the FITS standard. "
            "Did you mean pixel? "
            "If this is meant to be a custom unit, define it with 'u.def_unit'. "
            "To have it recognized inside a file "
            "reader or other code, enable it with 'u.add_enabled_units'. For details, "
            "see https://docs.astropy.org/en/latest/units/combining_and_defining.html",
        ),
        (UserWarning, "Ignoring 815 rows with NaN times"),
    }
    assert unique_warnings == expected, (
        f"Got some unexpected warnings\n{unique_warnings - expected}"
    )

    assert timeseries["time"].format == "isot"
    assert timeseries["time"].scale == "tdb"
    assert timeseries["sap_flux"].unit.to_string() == "electron / s"
    assert len(timeseries) == 19261
    assert len(timeseries.columns) == 20


def test_required_columns():
    # Test the machinery that makes sure that the required columns are present

    ts = TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=["a", "b"])

    # In the examples below, the operation (e.g. remove_column) is actually
    # carried out before the checks are made, so we need to use copy() so that
    # we don't change the main version of the time series.

    # Make sure copy works fine
    ts.copy()

    MESSAGE = (
        r"TimeSeries object is invalid - expected 'time' as the first column but found"
        r" '{}'"
    )

    with pytest.raises(ValueError, match=MESSAGE.format("c")):
        ts.copy().add_column(Column([3, 4, 5], name="c"), index=0)

    with pytest.raises(ValueError, match=MESSAGE.format("d")):
        ts.copy().add_columns(
            [Column([3, 4, 5], name="d"), Column([3, 4, 5], name="e")], indexes=[0, 1]
        )

    with pytest.raises(ValueError, match=MESSAGE.format("a")):
        ts.copy().keep_columns(["a", "b"])

    with pytest.raises(ValueError, match=MESSAGE.format("a")):
        ts.copy().remove_column("time")

    with pytest.raises(ValueError, match=MESSAGE.format("b")):
        ts.copy().remove_columns(["time", "a"])

    with pytest.raises(ValueError, match=MESSAGE.format("banana")):
        ts.copy().rename_column("time", "banana")

    # https://github.com/astropy/astropy/issues/13009
    MESSAGE = (
        r"TimeSeries object is invalid - expected \['time', 'a'\] as the first columns"
        r" but found \['time', 'b'\]"
    )
    ts_2cols_required = ts.copy()
    ts_2cols_required._required_columns = ["time", "a"]
    with pytest.raises(ValueError, match=MESSAGE):
        ts_2cols_required.remove_column("a")


@pytest.mark.parametrize("cls", [BoxLeastSquares, LombScargle])
def test_periodogram(cls):
    # Note that we don't need to check the actual results from the periodogram
    # classes here since these are tested extensively in
    # astropy.timeseries.periodograms.

    ts = TimeSeries(time=INPUT_TIME, data=[[10, 2, 3], [4, 5, 6]], names=["a", "b"])

    p1 = cls.from_timeseries(ts, "a")
    assert isinstance(p1, cls)
    assert_allclose(p1.t.jd, ts.time.jd)
    assert_equal(p1.y, ts["a"])
    assert p1.dy is None

    p2 = cls.from_timeseries(ts, "a", uncertainty="b")
    assert_quantity_allclose(p2.dy, ts["b"])

    p3 = cls.from_timeseries(ts, "a", uncertainty=0.1)
    assert_allclose(p3.dy, 0.1)
