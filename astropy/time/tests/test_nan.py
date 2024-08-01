import numpy as np
import pytest

import astropy.units as u
from astropy.time import Time, TimeDelta

FLOAT_FORMATS = [
    "byear",
    "cxcsec",
    "decimalyear",
    "gps",
    "jd",
    "jyear",
    "mjd",
    "plot_date",
    "stardate",
    "unix",
    "unix_tai",
]

STRING_FORMATS = {
    "byear_str": "B1950.0",
    "jyear_str": "J2000.0",
    "iso": "2020-01-01 20:00:00.5",
    "isot": "2020-01-01T20:00:00.5",
    "yday": "2020:300:20:00:05",
}


def test_time_creation():
    """Test a time scalar can be created as nan"""
    t = Time(np.nan, format="jd")  # for now, in principle, it shouldn't matter
    assert np.isnan(t)
    assert np.isnat(t)

    ts = Time([0, 1, np.nan, 2, 3], format="jd")
    assert np.all(np.isnan(ts) == [False, False, True, False, False])


def test_timedelta_creation():
    """Test a timedelta scalar can be created as nan"""
    t = TimeDelta(np.nan)

    t = TimeDelta(np.nan * u.s)
    assert np.isnan(t)
    assert np.isnat(t)

    assert np.isnan(t.to_value(u.h))
    assert np.isnan(t.to_value(u.s))

    ts = TimeDelta([0, 1, np.nan, 2, 3] * u.s)
    assert np.all(np.isnan(ts) == [False, False, True, False, False])


def test_nan_behaviour():
    """Test if times follow the usual nan behaviour"""

    valid_time = Time("2020-01-01T00:00")
    valid_delta = TimeDelta(1 * u.hour)
    nan_time = Time(np.nan, format="jd")
    nan_delta = TimeDelta(np.nan)
    valid_quantity = 5 * u.hour
    nan_quantity = np.nan * u.s

    assert np.isnan(valid_time + nan_delta)
    assert np.isnan(valid_time + nan_quantity)

    assert np.isnan(nan_time + valid_delta)
    assert np.isnan(nan_time + valid_quantity)

    assert np.isnan(nan_time + nan_quantity)
    assert np.isnan(nan_time + nan_delta)

    assert np.isnan(nan_time - valid_time)
    assert np.isnan(valid_time - nan_time)

    assert Time(np.nan, format="jd") != Time(np.nan, format="jd")
    assert TimeDelta(np.nan) != TimeDelta(np.nan)


@pytest.mark.parametrize("fmt", FLOAT_FORMATS)
def test_float_formats(fmt):
    """Test a time scalar can be created as nan"""
    t = Time(np.nan, format=fmt)
    assert np.isnat(t)

    t = Time([np.nan, 2000.0], format=fmt)
    assert np.all(np.isnat(t) == [True, False])


@pytest.mark.parametrize("fmt,example", STRING_FORMATS.items())
def test_string_formats(fmt, example):
    """Test a time scalar can be created as nan"""
    t = Time("NaT", format=fmt)
    assert np.isnat(t)

    t = Time(["NaT", example], format=fmt)
    assert np.all(np.isnat(t) == [True, False])


@pytest.mark.parametrize("func", (np.isnan, np.isnat))
def test_masked(func):
    """Test isnan/isnat behaviour on masked Time"""
    time = Time(np.ma.masked_array(np.nan, mask=True), format="mjd")
    assert func(time).mask

    mjd = np.ma.masked_array([np.nan, np.nan, 60e3], mask=[False, True, False])
    time = Time(mjd, format="mjd")
    expected = np.ma.masked_array([True, True, False], mask=time.mask)
    np.testing.assert_array_equal(func(time), expected)
