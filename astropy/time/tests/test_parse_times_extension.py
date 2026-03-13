# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for the _parse_times C extension in astropy.time."""

import numpy as np
import pytest

from astropy.time import TimeISO, TimeISOT, TimeYearDayTime, _parse_times


def _make_pars(fast_parser_pars):
    """Build a dt_pars structured array from a TimeString fast_parser_pars dict."""
    p = fast_parser_pars
    pars = np.zeros(7, dtype=_parse_times.dt_pars)
    pars["delim"] = np.asarray(p["delims"], dtype=np.uint8).view("S1")
    pars["start"] = p["starts"]
    pars["stop"] = p["stops"]
    pars["break_allowed"] = p["break_allowed"]
    return pars


# ---------------------------------------------------------------------------
# create_parser: invalid input
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("size", [0, 6, 8], ids=["empty", "too_small", "too_large"])
def test_create_parser_wrong_size_raises(size):
    """create_parser must raise ValueError for parameter arrays != 7 entries."""
    bad_pars = np.zeros(size, dtype=_parse_times.dt_pars)
    with pytest.raises(ValueError, match="Parameter array must have 7 entries"):
        _parse_times.create_parser(bad_pars)


# ---------------------------------------------------------------------------
# create_parser: valid input returns a ufunc
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "fmt,pars",
    [
        ("iso", TimeISO.fast_parser_pars),
        ("yday", TimeYearDayTime.fast_parser_pars),
        ("isot", TimeISOT.fast_parser_pars),
    ],
)
def test_create_parser_returns_ufunc(fmt, pars):
    """create_parser with valid pars returns a numpy ufunc."""
    parser = _parse_times.create_parser(_make_pars(pars))
    assert isinstance(parser, np.ufunc)


# ---------------------------------------------------------------------------
# Parser output: parse a known ISO string
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "time_str,expected",
    [
        ("2000-01-01 12:00:00.000", (2000, 1, 1, 12, 0, 0.0)),
        ("1981-12-31 12:13:14.000", (1981, 12, 31, 12, 13, 14.0)),
    ],
)
def test_iso_parser_parses_string(time_str, expected):
    """ISO parser correctly parses known date strings."""
    pars = _make_pars(TimeISO.fast_parser_pars)
    parser = _parse_times.create_parser(pars)
    val1 = np.array([time_str], dtype="S24")
    chars = val1.view((_parse_times.dt_u1, val1.dtype.itemsize))
    result = parser(chars)
    year, month, day, hour, minute, second = expected
    assert result["year"][0] == year
    assert result["month"][0] == month
    assert result["day"][0] == day
    assert result["hour"][0] == hour
    assert result["minute"][0] == minute
    assert result["second"][0] == pytest.approx(second)
