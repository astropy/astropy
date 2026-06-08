# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for the _parse_times C extension in astropy.time."""

import dataclasses
from collections.abc import Mapping
from typing import Any

import numpy as np
import pytest

from astropy.time import TimeISO, TimeISOT, TimeYearDayTime, _parse_times


@dataclasses.dataclass(kw_only=True, slots=True, frozen=True)
class ExpectedTime:
    """Strictly typed expected output for the C-level time struct."""

    year: int
    month: int
    day: int
    hour: int
    minute: int
    second: float


@dataclasses.dataclass(kw_only=True, slots=True, frozen=True)
class ParserTestCase:
    """Dataclass parameterization for standard string parsing."""

    time_str: str
    expected: ExpectedTime


def _make_pars(
    fast_parser_pars: Mapping[str, Any],
) -> np.ndarray[Any, np.dtype[np.void]]:
    """
    Build a dt_pars structured array from a TimeString fast_parser_pars dict.
    Returns a strict 1D numpy array with a void (structured) dtype.
    """
    p = fast_parser_pars
    pars: np.ndarray[Any, np.dtype[np.void]] = np.zeros(7, dtype=_parse_times.dt_pars)
    pars["delim"] = np.asarray(p["delims"], dtype=np.uint8).view("S1")
    pars["start"] = p["starts"]
    pars["stop"] = p["stops"]
    pars["break_allowed"] = p["break_allowed"]
    return pars


# ---------------------------------------------------------------------------
# create_parser: invalid input
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "size",
    [
        pytest.param(0, id="empty_array"),
        pytest.param(6, id="undersized_array"),
        pytest.param(8, id="oversized_array"),
    ],
)
def test_create_parser_wrong_size_raises(size: int) -> None:
    """create_parser must raise ValueError for parameter arrays != 7 entries."""
    bad_pars: np.ndarray[Any, np.dtype[np.void]] = np.zeros(
        size, dtype=_parse_times.dt_pars
    )
    with pytest.raises(ValueError, match="Parameter array must have 7 entries"):
        _parse_times.create_parser(bad_pars)


# ---------------------------------------------------------------------------
# create_parser: valid input returns a ufunc
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "fmt_name, pars_dict",
    [
        pytest.param("iso", TimeISO.fast_parser_pars, id="fmt_iso"),
        pytest.param("yday", TimeYearDayTime.fast_parser_pars, id="fmt_yday"),
        pytest.param("isot", TimeISOT.fast_parser_pars, id="fmt_isot"),
    ],
)
def test_create_parser_returns_ufunc(
    fmt_name: str, pars_dict: Mapping[str, Any]
) -> None:
    """create_parser with valid pars returns a strictly typed numpy ufunc."""
    parser = _parse_times.create_parser(_make_pars(pars_dict))

    # Strict subclass assertion to avoid sequence preservation bypasses
    assert type(parser) is np.ufunc


# ---------------------------------------------------------------------------
# Parser output: parse a known ISO string
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "case",
    [
        pytest.param(
            ParserTestCase(
                time_str="2000-01-01 12:00:00.000",
                expected=ExpectedTime(
                    year=2000, month=1, day=1, hour=12, minute=0, second=0.0
                ),
            ),
            id="y2k_midnight_standard",
        ),
        pytest.param(
            ParserTestCase(
                time_str="1981-12-31 12:13:14.000",
                expected=ExpectedTime(
                    year=1981, month=12, day=31, hour=12, minute=13, second=14.0
                ),
            ),
            id="end_of_year_arbitrary_time",
        ),
    ],
)
def test_iso_parser_parses_string(case: ParserTestCase) -> None:
    """ISO parser correctly parses known date strings directly at the C boundary."""
    pars = _make_pars(TimeISO.fast_parser_pars)
    parser = _parse_times.create_parser(pars)

    # Isolate the memoryview buffers and raw byte arrays
    val1: np.ndarray[Any, np.dtype[np.bytes_]] = np.array([case.time_str], dtype="S24")
    chars: np.ndarray[Any, np.dtype[np.uint8]] = val1.view(
        (_parse_times.dt_u1, val1.dtype.itemsize)
    )

    result = parser(chars)

    assert result["year"][0] == case.expected.year
    assert result["month"][0] == case.expected.month
    assert result["day"][0] == case.expected.day
    assert result["hour"][0] == case.expected.hour
    assert result["minute"][0] == case.expected.minute
    assert result["second"][0] == pytest.approx(case.expected.second)
