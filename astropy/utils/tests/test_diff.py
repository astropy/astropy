# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Some might be indirectly tested already in ``astropy.io.fits.tests``.
"""
import io

import numpy as np
import pytest

from astropy.table import Table
from astropy.utils.diff import diff_values, report_diff_values, where_not_allclose


@pytest.mark.parametrize("a", [np.nan, np.inf, 1.11, 1, "a"])
def test_diff_values_false(a):
    assert not diff_values(a, a)


@pytest.mark.parametrize(
    ("a", "b"), [(np.inf, np.nan), (1.11, 1.1), (1, 2), (1, "a"), ("a", "b")]
)
def test_diff_values_true(a, b):
    assert diff_values(a, b)


def test_float_comparison():
    """
    Regression test for https://github.com/spacetelescope/PyFITS/issues/21
    """
    f = io.StringIO()
    a = np.float32(0.029751372)
    b = np.float32(0.029751368)
    identical = report_diff_values(a, b, fileobj=f)
    assert not identical
    out = f.getvalue()

    # This test doesn't care about what the exact output is, just that it
    # did show a difference in their text representations
    assert "a>" in out
    assert "b>" in out


def test_diff_types():
    """
    Regression test for https://github.com/astropy/astropy/issues/4122
    """
    f = io.StringIO()
    a = 1.0
    b = "1.0"
    identical = report_diff_values(a, b, fileobj=f)
    assert not identical
    out = f.getvalue()
    # fmt: off
    assert out == (
        "  (float) a> 1.0\n"
        "    (str) b> '1.0'\n"
        "           ? +   +\n"
    )
    # fmt: on


def test_diff_numeric_scalar_types():
    """Test comparison of different numeric scalar types."""
    f = io.StringIO()
    assert not report_diff_values(1.0, 1, fileobj=f)
    out = f.getvalue()
    assert out == "  (float) a> 1.0\n    (int) b> 1\n"


def test_array_comparison():
    """
    Test diff-ing two arrays.
    """
    f = io.StringIO()
    a = np.arange(9).reshape(3, 3)
    b = a + 1
    identical = report_diff_values(a, b, fileobj=f)
    assert not identical
    out = f.getvalue()
    assert (
        out == "  at [0, 0]:\n"
        "    a> 0\n"
        "    b> 1\n"
        "  at [0, 1]:\n"
        "    a> 1\n"
        "    b> 2\n"
        "  at [0, 2]:\n"
        "    a> 2\n"
        "    b> 3\n"
        "  ...and at 6 more indices.\n"
    )


def test_diff_shaped_array_comparison():
    """
    Test diff-ing two differently shaped arrays.
    """
    f = io.StringIO()
    a = np.empty((1, 2, 3))
    identical = report_diff_values(a, a[0], fileobj=f)
    assert not identical
    out = f.getvalue()
    assert (
        out
        == "  Different array shapes:\n    a> (1, 2, 3)\n     ?  ---\n    b> (2, 3)\n"
    )


def test_tablediff():
    """
    Test diff-ing two simple Table objects.
    """
    a = Table.read(
        """name    obs_date    mag_b  mag_v
M31     2012-01-02  17.0   16.0
M82     2012-10-29  16.2   15.2
M101    2012-10-31  15.1   15.5""",
        format="ascii",
    )
    b = Table.read(
        """name    obs_date    mag_b  mag_v
M31     2012-01-02  17.0   16.5
M82     2012-10-29  16.2   15.2
M101    2012-10-30  15.1   15.5
NEW     2018-05-08   nan    9.0""",
        format="ascii",
    )
    f = io.StringIO()
    identical = report_diff_values(a, b, fileobj=f)
    assert not identical
    out = f.getvalue()
    assert (
        out == "     name  obs_date  mag_b mag_v\n"
        "     ---- ---------- ----- -----\n"
        "  a>  M31 2012-01-02  17.0  16.0\n"
        "   ?                           ^\n"
        "  b>  M31 2012-01-02  17.0  16.5\n"
        "   ?                           ^\n"
        "      M82 2012-10-29  16.2  15.2\n"
        "  a> M101 2012-10-31  15.1  15.5\n"
        "   ?               ^\n"
        "  b> M101 2012-10-30  15.1  15.5\n"
        "   ?               ^\n"
        "  b>  NEW 2018-05-08   nan   9.0\n"
    )

    # Identical
    assert report_diff_values(a, a, fileobj=f)


def test_large_table_diff():
    # see https://github.com/astropy/astropy/issues/14010
    colnames = [f"column{i}" for i in range(100)]
    t1 = Table(names=colnames)

    colnames.insert(50, "test")
    t2 = Table(names=colnames)

    assert not report_diff_values(t1, t2, fileobj=io.StringIO())


@pytest.mark.parametrize("kwargs", [{}, {"atol": 0, "rtol": 0}])
def test_where_not_allclose(kwargs):
    a = np.array([1, np.nan, np.inf, 4.5])
    b = np.array([1, np.inf, np.nan, 4.6])

    assert where_not_allclose(a, b, **kwargs) == ([3],)
    assert len(where_not_allclose(a, a, **kwargs)[0]) == 0
