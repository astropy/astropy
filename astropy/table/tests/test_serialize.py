# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for functions in astropy.table.serialize"""

import numpy as np
import pytest

from astropy.table import QTable, Table
from astropy.table.serialize import represent_nd_columns_as_1d_columns
from astropy.time import Time

# ---------------------------------------------------------------------------
# Tests of represent_nd_columns_as_1d_columns()
# ---------------------------------------------------------------------------


@pytest.fixture
def tbl_1d():
    """Table with only 1-D columns."""
    return Table({"a": [1, 2, 3], "b": [4.0, 5.0, 6.0]})


@pytest.fixture
def tbl_2d():
    """Table with a 2-D column (shape (3, 2)) and a 1-D column."""
    return Table({"data": np.arange(6).reshape(3, 2), "idx": [10, 20, 30]})


@pytest.fixture
def tbl_3d():
    """Table with a 3-D column (shape (3, 2, 1)) and a 1-D column."""
    return Table({"cube": np.arange(6).reshape(3, 2, 1), "idx": [1, 2, 3]})


def test_all_1d_copy_false_returns_same_object(tbl_1d):
    """copy=False with all 1-D columns → the original table object is returned."""
    out = represent_nd_columns_as_1d_columns(tbl_1d, copy=False)
    assert out is tbl_1d


def test_all_1d_copy_true_returns_copy(tbl_1d):
    """copy=True with all 1-D columns → a new table with equal contents."""
    out = represent_nd_columns_as_1d_columns(tbl_1d, copy=True)
    assert out is not tbl_1d
    assert out.colnames == tbl_1d.colnames
    np.testing.assert_array_equal(out["a"], tbl_1d["a"])
    np.testing.assert_array_equal(out["b"], tbl_1d["b"])


def test_3d_column_names(tbl_3d):
    """A (3, 2, 1) column 'cube' yields 'cube.0_0' and 'cube.1_0'."""
    out = represent_nd_columns_as_1d_columns(tbl_3d)
    assert out.colnames == ["cube.0_0", "cube.1_0", "idx"]


def test_3d_column_values(tbl_3d):
    """Flattened 3-D slices carry the correct values."""
    out = represent_nd_columns_as_1d_columns(tbl_3d)
    np.testing.assert_array_equal(out["cube.0_0"], tbl_3d["cube"][:, 0, 0])
    np.testing.assert_array_equal(out["cube.1_0"], tbl_3d["cube"][:, 1, 0])


def test_output_all_1d(tbl_3d):
    """Every column in the output must be 1-D."""
    out = represent_nd_columns_as_1d_columns(tbl_3d)
    assert all(len(col.shape) == 1 for col in out.itercols())


def test_collision_one_underscore():
    """
    'a' (2-D) would produce 'a.0' / 'a.1', but 'a.0' is already a 1-D column.
    Base name should escalate to 'a_', giving 'a_.0' and 'a_.1'.
    """
    tbl = Table({"a": np.arange(6).reshape(3, 2), "a.0": [10, 20, 30]})
    out = represent_nd_columns_as_1d_columns(tbl)
    # Original 1-D column is preserved unchanged
    assert "a.0" in out.colnames
    np.testing.assert_array_equal(out["a.0"], [10, 20, 30])
    # Flattened N-D sub-columns use the escalated base name
    assert "a_.0" in out.colnames
    assert "a_.1" in out.colnames
    np.testing.assert_array_equal(out["a_.0"], [0, 2, 4])
    np.testing.assert_array_equal(out["a_.1"], [1, 3, 5])


def test_collision_two_underscores():
    """
    Both 'a.0' and 'a_.0' pre-exist → base escalates twice to 'a__'.
    """
    tbl = Table(
        {
            "a": np.arange(6).reshape(3, 2),
            "a.0": [1, 2, 3],
            "a_.0": [4, 5, 6],
        }
    )
    out = represent_nd_columns_as_1d_columns(tbl)
    assert "a.0" in out.colnames
    assert "a_.0" in out.colnames
    assert "a__.0" in out.colnames
    assert "a__.1" in out.colnames
    # No duplicate names
    assert len(set(out.colnames)) == len(out.colnames)


def test_collision_between_two_nd_columns():
    """
    Two N-D columns where the second column's name ('a.0') would clash with
    the flattened output of the first → output has no duplicate names.
    """
    tbl = Table(
        {
            "a": np.arange(6).reshape(3, 2),
            "a.0": np.arange(6).reshape(3, 2),
        }
    )
    out = represent_nd_columns_as_1d_columns(tbl)
    assert len(set(out.colnames)) == len(out.colnames)


def test_time_mixin_3d():
    """
    A 3-D Time column of shape (3, 2, 1) should flatten to two 1-D Time
    columns named 'time.0_0' and 'time.1_0', with correct values.
    """
    time_data = np.arange(6).reshape(3, 2, 1)
    tbl = Table({"time": Time(time_data, format="cxcsec"), "idx": [1, 2, 3]})
    out = represent_nd_columns_as_1d_columns(tbl)

    assert out.colnames == ["time.0_0", "time.1_0", "idx"]

    # Values must match the corresponding slices of the original Time array
    np.testing.assert_array_equal(
        out["time.0_0"].value, Time(time_data[:, 0, 0], format="cxcsec").value
    )
    np.testing.assert_array_equal(
        out["time.1_0"].value, Time(time_data[:, 1, 0], format="cxcsec").value
    )
    # Result columns are 1-D
    assert len(out["time.0_0"].shape) == 1
    assert len(out["time.1_0"].shape) == 1

    assert isinstance(out["time.0_0"], Time)
    assert isinstance(out["time.1_0"], Time)


def test_table_subclass_preserved():
    """QTable input should produce QTable output."""
    tbl = QTable({"v": np.arange(6).reshape(3, 2), "n": [1, 2, 3]})
    out = represent_nd_columns_as_1d_columns(tbl)
    assert type(out) is QTable
