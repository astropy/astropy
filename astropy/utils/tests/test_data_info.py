# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.table.index import SlicedIndex
from astropy.time import Time
from astropy.utils.data_info import dtype_info_name

STRING_TYPE_NAMES = {(True, "S"): "bytes", (True, "U"): "str"}

DTYPE_TESTS = (
    (np.array(b"abcd").dtype, STRING_TYPE_NAMES[(True, "S")] + "4"),
    (np.array("abcd").dtype, STRING_TYPE_NAMES[(True, "U")] + "4"),
    ("S4", STRING_TYPE_NAMES[(True, "S")] + "4"),
    ("U4", STRING_TYPE_NAMES[(True, "U")] + "4"),
    (np.void, "void"),
    (np.int32, "int32"),
    (bool, "bool"),
    (float, "float64"),
    ("<f4", "float32"),
    ("u8", "uint64"),
    ("c16", "complex128"),
    ("object", "object"),
)


@pytest.mark.parametrize("input,output", DTYPE_TESTS)
def test_dtype_info_name(input, output):
    """
    Test that dtype_info_name is giving the expected output

    Here the available types::

      'b' boolean
      'i' (signed) integer
      'u' unsigned integer
      'f' floating-point
      'c' complex-floating point
      'O' (Python) objects
      'S', 'a' (byte-)string
      'U' Unicode
      'V' raw data (void)
    """
    assert dtype_info_name(input) == output


def test_info_no_copy_numpy():
    """Test that getting a single item from Table column object does not copy info.
    See #10889.
    """
    col = [1, 2]
    t = QTable([col], names=["col"])
    t.add_index("col")
    val = t["col"][0]
    # Returns a numpy scalar (e.g. np.float64) with no .info
    assert isinstance(val, np.number)
    with pytest.raises(AttributeError):
        val.info
    val = t["col"][:]
    assert val.info.indices == []


cols = [[1, 2] * u.m, Time([1, 2], format="cxcsec")]


@pytest.mark.parametrize("col", cols)
def test_info_no_copy_mixin_with_index(col):
    """Test that getting a single item from Table column object does not copy info.
    See #10889.
    """
    t = QTable([col], names=["col"])
    t.add_index("col")
    val = t["col"][0]
    assert "info" not in val.__dict__
    assert val.info.indices == []
    val = t["col"][:]
    assert "info" in val.__dict__
    assert val.info.indices == []
    val = t[:]["col"]
    assert "info" in val.__dict__
    assert isinstance(val.info.indices[0], SlicedIndex)


def test_info_no_copy_skycoord():
    """Test that getting a single item from Table SkyCoord column object does
    not copy info.  Cannot create an index on a SkyCoord currently.
    """
    col = (SkyCoord([1, 2], [1, 2], unit="deg"),)
    t = QTable([col], names=["col"])
    val = t["col"][0]
    assert "info" not in val.__dict__
    assert val.info.indices == []
    val = t["col"][:]
    assert val.info.indices == []
    val = t[:]["col"]
    assert val.info.indices == []
