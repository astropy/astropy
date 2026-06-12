# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Tests for the Cython extension module ``astropy/table/_column_mixins.pyx``."""

import numpy as np
import pytest
from numpy import ma

from astropy.table._column_mixins import _ColumnGetitemShim, _MaskedColumnGetitemShim

# ---------------------------------------------------------------------------
# Minimal subclasses used by the direct / isolated tests
# ---------------------------------------------------------------------------


class MinimalColumn(_ColumnGetitemShim, np.ndarray):
    """Wrapper around ``_ColumnGetitemShim`` which serves its isolated testing"""

    def __new__(cls, data, dtype=None):
        return np.asarray(data, dtype=dtype).view(cls)

    @property
    def data(self):
        """Return a plain ndarray view (required by ``base_getitem``)."""
        return self.view(np.ndarray)


class MinimalMaskedColumn(_MaskedColumnGetitemShim, ma.MaskedArray):
    """Wrapper around ``_ColumnGetitemShim`` which serves its isolated testing"""

    def __new__(cls, data, mask=None, dtype=None):
        return ma.array(data, mask=mask, dtype=dtype).view(cls)

    @property
    def data(self):
        return self.view(ma.MaskedArray)

    def _copy_attrs_slice(self, value):
        """No extra attributes to copy in the minimal version."""
        return value


# ---------------------------------------------------------------------------
# Direct / isolated tests - no Column or MaskedColumn
# ---------------------------------------------------------------------------


class TestColumnGetitemShimDirect:
    """Test ``_ColumnGetitemShim`` through ``MinimalColumn``, completely
    independent of ``astropy.table.column``.
    """

    def test_multint_ind_ndarray(self):
        """The core purpose of the shim: integer-indexing a multi-dimensional
        subclass must return a plain ``ndarray``, not the subclass.
        """

        col = MinimalColumn([[7, 8, 9], [10, 11, 12], [13, 14, 15]])
        result = col[1]

        assert isinstance(result, np.ndarray)
        assert not isinstance(result, MinimalColumn), (
            "Integer indexing a multi-dimensional column must NOT return the subclass."
        )
        np.testing.assert_array_equal(result, [10, 11, 12])

    @pytest.mark.parametrize(
        "array",
        [
            [[7, 8, 9, 10]],
            [[7, 8, 9], [10, 11, 12], [13, 14, 15]],
            [[], [], []],
        ],
    )
    def test_multdim_int_ind(self, array):
        # Makes check on basic multi-dimensional column indexing facilities

        data = np.array(array)
        col = MinimalColumn(data)

        assert (col[0] == data[0]).all()
        assert (col[-1] == data[-1]).all()

        n = np.random.randint(0, len(data))
        assert (col[n] == data[n]).all()

        if data.ndim > 1 and data.shape[1] > 0:
            assert col[0][2] == data[0][2]

    def test_multdim_slice_subclass(self):
        """Slice indexing on a multi-dimensional shim array must return
        the subclass - only integer indexing strips it.
        """

        col = MinimalColumn([[7, 8, 9], [10, 11, 12], [13, 14, 15]])
        result = col[1:]

        # Without subscripting we should get a Column object
        assert isinstance(result, MinimalColumn)

    @pytest.mark.parametrize(
        "array",
        [
            [[7, 8, 9, 10]],
            [[]],
            [[7, 8, 9], [10, 11, 12], [13, 14, 15]],
            [[], [], []],
        ],
    )
    def test_multdim_slicing(self, array):
        data = np.array(array)
        col = MinimalColumn(data)
        n = len(col)

        assert (col[0:] == data[0:]).all()
        assert (col[-n:] == data[-n:]).all()

        if col.data.shape[0] > 1:
            assert (col[1 : n - 1] == data[1 : n - 1]).all()
            assert (col[1 : n - 1 : 2] == data[1 : n - 1 : 2]).all()
            assert (col[-n + 1 : -1 : -2] == data[-n + 1 : -1 : -2]).all()
            assert (col[1, ...] == data[1, ...]).all()

        assert (col[::2] == data[::2]).all()
        assert (col[::-2] == data[::-2]).all()
        assert (col[:] == data[:]).all()
        assert (col[...] == data[...]).all()

    def test_multidim_fancy_indexing(self):
        """Fancy indexing on a multi-dimensional column must return a copy not view."""

        data = np.array([[7, 8, 9], [4, 5, 2], [1, 5, 56]])
        col = MinimalColumn(data)

        assert (col[[0, 2]] == data[[0, 2]]).all()
        assert (col[:, [0, 2]] == data[:, [0, 2]]).all()

    def test_void_dtype_string_key(self):
        dt = np.dtype([("name", "U10"), ("H0", "f4")])
        data = np.array([("Planck18", 67.66), ("Planck13", 67.77)], dtype=dt)

        col = MinimalColumn(data)

        np.testing.assert_array_equal(col["name"], ["Planck18", "Planck13"])
        np.testing.assert_array_almost_equal(col["H0"], [67.66, 67.77], decimal=2)

    def test_byte_decoding(self):
        data = [b"Planck18", b"Planck13"]
        col = MinimalColumn(data)

        result = col[1]
        assert isinstance(result, str)
        assert result == "Planck13"

        result = col[0:]
        assert isinstance(result, MinimalColumn)

    def test_sigdim_int_ind(self):
        arr = MinimalColumn([10, 20, 30, 40])

        assert arr[0] == 10
        assert arr[-1] == 40

    def test_sigdim_slice(self):
        data = np.array([10, 20, 30, 40])
        arr = MinimalColumn(data)

        np.testing.assert_array_equal(arr[1:3], data[1:3])
        np.testing.assert_array_equal(arr[::-1], data[::-1])
        np.testing.assert_array_equal(arr[::2], data[::2])
        np.testing.assert_array_equal(arr[:], data[:])
        np.testing.assert_array_equal(arr[...], data[...])

    def test_sigdim_boolean_mask(self):
        """For ``dtype=object`` which is neither ``INTEGER_TYPES`` nor ``STRING_TYPES``, so
        ``base_getitem`` falls through to ``ndarray.__getitem__``.
        """

        arr = MinimalColumn([10, 20, 30, 40])
        mask = np.array([True, False, True, False])

        np.testing.assert_array_equal(arr[mask], [10, 30])

    def test_sigdim_fancy_indexing(self):
        arr = MinimalColumn([10, 20, 30, 40])

        np.testing.assert_array_equal(arr[[0, 3]], [10, 40])

    def test_out_of_bounds_raises_index_error(self):
        arr = MinimalColumn([1, 2, 3])

        with pytest.raises(IndexError):
            _ = arr[10]

        empty = MinimalColumn([])

        with pytest.raises(IndexError):
            _ = empty[0]


class TestMaskedColumnGetitemShimDirect:
    """Test ``_MaskedColumnGetitemShim`` through ``MinimalMaskedColumn``,
    independent of ``astropy.table.column``.
    """

    def test_sigdim_masked_scalar_access(self):
        """Tests indexing out masked and unmasked values from a masked-column"""

        mcol = MinimalMaskedColumn(
            [7, 8, 9, 10, 11], mask=[True, False, False, True, True]
        )
        assert mcol[0] is ma.masked
        assert mcol[1] == 8

    def test_multdim_int_ind_masked_array(self):
        mcol = MinimalMaskedColumn(
            [[7, 8, 9], [4, 5, 6], [4, 2, 3]],
            mask=[[True, False, False], [True, True, True], [False, False, False]],
        )
        row = mcol[0]
        assert isinstance(row, ma.MaskedArray)
        assert isinstance(row[0], ma.core.MaskedConstant)
        assert mcol[2][0] == 4

    def test_slice_preserves_mask(self):
        arr = MinimalMaskedColumn([1, 2, 3, 4], mask=[False, True, False, True])
        sliced = arr[1:3]
        np.testing.assert_array_equal(sliced.mask, [True, False])
        assert sliced[0] is ma.masked
        assert sliced[1] == 3
