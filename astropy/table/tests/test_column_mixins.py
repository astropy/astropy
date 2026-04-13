# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Tests for the Cython extension module ``astropy/table/_column_mixins.pyx``."""

import numpy as np
import pytest
from numpy import ma

from astropy.table._column_mixins import _ColumnGetitemShim, _MaskedColumnGetitemShim
from astropy.table.column import Column, MaskedColumn

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

    def test_multidim_integer_index_returns_ndarray(self):
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

    def test_multidim_slice_preserves_subclass(self):
        """Slice indexing on a multi-dimensional shim array must still return
        the subclass - only integer indexing strips it.
        """

        col = MinimalColumn([[7, 8, 9], [10, 11, 12], [13, 14, 15]])
        result = col[1:]

        # Without subscripting we should get a Column object
        assert isinstance(result, MinimalColumn)

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

    def test_bytes_dtype_slice_not_decoded(self):
        arr = MinimalColumn([b"Planck18", b"Planck13"], dtype="S10")

        result = arr[0:]

        assert result.dtype.kind == "S"

    # --- Regression Tests ------------------------------------------

    def test_sigdim_integer_index(self):
        arr = MinimalColumn([10, 20, 30, 40])

        assert arr[0] == 10
        assert arr[-1] == 40

    def test_sigdim_slice(self):
        data = np.array([10, 20, 30, 40])
        arr = MinimalColumn(data)

        np.testing.assert_array_equal(arr[1:3], data[1:3])
        np.testing.assert_array_equal(arr[::-1], data[::-1])
        np.testing.assert_array_equal(arr[::2], data[::2])

    def test_sigdim_boolean_mask(self):
        """A boolean array is neither ``INTEGER_TYPES`` nor ``STRING_TYPES``, so
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

    def test_masked_scalar_returns_masked_constant(self):
        arr = MinimalMaskedColumn([1, 2, 3], mask=[True, False, False])

        assert arr[0] is ma.masked

    def test_unmasked_scalar_returns_value(self):
        arr = MinimalMaskedColumn([10, 20, 30], mask=[False, False, True])
        assert arr[1] == 20

    def test_multidim_integer_index_returns_masked_array(self):
        arr = MinimalMaskedColumn(
            [[1, 2, 3], [4, 5, 6]],
            mask=[[True, False, False], [False, False, True]],
        )
        row = arr[0]
        assert isinstance(row, ma.MaskedArray)
        assert row[0] is ma.masked
        assert row[1] == 2

    def test_slice_preserves_mask(self):
        arr = MinimalMaskedColumn([1, 2, 3, 4], mask=[False, True, False, True])
        sliced = arr[1:3]
        np.testing.assert_array_equal(sliced.mask, [True, False])
        assert sliced[0] is ma.masked
        assert sliced[1] == 3


# ---------------------------------------------------------------------------
# 2. Integration tests — through Column / MaskedColumn public API
# ---------------------------------------------------------------------------


class TestColumnGetitemShim:
    """Integration tests for ``_ColumnGetitemShim`` via ``Column``.

    These tests are similar to the tests above but more focused on interaction between shim and ```Column``` machinery
    """

    def test_multidim_integer_index_returns_ndarray_not_column(self):
        col = Column([[7, 8, 9], [10, 11, 12], [13, 14, 15]])

        result = col[1]
        assert isinstance(result, np.ndarray)
        assert not isinstance(result, Column)
        assert isinstance(col[1:], Column)

    def test_void_dtype_string_key(self):
        data = np.array(
            [("Planck18", 67.66), ("Planck13", 67.77)],
            dtype=[("name", "U10"), ("H0", "f4")],
        )
        col = Column(data)
        np.testing.assert_array_equal(col["name"], data["name"])

    def test_bytes_dtype_scalar_decoded(self):
        col = Column([b"Planck18", b"Planck13"], name="x")
        assert isinstance(col[0], str)
        assert col[0] == "Planck18"

    @pytest.mark.parametrize(
        "array",
        [
            [[7, 8, 9, 10]],
            [[7, 8, 9], [10, 11, 12], [13, 14, 15]],
            [[], [], []],
        ],
    )
    def test_multidim_integer_index(self, array):
        data = np.array(array)
        col = Column(data)

        assert (col[0] == data[0]).all()
        assert (col[-1] == data[-1]).all()

        if len(data) > 0:
            n = np.random.randint(0, len(data))
            assert (col[n] == data[n]).all()

        if data.ndim > 1 and data.shape[1] > 0:
            assert col[0][2] == data[0][2]

    @pytest.mark.parametrize(
        "array",
        [
            [[7, 8, 9, 10]],
            [[]],
            [[7, 8, 9], [10, 11, 12], [13, 14, 15]],
            [[], [], []],
        ],
    )
    def test_multidim_slicing(self, array):
        data = np.array(array)
        col = Column(data, name="X")
        n = len(col)

        assert (col[0:] == data[0:]).all()
        assert (col[-n:] == data[-n:]).all()

        if col.data.shape[0] > 1:
            assert (col[1 : n - 1] == data[1 : n - 1]).all()
            assert (col[1 : n - 1 : 2] == data[1 : n - 1 : 2]).all()
            assert (col[-n + 1 : -1 : -2] == data[-n + 1 : -1 : -2]).all()

        assert (col[::2] == data[::2]).all()
        assert (col[::-2] == data[::-2]).all()

    def test_multidim_fancy_indexing(self):
        """Fancy indexing on a multi-dimensional Column must return a copy not view."""
        data = np.array([[7, 8, 9], [4, 5, 2], [1, 5, 56]])
        col = Column(data, name="multidim")

        assert (col[[0, 2]] == data[[0, 2]]).all()
        assert (col[:, [0, 2]] == data[:, [0, 2]]).all()

    def test_sigdim_scalar_integer_index(self):
        data = np.array([7, 8, 9, 10])
        col = Column(data)

        assert col[0] == data[0]
        assert col[-1] == data[-1]

        n = np.random.randint(0, len(data))
        assert col[n] == data[n]

    def test_sigdim_slicing(self):
        data = np.array([7, 8, 9, 10])
        col = Column(data, name="X")
        n = len(col)

        assert (col[0:] == data[0:]).all()
        assert (col[1 : n - 1] == data[1 : n - 1]).all()
        assert (col[::2] == data[::2]).all()
        assert (col[::-2] == data[::-2]).all()

    def test_sigdim_boolean_mask(self):
        col = Column([10, 20, 30, 40])
        mask = np.array([True, False, True, False])
        assert list(col[mask]) == [10, 30]

    def test_sigdim_out_of_bounds_raises(self):
        col = Column([1, 2, 3])
        with pytest.raises(IndexError):
            _ = col[10]

        empty = Column([])
        with pytest.raises(IndexError):
            _ = empty[0]

    def test_sigdim_property_preservation_on_slice(self):
        col = Column(
            [7, 8, 9, 10, 11],
            name="velocity",
            unit="m/s",
            description="radial velocity",
            meta={"source": "Planck18"},
        )
        sliced = col[1:]

        assert sliced.name == col.name
        assert sliced.unit == col.unit
        assert sliced.description == col.description
        assert sliced.meta == col.meta


class TestMaskedColumnGetitemShim:
    """Integration tests for ``_MaskedColumnGetitemShim`` via ``MaskedColumn``."""

    def test_sigdim_masked_scalar_access(self):
        mcol = MaskedColumn([7, 8, 9, 10, 11], mask=[True, False, False, True, True])
        assert mcol[0] is ma.masked
        assert mcol[1] == 8

    def test_multidim_integer_index_returns_masked_array(self):
        mcol = MaskedColumn(
            [[7, 8, 9], [4, 5, 6], [4, 2, 3]],
            mask=[[True, False, False], [True, True, True], [False, False, False]],
        )
        row = mcol[0]
        assert isinstance(row, ma.MaskedArray)
        assert isinstance(row[0], ma.core.MaskedConstant)
        assert mcol[2][0] == 4

    def test_getitem_preserves_column_attributes(self):
        mcol = MaskedColumn(
            [7, 8, 9, 10, 11],
            mask=[False, True, False, True, False],
            name="flux",
            unit="Jy",
            description="flux density",
            meta={"instrument": "HST"},
        )
        sliced = mcol[1:]

        assert sliced.name == mcol.name
        assert sliced.unit == mcol.unit
        assert sliced.description == mcol.description
        assert sliced.meta == mcol.meta

    def test_mask_propagation_on_slice(self):
        mcol = MaskedColumn([1, 2, 3, 4, 5], mask=[False, True, False, True, False])
        sliced = mcol[1:4]

        np.testing.assert_array_equal(sliced.mask, [True, False, True])
        assert sliced[0] is ma.masked
        assert sliced[1] == 3
