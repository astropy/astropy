# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy.io.fits._logical_helpers import (
    _detect_legacy_logical_vla_heap,
    _logical_row_to_byte_storage,
    _logical_row_uses_byte_storage,
    _logical_to_fits_bytes,
    _logical_vla_heap_has_null,
)


class TestLogicalRowUsesByteStorage:
    @pytest.mark.parametrize(
        "row",
        [
            np.array([b"T", b"F"], dtype="S1"),
            np.ma.masked_array([True, False], mask=[True, False]),
            [True, None, False],
            (None, True),
        ],
    )
    def test_yes(self, row):
        assert _logical_row_uses_byte_storage(row) is True

    @pytest.mark.parametrize(
        "row",
        [
            [True, False],
            np.array([True, False]),
            np.array([1, 0], dtype=np.int8),
        ],
    )
    def test_no(self, row):
        assert _logical_row_uses_byte_storage(row) is False


class TestLogicalRowToByteStorage:
    def test_preserves_S1_input_verbatim(self):
        row = np.array([b"T", b"\x00", b"F"], dtype="S1")
        out = _logical_row_to_byte_storage(row)
        assert out.dtype == np.dtype("S1")
        assert out.tobytes() == b"T\x00F"

    def test_masked_array_writes_null_at_mask(self):
        row = np.ma.masked_array([True, False, True], mask=[False, True, False])
        out = _logical_row_to_byte_storage(row)
        assert out.tobytes() == b"T\x00T"

    def test_list_with_none_writes_null(self):
        out = _logical_row_to_byte_storage([True, None, False])
        assert out.tobytes() == b"T\x00F"

    def test_pure_bool_list(self):
        out = _logical_row_to_byte_storage([True, False, True])
        assert out.tobytes() == b"TFT"

    def test_2d_masked_array_preserves_shape(self):
        # Whole-array MaskedArray (used for fixed-length L columns where
        # the user passes an N-D MaskedArray rather than a list of rows).
        ma = np.ma.masked_array(
            [[False, False, True], [True, False, False]],
            mask=[[True, False, False], [False, True, False]],
        )
        out = _logical_row_to_byte_storage(ma)
        assert out.shape == (2, 3)
        assert out[0].tobytes() == b"\x00FT"
        assert out[1].tobytes() == b"T\x00F"


class TestLogicalToFitsBytes:
    def test_S1_viewed_as_int8(self):
        row = np.array([b"T", b"\x00", b"F"], dtype="S1")
        out = _logical_to_fits_bytes(row)
        assert out.dtype == np.int8
        assert_array_equal(out, [ord("T"), 0, ord("F")])

    def test_bool_input(self):
        out = _logical_to_fits_bytes(np.array([True, False, True]))
        assert_array_equal(out, [ord("T"), ord("F"), ord("T")])

    def test_numeric_nonzero_is_true(self):
        out = _logical_to_fits_bytes(np.array([0, 1, -3, 7]))
        assert_array_equal(out, [ord("F"), ord("T"), ord("T"), ord("T")])


class TestLogicalVlaHeapHasNull:
    @staticmethod
    def _heap(rows):
        """Build a (raw_data, field, heap_offset) triple from per-row bytes."""
        offsets, counts, blob = [], [], b""
        for r in rows:
            offsets.append(len(blob))
            counts.append(len(r))
            blob += r
        field = np.array(list(zip(counts, offsets)), dtype=np.int64)
        raw_data = np.frombuffer(blob, dtype=np.uint8)
        return raw_data, field, 0

    def test_mixed_null_and_tf(self):
        raw, field, off = self._heap([b"T\x00F", b"T"])
        assert _logical_vla_heap_has_null(raw, field, off) is True

    def test_only_nulls_is_ambiguous(self):
        # All-zero heap could be all-NULL or pre-fix all-False; not flagged.
        raw, field, off = self._heap([b"\x00\x00", b"\x00"])
        assert _logical_vla_heap_has_null(raw, field, off) is False

    def test_only_tf_no_null(self):
        raw, field, off = self._heap([b"TF", b"FT"])
        assert _logical_vla_heap_has_null(raw, field, off) is False

    def test_empty_rows_skipped(self):
        raw, field, off = self._heap([b"", b"T\x00"])
        assert _logical_vla_heap_has_null(raw, field, off) is True


class TestDetectLegacyLogicalVlaHeap:
    @staticmethod
    def _heap(blob):
        raw_data = np.frombuffer(blob, dtype=np.uint8)
        field = np.array([[len(blob), 0]], dtype=np.int64)
        return raw_data, field, 0

    def test_legacy_zeros_and_ones(self):
        raw, field, off = self._heap(b"\x00\x01\x01\x00")
        assert _detect_legacy_logical_vla_heap(raw, field, off) is True

    def test_modern_TF_bytes_not_legacy(self):
        raw, field, off = self._heap(b"TFT\x00")
        assert _detect_legacy_logical_vla_heap(raw, field, off) is False

    def test_all_zero_not_legacy_no_anchor(self):
        # Without a 0x01 anchor we cannot distinguish all-NULL from all-False;
        # the function returns False so the modern decoder runs.
        raw, field, off = self._heap(b"\x00\x00\x00")
        assert _detect_legacy_logical_vla_heap(raw, field, off) is False

    def test_no_rows_returns_false(self):
        raw = np.empty(0, dtype=np.uint8)
        field = np.empty((0, 2), dtype=np.int64)
        assert _detect_legacy_logical_vla_heap(raw, field, 0) is False
