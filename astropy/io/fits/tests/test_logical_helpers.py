# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy.io.fits._logical_helpers import (
    _detect_legacy_logical_vla_heap,
    _logical_to_fits_bytes,
    _logical_vla_heap_has_null,
    _validate_logical_input,
)
from astropy.utils.exceptions import AstropyDeprecationWarning


class TestValidateLogicalInput:
    @pytest.mark.parametrize(
        "row",
        [
            True,
            False,
            np.True_,
            np.array([True, False, True]),
            [True, False, True],
            (True, False),
            np.array([b"T", b"F", b"\x00"], dtype="S1"),
            np.array([], dtype=bool),
            np.array([], dtype="S1"),
            np.array([0, 1, 1], dtype=np.int8),
            [1, 0, 1, 1],
            (0, 1),
            np.array([0, 1, 1]),
            np.array([0, 1, 1], dtype=np.uint16),
            np.array([], dtype=np.int64),
        ],
    )
    def test_accepts(self, row):
        # No warning, no error.
        _validate_logical_input(row)

    @pytest.mark.parametrize(
        "row",
        [
            [2, 0, 1],
            (1, -1),
            np.array([1, 2, 3]),
            np.array([0.0, 1.0]),
            [True, None, False],
            # ``int8`` is the on-disk storage dtype but user-facing
            # int8 input is still subject to the 0/1 rule; internal
            # load paths view-cast i1 to |S1 before validation.
            np.array([5, 6, 7], dtype=np.int8),
            np.array([ord("T"), ord("F"), 0], dtype=np.int8),
        ],
    )
    def test_deprecates_non_bool_non_bytes(self, row):
        with pytest.warns(AstropyDeprecationWarning, match="deprecated"):
            _validate_logical_input(row)

    def test_rejects_invalid_byte_value(self):
        row = np.array([b"T", b"X", b"F"], dtype="S1")
        with pytest.raises(ValueError, match="invalid byte"):
            _validate_logical_input(row)

    def test_rejects_multibyte_string_dtype(self):
        row = np.array([b"TT", b"FF"], dtype="S2")
        with pytest.raises(ValueError, match="itemsize 1"):
            _validate_logical_input(row)


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
