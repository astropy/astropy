# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

from astropy.table import _np_utils


def _make_join_inputs(left_keys, right_keys):
    """
    Reproduces the pre-processing from `operations.py` that prepares
    index arrays for `join_inner()`.
    """
    combined = np.concatenate([left_keys, right_keys])
    idx_sort = combined.argsort(kind="stable")
    sorted_keys = combined[idx_sort]
    diffs = np.concatenate(([True], sorted_keys[1:] != sorted_keys[:-1], [True]))
    idxs = np.flatnonzero(diffs)

    # idxs and idx_sort must be passed as intp arrays to match Cython DTYPE_t
    idxs = np.asarray(idxs, dtype=np.intp)
    idx_sort = np.asarray(idx_sort, dtype=np.intp)

    return idxs, idx_sort, len(left_keys)


def test_returns_six_values():
    idxs, idx_sort, len_left = _make_join_inputs([1], [1])
    result = _np_utils.join_inner(idxs, idx_sort, len_left, 0)
    assert isinstance(result, tuple)
    assert len(result) == 6


class TestInnerJoin:
    join_type = 0

    def test_perfect_match(self):
        left, right = [1, 2, 3], [1, 2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = (
            _np_utils.join_inner(idxs, idx_sort, len_left, self.join_type)
        )
        assert n_out == len(left)
        assert not masked

        # Inner join of identical arrays -> indices map 1:1
        np.testing.assert_array_equal(left_out, [0, 1, 2])
        np.testing.assert_array_equal(right_out, [0, 1, 2])
        # Masks should be all False because technically all data is present
        np.testing.assert_array_equal(left_mask, [False, False, False])
        np.testing.assert_array_equal(right_mask, [False, False, False])

    def test_partial_overlap(self):
        left, right = [1, 2], [2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = (
            _np_utils.join_inner(idxs, idx_sort, len_left, self.join_type)
        )
        assert n_out == 1

        np.testing.assert_array_equal(left_out, [1])
        np.testing.assert_array_equal(right_out, [0])
        np.testing.assert_array_equal(left_mask, [False])
        np.testing.assert_array_equal(right_mask, [False])

    def test_no_overlap(self):
        left, right = [1, 2], [3, 4]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = (
            _np_utils.join_inner(idxs, idx_sort, len_left, self.join_type)
        )
        assert n_out == 0

    def test_duplicate_keys_cartesian(self):
        left, right = [1, 1], [1, 1]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = (
            _np_utils.join_inner(idxs, idx_sort, len_left, self.join_type)
        )
        assert n_out == 4

    def test_large_cartesian(self):
        left, right = [7] * 5, [7] * 5
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = (
            _np_utils.join_inner(idxs, idx_sort, len_left, self.join_type)
        )
        assert n_out == 25


class TestOuterJoin:
    join_type = 1

    def test_disjoint_all_rows_kept(self):
        left, right = [1, 2], [3, 4]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = (
            _np_utils.join_inner(idxs, idx_sort, len_left, self.join_type)
        )
        assert n_out == 4
        assert masked

        # All left rows kept, right missing
        # All right rows kept, left missing
        np.testing.assert_array_equal(left_out, [0, 1, 0, 0])
        np.testing.assert_array_equal(left_mask, [False, False, True, True])
        np.testing.assert_array_equal(right_out, [0, 0, 0, 1])
        np.testing.assert_array_equal(right_mask, [True, True, False, False])

    def test_partial_overlap(self):
        left, right = [1, 2], [2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = (
            _np_utils.join_inner(idxs, idx_sort, len_left, self.join_type)
        )
        assert n_out == 3
        assert masked

        # Keeps 1 (left only), 2 (both), 3 (right only)
        np.testing.assert_array_equal(left_out, [0, 1, 0])
        np.testing.assert_array_equal(left_mask, [False, False, True])
        np.testing.assert_array_equal(right_out, [0, 0, 1])
        np.testing.assert_array_equal(right_mask, [True, False, False])

    def test_perfect_match_no_masking(self):
        left, right = [1, 2, 3], [1, 2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = (
            _np_utils.join_inner(idxs, idx_sort, len_left, self.join_type)
        )
        assert not masked
        np.testing.assert_array_equal(left_out, [0, 1, 2])
        np.testing.assert_array_equal(right_out, [0, 1, 2])


class TestLeftJoin:
    join_type = 2

    def test_all_left_rows_kept(self):
        left, right = [1, 2, 3], [2]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = (
            _np_utils.join_inner(idxs, idx_sort, len_left, self.join_type)
        )
        assert n_out == 3
        assert masked

        np.testing.assert_array_equal(left_out, [0, 1, 2])
        np.testing.assert_array_equal(left_mask, [False, False, False])

        # Right only has a match for the middle element
        np.testing.assert_array_equal(right_out, [0, 0, 0])
        np.testing.assert_array_equal(right_mask, [True, False, True])

    def test_perfect_match_no_masking(self):
        left, right = [1, 2, 3], [1, 2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = (
            _np_utils.join_inner(idxs, idx_sort, len_left, self.join_type)
        )
        assert not masked
        np.testing.assert_array_equal(left_out, [0, 1, 2])
        np.testing.assert_array_equal(right_out, [0, 1, 2])

    def test_right_only_keys_ignored(self):
        left, right = [1], [1, 2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = (
            _np_utils.join_inner(idxs, idx_sort, len_left, self.join_type)
        )
        assert n_out == 1
        np.testing.assert_array_equal(left_out, [0])
        np.testing.assert_array_equal(right_out, [0])


class TestRightJoin:
    join_type = 3

    def test_all_right_rows_kept(self):
        left, right = [2], [1, 2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = (
            _np_utils.join_inner(idxs, idx_sort, len_left, self.join_type)
        )
        assert n_out == 3
        assert masked

        # Left only has a match for the middle element
        np.testing.assert_array_equal(left_out, [0, 0, 0])
        np.testing.assert_array_equal(left_mask, [True, False, True])

        np.testing.assert_array_equal(right_out, [0, 1, 2])
        np.testing.assert_array_equal(right_mask, [False, False, False])

    def test_perfect_match_no_masking(self):
        left, right = [1, 2], [1, 2]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = (
            _np_utils.join_inner(idxs, idx_sort, len_left, self.join_type)
        )
        assert not masked
        np.testing.assert_array_equal(left_out, [0, 1])
        np.testing.assert_array_equal(right_out, [0, 1])


@pytest.mark.parametrize("join_type", [0, 1, 2, 3])
def test_output_and_mask_array_lengths_match_n_out(join_type):
    """left_out, right_out, left_mask, right_mask must all have length n_out."""
    left, right = [1, 2], [2, 3]
    idxs, idx_sort, len_left = _make_join_inputs(left, right)
    masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
        idxs, idx_sort, len_left, join_type
    )
    assert len(left_out) == n_out
    assert len(right_out) == n_out
    assert len(left_mask) == n_out
    assert len(right_mask) == n_out


def test_single_matching_row():
    left, right = [42], [42]
    idxs, idx_sort, len_left = _make_join_inputs(left, right)
    masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
        idxs, idx_sort, len_left, 0
    )
    assert n_out == 1
    assert not masked


def test_single_non_matching_row_inner():
    left, right = [1], [2]
    idxs, idx_sort, len_left = _make_join_inputs(left, right)
    masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
        idxs, idx_sort, len_left, 0
    )
    assert n_out == 0


def test_many_unique_keys_inner():
    left = np.arange(100)
    right = np.arange(100)
    idxs, idx_sort, len_left = _make_join_inputs(left, right)
    masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
        idxs, idx_sort, len_left, 0
    )
    assert n_out == 100
    assert not masked
