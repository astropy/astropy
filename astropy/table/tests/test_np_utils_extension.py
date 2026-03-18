"""
Direct tests for the compiled Cython extension `astropy.table._np_utils`.

`_np_utils` contains a single function, `join_inner()`, which is the core
index-computation engine for all table join operations in astropy.
These tests bypass the public `astropy.table.join()` API to test the
compiled extension directly, fulfilling the goal in #19249 of building
standalone test coverage for performance-critical extensions.

Valid inputs to `join_inner()` are manually constructed here, reproducing
the pre-processing logic from `astropy/table/operations.py`.
The join_type integer encoding is:
0 = inner
1 = outer
2 = left
3 = right
"""

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

    # idxs and idx_sort must be passed as intp arrays
    idxs = np.asarray(idxs, dtype=np.intp)
    idx_sort = np.asarray(idx_sort, dtype=np.intp)

    return idxs, idx_sort, len(left_keys)

def test_module_importable():
    assert _np_utils is not None

def test_join_inner_exists_and_is_callable():
    assert hasattr(_np_utils, "join_inner")
    assert callable(_np_utils.join_inner)

def test_single_public_symbol():
    # Only 'join_inner' is expected to be public
    public_symbols = [s for s in dir(_np_utils) if not s.startswith("_")]
    assert public_symbols == ["join_inner"]

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
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, self.join_type
        )
        assert n_out == len(left)
        assert not masked

    def test_partial_overlap(self):
        left, right = [1, 2], [2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, self.join_type
        )
        assert n_out == 1

    def test_no_overlap(self):
        left, right = [1, 2], [3, 4]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, self.join_type
        )
        assert n_out == 0

    def test_duplicate_keys_cartesian(self):
        left, right = [1, 1], [1, 1]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, self.join_type
        )
        assert n_out == 4

    def test_large_cartesian(self):
        left, right = [7] * 5, [7] * 5
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, self.join_type
        )
        assert n_out == 25


class TestOuterJoin:
    join_type = 1

    def test_disjoint_all_rows_kept(self):
        left, right = [1, 2], [3, 4]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, self.join_type
        )
        assert n_out == 4
        assert masked

    def test_partial_overlap(self):
        left, right = [1, 2], [2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, self.join_type
        )
        assert n_out == 3
        assert masked

    def test_perfect_match_no_masking(self):
        left, right = [1, 2, 3], [1, 2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, self.join_type
        )
        assert not masked


class TestLeftJoin:
    join_type = 2

    def test_all_left_rows_kept(self):
        left, right = [1, 2, 3], [2]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, self.join_type
        )
        assert n_out == 3
        assert masked

    def test_perfect_match_no_masking(self):
        left, right = [1, 2, 3], [1, 2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, self.join_type
        )
        assert not masked

    def test_right_only_keys_ignored(self):
        left, right = [1], [1, 2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, self.join_type
        )
        assert n_out == 1


class TestRightJoin:
    join_type = 3

    def test_all_right_rows_kept(self):
        left, right = [2], [1, 2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, self.join_type
        )
        assert n_out == 3
        assert masked

    def test_perfect_match_no_masking(self):
        left, right = [1, 2], [1, 2]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, self.join_type
        )
        assert not masked


@pytest.mark.parametrize("join_type", [0, 1, 2, 3])
def test_output_array_lengths_match_n_out(join_type):
    left, right = [1, 2], [2, 3]
    idxs, idx_sort, len_left = _make_join_inputs(left, right)
    masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
        idxs, idx_sort, len_left, join_type
    )
    assert len(left_out) == n_out
    assert len(right_out) == n_out

@pytest.mark.parametrize("join_type", [0, 1, 2, 3])
def test_mask_array_lengths_match_n_out(join_type):
    left, right = [1, 2], [2, 3]
    idxs, idx_sort, len_left = _make_join_inputs(left, right)
    masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
        idxs, idx_sort, len_left, join_type
    )
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
    left = list(range(100))
    right = list(range(100))
    idxs, idx_sort, len_left = _make_join_inputs(left, right)
    masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
        idxs, idx_sort, len_left, 0
    )
    assert n_out == 100
    assert not masked
