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


class TestJoinInner:
    """
    Parametrized tests for `join_inner`. 
    Join types: 0 (Inner), 1 (Outer), 2 (Left), 3 (Right)
    """

    @pytest.mark.parametrize("join_type", [0, 1, 2, 3])
    def test_perfect_match(self, join_type):
        """A perfect match behaves identically across all join types."""
        left, right = [1, 2, 3], [1, 2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, join_type
        )

        assert n_out == 3
        assert not masked
        np.testing.assert_array_equal(left_out, [0, 1, 2])
        np.testing.assert_array_equal(right_out, [0, 1, 2])
        np.testing.assert_array_equal(left_mask, [False, False, False])
        np.testing.assert_array_equal(right_mask, [False, False, False])

    @pytest.mark.parametrize(
        "join_type, exp_n, exp_l_out, exp_l_mask, exp_r_out, exp_r_mask",
        [
            (0, 1, [1], [False], [0], [False]),  # Inner
            (1, 3, [0, 1, 0], [False, False, True], [0, 0, 1], [True, False, False]),  # Outer
            (2, 2, [0, 1], [False, False], [0, 0], [True, False]),  # Left
            (3, 2, [0, 0], [True, False], [0, 1], [False, False]),  # Right
        ],
    )
    def test_partial_overlap(self, join_type, exp_n, exp_l_out, exp_l_mask, exp_r_out, exp_r_mask):
        """Tests left=[1, 2] and right=[2, 3] across all join types."""
        left, right = [1, 2], [2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, join_type
        )

        assert n_out == exp_n
        np.testing.assert_array_equal(left_out, exp_l_out)
        np.testing.assert_array_equal(left_mask, exp_l_mask)
        np.testing.assert_array_equal(right_out, exp_r_out)
        np.testing.assert_array_equal(right_mask, exp_r_mask)

    @pytest.mark.parametrize(
        "join_type, exp_n, exp_l_out, exp_l_mask, exp_r_out, exp_r_mask",
        [
            (0, 0, [], [], [], []),  # Inner
            (1, 4, [0, 1, 0, 0], [False, False, True, True], [0, 0, 0, 1], [True, True, False, False]),  # Outer
            (2, 2, [0, 1], [False, False], [0, 0], [True, True]),  # Left
            (3, 2, [0, 0], [True, True], [0, 1], [False, False]),  # Right
        ],
    )
    def test_no_overlap(self, join_type, exp_n, exp_l_out, exp_l_mask, exp_r_out, exp_r_mask):
        """Tests disjoint arrays left=[1, 2] and right=[3, 4] across all join types."""
        left, right = [1, 2], [3, 4]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, join_type
        )

        assert n_out == exp_n
        np.testing.assert_array_equal(left_out, exp_l_out)
        np.testing.assert_array_equal(left_mask, exp_l_mask)
        np.testing.assert_array_equal(right_out, exp_r_out)
        np.testing.assert_array_equal(right_mask, exp_r_mask)

    @pytest.mark.parametrize("join_type", [0, 1, 2, 3])
    def test_duplicate_keys_cartesian(self, join_type):
        """Cartesian expansion behaves identically across all join types."""
        left, right = [1, 1], [1, 1]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, join_type
        )
        assert n_out == 4

        assert len(left_out) == n_out
        assert len(right_out) == n_out

    @pytest.mark.parametrize("join_type", [0, 1, 2, 3])
    def test_large_cartesian(self, join_type):
        """Scale test for cartesian expansion."""
        left, right = [7] * 5, [7] * 5
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, join_type
        )
        assert n_out == 25

    @pytest.mark.parametrize("join_type", [0, 1, 2, 3])
    def test_single_matching_row(self, join_type):
        """Edge case: Arrays of length 1."""
        left, right = [42], [42]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, join_type
        )
        assert n_out == 1
        assert not masked

    def test_many_unique_keys_inner(self):
        """Scale test for inner join specifically."""
        left = np.arange(100)
        right = np.arange(100)
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, 0
        )
        assert n_out == 100
        assert not masked