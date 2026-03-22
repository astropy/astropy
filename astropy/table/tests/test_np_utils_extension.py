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


JOIN_TYPES = [
    pytest.param(0, id="inner"),
    pytest.param(1, id="outer"),
    pytest.param(2, id="left"),
    pytest.param(3, id="right"),
]


def test_returns_six_values():
    left = [1]
    right = [1]
    idxs, idx_sort, len_left = _make_join_inputs(left, right)
    result = _np_utils.join_inner(idxs, idx_sort, len_left, 0)
    assert isinstance(result, tuple)
    assert len(result) == 6


class TestJoinInner:
    """
    Parametrized tests for `join_inner`.
    """

    @pytest.mark.parametrize("join_type", JOIN_TYPES)
    def test_perfect_match(self, join_type):
        """A perfect match behaves identically across all join types."""
        left = [1, 2, 3]
        right = [1, 2, 3]
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
            pytest.param(0, 1, [1], [False], [0], [False], id="inner"),
            pytest.param(1, 3, [0, 1, 0], [False, False, True], [0, 0, 1], [True, False, False], id="outer"),
            pytest.param(2, 2, [0, 1], [False, False], [0, 0], [True, False], id="left"),
            pytest.param(3, 2, [1, 0], [False, True], [0, 1], [False, False], id="right"),
        ],
    )
    def test_partial_overlap(self, join_type, exp_n, exp_l_out, exp_l_mask, exp_r_out, exp_r_mask):
        """Tests left=[1, 2] and right=[2, 3] across all join types."""
        left = [1, 2]
        right = [2, 3]
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
            pytest.param(0, 0, [], [], [], [], id="inner"),
            pytest.param(1, 4, [0, 1, 0, 0], [False, False, True, True], [0, 0, 0, 1], [True, True, False, False], id="outer"),
            pytest.param(2, 2, [0, 1], [False, False], [0, 0], [True, True], id="left"),
            pytest.param(3, 2, [0, 0], [True, True], [0, 1], [False, False], id="right"),
        ],
    )
    def test_no_overlap(self, join_type, exp_n, exp_l_out, exp_l_mask, exp_r_out, exp_r_mask):
        """Tests disjoint arrays left=[1, 2] and right=[3, 4] across all join types."""
        left = [1, 2]
        right = [3, 4]
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
            pytest.param(0, 2, [0, 1], [False, False], [1, 2], [False, False], id="inner"),
            pytest.param(1, 4, [0, 1, 0, 0], [False, False, True, True], [1, 2, 0, 3], [False, False, True, False], id="outer"),
            pytest.param(2, 3, [0, 1, 2], [False, False, False], [1, 2, 0], [False, False, True], id="left"),
            pytest.param(3, 3, [0, 0, 1], [True, False, False], [0, 1, 2], [False, False, False], id="right"),
        ],
    )
    def test_different_sizes(self, join_type, exp_n, exp_l_out, exp_l_mask, exp_r_out, exp_r_mask):
        """Tests asymmetric overlap where left=[1, 2, 3] and right=[0, 1, 2, 4]."""
        left = [1, 2, 3]
        right = [0, 1, 2, 4]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, join_type
        )

        assert n_out == exp_n
        np.testing.assert_array_equal(left_out, exp_l_out)
        np.testing.assert_array_equal(left_mask, exp_l_mask)
        np.testing.assert_array_equal(right_out, exp_r_out)
        np.testing.assert_array_equal(right_mask, exp_r_mask)

    @pytest.mark.parametrize("join_type", JOIN_TYPES)
    def test_duplicate_keys_cartesian(self, join_type):
        """Cartesian expansion behaves identically across all join types."""
        left = [1, 1]
        right = [1, 1]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, join_type
        )
        assert n_out == 4
        assert len(left_out) == n_out
        assert len(right_out) == n_out

    @pytest.mark.parametrize("join_type", JOIN_TYPES)
    def test_large_cartesian(self, join_type):
        """Scale test for cartesian expansion."""
        left = [7] * 5
        right = [7] * 5
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, join_type
        )
        assert n_out == 25

    @pytest.mark.parametrize("join_type", JOIN_TYPES)
    def test_single_matching_row(self, join_type):
        """Edge case: Arrays of length 1."""
        left = [42]
        right = [42]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, join_type
        )
        assert n_out == 1
        assert not masked

    def test_nan_handling(self):
        """
        Verify that NaNs are treated as unique, non-matching keys.
        Because np.nan != np.nan evaluates to True, `_make_join_inputs` correctly
        flags them as different keys, resulting in no inner join matches.
        """
        left = np.array([1.0, np.nan, np.nan])
        right = np.array([1.0, np.nan])
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = _np_utils.join_inner(
            idxs, idx_sort, len_left, 0  # Inner join
        )

        assert n_out == 1
        assert not masked
        np.testing.assert_array_equal(left_out, [0])
        np.testing.assert_array_equal(right_out, [0])

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