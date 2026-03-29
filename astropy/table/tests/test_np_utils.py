# Licensed under a 3-clause BSD style license - see LICENSE.rst

from dataclasses import dataclass

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy.table._np_utils import join_inner

# Strict type aliases for 1D index arrays and 1D boolean mask arrays
IndexArray = np.ndarray[tuple[int], np.dtype[np.intp]]
MaskArray = np.ndarray[tuple[int], np.dtype[np.bool_]]


@dataclass(kw_only=True, slots=True, frozen=True)
class ArrayMaskPair:
    array: IndexArray
    mask: MaskArray


@dataclass(kw_only=True, slots=True, frozen=True)
class ExpectedResults:
    left: ArrayMaskPair
    right: ArrayMaskPair


def _make_join_inputs(left_keys, right_keys):
    """
    Reproduces the pre-processing from `operations.py` that prepares
    index arrays for `join_inner()`.

    Note for future refactoring: This helper exists because `_get_join_sort_idxs`
    is currently tightly coupled to the high-level `Table` interface. Ideally,
    the internals should be re-wired so this low-level logic doesn't require
    full Table objects, which would remove the need for this mock function.
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


# Common parametrization for simple join type loops to keep code clean
JOIN_TYPES = [
    pytest.param(0, id="inner"),
    pytest.param(1, id="outer"),
    pytest.param(2, id="left"),
    pytest.param(3, id="right"),
]


class TestJoinInner:
    @pytest.mark.parametrize("join_type", JOIN_TYPES)
    def test_perfect_match(self, join_type):
        """A perfect match behaves identically across all join types."""
        left = [1, 2, 3]
        right = [1, 2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = join_inner(
            idxs, idx_sort, len_left, join_type
        )

        assert n_out == 3
        assert masked == 0
        assert_array_equal(left_out, [0, 1, 2])
        assert_array_equal(right_out, [0, 1, 2])
        assert_array_equal(left_mask, [False, False, False])
        assert_array_equal(right_mask, [False, False, False])

    @pytest.mark.parametrize(
        "join_type, expected",
        [
            pytest.param(
                0,
                ExpectedResults(
                    left=ArrayMaskPair(array=np.array([1]), mask=np.array([False])),
                    right=ArrayMaskPair(array=np.array([0]), mask=np.array([False])),
                ),
                id="inner",
            ),
            pytest.param(
                1,
                ExpectedResults(
                    left=ArrayMaskPair(
                        array=np.array([0, 1, 0]), mask=np.array([False, False, True])
                    ),
                    right=ArrayMaskPair(
                        array=np.array([0, 0, 1]), mask=np.array([True, False, False])
                    ),
                ),
                id="outer",
            ),
            pytest.param(
                2,
                ExpectedResults(
                    left=ArrayMaskPair(
                        array=np.array([0, 1]), mask=np.array([False, False])
                    ),
                    right=ArrayMaskPair(
                        array=np.array([0, 0]), mask=np.array([True, False])
                    ),
                ),
                id="left",
            ),
            pytest.param(
                3,
                ExpectedResults(
                    left=ArrayMaskPair(
                        array=np.array([1, 0]), mask=np.array([False, True])
                    ),
                    right=ArrayMaskPair(
                        array=np.array([0, 1]), mask=np.array([False, False])
                    ),
                ),
                id="right",
            ),
        ],
    )
    def test_partial_overlap(self, join_type, expected):
        """Tests left=[1, 2] and right=[2, 3] across all join types."""
        left = [1, 2]
        right = [2, 3]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = join_inner(
            idxs, idx_sort, len_left, join_type
        )

        assert n_out == len(expected.left.array)
        assert_array_equal(left_out, expected.left.array)
        assert_array_equal(left_mask, expected.left.mask)
        assert_array_equal(right_out, expected.right.array)
        assert_array_equal(right_mask, expected.right.mask)

    @pytest.mark.parametrize(
        "join_type, expected",
        [
            pytest.param(
                0,
                ExpectedResults(
                    left=ArrayMaskPair(array=np.array([]), mask=np.array([])),
                    right=ArrayMaskPair(array=np.array([]), mask=np.array([])),
                ),
                id="inner",
            ),
            pytest.param(
                1,
                ExpectedResults(
                    left=ArrayMaskPair(
                        array=np.array([0, 1, 0, 0]),
                        mask=np.array([False, False, True, True]),
                    ),
                    right=ArrayMaskPair(
                        array=np.array([0, 0, 0, 1]),
                        mask=np.array([True, True, False, False]),
                    ),
                ),
                id="outer",
            ),
            pytest.param(
                2,
                ExpectedResults(
                    left=ArrayMaskPair(
                        array=np.array([0, 1]), mask=np.array([False, False])
                    ),
                    right=ArrayMaskPair(
                        array=np.array([0, 0]), mask=np.array([True, True])
                    ),
                ),
                id="left",
            ),
            pytest.param(
                3,
                ExpectedResults(
                    left=ArrayMaskPair(
                        array=np.array([0, 0]), mask=np.array([True, True])
                    ),
                    right=ArrayMaskPair(
                        array=np.array([0, 1]), mask=np.array([False, False])
                    ),
                ),
                id="right",
            ),
        ],
    )
    def test_no_overlap(self, join_type, expected):
        """Tests disjoint arrays left=[1, 2] and right=[3, 4] across all join types."""
        left = [1, 2]
        right = [3, 4]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = join_inner(
            idxs, idx_sort, len_left, join_type
        )

        assert n_out == len(expected.left.array)
        assert_array_equal(left_out, expected.left.array)
        assert_array_equal(left_mask, expected.left.mask)
        assert_array_equal(right_out, expected.right.array)
        assert_array_equal(right_mask, expected.right.mask)

    @pytest.mark.parametrize(
        "join_type, expected",
        [
            pytest.param(
                0,
                ExpectedResults(
                    left=ArrayMaskPair(
                        array=np.array([0, 1]), mask=np.array([False, False])
                    ),
                    right=ArrayMaskPair(
                        array=np.array([1, 2]), mask=np.array([False, False])
                    ),
                ),
                id="inner",
            ),
            pytest.param(
                1,
                ExpectedResults(
                    left=ArrayMaskPair(
                        array=np.array([0, 0, 1, 2, 0]),
                        mask=np.array([True, False, False, False, True]),
                    ),
                    right=ArrayMaskPair(
                        array=np.array([0, 1, 2, 0, 3]),
                        mask=np.array([False, False, False, True, False]),
                    ),
                ),
                id="outer",
            ),
            pytest.param(
                2,
                ExpectedResults(
                    left=ArrayMaskPair(
                        array=np.array([0, 1, 2]), mask=np.array([False, False, False])
                    ),
                    right=ArrayMaskPair(
                        array=np.array([1, 2, 0]), mask=np.array([False, False, True])
                    ),
                ),
                id="left",
            ),
            pytest.param(
                3,
                ExpectedResults(
                    left=ArrayMaskPair(
                        array=np.array([0, 0, 1, 0]),
                        mask=np.array([True, False, False, True]),
                    ),
                    right=ArrayMaskPair(
                        array=np.array([0, 1, 2, 3]),
                        mask=np.array([False, False, False, False]),
                    ),
                ),
                id="right",
            ),
        ],
    )
    def test_different_sizes(self, join_type, expected):
        """Tests asymmetric overlap where left=[1, 2, 3] and right=[0, 1, 2, 4]."""
        left = [1, 2, 3]
        right = [0, 1, 2, 4]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = join_inner(
            idxs, idx_sort, len_left, join_type
        )

        assert n_out == len(expected.left.array)
        assert_array_equal(left_out, expected.left.array)
        assert_array_equal(left_mask, expected.left.mask)
        assert_array_equal(right_out, expected.right.array)
        assert_array_equal(right_mask, expected.right.mask)

    @pytest.mark.parametrize("join_type", JOIN_TYPES)
    def test_duplicate_keys_cartesian(self, join_type):
        """Cartesian expansion behaves identically across all join types."""
        left = [1, 1]
        right = [1, 1]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = join_inner(
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
        masked, n_out, left_out, left_mask, right_out, right_mask = join_inner(
            idxs, idx_sort, len_left, join_type
        )
        assert n_out == 25

    @pytest.mark.parametrize("join_type", JOIN_TYPES)
    def test_single_matching_row(self, join_type):
        """Edge case: Arrays of length 1."""
        left = [42]
        right = [42]
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = join_inner(
            idxs, idx_sort, len_left, join_type
        )
        assert n_out == 1
        assert masked == 0

    def test_nan_handling(self):
        """
        Verify that NaNs are treated as unique, non-matching keys.
        Because np.nan != np.nan evaluates to True, `_make_join_inputs` correctly
        flags them as different keys, resulting in no inner join matches.
        """
        left = np.array([1.0, np.nan, np.nan])
        right = np.array([1.0, np.nan])
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = join_inner(
            idxs,
            idx_sort,
            len_left,
            0,  # Inner join
        )

        assert n_out == 1
        assert masked == 0
        assert_array_equal(left_out, [0])
        assert_array_equal(right_out, [0])

    def test_many_unique_keys_inner(self):
        """Scale test for inner join specifically."""
        left = np.arange(100)
        right = np.arange(100)
        idxs, idx_sort, len_left = _make_join_inputs(left, right)
        masked, n_out, left_out, left_mask, right_out, right_mask = join_inner(
            idxs, idx_sort, len_left, 0
        )
        assert n_out == 100
        assert masked == 0
