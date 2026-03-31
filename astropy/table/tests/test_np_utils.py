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

    This helper exists because `_get_join_sort_idxs` is currently tightly
    coupled to the high-level `Table` interface. Ideally, the internals
    should be re-wired so this low-level logic doesn't require full Table
    objects, which would remove the need for this mock function.
    """
    combined = np.concatenate([left_keys, right_keys])
    idx_sort = combined.argsort(kind="stable")
    sorted_keys = combined[idx_sort]

    # Identify the boundaries between unique keys.
    # `diffs` is True wherever a key changes, allowing `idxs` to store
    # the index of the first occurrence of each unique key in the sorted array.
    diffs = np.concatenate(([True], sorted_keys[1:] != sorted_keys[:-1], [True]))
    idxs = np.flatnonzero(diffs)

    # idxs and idx_sort must be passed as intp arrays to match Cython DTYPE_t
    idxs = np.asarray(idxs, dtype=np.intp)
    idx_sort = np.asarray(idx_sort, dtype=np.intp)

    return idxs, idx_sort


def _check_n_out_unique(n_out, join_type, left_keys, right_keys, expected=None):
    """
    Helper function to verify the output length using explicit set theory.
    NOTE: This logic is only mathematically valid when all keys within both
    the left and right arrays are strictly unique.
    """
    if len(set(left_keys)) < len(left_keys) or len(set(right_keys)) < len(right_keys):
        raise AssertionError
    match join_type:
        case 0:  # inner
            n = len(set(left_keys).intersection(right_keys))
        case 1:  # outer
            n = len(set(left_keys).union(right_keys))
        case 2:  # left
            n = len(left_keys)
        case 3:  # right
            n = len(right_keys)
        case _:
            raise ValueError(f"Unknown join_type: {join_type}")

    assert n_out == n
    if expected is not None:
        assert len(expected.left.array) == n
        assert len(expected.right.array) == n


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
        keys = [1, 2, 3]
        idxs, idx_sort = _make_join_inputs(keys, keys)
        masked, n_out, left_indices, left_mask, right_indices, right_mask = join_inner(
            idxs, idx_sort, len(keys), join_type
        )

        _check_n_out_unique(n_out, join_type, keys, keys)
        assert masked == 0
        assert_array_equal(left_indices, [0, 1, 2])
        assert_array_equal(right_indices, [0, 1, 2])
        assert_array_equal(left_mask, np.full(3, False))
        assert_array_equal(right_mask, np.full(3, False))

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
                    left=ArrayMaskPair(array=np.array([0, 1]), mask=np.full(2, False)),
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
                    right=ArrayMaskPair(array=np.array([0, 1]), mask=np.full(2, False)),
                ),
                id="right",
            ),
        ],
    )
    def test_partial_overlap(self, join_type, expected):
        """Tests left=[1, 2] and right=[2, 3] across all join types."""
        left = [1, 2]
        right = [2, 3]
        idxs, idx_sort = _make_join_inputs(left, right)
        masked, n_out, left_indices, left_mask, right_indices, right_mask = join_inner(
            idxs, idx_sort, len(left), join_type
        )

        _check_n_out_unique(n_out, join_type, left, right, expected)
        assert_array_equal(left_indices, expected.left.array)
        assert_array_equal(left_mask, expected.left.mask)
        assert_array_equal(right_indices, expected.right.array)
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
                    left=ArrayMaskPair(array=np.array([0, 1]), mask=np.full(2, False)),
                    right=ArrayMaskPair(array=np.array([0, 0]), mask=np.full(2, True)),
                ),
                id="left",
            ),
            pytest.param(
                3,
                ExpectedResults(
                    left=ArrayMaskPair(array=np.array([0, 0]), mask=np.full(2, True)),
                    right=ArrayMaskPair(array=np.array([0, 1]), mask=np.full(2, False)),
                ),
                id="right",
            ),
        ],
    )
    def test_no_overlap(self, join_type, expected):
        """Tests disjoint arrays left=[1, 2] and right=[3, 4] across all join types."""
        left = [1, 2]
        right = [3, 4]
        idxs, idx_sort = _make_join_inputs(left, right)
        masked, n_out, left_indices, left_mask, right_indices, right_mask = join_inner(
            idxs, idx_sort, len(left), join_type
        )

        _check_n_out_unique(n_out, join_type, left, right, expected)
        assert_array_equal(left_indices, expected.left.array)
        assert_array_equal(left_mask, expected.left.mask)
        assert_array_equal(right_indices, expected.right.array)
        assert_array_equal(right_mask, expected.right.mask)

    @pytest.mark.parametrize(
        "join_type, expected",
        [
            pytest.param(
                0,
                ExpectedResults(
                    left=ArrayMaskPair(array=np.array([0, 1]), mask=np.full(2, False)),
                    right=ArrayMaskPair(array=np.array([1, 2]), mask=np.full(2, False)),
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
                        array=np.array([0, 1, 2]), mask=np.full(3, False)
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
                        array=np.array([0, 1, 2, 3]), mask=np.full(4, False)
                    ),
                ),
                id="right",
            ),
        ],
    )
    def test_different_sizes(self, join_type, expected):
        """Tests overlap with arrays of different lengths and partially disjoint keys."""
        left = [1, 2, 3]
        right = [0, 1, 2, 4]
        idxs, idx_sort = _make_join_inputs(left, right)
        masked, n_out, left_indices, left_mask, right_indices, right_mask = join_inner(
            idxs, idx_sort, len(left), join_type
        )

        _check_n_out_unique(n_out, join_type, left, right, expected)
        assert_array_equal(left_indices, expected.left.array)
        assert_array_equal(left_mask, expected.left.mask)
        assert_array_equal(right_indices, expected.right.array)
        assert_array_equal(right_mask, expected.right.mask)

    @pytest.mark.parametrize("join_type", JOIN_TYPES)
    def test_duplicate_keys_cartesian(self, join_type):
        """Cartesian expansion behaves identically across all join types."""
        keys = [1, 1]
        idxs, idx_sort = _make_join_inputs(keys, keys)
        masked, n_out, left_indices, left_mask, right_indices, right_mask = join_inner(
            idxs, idx_sort, len(keys), join_type
        )

        assert n_out == 4
        assert len(left_indices) == n_out
        assert len(right_indices) == n_out

    @pytest.mark.parametrize("join_type", JOIN_TYPES)
    def test_large_cartesian(self, join_type):
        """Test Cartesian expansion with larger arrays."""
        keys = [7] * 5
        idxs, idx_sort = _make_join_inputs(keys, keys)
        masked, n_out, left_indices, left_mask, right_indices, right_mask = join_inner(
            idxs, idx_sort, len(keys), join_type
        )

        assert n_out == 25

    @pytest.mark.parametrize("join_type", JOIN_TYPES)
    def test_single_matching_row(self, join_type):
        """Edge case: Arrays of length 1."""
        keys = [42]
        idxs, idx_sort = _make_join_inputs(keys, keys)
        masked, n_out, left_indices, left_mask, right_indices, right_mask = join_inner(
            idxs, idx_sort, len(keys), join_type
        )

        _check_n_out_unique(n_out, join_type, keys, keys)
        assert masked == 0

    def test_nan_handling(self):
        """
        Verify that NaNs are treated as unique, non-matching keys.
        Because np.nan != np.nan evaluates to True, `_make_join_inputs` correctly
        flags them as different keys, resulting in no inner join matches.
        """
        left = np.array([1.0, np.nan, np.nan])
        right = np.array([1.0, np.nan])
        idxs, idx_sort = _make_join_inputs(left, right)
        masked, n_out, left_indices, left_mask, right_indices, right_mask = join_inner(
            idxs,
            idx_sort,
            len(left),
            0,  # Inner join
        )

        assert n_out == 1
        assert masked == 0
        assert_array_equal(left_indices, [0])
        assert_array_equal(right_indices, [0])

    def test_many_unique_keys_inner(self):
        """Scale test for inner join specifically."""
        keys = np.arange(100)
        idxs, idx_sort = _make_join_inputs(keys, keys)
        masked, n_out, left_indices, left_mask, right_indices, right_mask = join_inner(
            idxs, idx_sort, len(keys), 0
        )

        _check_n_out_unique(n_out, 0, keys, keys)
        assert masked == 0
