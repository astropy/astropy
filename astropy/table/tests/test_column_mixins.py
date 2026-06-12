from dataclasses import dataclass
from typing import Generic, TypeVar

import numpy as np
import pytest
from numpy.typing import NDArray

from astropy.table._column_mixins import _ColumnGetitemShim, _MaskedColumnGetitemShim


class MinimalColumn(_ColumnGetitemShim, np.ndarray):
    """Minimal concrete class derived from _ColumnGetitemShim"""

    @property
    def data(self):
        return self.view(np.ndarray)


class MinimalMaskedColumn(_MaskedColumnGetitemShim, np.ma.MaskedArray):
    """Minimal concrete class derived from _MaskedColumnGetitemShim."""

    @property
    def data(self):
        return self.view(np.ndarray)

    def _copy_attrs_slice(self, value: np.ndarray) -> np.ndarray:
        return value.copy()


DType = TypeVar("DType", bound=np.generic)


@dataclass(kw_only=True, slots=True, frozen=True)
class MixinTestCase(Generic[DType]):
    input_array: NDArray[DType]
    index: int | str | slice | tuple[int, ...] | list[str]
    expected_output: NDArray[DType] | str


@dataclass(kw_only=True, slots=True, frozen=True)
class MaskedMixinTestCase(Generic[DType]):
    input_array: NDArray[DType]
    input_mask: NDArray[np.bool_]
    index: slice
    expected_data: NDArray[DType]
    expected_mask: NDArray[np.bool_]


@pytest.mark.parametrize(
    "case",
    [
        pytest.param(
            MixinTestCase(
                input_array=np.arange(10, dtype=np.intp).reshape(2, 5),
                index=0,
                expected_output=np.arange(5, dtype=np.intp),
            ),
            id="multidim_integer_indexing",
        ),
        pytest.param(
            MixinTestCase(
                input_array=np.array(
                    [(1, 2.5), (3, 4.5)], dtype=[("a", "i4"), ("b", "f8")]
                ),
                index="a",  # Pulling a single field
                expected_output=np.array([1, 3], dtype="i4"),
            ),
            id="structured_single_field_strips_wrapper",
        ),
        pytest.param(
            MixinTestCase(
                input_array=np.array(
                    [(1, 2.5), (3, 4.5)], dtype=[("a", "i4"), ("b", "f8")]
                ),
                index=["a"],  # Slicing with a list
                expected_output=np.array([(1,), (3,)], dtype=[("a", "i4")]),
            ),
            id="structured_list_slice_keeps_wrapper",
        ),
    ],
)
def test_minimal_column_getitem_standard(case: MixinTestCase):
    minimal_col = case.input_array.view(MinimalColumn)
    result = minimal_col[case.index]

    # 1. Pulling a single field (e.g., index='a') intentionally strips the wrapper
    #    and returns the raw numpy array for stability.
    # 2. Slicing with a list (e.g., index=['a']) intentionally keeps the wrapper
    #    so it acts like a mini-table.
    # We explicitly assert the exact type here to prevent future regressions.

    if isinstance(case.index, list):
        # When slicing with a list, we EXPECT the wrapper to remain.
        assert type(result) is MinimalColumn
    else:
        # For single indexing, we EXPECT the wrapper to be stripped.
        assert type(result) is np.ndarray

    np.testing.assert_array_equal(result, case.expected_output)


@pytest.mark.parametrize(
    "case",
    [
        pytest.param(
            MixinTestCase(
                input_array=np.array([b"test"], dtype="S4"),
                index=0,
                expected_output="test",
            ),
            id="byte_decode_indexing",
        )
    ],
)
def test_minimal_column_getitem_byte_decode(case: MixinTestCase):
    minimal_col = case.input_array.view(MinimalColumn)
    result = minimal_col[case.index]

    assert type(result) is str
    assert result == case.expected_output


@pytest.mark.parametrize(
    "case",
    [
        pytest.param(
            MaskedMixinTestCase(
                input_array=np.arange(6, dtype=np.float64).reshape(2, 3),
                input_mask=np.array(
                    [[False, True, False], [True, False, False]], dtype=np.bool_
                ),
                index=slice(1, 2),
                expected_data=np.array([[3.0, 4.0, 5.0]], dtype=np.float64),
                expected_mask=np.array([[True, False, False]], dtype=np.bool_),
            ),
            id="masked_multidim_slicing",
        )
    ],
)
def test_minimal_masked_column_getitem(case: MaskedMixinTestCase):
    raw_masked = np.ma.MaskedArray(data=case.input_array, mask=case.input_mask)
    minimal_masked = raw_masked.view(MinimalMaskedColumn)

    result = minimal_masked[case.index]

    assert type(result) is MinimalMaskedColumn
    np.testing.assert_array_equal(result.data, case.expected_data)
    np.testing.assert_array_equal(result.mask, case.expected_mask)
