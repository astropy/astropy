import numpy as np
import pytest

from astropy.table import np_utils

DTYPES = [
        ("int", int),
        ("uint8", np.uint8),
        ("float32", np.float32),
        ("float64", np.float64),
        ("str", "S2"),
        ("uni", "U2"),
        ("bool", bool),
        ("object", np.object_),
    ]

# Known bad combinations that must raise TableMergeError
BAD_PAIRS = [
    ("str", "int"), ("str", "bool"), ("uint8", "bool"), ("uint8", "str"),
    ("object", "float32"), ("bool", "object"), ("uni", "uint8"), ("int", "str"),
    ("bool", "str"), ("bool", "float64"), ("bool", "uni"), ("str", "float32"),
    ("uni", "float64"), ("uni", "object"), ("bool", "uint8"), ("object", "float64"),
    ("float32", "bool"), ("str", "uint8"), ("uni", "bool"), ("float64", "bool"),
    ("float64", "object"), ("int", "bool"), ("uni", "int"), ("uint8", "object"),
    ("int", "uni"), ("uint8", "uni"), ("float32", "uni"), ("object", "uni"),
    ("bool", "float32"), ("uni", "float32"), ("object", "str"), ("int", "object"),
    ("str", "float64"), ("object", "int"), ("float64", "uni"), ("bool", "int"),
    ("object", "bool"), ("object", "uint8"), ("float32", "object"), ("str", "object"),
    ("float64", "str"), ("float32", "str"),
]

GOOD_PAIRS = [
    (name1, name2)
    for name1, _ in DTYPES
    for name2, _ in DTYPES
    if (name1, name2) not in BAD_PAIRS
]


@pytest.mark.parametrize("name1, name2", BAD_PAIRS)
def test_common_dtype_incompatible(name1, name2):
    """Incompatible dtype combinations must raise TableMergeError."""
    arr = np.empty(1, dtype=DTYPES)
    with pytest.raises(np_utils.TableMergeError):
        np_utils.common_dtype([arr[name1], arr[name2]])


@pytest.mark.parametrize("name1, name2", GOOD_PAIRS)
def test_common_dtype_compatible(name1, name2):
    """Compatible dtype combinations must not raise."""
    arr = np.empty(1, dtype=DTYPES)
    np_utils.common_dtype([arr[name1], arr[name2]])


