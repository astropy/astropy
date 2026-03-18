import itertools

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
BAD = {
    "str int", "str bool", "uint8 bool", "uint8 str",
    "object float32", "bool object", "uni uint8", "int str",
    "bool str", "bool float64", "bool uni", "str float32",
    "uni float64", "uni object", "bool uint8", "object float64",
    "float32 bool", "str uint8", "uni bool", "float64 bool",
    "float64 object", "int bool", "uni int", "uint8 object",
    "int uni", "uint8 uni", "float32 uni", "object uni",
    "bool float32", "uni float32", "object str", "int object",
    "str float64", "object int", "float64 uni", "bool int",
    "object bool", "object uint8", "float32 object", "str object",
    "float64 str", "float32 str",
}


@pytest.mark.parametrize(
    "name1, name2",
    [
        (name1, name2)
        for (name1, _), (name2, _) in itertools.product(DTYPES, DTYPES)
    ],
)
def test_common_dtype(name1, name2):
    """
    Test common_dtype() for every combination of supported dtypes.
    Known-bad combinations must raise TableMergeError;
    all others must succeed.
    """
    arr = np.empty(1, dtype=DTYPES)
    pair = f"{name1} {name2}"
    if pair in BAD:
        with pytest.raises(np_utils.TableMergeError):
            np_utils.common_dtype([arr[name1], arr[name2]])
    else:
        # Should not raise
        np_utils.common_dtype([arr[name1], arr[name2]])