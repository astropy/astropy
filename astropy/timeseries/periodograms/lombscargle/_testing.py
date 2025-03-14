import numpy as np


def assert_not_strictly_equal(
    a, b, err="Output arrays are the same - dispatch test has failed"
):
    """
    Asserts that passed arguments are not bitwise equal.

    Failure means that arrays can be considered to be a copy of each other.

    Parameters
    ----------
    a : array-like
    b : array-like
    err : str (optional)
        Custom error message for the assertion failure.
    """
    assert np.bitwise_xor(
        np.ascontiguousarray(a).view(np.uint8),
        np.ascontiguousarray(b).view(np.uint8),
    ).any(), err
