# Licensed under a 3-clause BSD style license - see LICENSE.rst

from math import inf

import pytest

import numpy as np

from astropy.cosmology.utils import _float_or_none, inf_like, vectorize_if_needed


@pytest.mark.parametrize("x, digits, expected",
                         [(None, 1, "None"),
                          (10.1234, 3, "10.1"),
                          (10.1234, 7, "10.1234"),  # no more digits
                          # some edge cases I can think of
                          (10.0, 0, "1e+01"),  # weird
                          (10, 5, "10"),  # integer
                          # errors
                          (10, None, "missing precision"),
                          (10, 1.2, "Invalid format specifier"),
                          (10, -3, "missing precision")])
def test__float_or_none(x, digits, expected):
    """Test :func:`astropy.cosmology.utils._float_or_none`."""
    # handle errors
    if not isinstance(digits, int) or digits < 0:
        with pytest.raises(ValueError, match=expected):
            _float_or_none(x, digits=digits)
    # normal use cases
    else:
        assert _float_or_none(x, digits=digits) == expected


def test_vectorize_if_needed():
    """
    Test :func:`astropy.cosmology.utils.vectorize_if_needed`.
    There's no need to test 'veckw' because that is directly pasased to
    `numpy.vectorize` which thoroughly tests the various inputs.

    """
    func = lambda x: x ** 2

    # not vectorized
    assert vectorize_if_needed(func, 2) == 4

    # vectorized
    assert all(vectorize_if_needed(func, [2, 3]) == [4, 9])


@pytest.mark.parametrize("arr, expected",
                         [(0.0, inf),  # float scalar
                          (1, inf),  # integer scalar should give float output
                          ([0.0, 1.0, 2.0, 3.0], (inf, inf, inf, inf)),
                          ([0, 1, 2, 3], (inf, inf, inf, inf)),  # integer list
                         ])
def test_inf_like(arr, expected):
    """
    Test :func:`astropy.cosmology.utils.inf_like`.
    All inputs should give a float output.
    These tests are also in the docstring, but it's better to have them also
    in one consolidated location.
    """
    assert np.all(inf_like(arr) == expected)
