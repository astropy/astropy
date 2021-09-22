# Licensed under a 3-clause BSD style license - see LICENSE.rst

from math import inf

import pytest

import numpy as np

from astropy.cosmology.utils import inf_like, vectorize_if_needed


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
