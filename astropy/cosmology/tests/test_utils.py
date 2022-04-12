# Licensed under a 3-clause BSD style license - see LICENSE.rst

from math import inf

import numpy as np
import pytest

from astropy.cosmology.utils import aszarr, inf_like, vectorize_if_needed, vectorize_redshift_method
from astropy.utils.exceptions import AstropyDeprecationWarning

from .test_core import _zarr, invalid_zs, valid_zs


def test_vectorize_redshift_method():
    """Test :func:`astropy.cosmology.utils.vectorize_redshift_method`."""
    class Class:

        @vectorize_redshift_method
        def method(self, z):
            return z

    c = Class()

    assert hasattr(c.method, "__vectorized__")
    assert isinstance(c.method.__vectorized__, np.vectorize)

    # calling with Number
    assert c.method(1) == 1
    assert isinstance(c.method(1), int)

    # calling with a numpy scalar
    assert c.method(np.float64(1)) == np.float64(1)
    assert isinstance(c.method(np.float64(1)), np.float64)

    # numpy array
    assert all(c.method(np.array([1, 2])) == np.array([1, 2]))
    assert isinstance(c.method(np.array([1, 2])), np.ndarray)

    # non-scalar
    assert all(c.method([1, 2]) == np.array([1, 2]))
    assert isinstance(c.method([1, 2]), np.ndarray)


def test_vectorize_if_needed():
    """
    Test :func:`astropy.cosmology.utils.vectorize_if_needed`.
    There's no need to test 'veckw' because that is directly pasased to
    `numpy.vectorize` which thoroughly tests the various inputs.

    """
    def func(x):
        return x ** 2

    with pytest.warns(AstropyDeprecationWarning):
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
    with pytest.warns(AstropyDeprecationWarning):
        assert np.all(inf_like(arr) == expected)


# -------------------------------------------------------------------


class Test_aszarr:

    @pytest.mark.parametrize("z, expect", list(zip(valid_zs, [
        0, 1, 1100, np.float64(3300), 2.0, 3.0, _zarr, _zarr, _zarr, _zarr
    ])))
    def test_valid(self, z, expect):
        """Test :func:`astropy.cosmology.utils.aszarr`."""
        got = aszarr(z)
        assert np.array_equal(got, expect)

    @pytest.mark.parametrize("z, exc", invalid_zs)
    def test_invalid(self, z, exc):
        """Test :func:`astropy.cosmology.utils.aszarr`."""
        with pytest.raises(exc):
            aszarr(z)
