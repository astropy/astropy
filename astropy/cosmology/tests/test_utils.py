# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

from astropy.cosmology.utils import aszarr, vectorize_redshift_method

from .test_core import invalid_zs, valid_zs, z_arr


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


# -------------------------------------------------------------------


class Test_aszarr:
    @pytest.mark.parametrize(
        "z, expect",
        list(
            zip(
                valid_zs,
                [0, 1, 1100, np.float64(3300), 2.0, 3.0, z_arr, z_arr, z_arr, z_arr],
            )
        ),
    )
    def test_valid(self, z, expect):
        """Test :func:`astropy.cosmology.utils.aszarr`."""
        got = aszarr(z)
        assert np.array_equal(got, expect)

    @pytest.mark.parametrize("z, exc", invalid_zs)
    def test_invalid(self, z, exc):
        """Test :func:`astropy.cosmology.utils.aszarr`."""
        with pytest.raises(exc):
            aszarr(z)
