# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

import astropy.units as u
from astropy.cosmology import utils
from astropy.cosmology._utils import all_cls_vars, aszarr, vectorize_redshift_method
from astropy.utils.compat.optional_deps import HAS_PANDAS
from astropy.utils.exceptions import AstropyDeprecationWarning

from .test_core import invalid_zs, valid_zs, z_arr


def test_deprecated():
    match = "this private function has been moved to the private module `astropy.cosmology._utils`"

    with pytest.warns(AstropyDeprecationWarning, match=match):
        utils.vectorize_redshift_method()

    with pytest.warns(AstropyDeprecationWarning, match=match):
        utils.aszarr(1)


def test_vectorize_redshift_method():
    """Test :func:`astropy.cosmology._utils.vectorize_redshift_method`."""

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
        """Test :func:`astropy.cosmology._utils.aszarr`."""
        got = aszarr(z)
        assert np.array_equal(got, expect)

    @pytest.mark.parametrize("z, exc", invalid_zs)
    def test_invalid(self, z, exc):
        """Test :func:`astropy.cosmology._utils.aszarr`."""
        with pytest.raises(exc):
            aszarr(z)

    @pytest.mark.skipif(not HAS_PANDAS, reason="requires pandas")
    def test_pandas(self):
        import pandas as pd

        x = pd.Series([1, 2, 3, 4, 5])

        # Demonstrate Pandas doesn't work with units
        assert not isinstance(x * u.km, u.Quantity)

        # Test aszarr works with Pandas
        assert isinstance(aszarr(x), np.ndarray)
        np.testing.assert_array_equal(aszarr(x), x.values)


# -------------------------------------------------------------------


def test_all_cls_vars():
    """Test :func:`astropy.cosmology._utils.all_cls_vars`."""

    class ClassA:
        a = 1
        b = 2

    all_vars = all_cls_vars(ClassA)
    public_all_vars = {k: v for k, v in all_vars.items() if not k.startswith("_")}
    assert public_all_vars == {"a": 1, "b": 2}

    class ClassB(ClassA):
        c = 3

    all_vars = all_cls_vars(ClassB)
    public_all_vars = {k: v for k, v in all_vars.items() if not k.startswith("_")}
    assert public_all_vars == {"a": 1, "b": 2, "c": 3}
    assert "a" not in vars(ClassB)
    assert "b" not in vars(ClassB)
