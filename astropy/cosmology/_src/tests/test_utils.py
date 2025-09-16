# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

import astropy.units as u
from astropy.cosmology._src.utils import (
    all_cls_vars,
    aszarr,
    deprecated_keywords,
    vectorize_redshift_method,
)
from astropy.utils.compat.optional_deps import HAS_PANDAS

from .test_core import invalid_zs, valid_zs, z_arr


def test_vectorize_redshift_method():
    """Test :func:`astropy.cosmology._src.utils.vectorize_redshift_method`."""

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
        """Test :func:`astropy.cosmology._src.utils.aszarr`."""
        got = aszarr(z)
        assert np.array_equal(got, expect)

    @pytest.mark.parametrize("z, exc", invalid_zs)
    def test_invalid(self, z, exc):
        """Test :func:`astropy.cosmology._src.utils.aszarr`."""
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
    """Test :func:`astropy.cosmology._src.utils.all_cls_vars`."""

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


class TestDeprecatedKeywords:
    @classmethod
    def setup_class(cls):
        def noop(a, b, c, d):
            # a minimal function that does nothing,
            # with multiple positional-or-keywords arguments
            return

        cls.base_func = noop
        cls.depr_funcs = {
            1: deprecated_keywords("a", since="999.999.999")(noop),
            2: deprecated_keywords("a", "b", since="999.999.999")(noop),
            4: deprecated_keywords("a", "b", "c", "d", since="999.999.999")(noop),
        }

    def test_type_safety(self):
        dec = deprecated_keywords(b"a", since="999.999.999")
        with pytest.raises(TypeError, match=r"names\[0\] must be a string"):
            dec(self.base_func)

        dec = deprecated_keywords("a", since=b"999.999.999")
        with pytest.raises(TypeError, match=r"since must be a string"):
            dec(self.base_func)

    @pytest.mark.parametrize("n_deprecated_keywords", [1, 2, 4])
    def test_no_warn(self, n_deprecated_keywords):
        func = self.depr_funcs[n_deprecated_keywords]
        func(1, 2, 3, 4)

    @pytest.mark.parametrize(
        "n_deprecated_keywords, args, kwargs, match",
        [
            pytest.param(
                1,
                (),
                {"a": 1, "b": 2, "c": 3, "d": 4},
                r"Passing 'a' as keyword is deprecated since",
                id="1 deprecation, 1 warn",
            ),
            pytest.param(
                2,
                (1,),
                {"b": 2, "c": 3, "d": 4},
                r"Passing 'b' as keyword is deprecated since",
                id="2 deprecation, 1 warn",
            ),
            pytest.param(
                2,
                (),
                {"a": 1, "b": 2, "c": 3, "d": 4},
                r"Passing \['a', 'b'\] arguments as keywords is deprecated since",
                id="2 deprecations, 2 warns",
            ),
            pytest.param(
                4,
                (),
                {"a": 1, "b": 2, "c": 3, "d": 4},
                (
                    r"Passing \['a', 'b', 'c', 'd'\] arguments as keywords "
                    "is deprecated since"
                ),
                id="4 deprecations, 4 warns",
            ),
            pytest.param(
                4,
                (1,),
                {"b": 2, "c": 3, "d": 4},
                r"Passing \['b', 'c', 'd'\] arguments as keywords is deprecated since",
                id="4 deprecations, 3 warns",
            ),
            pytest.param(
                4,
                (1, 2),
                {"c": 3, "d": 4},
                r"Passing \['c', 'd'\] arguments as keywords is deprecated since",
                id="4 deprecations, 2 warns",
            ),
            pytest.param(
                4,
                (1, 2, 3),
                {"d": 4},
                r"Passing 'd' as keyword is deprecated since",
                id="4 deprecations, 1 warn",
            ),
        ],
    )
    def test_warn(self, n_deprecated_keywords, args, kwargs, match):
        func = self.depr_funcs[n_deprecated_keywords]
        with pytest.warns(FutureWarning, match=match):
            func(*args, **kwargs)
