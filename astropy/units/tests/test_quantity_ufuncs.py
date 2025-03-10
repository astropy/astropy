# The purpose of these tests are to ensure that calling ufuncs with quantities
# returns quantities with the right units, or raises exceptions.

from __future__ import annotations

import concurrent.futures
import dataclasses
import warnings
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import pytest
from erfa import ufunc as erfa_ufunc
from numpy.testing import assert_allclose, assert_array_equal

from astropy import units as u
from astropy.units import quantity_helper as qh
from astropy.units.quantity_helper.converters import UfuncHelpers
from astropy.units.quantity_helper.helpers import helper_sqrt
from astropy.utils.compat.numpycompat import NUMPY_LT_1_25, NUMPY_LT_2_0, NUMPY_LT_2_3
from astropy.utils.compat.optional_deps import HAS_SCIPY

if TYPE_CHECKING:
    from collections.abc import Callable

if NUMPY_LT_2_0:
    from numpy.core import umath as np_umath
else:
    from numpy._core import umath as np_umath


class testcase(NamedTuple):
    """A test case for a ufunc."""

    f: Callable
    """The ufunc to test."""

    q_in: tuple[u.Quantity]
    """The input quantities."""

    q_out: tuple[u.Quantity]
    """The expected output quantities."""


class testexc(NamedTuple):
    """A test case for a ufunc that should raise an exception."""

    f: Callable
    """The ufunc to test."""

    q_in: tuple[u.Quantity]
    """The input quantities."""

    exc: type
    """The expected exception type."""

    msg: str | None
    """The expected exception message."""


class testwarn(NamedTuple):
    """A test case for a ufunc that should raise a warning."""

    f: Callable
    """The ufunc to test."""

    q_in: tuple[u.Quantity]
    """The input quantities."""

    wfilter: str
    """The expected warning filter."""


@pytest.mark.skip
def test_testcase(tc):
    results = tc.f(*tc.q_in)
    # careful of the following line, would break on a function returning
    # a single tuple (as opposed to tuple of return values)
    results = (results,) if not isinstance(results, tuple) else results
    for result, expected in zip(results, tc.q_out):
        assert result.unit == expected.unit
        assert_allclose(result.value, expected.value, atol=1.0e-15)


@pytest.mark.skip
def test_testexc(te):
    with pytest.raises(te.exc) as exc:
        te.f(*te.q_in)
    if te.msg is not None:
        assert te.msg in exc.value.args[0]


@pytest.mark.skip
def test_testwarn(tw):
    with warnings.catch_warnings():
        warnings.filterwarnings(tw.wfilter)
        tw.f(*tw.q_in)


class TestUfuncHelpers:
    # Note that this test may fail if scipy is present, although the
    # scipy.special ufuncs are only loaded on demand. This is because
    # if a prior test has already imported scipy.special, then this test will be
    # disrupted.
    # The test passes independently of whether erfa is already loaded
    # (which will be the case for a full test, since coordinates uses it).
    @pytest.mark.skipif(HAS_SCIPY, reason="scipy coverage is known to be incomplete")
    def test_coverage(self):
        """Test that we cover all ufunc's"""
        all_np_ufuncs = {
            ufunc for ufunc in np_umath.__dict__.values() if isinstance(ufunc, np.ufunc)
        }

        all_q_ufuncs = qh.UNSUPPORTED_UFUNCS | set(qh.UFUNC_HELPERS.keys())
        # Check that every numpy ufunc is covered.
        assert all_np_ufuncs - all_q_ufuncs == set()
        # Check that all ufuncs we cover come from numpy or erfa.
        # (Since coverage for erfa is incomplete, we do not check
        # this the other way).
        all_erfa_ufuncs = {
            ufunc
            for ufunc in erfa_ufunc.__dict__.values()
            if isinstance(ufunc, np.ufunc)
        }
        assert all_q_ufuncs - all_np_ufuncs - all_erfa_ufuncs == set()

    @pytest.mark.skipif(
        HAS_SCIPY,
        reason=(
            "UFUNC_HELPERS.modules might be in a different state "
            "(by design) if scipy.special already registered"
        ),
    )
    def test_scipy_registered(self):
        # Should be registered as existing even if scipy is not available.
        assert "scipy.special" in qh.UFUNC_HELPERS.modules

    def test_removal_addition(self):
        assert np.add in qh.UFUNC_HELPERS
        assert np.add not in qh.UNSUPPORTED_UFUNCS
        qh.UFUNC_HELPERS[np.add] = None
        assert np.add not in qh.UFUNC_HELPERS
        assert np.add in qh.UNSUPPORTED_UFUNCS
        qh.UFUNC_HELPERS[np.add] = qh.UFUNC_HELPERS[np.subtract]
        assert np.add in qh.UFUNC_HELPERS
        assert np.add not in qh.UNSUPPORTED_UFUNCS

    @pytest.mark.slow
    def test_thread_safety(self, fast_thread_switching):
        def dummy_ufunc(*args, **kwargs):
            return np.sqrt(*args, **kwargs)

        def register():
            return {dummy_ufunc: helper_sqrt}

        workers = 8
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            for _ in range(10000):
                helpers = UfuncHelpers()
                helpers.register_module(
                    "astropy.units.tests.test_quantity_ufuncs",
                    ["dummy_ufunc"],
                    register,
                )
                futures = [
                    executor.submit(lambda: helpers[dummy_ufunc])
                    for i in range(workers)
                ]
                values = [future.result() for future in futures]
                assert values == [helper_sqrt] * workers


class TestQuantityTrigonometricFuncs:
    """
    Test trigonometric functions
    """

    @pytest.mark.parametrize(
        "tc",
        (
            testcase(
                f=np.sin,
                q_in=(30.0 * u.degree,),
                q_out=(0.5 * u.dimensionless_unscaled,),
            ),
            testcase(
                f=np.sin,
                q_in=(np.array([0.0, np.pi / 4.0, np.pi / 2.0]) * u.radian,),
                q_out=(np.array([0.0, 1.0 / np.sqrt(2.0), 1.0]) * u.one,),
            ),
            testcase(
                f=np.arcsin,
                q_in=(np.sin(30.0 * u.degree),),
                q_out=(np.radians(30.0) * u.radian,),
            ),
            testcase(
                f=np.arcsin,
                q_in=(np.sin(np.array([0.0, np.pi / 4.0, np.pi / 2.0]) * u.radian),),
                q_out=(np.array([0.0, np.pi / 4.0, np.pi / 2.0]) * u.radian,),
            ),
            testcase(
                f=np.cos,
                q_in=(np.pi / 3.0 * u.radian,),
                q_out=(0.5 * u.dimensionless_unscaled,),
            ),
            testcase(
                f=np.cos,
                q_in=(np.array([0.0, np.pi / 4.0, np.pi / 2.0]) * u.radian,),
                q_out=(np.array([1.0, 1.0 / np.sqrt(2.0), 0.0]) * u.one,),
            ),
            testcase(
                f=np.arccos,
                q_in=(np.cos(np.pi / 3.0 * u.radian),),
                q_out=(np.pi / 3.0 * u.radian,),
            ),
            testcase(
                f=np.arccos,
                q_in=(np.cos(np.array([0.0, np.pi / 4.0, np.pi / 2.0]) * u.radian),),
                q_out=(np.array([0.0, np.pi / 4.0, np.pi / 2.0]) * u.radian,),
            ),
            testcase(
                f=np.tan,
                q_in=(np.pi / 3.0 * u.radian,),
                q_out=(np.sqrt(3.0) * u.dimensionless_unscaled,),
            ),
            testcase(
                f=np.tan,
                q_in=(np.array([0.0, 45.0, 135.0, 180.0]) * u.degree,),
                q_out=(np.array([0.0, 1.0, -1.0, 0.0]) * u.dimensionless_unscaled,),
            ),
            testcase(
                f=np.arctan,
                q_in=(np.tan(np.pi / 3.0 * u.radian),),
                q_out=(np.pi / 3.0 * u.radian,),
            ),
            testcase(
                f=np.arctan,
                q_in=(np.tan(np.array([10.0, 30.0, 70.0, 80.0]) * u.degree),),
                q_out=(np.radians(np.array([10.0, 30.0, 70.0, 80.0]) * u.degree),),
            ),
            testcase(
                f=np.arctan2,
                q_in=(np.array([10.0, 30.0, 70.0, 80.0]) * u.m, 2.0 * u.km),
                q_out=(
                    np.arctan2(np.array([10.0, 30.0, 70.0, 80.0]), 2000.0) * u.radian,
                ),
            ),
            testcase(
                f=np.arctan2,
                q_in=((np.array([10.0, 80.0]) * u.m / (2.0 * u.km)).to(u.one), 1.0),
                q_out=(np.arctan2(np.array([10.0, 80.0]) / 2000.0, 1.0) * u.radian,),
            ),
            testcase(f=np.deg2rad, q_in=(180.0 * u.degree,), q_out=(np.pi * u.radian,)),
            testcase(f=np.radians, q_in=(180.0 * u.degree,), q_out=(np.pi * u.radian,)),
            testcase(f=np.deg2rad, q_in=(3.0 * u.radian,), q_out=(3.0 * u.radian,)),
            testcase(f=np.radians, q_in=(3.0 * u.radian,), q_out=(3.0 * u.radian,)),
            testcase(f=np.rad2deg, q_in=(60.0 * u.degree,), q_out=(60.0 * u.degree,)),
            testcase(f=np.degrees, q_in=(60.0 * u.degree,), q_out=(60.0 * u.degree,)),
            testcase(f=np.rad2deg, q_in=(np.pi * u.radian,), q_out=(180.0 * u.degree,)),
            testcase(f=np.degrees, q_in=(np.pi * u.radian,), q_out=(180.0 * u.degree,)),
        ),
    )
    def test_testcases(self, tc):
        return test_testcase(tc)

    @pytest.mark.parametrize(
        "te",
        (
            testexc(f=np.deg2rad, q_in=(3.0 * u.m,), exc=TypeError, msg=None),
            testexc(f=np.radians, q_in=(3.0 * u.m,), exc=TypeError, msg=None),
            testexc(f=np.rad2deg, q_in=(3.0 * u.m), exc=TypeError, msg=None),
            testexc(f=np.degrees, q_in=(3.0 * u.m), exc=TypeError, msg=None),
            testexc(
                f=np.sin,
                q_in=(3.0 * u.m,),
                exc=TypeError,
                msg="Can only apply 'sin' function to quantities with angle units",
            ),
            testexc(
                f=np.arcsin,
                q_in=(3.0 * u.m,),
                exc=TypeError,
                msg="Can only apply 'arcsin' function to dimensionless quantities",
            ),
            testexc(
                f=np.cos,
                q_in=(3.0 * u.s,),
                exc=TypeError,
                msg="Can only apply 'cos' function to quantities with angle units",
            ),
            testexc(
                f=np.arccos,
                q_in=(3.0 * u.s,),
                exc=TypeError,
                msg="Can only apply 'arccos' function to dimensionless quantities",
            ),
            testexc(
                f=np.tan,
                q_in=(np.array([1, 2, 3]) * u.N,),
                exc=TypeError,
                msg="Can only apply 'tan' function to quantities with angle units",
            ),
            testexc(
                f=np.arctan,
                q_in=(np.array([1, 2, 3]) * u.N,),
                exc=TypeError,
                msg="Can only apply 'arctan' function to dimensionless quantities",
            ),
            testexc(
                f=np.arctan2,
                q_in=(np.array([1, 2, 3]) * u.N, 1.0 * u.s),
                exc=u.UnitsError,
                msg="compatible dimensions",
            ),
            testexc(
                f=np.arctan2,
                q_in=(np.array([1, 2, 3]) * u.N, 1.0),
                exc=u.UnitsError,
                msg="dimensionless quantities when other arg",
            ),
        ),
    )
    def test_testexcs(self, te):
        return test_testexc(te)

    @pytest.mark.parametrize(
        "tw",
        (testwarn(f=np.arcsin, q_in=(27.0 * u.pc / (15 * u.kpc),), wfilter="error"),),
    )
    def test_testwarns(self, tw):
        return test_testwarn(tw)

    def test_sin_with_quantity_out(self):
        # Test for a useful error message - see gh-16873.
        # Non-quantity input should be treated as dimensionless and thus cannot
        # be converted to radians.
        out = u.Quantity(0)
        with pytest.raises(
            AttributeError,
            match=(
                "'NoneType' object has no attribute 'get_converter'"
                ".*\n.*treated as dimensionless"
            ),
        ):
            np.sin(0.5, out=out)

        # Except if we have the right equivalency in place.
        with u.add_enabled_equivalencies(u.dimensionless_angles()):
            result = np.sin(0.5, out=out)

        assert result is out
        assert result == np.sin(0.5) * u.dimensionless_unscaled


class TestQuantityMathFuncs:
    """
    Test other mathematical functions
    """

    def test_multiply_scalar(self):
        assert np.multiply(4.0 * u.m, 2.0 / u.s) == 8.0 * u.m / u.s
        assert np.multiply(4.0 * u.m, 2.0) == 8.0 * u.m
        assert np.multiply(4.0, 2.0 / u.s) == 8.0 / u.s

    def test_multiply_array(self):
        assert np.all(
            np.multiply(np.arange(3.0) * u.m, 2.0 / u.s)
            == np.arange(0, 6.0, 2.0) * u.m / u.s
        )

    def test_matmul(self):
        q = np.arange(3.0) * u.m
        r = np.matmul(q, q)
        assert r == 5.0 * u.m**2
        # less trivial case.
        q1 = np.eye(3) * u.m
        q2 = np.array(
            [[[1., 0., 0.],
              [0., 1., 0.],
              [0., 0., 1.]],
             [[0., 1., 0.],
              [0., 0., 1.],
              [1., 0., 0.]],
             [[0., 0., 1.],
              [1., 0., 0.],
              [0., 1., 0.]]]
        ) / u.s  # fmt: skip
        r2 = np.matmul(q1, q2)
        assert np.all(r2 == np.matmul(q1.value, q2.value) * q1.unit * q2.unit)

    @pytest.mark.skipif(NUMPY_LT_2_0, reason="vecdot only added in numpy 2.0")
    def test_vecdot(self):
        q1 = np.array([1j, 2j, 3j]) * u.m
        q2 = np.array([4j, 5j, 6j]) / u.s
        o = np.vecdot(q1, q2)
        assert o == (32.0 + 0j) * u.m / u.s

    @pytest.mark.skipif(
        NUMPY_LT_2_3, reason="np.matvec and np.vecmat are new in NumPy 2.3"
    )
    def test_matvec(self):
        vec = np.arange(3) << u.s
        mat = (
            np.array(
                [
                    [1.0, -1.0, 2.0],
                    [0.0, 3.0, -1.0],
                    [-1.0, -1.0, 1.0],
                ]
            )
            << u.m
        )
        ref_matvec = (vec * mat).sum(-1)
        res_matvec = np.matvec(mat, vec)
        assert_array_equal(res_matvec, ref_matvec)

        ref_vecmat = (vec * mat.T).sum(-1)
        res_vecmat = np.vecmat(vec, mat)
        assert_array_equal(res_vecmat, ref_vecmat)

    @pytest.mark.parametrize("function", (np.divide, np.true_divide))
    def test_divide_scalar(self, function):
        assert function(4.0 * u.m, 2.0 * u.s) == function(4.0, 2.0) * u.m / u.s
        assert function(4.0 * u.m, 2.0) == function(4.0, 2.0) * u.m
        assert function(4.0, 2.0 * u.s) == function(4.0, 2.0) / u.s

    @pytest.mark.parametrize("function", (np.divide, np.true_divide))
    def test_divide_array(self, function):
        assert np.all(
            function(np.arange(3.0) * u.m, 2.0 * u.s)
            == function(np.arange(3.0), 2.0) * u.m / u.s
        )

    def test_floor_divide_remainder_and_divmod(self):
        inch = u.Unit(0.0254 * u.m)
        dividend = np.array([1.0, 2.0, 3.0]) * u.m
        divisor = np.array([3.0, 4.0, 5.0]) * inch
        quotient = dividend // divisor
        remainder = dividend % divisor
        assert_allclose(quotient.value, [13.0, 19.0, 23.0])
        assert quotient.unit == u.dimensionless_unscaled
        assert_allclose(remainder.value, [0.0094, 0.0696, 0.079])
        assert remainder.unit == dividend.unit
        quotient2 = np.floor_divide(dividend, divisor)
        remainder2 = np.remainder(dividend, divisor)
        assert np.all(quotient2 == quotient)
        assert np.all(remainder2 == remainder)
        quotient3, remainder3 = divmod(dividend, divisor)
        assert np.all(quotient3 == quotient)
        assert np.all(remainder3 == remainder)

        with pytest.raises(TypeError):
            divmod(dividend, u.km)

        with pytest.raises(TypeError):
            dividend // u.km

        with pytest.raises(TypeError):
            dividend % u.km

        quotient4, remainder4 = np.divmod(dividend, divisor)
        assert np.all(quotient4 == quotient)
        assert np.all(remainder4 == remainder)
        with pytest.raises(TypeError):
            np.divmod(dividend, u.km)

    def test_sqrt_scalar(self):
        assert np.sqrt(4.0 * u.m) == 2.0 * u.m**0.5

    def test_sqrt_array(self):
        assert np.all(
            np.sqrt(np.array([1.0, 4.0, 9.0]) * u.m)
            == np.array([1.0, 2.0, 3.0]) * u.m**0.5
        )

    def test_square_scalar(self):
        assert np.square(4.0 * u.m) == 16.0 * u.m**2

    def test_square_array(self):
        assert np.all(
            np.square(np.array([1.0, 2.0, 3.0]) * u.m)
            == np.array([1.0, 4.0, 9.0]) * u.m**2
        )

    def test_reciprocal_scalar(self):
        assert np.reciprocal(4.0 * u.m) == 0.25 / u.m

    def test_reciprocal_array(self):
        assert np.all(
            np.reciprocal(np.array([1.0, 2.0, 4.0]) * u.m)
            == np.array([1.0, 0.5, 0.25]) / u.m
        )

    def test_heaviside_scalar(self):
        assert np.heaviside(0.0 * u.m, 0.5) == 0.5 * u.dimensionless_unscaled
        assert (
            np.heaviside(0.0 * u.s, 25 * u.percent) == 0.25 * u.dimensionless_unscaled
        )
        assert np.heaviside(2.0 * u.J, 0.25) == 1.0 * u.dimensionless_unscaled

    def test_heaviside_array(self):
        values = np.array([-1.0, 0.0, 0.0, +1.0])
        halfway = np.array([0.75, 0.25, 0.75, 0.25]) * u.dimensionless_unscaled
        assert np.all(
            np.heaviside(values * u.m, halfway * u.dimensionless_unscaled)
            == [0, 0.25, 0.75, +1.0] * u.dimensionless_unscaled
        )

    @pytest.mark.parametrize("function", (np.cbrt,))
    def test_cbrt_scalar(self, function):
        assert function(8.0 * u.m**3) == 2.0 * u.m

    @pytest.mark.parametrize("function", (np.cbrt,))
    def test_cbrt_array(self, function):
        # Calculate cbrt on both sides since on Windows the cube root of 64
        # does not exactly equal 4.  See 4388.
        values = np.array([1.0, 8.0, 64.0])
        assert np.all(function(values * u.m**3) == function(values) * u.m)

    def test_power_scalar(self):
        assert np.power(4.0 * u.m, 2.0) == 16.0 * u.m**2
        assert np.power(4.0, 200.0 * u.cm / u.m) == u.Quantity(
            16.0, u.dimensionless_unscaled
        )
        # regression check on #1696
        assert np.power(4.0 * u.m, 0.0) == 1.0 * u.dimensionless_unscaled

    def test_power_scalar_filledarray(self):
        result = np.power(4.0 * u.m, np.array([2.0, 2.0]))
        assert np.all(result == 16.0 * u.m**2)

    def test_power_scalar_strarray(self):
        with pytest.raises(
            expected_exception=ValueError,
            match="could not convert string to float",
        ):
            np.power(4.0 * u.m, np.array(["foo"]))

    def test_power_array(self):
        assert np.all(
            np.power(np.array([1.0, 2.0, 3.0]) * u.m, 3.0)
            == np.array([1.0, 8.0, 27.0]) * u.m**3
        )
        # regression check on #1696
        assert np.all(
            np.power(np.arange(4.0) * u.m, 0.0) == 1.0 * u.dimensionless_unscaled
        )

    def test_float_power_array(self):
        assert np.all(
            np.float_power(np.array([1.0, 2.0, 3.0]) * u.m, 3.0)
            == np.array([1.0, 8.0, 27.0]) * u.m**3
        )
        # regression check on #1696
        assert np.all(
            np.float_power(np.arange(4.0) * u.m, 0.0) == 1.0 * u.dimensionless_unscaled
        )

    def test_power_array_array(self):
        with pytest.raises(ValueError):
            np.power(4.0 * u.m, [2.0, 4.0])

    def test_power_array_array2(self):
        with pytest.raises(ValueError):
            np.power([2.0, 4.0] * u.m, [2.0, 4.0])

    def test_power_array_array3(self):
        # Identical unit fractions are converted automatically to dimensionless
        # and should be allowed as base for np.power: #4764
        q = [2.0, 4.0] * u.m / u.m
        powers = [2.0, 4.0]
        res = np.power(q, powers)
        assert np.all(res.value == q.value**powers)
        assert res.unit == u.dimensionless_unscaled
        # The same holds for unit fractions that are scaled dimensionless.
        q2 = [2.0, 4.0] * u.m / u.cm
        # Test also against different types of exponent
        for cls in (list, tuple, np.array, np.ma.array, u.Quantity):
            res2 = np.power(q2, cls(powers))
            assert np.all(res2.value == q2.to_value(1) ** powers)
            assert res2.unit == u.dimensionless_unscaled
        # Though for single powers, we keep the composite unit.
        res3 = q2**2
        assert np.all(res3.value == q2.value**2)
        assert res3.unit == q2.unit**2
        assert np.all(res3 == q2 ** [2, 2])

    def test_power_invalid(self):
        with pytest.raises(TypeError, match="raise something to a dimensionless"):
            np.power(3.0, 4.0 * u.m)

    def test_copysign_scalar(self):
        assert np.copysign(3 * u.m, 1.0) == 3.0 * u.m
        assert np.copysign(3 * u.m, 1.0 * u.s) == 3.0 * u.m
        assert np.copysign(3 * u.m, -1.0) == -3.0 * u.m
        assert np.copysign(3 * u.m, -1.0 * u.s) == -3.0 * u.m

    def test_copysign_array(self):
        assert np.all(
            np.copysign(np.array([1.0, 2.0, 3.0]) * u.s, -1.0)
            == -np.array([1.0, 2.0, 3.0]) * u.s
        )
        assert np.all(
            np.copysign(np.array([1.0, 2.0, 3.0]) * u.s, -1.0 * u.m)
            == -np.array([1.0, 2.0, 3.0]) * u.s
        )
        assert np.all(
            np.copysign(
                np.array([1.0, 2.0, 3.0]) * u.s, np.array([-2.0, 2.0, -4.0]) * u.m
            )
            == np.array([-1.0, 2.0, -3.0]) * u.s
        )

        q = np.copysign(np.array([1.0, 2.0, 3.0]), -3 * u.m)
        assert np.all(q == np.array([-1.0, -2.0, -3.0]))
        assert not isinstance(q, u.Quantity)

    def test_ldexp_scalar(self):
        assert np.ldexp(4.0 * u.m, 2) == 16.0 * u.m

    def test_ldexp_array(self):
        assert np.all(
            np.ldexp(np.array([1.0, 2.0, 3.0]) * u.m, [3, 2, 1])
            == np.array([8.0, 8.0, 6.0]) * u.m
        )

    def test_ldexp_invalid(self):
        with pytest.raises(TypeError):
            np.ldexp(3.0 * u.m, 4.0)

        with pytest.raises(TypeError):
            np.ldexp(3.0, u.Quantity(4, u.m, dtype=int))

    @pytest.mark.parametrize(
        "function", (np.exp, np.expm1, np.exp2, np.log, np.log2, np.log10, np.log1p)
    )
    def test_exp_scalar(self, function):
        q = function(3.0 * u.m / (6.0 * u.m))
        assert q.unit == u.dimensionless_unscaled
        assert q.value == function(0.5)

    @pytest.mark.parametrize(
        "function", (np.exp, np.expm1, np.exp2, np.log, np.log2, np.log10, np.log1p)
    )
    def test_exp_array(self, function):
        q = function(np.array([2.0, 3.0, 6.0]) * u.m / (6.0 * u.m))
        assert q.unit == u.dimensionless_unscaled
        assert np.all(q.value == function(np.array([1.0 / 3.0, 1.0 / 2.0, 1.0])))
        # should also work on quantities that can be made dimensionless
        q2 = function(np.array([2.0, 3.0, 6.0]) * u.m / (6.0 * u.cm))
        assert q2.unit == u.dimensionless_unscaled
        assert_allclose(q2.value, function(np.array([100.0 / 3.0, 100.0 / 2.0, 100.0])))

    @pytest.mark.parametrize(
        "function", (np.exp, np.expm1, np.exp2, np.log, np.log2, np.log10, np.log1p)
    )
    def test_exp_invalid_units(self, function):
        # Can't use exp() with non-dimensionless quantities
        with pytest.raises(
            TypeError,
            match=(
                f"Can only apply '{function.__name__}' function "
                "to dimensionless quantities"
            ),
        ):
            function(3.0 * u.m / u.s)

    def test_modf_scalar(self):
        q = np.modf(9.0 * u.m / (600.0 * u.cm))
        assert q == (0.5 * u.dimensionless_unscaled, 1.0 * u.dimensionless_unscaled)

    def test_modf_array(self):
        v = np.arange(10.0) * u.m / (500.0 * u.cm)
        q = np.modf(v)
        n = np.modf(v.to_value(u.dimensionless_unscaled))
        assert q[0].unit == u.dimensionless_unscaled
        assert q[1].unit == u.dimensionless_unscaled
        assert all(q[0].value == n[0])
        assert all(q[1].value == n[1])

    def test_frexp_scalar(self):
        q = np.frexp(3.0 * u.m / (6.0 * u.m))
        assert q == (np.array(0.5), np.array(0.0))

    def test_frexp_array(self):
        q = np.frexp(np.array([2.0, 3.0, 6.0]) * u.m / (6.0 * u.m))
        assert all(
            (_q0, _q1) == np.frexp(_d)
            for _q0, _q1, _d in zip(q[0], q[1], [1.0 / 3.0, 1.0 / 2.0, 1.0])
        )

    def test_frexp_invalid_units(self):
        # Can't use prod() with non-dimensionless quantities
        with pytest.raises(
            TypeError,
            match=(
                "Can only apply 'frexp' function to unscaled dimensionless quantities"
            ),
        ):
            np.frexp(3.0 * u.m / u.s)

        # also does not work on quantities that can be made dimensionless
        with pytest.raises(
            TypeError,
            match=(
                "Can only apply 'frexp' function to unscaled dimensionless quantities"
            ),
        ):
            np.frexp(np.array([2.0, 3.0, 6.0]) * u.m / (6.0 * u.cm))

    @pytest.mark.parametrize("function", (np.logaddexp, np.logaddexp2))
    def test_dimensionless_twoarg_array(self, function):
        q = function(np.array([2.0, 3.0, 6.0]) * u.m / (6.0 * u.cm), 1.0)
        assert q.unit == u.dimensionless_unscaled
        assert_allclose(
            q.value, function(np.array([100.0 / 3.0, 100.0 / 2.0, 100.0]), 1.0)
        )

    @pytest.mark.parametrize("function", (np.logaddexp, np.logaddexp2))
    def test_dimensionless_twoarg_invalid_units(self, function):
        with pytest.raises(
            TypeError,
            match=(
                f"Can only apply '{function.__name__}' function to dimensionless"
                " quantities"
            ),
        ):
            function(1.0 * u.km / u.s, 3.0 * u.m / u.s)


class TestInvariantUfuncs:
    @pytest.mark.parametrize(
        "ufunc",
        [
            np.absolute,
            np.fabs,
            np.conj,
            np.conjugate,
            np.negative,
            np.spacing,
            np.rint,
            np.floor,
            np.ceil,
            np.positive,
        ],
    )
    def test_invariant_scalar(self, ufunc):
        q_i = 4.7 * u.m
        q_o = ufunc(q_i)
        assert isinstance(q_o, u.Quantity)
        assert q_o.unit == q_i.unit
        assert q_o.value == ufunc(q_i.value)

    @pytest.mark.parametrize(
        "ufunc", [np.absolute, np.conjugate, np.negative, np.rint, np.floor, np.ceil]
    )
    def test_invariant_array(self, ufunc):
        q_i = np.array([-3.3, 2.1, 10.2]) * u.kg / u.s
        q_o = ufunc(q_i)
        assert isinstance(q_o, u.Quantity)
        assert q_o.unit == q_i.unit
        assert np.all(q_o.value == ufunc(q_i.value))

    @pytest.mark.parametrize(
        "ufunc",
        [
            np.add,
            np.subtract,
            np.hypot,
            np.maximum,
            np.minimum,
            np.nextafter,
            np.remainder,
            np.mod,
            np.fmod,
        ],
    )
    def test_invariant_twoarg_scalar(self, ufunc):
        q_i1 = 4.7 * u.m
        q_i2 = 9.4 * u.km
        q_o = ufunc(q_i1, q_i2)
        assert isinstance(q_o, u.Quantity)
        assert q_o.unit == q_i1.unit
        assert_allclose(q_o.value, ufunc(q_i1.value, q_i2.to_value(q_i1.unit)))

    @pytest.mark.parametrize(
        "ufunc",
        [
            np.add,
            np.subtract,
            np.hypot,
            np.maximum,
            np.minimum,
            np.nextafter,
            np.remainder,
            np.mod,
            np.fmod,
        ],
    )
    def test_invariant_twoarg_array(self, ufunc):
        q_i1 = np.array([-3.3, 2.1, 10.2]) * u.kg / u.s
        q_i2 = np.array([10.0, -5.0, 1.0e6]) * u.g / u.us
        q_o = ufunc(q_i1, q_i2)
        assert isinstance(q_o, u.Quantity)
        assert q_o.unit == q_i1.unit
        assert_allclose(q_o.value, ufunc(q_i1.value, q_i2.to_value(q_i1.unit)))

    @pytest.mark.parametrize(
        ("ufunc", "arbitrary"),
        [
            (np.add, 0.0),
            (np.subtract, 0.0),
            (np.hypot, 0.0),
            (np.maximum, 0.0),
            (np.minimum, 0.0),
            (np.nextafter, 0.0),
            (np.remainder, np.inf),
            (np.mod, np.inf),
            (np.fmod, np.inf),
        ],
    )
    def test_invariant_twoarg_one_arbitrary(self, ufunc, arbitrary):
        q_i1 = np.array([-3.3, 2.1, 10.2]) * u.kg / u.s
        q_o = ufunc(q_i1, arbitrary)
        assert isinstance(q_o, u.Quantity)
        assert q_o.unit == q_i1.unit
        assert_allclose(q_o.value, ufunc(q_i1.value, arbitrary))

    @pytest.mark.parametrize(
        "ufunc",
        [
            np.add,
            np.subtract,
            np.hypot,
            np.maximum,
            np.minimum,
            np.nextafter,
            np.remainder,
            np.mod,
            np.fmod,
        ],
    )
    def test_invariant_twoarg_invalid_units(self, ufunc):
        q_i1 = 4.7 * u.m
        q_i2 = 9.4 * u.s
        with pytest.raises(u.UnitsError, match="compatible dimensions"):
            ufunc(q_i1, q_i2)


class TestComparisonUfuncs:
    @pytest.mark.parametrize(
        "ufunc",
        [np.greater, np.greater_equal, np.less, np.less_equal, np.not_equal, np.equal],
    )
    def test_comparison_valid_units(self, ufunc):
        q_i1 = np.array([-3.3, 2.1, 10.2]) * u.kg / u.s
        q_i2 = np.array([10.0, -5.0, 1.0e6]) * u.g / u.Ms
        q_o = ufunc(q_i1, q_i2)
        assert not isinstance(q_o, u.Quantity)
        assert q_o.dtype == bool
        assert np.all(q_o == ufunc(q_i1.value, q_i2.to_value(q_i1.unit)))
        q_o2 = ufunc(q_i1 / q_i2, 2.0)
        assert not isinstance(q_o2, u.Quantity)
        assert q_o2.dtype == bool
        assert np.all(
            q_o2 == ufunc((q_i1 / q_i2).to_value(u.dimensionless_unscaled), 2.0)
        )
        # comparison with 0., inf, nan is OK even for dimensional quantities
        # (though ignore numpy runtime warnings for comparisons with nan).
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            for arbitrary_unit_value in (0.0, np.inf, np.nan):
                ufunc(q_i1, arbitrary_unit_value)
                ufunc(q_i1, arbitrary_unit_value * np.ones(len(q_i1)))
            # and just for completeness
            ufunc(q_i1, np.array([0.0, np.inf, np.nan]))

    @pytest.mark.parametrize(
        "ufunc",
        [np.greater, np.greater_equal, np.less, np.less_equal, np.not_equal, np.equal],
    )
    def test_comparison_invalid_units(self, ufunc):
        q_i1 = 4.7 * u.m
        q_i2 = 9.4 * u.s
        with pytest.raises(u.UnitsError, match="compatible dimensions"):
            ufunc(q_i1, q_i2)

    @pytest.mark.parametrize("ufunc", (np.isfinite, np.isinf, np.isnan, np.signbit))
    def test_onearg_test_ufuncs(self, ufunc):
        q = [1.0, np.inf, -np.inf, np.nan, -1.0, 0.0] * u.m
        out = ufunc(q)
        assert not isinstance(out, u.Quantity)
        assert out.dtype == bool
        assert np.all(out == ufunc(q.value))

    # Ignore RuntimeWarning raised on Windows and s390.
    @pytest.mark.filterwarnings("ignore:.*invalid value encountered in sign")
    def test_sign(self):
        q = [1.0, np.inf, -np.inf, np.nan, -1.0, 0.0] * u.m

        out = np.sign(q)
        assert not isinstance(out, u.Quantity)
        assert out.dtype == q.dtype
        assert np.all((out == np.sign(q.value)) | (np.isnan(out) & np.isnan(q.value)))


class TestInplaceUfuncs:
    @pytest.mark.parametrize("value", [1.0, np.arange(10.0)])
    def test_one_argument_ufunc_inplace(self, value):
        # without scaling
        s = value * u.rad
        check = s
        np.sin(s, out=s)
        assert check is s
        assert check.unit == u.dimensionless_unscaled
        # with scaling
        s2 = (value * u.rad).to(u.deg)
        check2 = s2
        np.sin(s2, out=s2)
        assert check2 is s2
        assert check2.unit == u.dimensionless_unscaled
        assert_allclose(s.value, s2.value)

    @pytest.mark.parametrize("value", [1.0, np.arange(10.0)])
    def test_one_argument_ufunc_inplace_2(self, value):
        """Check inplace works with non-quantity input and quantity output"""
        s = value * u.m
        check = s
        np.absolute(value, out=s)
        assert check is s
        assert np.all(check.value == np.absolute(value))
        assert check.unit is u.dimensionless_unscaled
        np.sqrt(value, out=s)
        assert check is s
        assert np.all(check.value == np.sqrt(value))
        assert check.unit is u.dimensionless_unscaled
        np.exp(value, out=s)
        assert check is s
        assert np.all(check.value == np.exp(value))
        assert check.unit is u.dimensionless_unscaled
        np.arcsin(value / 10.0, out=s)
        assert check is s
        assert np.all(check.value == np.arcsin(value / 10.0))
        assert check.unit is u.radian

    @pytest.mark.parametrize("value", [1.0, np.arange(10.0)])
    def test_one_argument_two_output_ufunc_inplace(self, value):
        v = 100.0 * value * u.cm / u.m
        v_copy = v.copy()
        tmp = v.copy()
        check = v
        np.modf(v, tmp, v)
        assert check is v
        assert check.unit == u.dimensionless_unscaled
        v2 = v_copy.to(u.dimensionless_unscaled)
        check2 = v2
        np.modf(v2, tmp, v2)
        assert check2 is v2
        assert check2.unit == u.dimensionless_unscaled
        # can also replace in last position if no scaling is needed
        v3 = v_copy.to(u.dimensionless_unscaled)
        check3 = v3
        np.modf(v3, v3, tmp)
        assert check3 is v3
        assert check3.unit == u.dimensionless_unscaled
        # can also replace input with first output when scaling
        v4 = v_copy.copy()
        check4 = v4
        np.modf(v4, v4, tmp)
        assert check4 is v4
        assert check4.unit == u.dimensionless_unscaled

    @pytest.mark.parametrize("value", [1.0, np.arange(10.0)])
    def test_two_argument_ufunc_inplace_1(self, value):
        s = value * u.cycle
        check = s
        s /= 2.0
        assert check is s
        assert np.all(check.value == value / 2.0)
        s /= u.s
        assert check is s
        assert check.unit == u.cycle / u.s
        s *= 2.0 * u.s
        assert check is s
        assert np.all(check == value * u.cycle)

    @pytest.mark.parametrize("value", [1.0, np.arange(10.0)])
    def test_two_argument_ufunc_inplace_2(self, value):
        s = value * u.cycle
        check = s
        np.arctan2(s, s, out=s)
        assert check is s
        assert check.unit == u.radian
        with pytest.raises(u.UnitsError):
            s += 1.0 * u.m
        assert check is s
        assert check.unit == u.radian
        np.arctan2(1.0 * u.deg, s, out=s)
        assert check is s
        assert check.unit == u.radian
        np.add(1.0 * u.deg, s, out=s)
        assert check is s
        assert check.unit == u.deg
        np.multiply(2.0 / u.s, s, out=s)
        assert check is s
        assert check.unit == u.deg / u.s

    def test_two_argument_ufunc_inplace_3(self):
        s = np.array([1.0, 2.0, 3.0]) * u.dimensionless_unscaled
        np.add(np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0]) * 2.0, out=s)
        assert np.all(s.value == np.array([3.0, 6.0, 9.0]))
        assert s.unit is u.dimensionless_unscaled
        np.arctan2(np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0]) * 2.0, out=s)
        assert_allclose(s.value, np.arctan2(1.0, 2.0))
        assert s.unit is u.radian

    @pytest.mark.parametrize("value", [1.0, np.arange(10.0)])
    def test_two_argument_two_output_ufunc_inplace(self, value):
        v = value * u.m
        divisor = 70.0 * u.cm
        v1 = v.copy()
        tmp = v.copy()
        check = np.divmod(v1, divisor, out=(tmp, v1))
        assert check[0] is tmp and check[1] is v1
        assert tmp.unit == u.dimensionless_unscaled
        assert v1.unit == v.unit
        v2 = v.copy()
        check2 = np.divmod(v2, divisor, out=(v2, tmp))
        assert check2[0] is v2 and check2[1] is tmp
        assert v2.unit == u.dimensionless_unscaled
        assert tmp.unit == v.unit
        v3a = v.copy()
        v3b = v.copy()
        check3 = np.divmod(v3a, divisor, out=(v3a, v3b))
        assert check3[0] is v3a and check3[1] is v3b
        assert v3a.unit == u.dimensionless_unscaled
        assert v3b.unit == v.unit

    def test_ufunc_inplace_non_contiguous_data(self):
        # ensure inplace works also for non-contiguous data (closes #1834)
        s = np.arange(10.0) * u.m
        s_copy = s.copy()
        s2 = s[::2]
        s2 += 1.0 * u.cm
        assert np.all(s[::2] > s_copy[::2])
        assert np.all(s[1::2] == s_copy[1::2])

    def test_ufunc_inplace_non_standard_dtype(self):
        """Check that inplace operations check properly for casting.

        First two tests that check that float32 is kept close #3976.
        """
        a1 = u.Quantity([1, 2, 3, 4], u.m, dtype=np.float32)
        a1 *= np.float32(10)
        assert a1.unit is u.m
        assert a1.dtype == np.float32
        a2 = u.Quantity([1, 2, 3, 4], u.m, dtype=np.float32)
        a2 += 20.0 * u.km
        assert a2.unit is u.m
        assert a2.dtype == np.float32
        # For integer, in-place only works if no conversion is done.
        a3 = u.Quantity([1, 2, 3, 4], u.m, dtype=np.int32)
        a3 += u.Quantity(10, u.m, dtype=np.int64)
        assert a3.unit is u.m
        assert a3.dtype == np.int32
        a4 = u.Quantity([1, 2, 3, 4], u.m, dtype=np.int32)
        with pytest.raises(TypeError):
            a4 += u.Quantity(10, u.mm, dtype=np.int64)

    @pytest.mark.parametrize("ufunc", (np.equal, np.greater))
    def test_comparison_ufuncs_inplace(self, ufunc):
        q_i1 = np.array([-3.3, 2.1, 10.2]) * u.kg / u.s
        q_i2 = np.array([10.0, -5.0, 1.0e6]) * u.g / u.Ms
        check = np.empty(q_i1.shape, bool)
        ufunc(q_i1.value, q_i2.to_value(q_i1.unit), out=check)

        result = np.empty(q_i1.shape, bool)
        q_o = ufunc(q_i1, q_i2, out=result)
        assert q_o is result
        assert type(q_o) is np.ndarray
        assert q_o.dtype == bool
        assert np.all(q_o == check)

    @pytest.mark.parametrize("ufunc", (np.isfinite, np.signbit))
    def test_onearg_test_ufuncs_inplace(self, ufunc):
        q = [1.0, np.inf, -np.inf, np.nan, -1.0, 0.0] * u.m
        check = np.empty(q.shape, bool)
        ufunc(q.value, out=check)

        result = np.empty(q.shape, bool)
        out = ufunc(q, out=result)
        assert out is result
        assert type(out) is np.ndarray
        assert out.dtype == bool
        assert np.all(out == ufunc(q.value))

    # Ignore RuntimeWarning raised on Windows and s390.
    @pytest.mark.filterwarnings("ignore:.*invalid value encountered in sign")
    def test_sign_inplace(self):
        q = [1.0, np.inf, -np.inf, np.nan, -1.0, 0.0] * u.m
        check = np.empty(q.shape, q.dtype)

        np.sign(q.value, out=check)

        result = np.empty(q.shape, q.dtype)
        out = np.sign(q, out=result)
        assert out is result
        assert type(out) is np.ndarray
        assert out.dtype == q.dtype
        assert np.all((out == np.sign(q.value)) | (np.isnan(out) & np.isnan(q.value)))

    def test_ndarray_inplace_op_with_quantity(self):
        """Regression test for gh-13911."""
        a = np.arange(3.0)
        q = u.Quantity([12.5, 25.0], u.percent)
        a[:2] += q  # This used to fail
        assert_array_equal(a, np.array([0.125, 1.25, 2.0]))


class TestWhere:
    """Test the where argument in ufuncs."""

    def test_where(self):
        q = np.arange(4.0) << u.m
        out = np.zeros(4) << u.m
        result = np.add(q, 1 * u.km, out=out, where=[True, True, True, False])
        assert result is out
        assert_array_equal(result, [1000.0, 1001.0, 1002.0, 0.0] << u.m)

    @pytest.mark.xfail(
        NUMPY_LT_1_25, reason="where array_ufunc support introduced in numpy 1.25"
    )
    def test_exception_with_where_quantity(self):
        a = np.ones(2)
        where = np.ones(2, bool) << u.m
        with pytest.raises(TypeError, match="all returned NotImplemented"):
            np.add(a, a, out=a, where=where)


@pytest.mark.skipif(not hasattr(np_umath, "clip"), reason="no clip ufunc available")
class TestClip:
    """Test the clip ufunc.

    In numpy, this is hidden behind a function that does not backwards
    compatibility checks.  We explicitly test the ufunc here.
    """

    def setup_method(self):
        self.clip = np_umath.clip

    def test_clip_simple(self):
        q = np.arange(-1.0, 10.0) * u.m
        q_min = 125 * u.cm
        q_max = 0.0055 * u.km
        result = self.clip(q, q_min, q_max)
        assert result.unit == q.unit
        expected = (
            self.clip(q.value, q_min.to_value(q.unit), q_max.to_value(q.unit)) * q.unit
        )
        assert np.all(result == expected)

    def test_clip_unitless_parts(self):
        q = np.arange(-1.0, 10.0) * u.m
        qlim = 0.0055 * u.km
        # one-sided
        result1 = self.clip(q, -np.inf, qlim)
        expected1 = self.clip(q.value, -np.inf, qlim.to_value(q.unit)) * q.unit
        assert np.all(result1 == expected1)
        result2 = self.clip(q, qlim, np.inf)
        expected2 = self.clip(q.value, qlim.to_value(q.unit), np.inf) * q.unit
        assert np.all(result2 == expected2)
        # Zero
        result3 = self.clip(q, np.zeros(q.shape), qlim)
        expected3 = self.clip(q.value, 0, qlim.to_value(q.unit)) * q.unit
        assert np.all(result3 == expected3)
        # Two unitless parts, array-shaped.
        result4 = self.clip(q, np.zeros(q.shape), np.full(q.shape, np.inf))
        expected4 = self.clip(q.value, 0, np.inf) * q.unit
        assert np.all(result4 == expected4)

    def test_clip_dimensionless(self):
        q = np.arange(-1.0, 10.0) * u.dimensionless_unscaled
        result = self.clip(q, 200 * u.percent, 5.0)
        expected = self.clip(q, 2.0, 5.0)
        assert result.unit == u.dimensionless_unscaled
        assert np.all(result == expected)

    def test_clip_ndarray(self):
        a = np.arange(-1.0, 10.0)
        result = self.clip(a, 200 * u.percent, 5.0 * u.dimensionless_unscaled)
        assert isinstance(result, u.Quantity)
        expected = self.clip(a, 2.0, 5.0) * u.dimensionless_unscaled
        assert np.all(result == expected)

    def test_clip_quantity_inplace(self):
        q = np.arange(-1.0, 10.0) * u.m
        q_min = 125 * u.cm
        q_max = 0.0055 * u.km
        expected = (
            self.clip(q.value, q_min.to_value(q.unit), q_max.to_value(q.unit)) * q.unit
        )
        result = self.clip(q, q_min, q_max, out=q)
        assert result is q
        assert np.all(result == expected)

    def test_clip_ndarray_dimensionless_output(self):
        a = np.arange(-1.0, 10.0)
        q = np.zeros_like(a) * u.m
        expected = self.clip(a, 2.0, 5.0) * u.dimensionless_unscaled
        result = self.clip(a, 200 * u.percent, 5.0 * u.dimensionless_unscaled, out=q)
        assert result is q
        assert result.unit == u.dimensionless_unscaled
        assert np.all(result == expected)

    def test_clip_errors(self):
        q = np.arange(-1.0, 10.0) * u.m
        with pytest.raises(u.UnitsError):
            self.clip(q, 0, 1 * u.s)
        with pytest.raises(u.UnitsError):
            self.clip(q.value, 0, 1 * u.s)
        with pytest.raises(u.UnitsError):
            self.clip(q, -1, 0.0)
        with pytest.raises(u.UnitsError):
            self.clip(q, 0.0, 1.0)


class TestUfuncAt:
    """Test that 'at' method for ufuncs (calculates in-place at given indices)

    For Quantities, since calculations are in-place, it makes sense only
    if the result is still a quantity, and if the unit does not have to change
    """

    def test_one_argument_ufunc_at(self):
        q = np.arange(10.0) * u.m
        i = np.array([1, 2])
        qv = q.value.copy()
        np.negative.at(q, i)
        np.negative.at(qv, i)
        assert np.all(q.value == qv)
        assert q.unit is u.m

        # cannot change from quantity to bool array
        with pytest.raises(TypeError):
            np.isfinite.at(q, i)

        # for selective in-place, cannot change the unit
        with pytest.raises(u.UnitsError):
            np.square.at(q, i)

        # except if the unit does not change (i.e., dimensionless)
        d = np.arange(10.0) * u.dimensionless_unscaled
        dv = d.value.copy()
        np.square.at(d, i)
        np.square.at(dv, i)
        assert np.all(d.value == dv)
        assert d.unit is u.dimensionless_unscaled

        d = np.arange(10.0) * u.dimensionless_unscaled
        dv = d.value.copy()
        np.log.at(d, i)
        np.log.at(dv, i)
        assert np.all(d.value == dv)
        assert d.unit is u.dimensionless_unscaled

        # also for sine it doesn't work, even if given an angle
        a = np.arange(10.0) * u.radian
        with pytest.raises(u.UnitsError):
            np.sin.at(a, i)

        # except, for consistency, if we have made radian equivalent to
        # dimensionless (though hopefully it will never be needed)
        av = a.value.copy()
        with u.add_enabled_equivalencies(u.dimensionless_angles()):
            np.sin.at(a, i)
            np.sin.at(av, i)
            assert_allclose(a.value, av)

            # but we won't do double conversion
            ad = np.arange(10.0) * u.degree
            with pytest.raises(u.UnitsError):
                np.sin.at(ad, i)

    def test_two_argument_ufunc_at(self):
        s = np.arange(10.0) * u.m
        i = np.array([1, 2])
        check = s.value.copy()
        np.add.at(s, i, 1.0 * u.km)
        np.add.at(check, i, 1000.0)
        assert np.all(s.value == check)
        assert s.unit is u.m

        with pytest.raises(u.UnitsError):
            np.add.at(s, i, 1.0 * u.s)

        # also raise UnitsError if unit would have to be changed
        with pytest.raises(u.UnitsError):
            np.multiply.at(s, i, 1 * u.s)

        # but be fine if it does not
        s = np.arange(10.0) * u.m
        check = s.value.copy()
        np.multiply.at(s, i, 2.0 * u.dimensionless_unscaled)
        np.multiply.at(check, i, 2)
        assert np.all(s.value == check)
        s = np.arange(10.0) * u.m
        np.multiply.at(s, i, 2.0)
        assert np.all(s.value == check)

        # of course cannot change class of data either
        with pytest.raises(TypeError):
            np.greater.at(s, i, 1.0 * u.km)


class TestUfuncReduceReduceatAccumulate:
    """Test 'reduce', 'reduceat' and 'accumulate' methods for ufuncs

    For Quantities, it makes sense only if the unit does not have to change
    """

    def test_one_argument_ufunc_reduce_accumulate(self):
        # one argument cannot be used
        s = np.arange(10.0) * u.radian
        i = np.array([0, 5, 1, 6])
        with pytest.raises(ValueError):
            np.sin.reduce(s)
        with pytest.raises(ValueError):
            np.sin.accumulate(s)
        with pytest.raises(ValueError):
            np.sin.reduceat(s, i)

    def test_two_argument_ufunc_reduce_accumulate(self):
        s = np.arange(10.0) * u.m
        i = np.array([0, 5, 1, 6])
        check = s.value.copy()
        s_add_reduce = np.add.reduce(s)
        check_add_reduce = np.add.reduce(check)
        assert s_add_reduce.value == check_add_reduce
        assert s_add_reduce.unit is u.m

        s_add_accumulate = np.add.accumulate(s)
        check_add_accumulate = np.add.accumulate(check)
        assert np.all(s_add_accumulate.value == check_add_accumulate)
        assert s_add_accumulate.unit is u.m

        s_add_reduceat = np.add.reduceat(s, i)
        check_add_reduceat = np.add.reduceat(check, i)
        assert np.all(s_add_reduceat.value == check_add_reduceat)
        assert s_add_reduceat.unit is u.m

        # reduce(at) or accumulate on comparisons makes no sense,
        # as intermediate result is not even a Quantity
        with pytest.raises(TypeError):
            np.greater.reduce(s)

        with pytest.raises(TypeError):
            np.greater.accumulate(s)

        with pytest.raises(TypeError):
            np.greater.reduceat(s, i)

        # raise UnitsError if unit would have to be changed
        with pytest.raises(u.UnitsError):
            np.multiply.reduce(s)

        with pytest.raises(u.UnitsError):
            np.multiply.accumulate(s)

        with pytest.raises(u.UnitsError):
            np.multiply.reduceat(s, i)

        # but be fine if it does not
        s = np.arange(10.0) * u.dimensionless_unscaled
        check = s.value.copy()
        s_multiply_reduce = np.multiply.reduce(s)
        check_multiply_reduce = np.multiply.reduce(check)
        assert s_multiply_reduce.value == check_multiply_reduce
        assert s_multiply_reduce.unit is u.dimensionless_unscaled
        s_multiply_accumulate = np.multiply.accumulate(s)
        check_multiply_accumulate = np.multiply.accumulate(check)
        assert np.all(s_multiply_accumulate.value == check_multiply_accumulate)
        assert s_multiply_accumulate.unit is u.dimensionless_unscaled
        s_multiply_reduceat = np.multiply.reduceat(s, i)
        check_multiply_reduceat = np.multiply.reduceat(check, i)
        assert np.all(s_multiply_reduceat.value == check_multiply_reduceat)
        assert s_multiply_reduceat.unit is u.dimensionless_unscaled


class TestUfuncOuter:
    """Test 'outer' methods for ufuncs

    Just a few spot checks, since it uses the same code as the regular
    ufunc call
    """

    def test_one_argument_ufunc_outer(self):
        # one argument cannot be used
        s = np.arange(10.0) * u.radian
        with pytest.raises(ValueError):
            np.sin.outer(s)

    def test_two_argument_ufunc_outer(self):
        s1 = np.arange(10.0) * u.m
        s2 = np.arange(2.0) * u.s
        check1 = s1.value
        check2 = s2.value
        s12_multiply_outer = np.multiply.outer(s1, s2)
        check12_multiply_outer = np.multiply.outer(check1, check2)
        assert np.all(s12_multiply_outer.value == check12_multiply_outer)
        assert s12_multiply_outer.unit == s1.unit * s2.unit

        # raise UnitsError if appropriate
        with pytest.raises(u.UnitsError):
            np.add.outer(s1, s2)

        # but be fine if it does not
        s3 = np.arange(2.0) * s1.unit
        check3 = s3.value
        s13_add_outer = np.add.outer(s1, s3)
        check13_add_outer = np.add.outer(check1, check3)
        assert np.all(s13_add_outer.value == check13_add_outer)
        assert s13_add_outer.unit is s1.unit

        s13_greater_outer = np.greater.outer(s1, s3)
        check13_greater_outer = np.greater.outer(check1, check3)
        assert type(s13_greater_outer) is np.ndarray
        assert np.all(s13_greater_outer == check13_greater_outer)


@dataclasses.dataclass
class DuckQuantity1:
    data: u.Quantity


@dataclasses.dataclass
class DuckQuantity2(DuckQuantity1):
    @property
    def unit(self) -> u.UnitBase:
        return self.data.unit


@dataclasses.dataclass(eq=False)
class DuckQuantity3(DuckQuantity2):
    def __array_ufunc__(self, function, method, *inputs, **kwargs):
        inputs = [inp.data if isinstance(inp, type(self)) else inp for inp in inputs]

        out = kwargs.get("out")

        kwargs_copy = {}
        for k, kwarg in kwargs.items():
            if isinstance(kwarg, type(self)):
                kwargs_copy[k] = kwarg.data
            elif isinstance(kwarg, (list, tuple)):
                kwargs_copy[k] = type(kwarg)(
                    item.data if isinstance(item, type(self)) else item
                    for item in kwarg
                )
            else:
                kwargs_copy[k] = kwarg
        kwargs = kwargs_copy

        for inp in inputs:
            if isinstance(inp, np.ndarray):
                result = inp.__array_ufunc__(function, method, *inputs, **kwargs)
                if result is not NotImplemented:
                    if out is None:
                        return type(self)(result)
                    else:
                        if function.nout == 1:
                            return out[0]
                        else:
                            return out

        return NotImplemented


class DuckQuantity4(DuckQuantity3):
    @property
    def unit(self):
        return DuckQuantity1(1 * self.data.unit)


class TestUfuncReturnsNotImplemented:
    @pytest.mark.parametrize("ufunc", (np.negative, np.abs))
    class TestUnaryUfuncs:
        @pytest.mark.parametrize(
            "duck_quantity",
            [DuckQuantity1(1 * u.mm), DuckQuantity2(1 * u.mm)],
        )
        def test_basic(self, ufunc, duck_quantity):
            with pytest.raises(TypeError, match="bad operand type for .*"):
                ufunc(duck_quantity)

        @pytest.mark.parametrize(
            "duck_quantity",
            [
                DuckQuantity3(1 * u.mm),
                DuckQuantity3([1, 2] * u.mm),
                DuckQuantity4(1 * u.mm),
            ],
        )
        @pytest.mark.parametrize("out", [None, "empty"])
        def test_full(self, ufunc, duck_quantity, out):
            out_expected = out
            if out == "empty":
                out = type(duck_quantity)(np.empty_like(ufunc(duck_quantity.data)))
                out_expected = np.empty_like(ufunc(duck_quantity.data))

            result = ufunc(duck_quantity, out=out)
            if out is not None:
                assert result is out

            result_expected = ufunc(duck_quantity.data, out=out_expected)
            assert np.all(result.data == result_expected)

    @pytest.mark.parametrize("ufunc", (np.add, np.multiply, np.less))
    @pytest.mark.parametrize("quantity", (1 * u.m, [1, 2] * u.m))
    class TestBinaryUfuncs:
        @pytest.mark.parametrize(
            "duck_quantity",
            [DuckQuantity1(1 * u.mm), DuckQuantity2(1 * u.mm)],
        )
        def test_basic(self, ufunc, quantity, duck_quantity):
            with pytest.raises(
                (TypeError, ValueError),
                match=(
                    r"(Unsupported operand type\(s\) for ufunc .*)|"
                    r"(unsupported operand type\(s\) for .*)|"
                    r"(Value not scalar compatible or convertible to an int, float, or complex array)"
                ),
            ):
                ufunc(quantity, duck_quantity)

        @pytest.mark.parametrize(
            "duck_quantity",
            [
                DuckQuantity3(1 * u.mm),
                DuckQuantity3([1, 2] * u.mm),
                DuckQuantity4(1 * u.mm),
            ],
        )
        @pytest.mark.parametrize("out", [None, "empty"])
        def test_full(self, ufunc, quantity, duck_quantity, out):
            out_expected = out
            if out == "empty":
                out = type(duck_quantity)(
                    np.empty_like(ufunc(quantity, duck_quantity.data))
                )
                out_expected = np.empty_like(ufunc(quantity, duck_quantity.data))

            result = ufunc(quantity, duck_quantity, out=out)
            if out is not None:
                assert result is out

            result_expected = ufunc(quantity, duck_quantity.data, out=out_expected)
            assert np.all(result.data == result_expected)


if HAS_SCIPY:
    from scipy import special as sps

    erf_like_ufuncs = (
        sps.erf, sps.erfc, sps.erfcx, sps.erfi,
        sps.gamma, sps.gammaln, sps.loggamma, sps.gammasgn, sps.psi,
        sps.rgamma, sps.digamma, sps.wofz, sps.dawsn,
        sps.entr, sps.exprel, sps.expm1, sps.log1p, sps.exp2, sps.exp10,
    )  # fmt: skip

    if isinstance(sps.erfinv, np.ufunc):
        erf_like_ufuncs += (sps.erfinv, sps.erfcinv)

    def test_scipy_registration():
        """Check that scipy gets loaded upon first use."""
        if sps.erf in qh.UFUNC_HELPERS:
            # Generally, scipy will not be loaded here, but in a double run it might.
            pytest.skip()
        sps.erf(1.0 * u.percent)
        assert sps.erf in qh.UFUNC_HELPERS

    class TestScipySpecialUfuncs:
        @pytest.mark.parametrize("function", erf_like_ufuncs)
        def test_erf_scalar(self, function):
            TestQuantityMathFuncs.test_exp_scalar(None, function)

        @pytest.mark.parametrize("function", erf_like_ufuncs)
        def test_erf_array(self, function):
            TestQuantityMathFuncs.test_exp_array(None, function)

        @pytest.mark.parametrize("function", erf_like_ufuncs)
        def test_erf_invalid_units(self, function):
            TestQuantityMathFuncs.test_exp_invalid_units(None, function)

        @pytest.mark.parametrize("function", (sps.cbrt,))
        def test_cbrt_scalar(self, function):
            TestQuantityMathFuncs.test_cbrt_scalar(None, function)

        @pytest.mark.parametrize("function", (sps.cbrt,))
        def test_cbrt_array(self, function):
            TestQuantityMathFuncs.test_cbrt_array(None, function)

        @pytest.mark.parametrize("function", (sps.radian,))
        def test_radian(self, function):
            q1 = function(180.0 * u.degree, 0.0 * u.arcmin, 0.0 * u.arcsec)
            assert_allclose(q1.value, np.pi)
            assert q1.unit == u.radian

            q2 = function(0.0 * u.degree, 30.0 * u.arcmin, 0.0 * u.arcsec)
            assert_allclose(q2.value, (30.0 * u.arcmin).to(u.radian).value)
            assert q2.unit == u.radian

            q3 = function(0.0 * u.degree, 0.0 * u.arcmin, 30.0 * u.arcsec)
            assert_allclose(q3.value, (30.0 * u.arcsec).to(u.radian).value)

            # the following doesn't make much sense in terms of the name of the
            # routine, but we check it gives the correct result.
            q4 = function(3.0 * u.radian, 0.0 * u.arcmin, 0.0 * u.arcsec)
            assert_allclose(q4.value, 3.0)
            assert q4.unit == u.radian

            with pytest.raises(TypeError):
                function(3.0 * u.m, 2.0 * u.s, 1.0 * u.kg)

        jv_like_ufuncs = (
            sps.jv, sps.jn, sps.jve, sps.yn, sps.yv, sps.yve, sps.kn, sps.kv,
            sps.kve, sps.iv, sps.ive, sps.hankel1, sps.hankel1e, sps.hankel2,
            sps.hankel2e,
        )  # fmt: skip

        @pytest.mark.parametrize("function", jv_like_ufuncs)
        def test_jv_scalar(self, function):
            q = function(2.0 * u.m / (2.0 * u.m), 3.0 * u.m / (6.0 * u.m))
            assert q.unit == u.dimensionless_unscaled
            assert q.value == function(1.0, 0.5)

        @pytest.mark.parametrize("function", jv_like_ufuncs)
        def test_jv_array(self, function):
            q = function(
                np.ones(3) * u.m / (1.0 * u.m),
                np.array([2.0, 3.0, 6.0]) * u.m / (6.0 * u.m),
            )
            assert q.unit == u.dimensionless_unscaled
            assert np.all(
                q.value == function(np.ones(3), np.array([1.0 / 3.0, 1.0 / 2.0, 1.0]))
            )
            # should also work on quantities that can be made dimensionless
            q2 = function(
                np.ones(3) * u.m / (1.0 * u.m),
                np.array([2.0, 3.0, 6.0]) * u.m / (6.0 * u.cm),
            )
            assert q2.unit == u.dimensionless_unscaled
            assert_allclose(
                q2.value,
                function(np.ones(3), np.array([100.0 / 3.0, 100.0 / 2.0, 100.0])),
            )

        @pytest.mark.parametrize("function", jv_like_ufuncs)
        def test_jv_invalid_units(self, function):
            # Can't use jv() with non-dimensionless quantities
            with pytest.raises(
                TypeError,
                match=(
                    f"Can only apply '{function.__name__}' function to dimensionless"
                    " quantities"
                ),
            ):
                function(1.0 * u.kg, 3.0 * u.m / u.s)
