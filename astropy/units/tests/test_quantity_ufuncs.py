# The purpose of these tests are to ensure that calling ufuncs with quantities
# returns quantities with the right units, or raises exceptions.

import warnings
from collections import namedtuple

import pytest
import numpy as np
from numpy.testing import assert_allclose

from ... import units as u
from .. import quantity_helper as qh
from ..._erfa import ufunc as erfa_ufunc
from ...tests.helper import raises

try:
    import scipy  # pylint: disable=W0611
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True


testcase = namedtuple('testcase', ['f', 'q_in', 'q_out'])
testexc = namedtuple('testexc', ['f', 'q_in', 'exc', 'msg'])
testwarn = namedtuple('testwarn', ['f', 'q_in', 'wfilter'])


@pytest.mark.skip
def test_testcase(tc):
        results = tc.f(*tc.q_in)
        # careful of the following line, would break on a function returning
        # a single tuple (as opposed to tuple of return values)
        results = (results, ) if type(results) != tuple else results
        for result, expected in zip(results, tc.q_out):
            assert result.unit == expected.unit
            assert_allclose(result.value, expected.value, atol=1.E-15)


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
    # Note that this test should work even if scipy is present, since
    # the scipy.special ufuncs are only loaded on demand.
    # The test passes independently of whether erfa is already loaded
    # (which will be the case for a full test, since coordinates uses it).
    def test_coverage(self):
        """Test that we cover all ufunc's"""

        all_np_ufuncs = set([ufunc for ufunc in np.core.umath.__dict__.values()
                             if isinstance(ufunc, np.ufunc)])

        all_q_ufuncs = (qh.UNSUPPORTED_UFUNCS |
                        set(qh.UFUNC_HELPERS.keys()))
        # Check that every numpy ufunc is covered.
        assert all_np_ufuncs - all_q_ufuncs == set()
        # Check that all ufuncs we cover come from numpy or erfa.
        # (Since coverage for erfa is incomplete, we do not check
        # this the other way).
        all_erfa_ufuncs = set([ufunc for ufunc in erfa_ufunc.__dict__.values()
                               if isinstance(ufunc, np.ufunc)])
        assert (all_q_ufuncs - all_np_ufuncs - all_erfa_ufuncs == set())

    def test_scipy_registered(self):
        # Should be registered as existing even if scipy is not available.
        assert 'scipy.special' in qh.UFUNC_HELPERS.modules

    def test_removal_addition(self):
        assert np.add in qh.UFUNC_HELPERS
        assert np.add not in qh.UNSUPPORTED_UFUNCS
        qh.UFUNC_HELPERS[np.add] = None
        assert np.add not in qh.UFUNC_HELPERS
        assert np.add in qh.UNSUPPORTED_UFUNCS
        qh.UFUNC_HELPERS[np.add] = qh.UFUNC_HELPERS[np.subtract]
        assert np.add in qh.UFUNC_HELPERS
        assert np.add not in qh.UNSUPPORTED_UFUNCS


class TestQuantityTrigonometricFuncs:
    """
    Test trigonometric functions
    """
    @pytest.mark.parametrize('tc', (
        testcase(
            f=np.sin,
            q_in=(30. * u.degree, ),
            q_out=(0.5*u.dimensionless_unscaled, )
        ),
        testcase(
            f=np.sin,
            q_in=(np.array([0., np.pi / 4., np.pi / 2.]) * u.radian, ),
            q_out=(np.array([0., 1. / np.sqrt(2.), 1.]) * u.one, )
        ),
        testcase(
            f=np.arcsin,
            q_in=(np.sin(30. * u.degree), ),
            q_out=(np.radians(30.) * u.radian, )
        ),
        testcase(
            f=np.arcsin,
            q_in=(np.sin(np.array([0., np.pi / 4., np.pi / 2.]) * u.radian), ),
            q_out=(np.array([0., np.pi / 4., np.pi / 2.]) * u.radian, )
        ),
        testcase(
            f=np.cos,
            q_in=(np.pi / 3. * u.radian, ),
            q_out=(0.5 * u.dimensionless_unscaled, )
        ),
        testcase(
            f=np.cos,
            q_in=(np.array([0., np.pi / 4., np.pi / 2.]) * u.radian, ),
            q_out=(np.array([1., 1. / np.sqrt(2.), 0.]) * u.one, )
        ),
        testcase(
            f=np.arccos,
            q_in=(np.cos(np.pi / 3. * u.radian), ),
            q_out=(np.pi / 3. * u.radian, )
        ),
        testcase(
            f=np.arccos,
            q_in=(np.cos(np.array([0., np.pi / 4., np.pi / 2.]) * u.radian), ),
            q_out=(np.array([0., np.pi / 4., np.pi / 2.]) * u.radian, ),
        ),
        testcase(
            f=np.tan,
            q_in=(np.pi / 3. * u.radian, ),
            q_out=(np.sqrt(3.) * u.dimensionless_unscaled, )
        ),
        testcase(
            f=np.tan,
            q_in=(np.array([0., 45., 135., 180.]) * u.degree, ),
            q_out=(np.array([0., 1., -1., 0.]) * u.dimensionless_unscaled, )
        ),
        testcase(
            f=np.arctan,
            q_in=(np.tan(np.pi / 3. * u.radian), ),
            q_out=(np.pi / 3. * u.radian, )
        ),
        testcase(
            f=np.arctan,
            q_in=(np.tan(np.array([10., 30., 70., 80.]) * u.degree), ),
            q_out=(np.radians(np.array([10., 30., 70., 80.]) * u.degree), )
        ),
        testcase(
            f=np.arctan2,
            q_in=(np.array([10., 30., 70., 80.]) * u.m, 2.0 * u.km),
            q_out=(np.arctan2(np.array([10., 30., 70., 80.]),
                              2000.) * u.radian, )
        ),
        testcase(
            f=np.arctan2,
            q_in=((np.array([10., 80.]) * u.m / (2.0 * u.km)).to(u.one), 1.),
            q_out=(np.arctan2(np.array([10., 80.]) / 2000., 1.) * u.radian, )
        ),
        testcase(
            f=np.deg2rad,
            q_in=(180. * u.degree, ),
            q_out=(np.pi * u.radian, )
        ),
        testcase(
            f=np.radians,
            q_in=(180. * u.degree, ),
            q_out=(np.pi * u.radian, )
        ),
        testcase(
            f=np.deg2rad,
            q_in=(3. * u.radian, ),
            q_out=(3. * u.radian, )
        ),
        testcase(
            f=np.radians,
            q_in=(3. * u.radian, ),
            q_out=(3. * u.radian, )
        ),
        testcase(
            f=np.rad2deg,
            q_in=(60. * u.degree, ),
            q_out=(60. * u.degree, )
        ),
        testcase(
            f=np.degrees,
            q_in=(60. * u.degree, ),
            q_out=(60. * u.degree, )
        ),
        testcase(
            f=np.rad2deg,
            q_in=(np.pi * u.radian, ),
            q_out=(180. * u.degree, )
        ),
        testcase(
            f=np.degrees,
            q_in=(np.pi * u.radian, ),
            q_out=(180. * u.degree, )
        )
    ))
    def test_testcases(self, tc):
        return test_testcase(tc)

    @pytest.mark.parametrize('te', (
        testexc(
            f=np.deg2rad,
            q_in=(3. * u.m, ),
            exc=TypeError,
            msg=None
        ),
        testexc(
            f=np.radians,
            q_in=(3. * u.m, ),
            exc=TypeError,
            msg=None
        ),
        testexc(
            f=np.rad2deg,
            q_in=(3. * u.m),
            exc=TypeError,
            msg=None
        ),
        testexc(
            f=np.degrees,
            q_in=(3. * u.m),
            exc=TypeError,
            msg=None
        ),
        testexc(
            f=np.sin,
            q_in=(3. * u.m, ),
            exc=TypeError,
            msg="Can only apply 'sin' function to quantities with angle units"
        ),
        testexc(
            f=np.arcsin,
            q_in=(3. * u.m, ),
            exc=TypeError,
            msg="Can only apply 'arcsin' function to dimensionless quantities"
        ),
        testexc(
            f=np.cos,
            q_in=(3. * u.s, ),
            exc=TypeError,
            msg="Can only apply 'cos' function to quantities with angle units"
        ),
        testexc(
            f=np.arccos,
            q_in=(3. * u.s, ),
            exc=TypeError,
            msg="Can only apply 'arccos' function to dimensionless quantities"
        ),
        testexc(
            f=np.tan,
            q_in=(np.array([1, 2, 3]) * u.N, ),
            exc=TypeError,
            msg="Can only apply 'tan' function to quantities with angle units"
        ),
        testexc(
            f=np.arctan,
            q_in=(np.array([1, 2, 3]) * u.N, ),
            exc=TypeError,
            msg="Can only apply 'arctan' function to dimensionless quantities"
        ),
        testexc(
            f=np.arctan2,
            q_in=(np.array([1, 2, 3]) * u.N, 1. * u.s),
            exc=u.UnitsError,
            msg="compatible dimensions"
        ),
        testexc(
            f=np.arctan2,
            q_in=(np.array([1, 2, 3]) * u.N, 1.),
            exc=u.UnitsError,
            msg="dimensionless quantities when other arg"
        )
    ))
    def test_testexcs(self, te):
        return test_testexc(te)

    @pytest.mark.parametrize('tw', (
        testwarn(
            f=np.arcsin,
            q_in=(27. * u.pc / (15 * u.kpc), ),
            wfilter='error'
        ),
    ))
    def test_testwarns(self, tw):
        return test_testwarn(tw)


class TestQuantityMathFuncs:
    """
    Test other mathematical functions
    """

    def test_multiply_scalar(self):
        assert np.multiply(4. * u.m, 2. / u.s) == 8. * u.m / u.s
        assert np.multiply(4. * u.m, 2.) == 8. * u.m
        assert np.multiply(4., 2. / u.s) == 8. / u.s

    def test_multiply_array(self):
        assert np.all(np.multiply(np.arange(3.) * u.m, 2. / u.s) ==
                      np.arange(0, 6., 2.) * u.m / u.s)

    @pytest.mark.parametrize('function', (np.divide, np.true_divide))
    def test_divide_scalar(self, function):
        assert function(4. * u.m, 2. * u.s) == function(4., 2.) * u.m / u.s
        assert function(4. * u.m, 2.) == function(4., 2.) * u.m
        assert function(4., 2. * u.s) == function(4., 2.) / u.s

    @pytest.mark.parametrize('function', (np.divide, np.true_divide))
    def test_divide_array(self, function):
        assert np.all(function(np.arange(3.) * u.m, 2. * u.s) ==
                      function(np.arange(3.), 2.) * u.m / u.s)

    def test_floor_divide_remainder_and_divmod(self):
        inch = u.Unit(0.0254 * u.m)
        dividend = np.array([1., 2., 3.]) * u.m
        divisor = np.array([3., 4., 5.]) * inch
        quotient = dividend // divisor
        remainder = dividend % divisor
        assert_allclose(quotient.value, [13., 19., 23.])
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
        assert np.sqrt(4. * u.m) == 2. * u.m ** 0.5

    def test_sqrt_array(self):
        assert np.all(np.sqrt(np.array([1., 4., 9.]) * u.m)
                      == np.array([1., 2., 3.]) * u.m ** 0.5)

    def test_square_scalar(self):
        assert np.square(4. * u.m) == 16. * u.m ** 2

    def test_square_array(self):
        assert np.all(np.square(np.array([1., 2., 3.]) * u.m)
                      == np.array([1., 4., 9.]) * u.m ** 2)

    def test_reciprocal_scalar(self):
        assert np.reciprocal(4. * u.m) == 0.25 / u.m

    def test_reciprocal_array(self):
        assert np.all(np.reciprocal(np.array([1., 2., 4.]) * u.m)
                      == np.array([1., 0.5, 0.25]) / u.m)

    def test_heaviside_scalar(self):
        assert np.heaviside(0. * u.m, 0.5) == 0.5 * u.dimensionless_unscaled
        assert np.heaviside(0. * u.s,
                            25 * u.percent) == 0.25 * u.dimensionless_unscaled
        assert np.heaviside(2. * u.J, 0.25) == 1. * u.dimensionless_unscaled

    def test_heaviside_array(self):
        values = np.array([-1., 0., 0., +1.])
        halfway = np.array([0.75, 0.25, 0.75, 0.25]) * u.dimensionless_unscaled
        assert np.all(np.heaviside(values * u.m,
                                   halfway * u.dimensionless_unscaled) ==
                      [0, 0.25, 0.75, +1.] * u.dimensionless_unscaled)

    @pytest.mark.parametrize('function', (np.cbrt, ))
    def test_cbrt_scalar(self, function):
        assert function(8. * u.m**3) == 2. * u.m

    @pytest.mark.parametrize('function', (np.cbrt, ))
    def test_cbrt_array(self, function):
        # Calculate cbrt on both sides since on Windows the cube root of 64
        # does not exactly equal 4.  See 4388.
        values = np.array([1., 8., 64.])
        assert np.all(function(values * u.m**3) ==
                      function(values) * u.m)

    def test_power_scalar(self):
        assert np.power(4. * u.m, 2.) == 16. * u.m ** 2
        assert np.power(4., 200. * u.cm / u.m) == \
            u.Quantity(16., u.dimensionless_unscaled)
        # regression check on #1696
        assert np.power(4. * u.m, 0.) == 1. * u.dimensionless_unscaled

    def test_power_array(self):
        assert np.all(np.power(np.array([1., 2., 3.]) * u.m, 3.)
                      == np.array([1., 8., 27.]) * u.m ** 3)
        # regression check on #1696
        assert np.all(np.power(np.arange(4.) * u.m, 0.) ==
                      1. * u.dimensionless_unscaled)

    # float_power only introduced in numpy 1.12
    @pytest.mark.skipif("not hasattr(np, 'float_power')")
    def test_float_power_array(self):
        assert np.all(np.float_power(np.array([1., 2., 3.]) * u.m, 3.)
                      == np.array([1., 8., 27.]) * u.m ** 3)
        # regression check on #1696
        assert np.all(np.float_power(np.arange(4.) * u.m, 0.) ==
                      1. * u.dimensionless_unscaled)

    @raises(ValueError)
    def test_power_array_array(self):
        np.power(4. * u.m, [2., 4.])

    @raises(ValueError)
    def test_power_array_array2(self):
        np.power([2., 4.] * u.m, [2., 4.])

    def test_power_array_array3(self):
        # Identical unit fractions are converted automatically to dimensionless
        # and should be allowed as base for np.power: #4764
        q = [2., 4.] * u.m / u.m
        powers = [2., 4.]
        res = np.power(q, powers)
        assert np.all(res.value == q.value ** powers)
        assert res.unit == u.dimensionless_unscaled
        # The same holds for unit fractions that are scaled dimensionless.
        q2 = [2., 4.] * u.m / u.cm
        # Test also against different types of exponent
        for cls in (list, tuple, np.array, np.ma.array, u.Quantity):
            res2 = np.power(q2, cls(powers))
            assert np.all(res2.value == q2.to_value(1) ** powers)
            assert res2.unit == u.dimensionless_unscaled
        # Though for single powers, we keep the composite unit.
        res3 = q2 ** 2
        assert np.all(res3.value == q2.value ** 2)
        assert res3.unit == q2.unit ** 2
        assert np.all(res3 == q2 ** [2, 2])

    def test_power_invalid(self):
        with pytest.raises(TypeError) as exc:
            np.power(3., 4. * u.m)
        assert "raise something to a dimensionless" in exc.value.args[0]

    def test_copysign_scalar(self):
        assert np.copysign(3 * u.m, 1.) == 3. * u.m
        assert np.copysign(3 * u.m, 1. * u.s) == 3. * u.m
        assert np.copysign(3 * u.m, -1.) == -3. * u.m
        assert np.copysign(3 * u.m, -1. * u.s) == -3. * u.m

    def test_copysign_array(self):
        assert np.all(np.copysign(np.array([1., 2., 3.]) * u.s, -1.) ==
                      -np.array([1., 2., 3.]) * u.s)
        assert np.all(np.copysign(np.array([1., 2., 3.]) * u.s, -1. * u.m) ==
                      -np.array([1., 2., 3.]) * u.s)
        assert np.all(np.copysign(np.array([1., 2., 3.]) * u.s,
                                  np.array([-2., 2., -4.]) * u.m) ==
                      np.array([-1., 2., -3.]) * u.s)

        q = np.copysign(np.array([1., 2., 3.]), -3 * u.m)
        assert np.all(q == np.array([-1., -2., -3.]))
        assert not isinstance(q, u.Quantity)

    def test_ldexp_scalar(self):
        assert np.ldexp(4. * u.m, 2) == 16. * u.m

    def test_ldexp_array(self):
        assert np.all(np.ldexp(np.array([1., 2., 3.]) * u.m, [3, 2, 1])
                      == np.array([8., 8., 6.]) * u.m)

    def test_ldexp_invalid(self):
        with pytest.raises(TypeError):
            np.ldexp(3. * u.m, 4.)

        with pytest.raises(TypeError):
            np.ldexp(3., u.Quantity(4, u.m, dtype=int))

    @pytest.mark.parametrize('function', (np.exp, np.expm1, np.exp2,
                                          np.log, np.log2, np.log10, np.log1p))
    def test_exp_scalar(self, function):
        q = function(3. * u.m / (6. * u.m))
        assert q.unit == u.dimensionless_unscaled
        assert q.value == function(0.5)

    @pytest.mark.parametrize('function', (np.exp, np.expm1, np.exp2,
                                          np.log, np.log2, np.log10, np.log1p))
    def test_exp_array(self, function):
        q = function(np.array([2., 3., 6.]) * u.m / (6. * u.m))
        assert q.unit == u.dimensionless_unscaled
        assert np.all(q.value
                      == function(np.array([1. / 3., 1. / 2., 1.])))
        # should also work on quantities that can be made dimensionless
        q2 = function(np.array([2., 3., 6.]) * u.m / (6. * u.cm))
        assert q2.unit == u.dimensionless_unscaled
        assert_allclose(q2.value,
                        function(np.array([100. / 3., 100. / 2., 100.])))

    @pytest.mark.parametrize('function', (np.exp, np.expm1, np.exp2,
                                          np.log, np.log2, np.log10, np.log1p))
    def test_exp_invalid_units(self, function):
        # Can't use exp() with non-dimensionless quantities
        with pytest.raises(TypeError) as exc:
            function(3. * u.m / u.s)
        assert exc.value.args[0] == ("Can only apply '{0}' function to "
                                     "dimensionless quantities"
                                     .format(function.__name__))

    def test_modf_scalar(self):
        q = np.modf(9. * u.m / (600. * u.cm))
        assert q == (0.5 * u.dimensionless_unscaled,
                     1. * u.dimensionless_unscaled)

    def test_modf_array(self):
        v = np.arange(10.) * u.m / (500. * u.cm)
        q = np.modf(v)
        n = np.modf(v.to_value(u.dimensionless_unscaled))
        assert q[0].unit == u.dimensionless_unscaled
        assert q[1].unit == u.dimensionless_unscaled
        assert all(q[0].value == n[0])
        assert all(q[1].value == n[1])

    def test_frexp_scalar(self):
        q = np.frexp(3. * u.m / (6. * u.m))
        assert q == (np.array(0.5), np.array(0.0))

    def test_frexp_array(self):
        q = np.frexp(np.array([2., 3., 6.]) * u.m / (6. * u.m))
        assert all((_q0, _q1) == np.frexp(_d) for _q0, _q1, _d
                   in zip(q[0], q[1], [1. / 3., 1. / 2., 1.]))

    def test_frexp_invalid_units(self):
        # Can't use prod() with non-dimensionless quantities
        with pytest.raises(TypeError) as exc:
            np.frexp(3. * u.m / u.s)
        assert exc.value.args[0] == ("Can only apply 'frexp' function to "
                                     "unscaled dimensionless quantities")

        # also does not work on quantities that can be made dimensionless
        with pytest.raises(TypeError) as exc:
            np.frexp(np.array([2., 3., 6.]) * u.m / (6. * u.cm))
        assert exc.value.args[0] == ("Can only apply 'frexp' function to "
                                     "unscaled dimensionless quantities")

    @pytest.mark.parametrize('function', (np.logaddexp, np.logaddexp2))
    def test_dimensionless_twoarg_array(self, function):
        q = function(np.array([2., 3., 6.]) * u.m / (6. * u.cm), 1.)
        assert q.unit == u.dimensionless_unscaled
        assert_allclose(q.value,
                        function(np.array([100. / 3., 100. / 2., 100.]), 1.))

    @pytest.mark.parametrize('function', (np.logaddexp, np.logaddexp2))
    def test_dimensionless_twoarg_invalid_units(self, function):

        with pytest.raises(TypeError) as exc:
            function(1. * u.km / u.s, 3. * u.m / u.s)
        assert exc.value.args[0] == ("Can only apply '{0}' function to "
                                     "dimensionless quantities"
                                     .format(function.__name__))


class TestInvariantUfuncs:

    @pytest.mark.parametrize(('ufunc'), [np.absolute, np.fabs,
                                         np.conj, np.conjugate,
                                         np.negative, np.spacing, np.rint,
                                         np.floor, np.ceil, np.positive])
    def test_invariant_scalar(self, ufunc):

        q_i = 4.7 * u.m
        q_o = ufunc(q_i)
        assert isinstance(q_o, u.Quantity)
        assert q_o.unit == q_i.unit
        assert q_o.value == ufunc(q_i.value)

    @pytest.mark.parametrize(('ufunc'), [np.absolute, np.conjugate,
                                         np.negative, np.rint,
                                         np.floor, np.ceil])
    def test_invariant_array(self, ufunc):

        q_i = np.array([-3.3, 2.1, 10.2]) * u.kg / u.s
        q_o = ufunc(q_i)
        assert isinstance(q_o, u.Quantity)
        assert q_o.unit == q_i.unit
        assert np.all(q_o.value == ufunc(q_i.value))

    @pytest.mark.parametrize(('ufunc'), [np.add, np.subtract, np.hypot,
                                         np.maximum, np.minimum, np.nextafter,
                                         np.remainder, np.mod, np.fmod])
    def test_invariant_twoarg_scalar(self, ufunc):

        q_i1 = 4.7 * u.m
        q_i2 = 9.4 * u.km
        q_o = ufunc(q_i1, q_i2)
        assert isinstance(q_o, u.Quantity)
        assert q_o.unit == q_i1.unit
        assert_allclose(q_o.value, ufunc(q_i1.value, q_i2.to_value(q_i1.unit)))

    @pytest.mark.parametrize(('ufunc'), [np.add, np.subtract, np.hypot,
                                         np.maximum, np.minimum, np.nextafter,
                                         np.remainder, np.mod, np.fmod])
    def test_invariant_twoarg_array(self, ufunc):

        q_i1 = np.array([-3.3, 2.1, 10.2]) * u.kg / u.s
        q_i2 = np.array([10., -5., 1.e6]) * u.g / u.us
        q_o = ufunc(q_i1, q_i2)
        assert isinstance(q_o, u.Quantity)
        assert q_o.unit == q_i1.unit
        assert_allclose(q_o.value, ufunc(q_i1.value, q_i2.to_value(q_i1.unit)))

    @pytest.mark.parametrize(('ufunc'), [np.add, np.subtract, np.hypot,
                                         np.maximum, np.minimum, np.nextafter,
                                         np.remainder, np.mod, np.fmod])
    def test_invariant_twoarg_one_arbitrary(self, ufunc):

        q_i1 = np.array([-3.3, 2.1, 10.2]) * u.kg / u.s
        arbitrary_unit_value = np.array([0.])
        q_o = ufunc(q_i1, arbitrary_unit_value)
        assert isinstance(q_o, u.Quantity)
        assert q_o.unit == q_i1.unit
        assert_allclose(q_o.value, ufunc(q_i1.value, arbitrary_unit_value))

    @pytest.mark.parametrize(('ufunc'), [np.add, np.subtract, np.hypot,
                                         np.maximum, np.minimum, np.nextafter,
                                         np.remainder, np.mod, np.fmod])
    def test_invariant_twoarg_invalid_units(self, ufunc):

        q_i1 = 4.7 * u.m
        q_i2 = 9.4 * u.s
        with pytest.raises(u.UnitsError) as exc:
            ufunc(q_i1, q_i2)
        assert "compatible dimensions" in exc.value.args[0]


class TestComparisonUfuncs:

    @pytest.mark.parametrize(('ufunc'), [np.greater, np.greater_equal,
                                         np.less, np.less_equal,
                                         np.not_equal, np.equal])
    def test_comparison_valid_units(self, ufunc):
        q_i1 = np.array([-3.3, 2.1, 10.2]) * u.kg / u.s
        q_i2 = np.array([10., -5., 1.e6]) * u.g / u.Ms
        q_o = ufunc(q_i1, q_i2)
        assert not isinstance(q_o, u.Quantity)
        assert q_o.dtype == bool
        assert np.all(q_o == ufunc(q_i1.value, q_i2.to_value(q_i1.unit)))
        q_o2 = ufunc(q_i1 / q_i2, 2.)
        assert not isinstance(q_o2, u.Quantity)
        assert q_o2.dtype == bool
        assert np.all(q_o2 == ufunc((q_i1 / q_i2)
                                    .to_value(u.dimensionless_unscaled), 2.))
        # comparison with 0., inf, nan is OK even for dimensional quantities
        for arbitrary_unit_value in (0., np.inf, np.nan):
            ufunc(q_i1, arbitrary_unit_value)
            ufunc(q_i1, arbitrary_unit_value*np.ones(len(q_i1)))
        # and just for completeness
        ufunc(q_i1, np.array([0., np.inf, np.nan]))

    @pytest.mark.parametrize(('ufunc'), [np.greater, np.greater_equal,
                                         np.less, np.less_equal,
                                         np.not_equal, np.equal])
    def test_comparison_invalid_units(self, ufunc):
        q_i1 = 4.7 * u.m
        q_i2 = 9.4 * u.s
        with pytest.raises(u.UnitsError) as exc:
            ufunc(q_i1, q_i2)
        assert "compatible dimensions" in exc.value.args[0]


class TestInplaceUfuncs:

    @pytest.mark.parametrize(('value'), [1., np.arange(10.)])
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

    @pytest.mark.parametrize(('value'), [1., np.arange(10.)])
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
        np.arcsin(value/10., out=s)
        assert check is s
        assert np.all(check.value == np.arcsin(value/10.))
        assert check.unit is u.radian

    @pytest.mark.parametrize(('value'), [1., np.arange(10.)])
    def test_one_argument_two_output_ufunc_inplace(self, value):
        v = 100. * value * u.cm / u.m
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
        # And now, with numpy >= 1.13, one can also replace input with
        # first output when scaling
        v4 = v_copy.copy()
        check4 = v4
        np.modf(v4, v4, tmp)
        assert check4 is v4
        assert check4.unit == u.dimensionless_unscaled

    @pytest.mark.parametrize(('value'), [1., np.arange(10.)])
    def test_two_argument_ufunc_inplace_1(self, value):
        s = value * u.cycle
        check = s
        s /= 2.
        assert check is s
        assert np.all(check.value == value / 2.)
        s /= u.s
        assert check is s
        assert check.unit == u.cycle / u.s
        s *= 2. * u.s
        assert check is s
        assert np.all(check == value * u.cycle)

    @pytest.mark.parametrize(('value'), [1., np.arange(10.)])
    def test_two_argument_ufunc_inplace_2(self, value):
        s = value * u.cycle
        check = s
        np.arctan2(s, s, out=s)
        assert check is s
        assert check.unit == u.radian
        with pytest.raises(u.UnitsError):
            s += 1. * u.m
        assert check is s
        assert check.unit == u.radian
        np.arctan2(1. * u.deg, s, out=s)
        assert check is s
        assert check.unit == u.radian
        np.add(1. * u.deg, s, out=s)
        assert check is s
        assert check.unit == u.deg
        np.multiply(2. / u.s, s, out=s)
        assert check is s
        assert check.unit == u.deg / u.s

    def test_two_argument_ufunc_inplace_3(self):
        s = np.array([1., 2., 3.]) * u.dimensionless_unscaled
        np.add(np.array([1., 2., 3.]), np.array([1., 2., 3.]) * 2., out=s)
        assert np.all(s.value == np.array([3., 6., 9.]))
        assert s.unit is u.dimensionless_unscaled
        np.arctan2(np.array([1., 2., 3.]), np.array([1., 2., 3.]) * 2., out=s)
        assert_allclose(s.value, np.arctan2(1., 2.))
        assert s.unit is u.radian

    @pytest.mark.parametrize(('value'), [1., np.arange(10.)])
    def test_two_argument_two_output_ufunc_inplace(self, value):
        v = value * u.m
        divisor = 70.*u.cm
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
        s = np.arange(10.) * u.m
        s_copy = s.copy()
        s2 = s[::2]
        s2 += 1. * u.cm
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
        a2 += (20.*u.km)
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


class TestUfuncAt:
    """Test that 'at' method for ufuncs (calculates in-place at given indices)

    For Quantities, since calculations are in-place, it makes sense only
    if the result is still a quantity, and if the unit does not have to change
    """

    def test_one_argument_ufunc_at(self):
        q = np.arange(10.) * u.m
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
        d = np.arange(10.) * u.dimensionless_unscaled
        dv = d.value.copy()
        np.square.at(d, i)
        np.square.at(dv, i)
        assert np.all(d.value == dv)
        assert d.unit is u.dimensionless_unscaled

        d = np.arange(10.) * u.dimensionless_unscaled
        dv = d.value.copy()
        np.log.at(d, i)
        np.log.at(dv, i)
        assert np.all(d.value == dv)
        assert d.unit is u.dimensionless_unscaled

        # also for sine it doesn't work, even if given an angle
        a = np.arange(10.) * u.radian
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
            ad = np.arange(10.) * u.degree
            with pytest.raises(u.UnitsError):
                np.sin.at(ad, i)

    def test_two_argument_ufunc_at(self):
        s = np.arange(10.) * u.m
        i = np.array([1, 2])
        check = s.value.copy()
        np.add.at(s, i, 1.*u.km)
        np.add.at(check, i, 1000.)
        assert np.all(s.value == check)
        assert s.unit is u.m

        with pytest.raises(u.UnitsError):
            np.add.at(s, i, 1.*u.s)

        # also raise UnitsError if unit would have to be changed
        with pytest.raises(u.UnitsError):
            np.multiply.at(s, i, 1*u.s)

        # but be fine if it does not
        s = np.arange(10.) * u.m
        check = s.value.copy()
        np.multiply.at(s, i, 2.*u.dimensionless_unscaled)
        np.multiply.at(check, i, 2)
        assert np.all(s.value == check)
        s = np.arange(10.) * u.m
        np.multiply.at(s, i, 2.)
        assert np.all(s.value == check)

        # of course cannot change class of data either
        with pytest.raises(TypeError):
            np.greater.at(s, i, 1.*u.km)


class TestUfuncReduceReduceatAccumulate:
    """Test 'reduce', 'reduceat' and 'accumulate' methods for ufuncs

    For Quantities, it makes sense only if the unit does not have to change
    """

    def test_one_argument_ufunc_reduce_accumulate(self):
        # one argument cannot be used
        s = np.arange(10.) * u.radian
        i = np.array([0, 5, 1, 6])
        with pytest.raises(ValueError):
            np.sin.reduce(s)
        with pytest.raises(ValueError):
            np.sin.accumulate(s)
        with pytest.raises(ValueError):
            np.sin.reduceat(s, i)

    def test_two_argument_ufunc_reduce_accumulate(self):
        s = np.arange(10.) * u.m
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
        s = np.arange(10.) * u.dimensionless_unscaled
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
        s = np.arange(10.) * u.radian
        with pytest.raises(ValueError):
            np.sin.outer(s)

    def test_two_argument_ufunc_outer(self):
        s1 = np.arange(10.) * u.m
        s2 = np.arange(2.) * u.s
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
        s3 = np.arange(2.) * s1.unit
        check3 = s3.value
        s13_add_outer = np.add.outer(s1, s3)
        check13_add_outer = np.add.outer(check1, check3)
        assert np.all(s13_add_outer.value == check13_add_outer)
        assert s13_add_outer.unit is s1.unit

        s13_greater_outer = np.greater.outer(s1, s3)
        check13_greater_outer = np.greater.outer(check1, check3)
        assert type(s13_greater_outer) is np.ndarray
        assert np.all(s13_greater_outer == check13_greater_outer)


if HAS_SCIPY:
    from scipy import special as sps

    def test_scipy_registration():
        """Check that scipy gets loaded upon first use."""
        assert sps.erf not in qh.UFUNC_HELPERS
        sps.erf(1. * u.percent)
        assert sps.erf in qh.UFUNC_HELPERS

    class TestScipySpecialUfuncs:

        erf_like_ufuncs = (
            sps.erf, sps.gamma, sps.loggamma, sps.gammasgn, sps.psi,
            sps.rgamma, sps.erfc, sps.erfcx, sps.erfi, sps.wofz, sps.dawsn,
            sps.entr, sps.exprel, sps.expm1, sps.log1p, sps.exp2, sps.exp10)

        @pytest.mark.parametrize('function', erf_like_ufuncs)
        def test_erf_scalar(self, function):
            TestQuantityMathFuncs.test_exp_scalar(None, function)

        @pytest.mark.parametrize('function', erf_like_ufuncs)
        def test_erf_array(self, function):
            TestQuantityMathFuncs.test_exp_array(None, function)

        @pytest.mark.parametrize('function', erf_like_ufuncs)
        def test_erf_invalid_units(self, function):
            TestQuantityMathFuncs.test_exp_invalid_units(None, function)

        @pytest.mark.parametrize('function', (sps.cbrt, ))
        def test_cbrt_scalar(self, function):
            TestQuantityMathFuncs.test_cbrt_scalar(None, function)

        @pytest.mark.parametrize('function', (sps.cbrt, ))
        def test_cbrt_array(self, function):
            TestQuantityMathFuncs.test_cbrt_array(None, function)

        @pytest.mark.parametrize('function', (sps.radian, ))
        def test_radian(self, function):
            q1 = function(180. * u.degree, 0. * u.arcmin, 0. * u.arcsec)
            assert_allclose(q1.value, np.pi)
            assert q1.unit == u.radian

            q2 = function(0. * u.degree, 30. * u.arcmin, 0. * u.arcsec)
            assert_allclose(q2.value, (30. * u.arcmin).to(u.radian).value)
            assert q2.unit == u.radian

            q3 = function(0. * u.degree, 0. * u.arcmin, 30. * u.arcsec)
            assert_allclose(q3.value, (30. * u.arcsec).to(u.radian).value)

            # the following doesn't make much sense in terms of the name of the
            # routine, but we check it gives the correct result.
            q4 = function(3. * u.radian, 0. * u.arcmin, 0. * u.arcsec)
            assert_allclose(q4.value, 3.)
            assert q4.unit == u.radian

            with pytest.raises(TypeError):
                function(3. * u.m, 2. * u.s, 1. * u.kg)

        jv_like_ufuncs = (
            sps.jv, sps.jn, sps.jve, sps.yn, sps.yv, sps.yve, sps.kn, sps.kv,
            sps.kve, sps.iv, sps.ive, sps.hankel1, sps.hankel1e, sps.hankel2,
            sps.hankel2e)

        @pytest.mark.parametrize('function', jv_like_ufuncs)
        def test_jv_scalar(self, function):
            q = function(2. * u.m / (2. * u.m), 3. * u.m / (6. * u.m))
            assert q.unit == u.dimensionless_unscaled
            assert q.value == function(1.0, 0.5)

        @pytest.mark.parametrize('function', jv_like_ufuncs)
        def test_jv_array(self, function):
            q = function(np.ones(3) * u.m / (1. * u.m),
                         np.array([2., 3., 6.]) * u.m / (6. * u.m))
            assert q.unit == u.dimensionless_unscaled
            assert np.all(q.value == function(
                np.ones(3),
                np.array([1. / 3., 1. / 2., 1.]))
            )
            # should also work on quantities that can be made dimensionless
            q2 = function(np.ones(3) * u.m / (1. * u.m),
                          np.array([2., 3., 6.]) * u.m / (6. * u.cm))
            assert q2.unit == u.dimensionless_unscaled
            assert_allclose(q2.value,
                            function(np.ones(3),
                                     np.array([100. / 3., 100. / 2., 100.])))

        @pytest.mark.parametrize('function', jv_like_ufuncs)
        def test_jv_invalid_units(self, function):
            # Can't use jv() with non-dimensionless quantities
            with pytest.raises(TypeError) as exc:
                function(1. * u.kg, 3. * u.m / u.s)
            assert exc.value.args[0] == ("Can only apply '{0}' function to "
                                         "dimensionless quantities"
                                         .format(function.__name__))
