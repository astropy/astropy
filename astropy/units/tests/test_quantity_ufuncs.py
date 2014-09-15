# The purpose of these tests are to ensure that calling ufuncs with quantities
# returns quantities with the right units, or raises exceptions.

import numpy as np
from numpy.testing.utils import assert_allclose

from ... import units as u
from ...tests.helper import pytest, raises
from ...utils.compat import NUMPY_LT_1_10


class TestUfuncCoverage(object):
    """Test that we cover all ufunc's"""
    def test_coverage(self):
        all_np_ufuncs = set([ufunc for ufunc in np.core.umath.__dict__.values()
                             if type(ufunc) == np.ufunc])

        # in numpy >=1.10, with __numpy_ufunc__, np.dot behaves like a ufunc.
        if not NUMPY_LT_1_10:
            all_np_ufuncs |= set([np.dot])

        from .. import quantity_helper as qh

        all_q_ufuncs = (qh.UNSUPPORTED_UFUNCS |
                        set(qh.UFUNC_HELPERS.keys()))

        assert all_np_ufuncs - all_q_ufuncs == set([])
        assert all_q_ufuncs - all_np_ufuncs == set([])


class TestQuantityTrigonometricFuncs(object):
    """
    Test trigonometric functions
    """

    def test_sin_scalar(self):
        q = np.sin(30. * u.degree)
        assert q.unit == u.dimensionless_unscaled
        assert_allclose(q.value, 0.5)

    def test_sin_array(self):
        q = np.sin(np.array([0., np.pi / 4., np.pi / 2.]) * u.radian)
        assert q.unit == u.dimensionless_unscaled
        assert_allclose(q.value,
                        np.array([0., 1. / np.sqrt(2.), 1.]), atol=1.e-15)

    def test_arcsin_scalar(self):
        q1 = 30. * u.degree
        q2 = np.arcsin(np.sin(q1)).to(q1.unit)
        assert_allclose(q1.value, q2.value)

    def test_arcsin_array(self):
        q1 = np.array([0., np.pi / 4., np.pi / 2.]) * u.radian
        q2 = np.arcsin(np.sin(q1)).to(q1.unit)
        assert_allclose(q1.value, q2.value)

    def test_sin_invalid_units(self):
        with pytest.raises(TypeError) as exc:
            np.sin(3. * u.m)
        assert exc.value.args[0] == ("Can only apply 'sin' function "
                                     "to quantities with angle units")

    def test_arcsin_invalid_units(self):
        with pytest.raises(TypeError) as exc:
            np.arcsin(3. * u.m)
        assert exc.value.args[0] == ("Can only apply 'arcsin' function to "
                                     "dimensionless quantities")

    def test_cos_scalar(self):
        q = np.cos(np.pi / 3. * u.radian)
        assert q.unit == u.dimensionless_unscaled
        assert_allclose(q.value, 0.5)

    def test_cos_array(self):
        q = np.cos(np.array([0., np.pi / 4., np.pi / 2.]) * u.radian)
        assert q.unit == u.dimensionless_unscaled
        assert_allclose(q.value,
                        np.array([1., 1. / np.sqrt(2.), 0.]), atol=1.e-15)

    def test_arccos_scalar(self):
        q1 = np.pi / 3. * u.radian
        q2 = np.arccos(np.cos(q1)).to(q1.unit)
        assert_allclose(q1.value, q2.value)

    def test_arccos_array(self):
        q1 = np.array([0., np.pi / 4., np.pi / 2.]) * u.radian
        q2 = np.arccos(np.cos(q1)).to(q1.unit)
        assert_allclose(q1.value, q2.value)

    def test_cos_invalid_units(self):
        with pytest.raises(TypeError) as exc:
            np.cos(3. * u.s)
        assert exc.value.args[0] == ("Can only apply 'cos' function "
                                     "to quantities with angle units")

    def test_arccos_invalid_units(self):
        with pytest.raises(TypeError) as exc:
            np.arccos(3. * u.s)
        assert exc.value.args[0] == ("Can only apply 'arccos' function to "
                                     "dimensionless quantities")

    def test_tan_scalar(self):
        q = np.tan(np.pi / 3. * u.radian)
        assert q.unit == u.dimensionless_unscaled
        assert_allclose(q.value, np.sqrt(3.))

    def test_tan_array(self):
        q = np.tan(np.array([0., 45., 135., 180.]) * u.degree)
        assert q.unit == u.dimensionless_unscaled
        assert_allclose(q.value,
                        np.array([0., 1., -1., 0.]), atol=1.e-15)

    def test_arctan_scalar(self):
        q = np.pi / 3. * u.radian
        assert np.arctan(np.tan(q))

    def test_arctan_array(self):
        q = np.array([10., 30., 70., 80.]) * u.degree
        assert_allclose(np.arctan(np.tan(q)).to(q.unit).value, q.value)

    def test_tan_invalid_units(self):
        with pytest.raises(TypeError) as exc:
            np.tan(np.array([1, 2, 3]) * u.N)
        assert exc.value.args[0] == ("Can only apply 'tan' function "
                                     "to quantities with angle units")

    def test_arctan_invalid_units(self):
        with pytest.raises(TypeError) as exc:
            np.arctan(np.array([1, 2, 3]) * u.N)
        assert exc.value.args[0] == ("Can only apply 'arctan' function to "
                                     "dimensionless quantities")

    def test_arctan2_valid(self):
        q1 = np.array([10., 30., 70., 80.]) * u.m
        q2 = 2.0 * u.km
        assert np.arctan2(q1, q2).unit == u.radian
        assert_allclose(np.arctan2(q1, q2).value,
                        np.arctan2(q1.value, q2.to(q1.unit).value))
        q3 = q1 / q2
        q4 = 1.
        at2 = np.arctan2(q3, q4)
        assert_allclose(at2.value, np.arctan2(q3.to(1).value, q4))

    def test_arctan2_invalid(self):
        with pytest.raises(u.UnitsError) as exc:
            np.arctan2(np.array([1, 2, 3]) * u.N, 1. * u.s)
        assert "compatible dimensions" in exc.value.args[0]
        with pytest.raises(u.UnitsError) as exc:
            np.arctan2(np.array([1, 2, 3]) * u.N, 1.)
        assert "dimensionless quantities when other arg" in exc.value.args[0]

    def test_radians(self):

        q1 = np.deg2rad(180. * u.degree)
        assert_allclose(q1.value, np.pi)
        assert q1.unit == u.radian

        q2 = np.radians(180. * u.degree)
        assert_allclose(q2.value, np.pi)
        assert q2.unit == u.radian

        # the following doesn't make much sense in terms of the name of the
        # routine, but we check it gives the correct result.
        q3 = np.deg2rad(3. * u.radian)
        assert_allclose(q3.value, 3.)
        assert q3.unit == u.radian

        q4 = np.radians(3. * u.radian)
        assert_allclose(q4.value, 3.)
        assert q4.unit == u.radian

        with pytest.raises(TypeError):
            np.deg2rad(3. * u.m)

        with pytest.raises(TypeError):
            np.radians(3. * u.m)

    def test_degrees(self):

        # the following doesn't make much sense in terms of the name of the
        # routine, but we check it gives the correct result.
        q1 = np.rad2deg(60. * u.degree)
        assert_allclose(q1.value, 60.)
        assert q1.unit == u.degree

        q2 = np.degrees(60. * u.degree)
        assert_allclose(q2.value, 60.)
        assert q2.unit == u.degree

        q3 = np.rad2deg(np.pi * u.radian)
        assert_allclose(q3.value, 180.)
        assert q3.unit == u.degree

        q4 = np.degrees(np.pi * u.radian)
        assert_allclose(q4.value, 180.)
        assert q4.unit == u.degree

        with pytest.raises(TypeError):
            np.rad2deg(3. * u.m)

        with pytest.raises(TypeError):
            np.degrees(3. * u.m)


class TestQuantityMathFuncs(object):
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

    @pytest.mark.parametrize('function', (np.divide, np.true_divide,
                                          np.floor_divide))
    def test_divide_scalar(self, function):
        assert function(4. * u.m, 2. * u.s) == function(4., 2.) * u.m / u.s
        assert function(4. * u.m, 2.) == function(4., 2.) * u.m
        assert function(4., 2. * u.s) == function(4., 2.) / u.s

    @pytest.mark.parametrize('function', (np.divide, np.true_divide,
                                          np.floor_divide))
    def test_divide_array(self, function):
        assert np.all(function(np.arange(3.) * u.m, 2. * u.s) ==
                      function(np.arange(3.), 2.) * u.m / u.s)

    def test_divmod(self):
        inch = u.Unit(0.0254 * u.m)
        quotient, remainder = divmod(
            np.array([1., 2., 3.]) * u.m,
            np.array([3., 4., 5.]) * inch)
        assert_allclose(quotient.value, [13., 19., 23.])
        assert quotient.unit == u.dimensionless_unscaled
        assert_allclose(remainder.value, [0.0094, 0.0696, 0.079])
        assert remainder.unit == u.m

        quotient, remainder = divmod(
            np.array([1., 2., 3.]) * u.m, u.km)
        assert_allclose(quotient.value, [1., 2., 3.])
        assert quotient.unit == u.m / u.km
        assert remainder.value == 0.
        assert remainder.unit == u.dimensionless_unscaled

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

    # cbrt only introduced in numpy 1.10
    @pytest.mark.skipif("not hasattr(np, 'cbrt')")
    def test_cbrt_scalar(self):
        assert np.cbrt(8. * u.m**3) == 2. * u.m

    @pytest.mark.skipif("not hasattr(np, 'cbrt')")
    def test_cbrt_array(self):
        assert np.all(np.cbrt(np.array([1., 8., 64.]) * u.m**3)
                      == np.array([1., 2., 4.]) * u.m)

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

    @raises(ValueError)
    def test_power_array_array(self):
        np.power(4. * u.m, [2., 4.])

    @raises(ValueError)
    def test_power_array_array2(self):
        np.power([2., 4.] * u.m, [2., 4.])

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
        assert np.all(np.copysign(np.array([1., 2., 3.]) * u.s, -1.) == -np.array([1., 2., 3.]) * u.s)
        assert np.all(np.copysign(np.array([1., 2., 3.]) * u.s, -1. * u.m) == -np.array([1., 2., 3.]) * u.s)
        assert np.all(np.copysign(np.array([1., 2., 3.]) * u.s, np.array([-2.,2.,-4.]) * u.m) == np.array([-1., 2., -3.]) * u.s)

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
        n = np.modf(v.to(1).value)
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


class TestInvariantUfuncs(object):

    @pytest.mark.parametrize(('ufunc'), [np.absolute, np.fabs,
                                         np.conj, np.conjugate,
                                         np.negative, np.spacing, np.rint,
                                         np.floor, np.ceil])
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
        assert_allclose(q_o.value, ufunc(q_i1.value, q_i2.to(q_i1.unit).value))

    @pytest.mark.parametrize(('ufunc'), [np.add, np.subtract, np.hypot,
                                         np.maximum, np.minimum, np.nextafter,
                                         np.remainder, np.mod, np.fmod])
    def test_invariant_twoarg_array(self, ufunc):

        q_i1 = np.array([-3.3, 2.1, 10.2]) * u.kg / u.s
        q_i2 = np.array([10., -5., 1.e6]) * u.g / u.us
        q_o = ufunc(q_i1, q_i2)
        assert isinstance(q_o, u.Quantity)
        assert q_o.unit == q_i1.unit
        assert_allclose(q_o.value, ufunc(q_i1.value, q_i2.to(q_i1.unit).value))

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


class TestComparisonUfuncs(object):

    @pytest.mark.parametrize(('ufunc'), [np.greater, np.greater_equal,
                                         np.less, np.less_equal,
                                         np.not_equal, np.equal])
    def test_comparison_valid_units(self, ufunc):
        q_i1 = np.array([-3.3, 2.1, 10.2]) * u.kg / u.s
        q_i2 = np.array([10., -5., 1.e6]) * u.g / u.Ms
        q_o = ufunc(q_i1, q_i2)
        assert not isinstance(q_o, u.Quantity)
        assert q_o.dtype == np.bool
        assert np.all(q_o == ufunc(q_i1.value, q_i2.to(q_i1.unit).value))
        q_o2 = ufunc(q_i1 / q_i2, 2.)
        assert not isinstance(q_o2, u.Quantity)
        assert q_o2.dtype == np.bool
        assert np.all(q_o2 == ufunc((q_i1 / q_i2).to(1).value, 2.))
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


class TestInplaceUfuncs(object):

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
        np.modf(v, tmp, v)  # cannot use out1,out2 keywords with numpy 1.7
        assert check is v
        assert check.unit == u.dimensionless_unscaled
        v2 = v_copy.to(1)
        check2 = v2
        np.modf(v2, tmp, v2)
        assert check2 is v2
        assert check2.unit == u.dimensionless_unscaled
        # can also replace in last position if no scaling is needed
        v3 = v_copy.to(1)
        check3 = v3
        np.modf(v3, v3, tmp)
        assert check3 is v3
        assert check3.unit == u.dimensionless_unscaled
        # in np<1.10, cannot replace input with first output when scaling
        v4 = v_copy.copy()
        if NUMPY_LT_1_10:
            with pytest.raises(TypeError):
                np.modf(v4, v4, tmp)
        else:
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

    def test_ufunc_inplace_non_contiguous_data(self):
        # ensure inplace works also for non-contiguous data (closes #1834)
        s = np.arange(10.) * u.m
        s_copy = s.copy()
        s2 = s[::2]
        s2 += 1. * u.cm
        assert np.all(s[::2] > s_copy[::2])
        assert np.all(s[1::2] == s_copy[1::2])
