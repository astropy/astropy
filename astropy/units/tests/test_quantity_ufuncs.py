# The purpose of these tests are to ensure that calling ufuncs with quantities
# returns quantities with the right units, or raises exceptions.

import numpy as np

from ... import units as u
from ...tests.helper import pytest
from ...tests.compat import assert_allclose


class TestUfuncCoverage(object):
    """Test that we cover all ufunc's"""
    def test_coverage(self):
        q = u.quantity
        all_np_ufuncs = set([np.__dict__[i] for i in np.__dict__
                             if getattr(np.__dict__[i],
                                        '__class__', None) == np.ufunc])
        all_q_ufuncs = (q.UNSUPPORTED_UFUNCS |
                        q.SPECIAL_UFUNCS |
                        q.ANGLE_UFUNCS |
                        q.DIMENSIONLESS_UFUNCS |
                        q.IS_UNITY_UFUNCS |
                        q.INVARIANT_UFUNCS |
                        q.TEST_UFUNCS |
                        q.SPECIAL_TWOARG_UFUNCS |
                        q.INVTRIG_TWOARG_UFUNCS |
                        q.INVARIANT_TWOARG_UFUNCS |
                        q.DIMENSIONLESS_TWOARG_UFUNCS |
                        q.COMPARISON_UFUNCS)
        assert all_np_ufuncs - all_q_ufuncs == set([])
        assert all_q_ufuncs - all_np_ufuncs == set([])


class TestQuantityStatsFuncs(object):
    """
    Test statistical functions
    """

    def test_mean(self):
        q1 = np.array([1.,2.,4.,5.,6.]) * u.m
        assert np.mean(q1) == 3.6 * u.m
        assert np.mean(q1).unit == u.m
        assert np.mean(q1).value == 3.6

    def test_std(self):
        q1 = np.array([1.,2.]) * u.m
        assert np.std(q1) == 0.5 * u.m
        assert np.std(q1).unit == u.m
        assert np.std(q1).value == 0.5

    def test_var(self):
        q1 = np.array([1.,2.]) * u.m
        assert np.var(q1) == 0.25 * u.m ** 2
        assert np.var(q1).unit == u.m ** 2
        assert np.var(q1).value == 0.25

    def test_median(self):
        q1 = np.array([1.,2.,4.,5.,6.]) * u.m
        assert np.median(q1) == 4. * u.m
        assert np.median(q1).unit == u.m
        assert np.median(q1).value == 4.

    def test_min(self):
        q1 = np.array([1.,2.,4.,5.,6.]) * u.m
        assert np.min(q1) == 1. * u.m
        assert np.min(q1).unit == u.m
        assert np.min(q1).value == 1.

    def test_max(self):
        q1 = np.array([1.,2.,4.,5.,6.]) * u.m
        assert np.max(q1) == 6. * u.m
        assert np.max(q1).unit == u.m
        assert np.max(q1).value == 6.

    def test_ptp(self):
        q1 = np.array([1.,2.,4.,5.,6.]) * u.m
        assert np.ptp(q1) == 5. * u.m
        assert np.ptp(q1).unit == u.m
        assert np.ptp(q1).value == 5.

    def test_round(self):
        q1 = np.array([1.2, 2.2, 3.2]) * u.kg
        assert np.all(np.round(q1) == np.array([1, 2, 3]) * u.kg)

    def test_sum(self):

        q1 = np.array([1., 2., 6.]) * u.m
        assert np.all(q1.sum() == 9. * u.m)
        assert np.all(np.sum(q1) == 9. * u.m)

        q2 = np.array([[4., 5., 9.], [1.,1.,1.]]) * u.s
        assert np.all(q2.sum(0) == np.array([5., 6., 10.]) * u.s)
        assert np.all(np.sum(q2, 0) == np.array([5., 6., 10.]) * u.s)

    def test_cumsum(self):

        q1 = np.array([1, 2, 6]) * u.m
        assert np.all(q1.cumsum() == np.array([1, 3, 9]) * u.m)
        assert np.all(np.cumsum(q1) == np.array([1, 3, 9]) * u.m)

        q2 = np.array([4, 5, 9]) * u.s
        assert np.all(q2.cumsum() == np.array([4, 9, 18]) * u.s)
        assert np.all(np.cumsum(q2) == np.array([4, 9, 18]) * u.s)

    def test_prod(self):

        q1 = np.array([1, 2, 6]) * u.m
        with pytest.raises(ValueError) as exc:
            q1.prod()
        assert 'cannot use prod' in exc.value.args[0]

        with pytest.raises(ValueError) as exc:
            np.prod(q1)
        assert 'cannot use prod' in exc.value.args[0]

        q2 = np.array([3., 4., 5.]) * u.Unit(1)
        assert q2.prod() == 60. * u.Unit(1)
        assert np.prod(q2) == 60. * u.Unit(1)

    def test_cumprod(self):

        q1 = np.array([1, 2, 6]) * u.m
        with pytest.raises(ValueError) as exc:
            q1.cumprod()
        assert 'cannot use cumprod' in exc.value.args[0]

        with pytest.raises(ValueError) as exc:
            np.cumprod(q1)
        assert 'cannot use cumprod' in exc.value.args[0]

        q2 = np.array([3, 4, 5]) * u.Unit(1)
        assert np.all(q2.cumprod() == np.array([3, 12, 60]) * u.Unit(1))
        assert np.all(np.cumprod(q2) == np.array([3, 12, 60]) * u.Unit(1))


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
            np.tan(np.array([1,2,3]) * u.N)
        assert exc.value.args[0] == ("Can only apply 'tan' function "
                                     "to quantities with angle units")

    def test_arctan_invalid_units(self):
        with pytest.raises(TypeError) as exc:
            np.arctan(np.array([1,2,3]) * u.N)
        assert exc.value.args[0] == ("Can only apply 'arctan' function to "
                                     "dimensionless quantities")

    def test_arctan2_valid(self):
        q1 = np.array([10., 30., 70., 80.]) * u.m
        q2 = 2.0 * u.km
        assert np.arctan2(q1, q2).unit == u.radian
        assert_allclose(np.arctan2(q1, q2).value,
                        np.arctan2(q1.value, q2.to(q1.unit).value))
        q3 = q1/q2
        q4 = 1.
        at2 = np.arctan2(q3, q4)
        assert_allclose(at2, np.arctan2(q3.to(1).value, q4))

    def test_arctan2_invalid(self):
        with pytest.raises(u.UnitsException) as exc:
            np.arctan2(np.array([1,2,3]) * u.N, 1. * u.s)
        assert "compatible dimensions" in exc.value.args[0]
        with pytest.raises(u.UnitsException) as exc:
            np.arctan2(np.array([1,2,3]) * u.N, 1.)
        assert "dimensionless quantities when other arg" in exc.value.args[0]


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
                      np.arange(0,6.,2.) * u.m / u.s)

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

    def test_power_scalar(self):
        assert np.power(4. * u.m, 2.) == 16. * u.m ** 2
        assert np.power(4., 200.*u.cm/u.m) == \
            u.Quantity(16., u.dimensionless_unscaled)

    def test_power_array(self):
        assert np.all(np.power(np.array([1., 2., 3.]) * u.m, 3.)
                      == np.array([1., 8., 27.]) * u.m ** 3)

    def test_power_invalid(self):
        with pytest.raises(TypeError) as exc:
            np.power(3., 4.*u.m)
        assert "raise something to a scalar dimensionless" in exc.value.args[0]
        with pytest.raises(ValueError) as exc:
            np.power(2.*u.m, 4.6)
        assert "must be integers" in exc.value.args[0]

    def test_ldexp_scalar(self):
        assert np.ldexp(4. * u.m, 2) == 16. * u.m

    def test_ldexp_array(self):
        assert np.all(np.ldexp(np.array([1., 2., 3.]) * u.m, [3,2,1])
                      == np.array([8., 8., 6.]) * u.m)

    def test_ldexp_invalid(self):
        with pytest.raises(TypeError) as exc:
            np.ldexp(3.*u.m, 4.)
        assert "not supported" in exc.value.args[0]
        with pytest.raises(TypeError) as exc:
            np.ldexp(3., 4*u.m)
        assert "Cannot use ldexp" in exc.value.args[0]

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

    @pytest.mark.parametrize('function', (np.modf, np.frexp))
    def test_modf_frexp_scalar(self, function):
        q = function(3. * u.m / (6. * u.m))
        assert q == (np.array(0.5), np.array(0.0))

    @pytest.mark.parametrize('function', (np.modf, np.frexp))
    def test_modf_frexp_array(self, function):
        q = function(np.array([2., 3., 6.]) * u.m / (6. * u.m))
        assert all((_q0,_q1) == function(_d) for _q0, _q1, _d
                   in zip(q[0], q[1], [1. / 3., 1. / 2., 1.]))

    @pytest.mark.parametrize('function', (np.modf, np.frexp))
    def test_mod_frexp_invalid_units(self, function):
        # Can't use prod() with non-dimensionless quantities
        with pytest.raises(TypeError) as exc:
            function(3. * u.m / u.s)
        assert exc.value.args[0] == ("Can only apply '{0}' function to "
                                     "unscaled dimensionless quantities"
                                     .format(function.__name__))

        # also does not work on quantities that can be made dimensionless
        with pytest.raises(TypeError) as exc:
            function(np.array([2., 3., 6.]) * u.m / (6. * u.cm))
        assert exc.value.args[0] == ("Can only apply '{0}' function to "
                                     "unscaled dimensionless quantities"
                                     .format(function.__name__))

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
                                         np.negative, np.ones_like, np.rint,
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
    def test_invariant_twoarg_invalid_units(self, ufunc):

        q_i1 = 4.7 * u.m
        q_i2 = 9.4 * u.s
        with pytest.raises(u.UnitsException) as exc:
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
        q_o2 = ufunc(q_i1/q_i2, 2.)
        assert not isinstance(q_o2, u.Quantity)
        assert q_o2.dtype == np.bool
        assert np.all(q_o2 == ufunc((q_i1/q_i2).to(1).value, 2.))

    @pytest.mark.parametrize(('ufunc'), [np.greater, np.greater_equal,
                                         np.less, np.less_equal,
                                         np.not_equal, np.equal])
    def test_comparison_invalid_units(self, ufunc):
        q_i1 = 4.7 * u.m
        q_i2 = 9.4 * u.s
        with pytest.raises(u.UnitsException) as exc:
            ufunc(q_i1, q_i2)
        assert "compatible dimensions" in exc.value.args[0]


class TestInplaceUfuncs(object):
    @pytest.mark.parametrize(('value'), [1., np.arange(10.)])
    def test_one_argument_ufunc_inplace(self, value):
        s = value*u.rad
        check = s
        np.sin(s, out=s)
        assert check is s
        assert check.unit == u.dimensionless_unscaled

    @pytest.mark.parametrize(('value'), [1., np.arange(10.)])
    def test_two_argument_ufunc_inplace(self, value):
        s = value*u.cycle
        check = s
        s /= 2.
        assert check is s
        assert np.all(check.value == value/2.)
        s /= u.s
        assert check is s
        assert check.unit == u.cycle/u.s
        s *= 2. * u.s
        assert check is s
        assert np.all(check == value * u.cycle)
        np.arctan2(s,s,out=s)
        assert check is s
        assert check.unit == u.radian
        with pytest.raises(u.UnitsException):
            s += 1.*u.m
        assert check is s
        assert check.unit == u.radian
