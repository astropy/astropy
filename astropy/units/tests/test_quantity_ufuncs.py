# The purpose of these tests are to ensure that calling ufuncs with quantities
# returns quantities with the right units, or raises exceptions.

import numpy as np

from ... import units as u
from ...tests.helper import pytest
from ...tests.compat import assert_allclose


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

    def test_cumsum(self):

        q1 = np.array([1, 2, 6]) * u.m
        assert np.all(q1.cumsum() == np.array([1, 3, 9]) * u.m)
        assert np.all(np.cumsum(q1) == np.array([1, 3, 9]) * u.m)

        q2 = np.array([4, 5, 9]) * u.s
        assert np.all(q2.cumsum() == np.array([4, 9, 18]) * u.s)
        assert np.all(np.cumsum(q2) == np.array([4, 9, 18]) * u.s)

    def test_cumprod(self):

        q1 = np.array([1, 2, 6]) * u.m
        with pytest.raises(ValueError) as exc:
            q1.cumprod()
        assert exc.value.args[0] == 'cannot use cumprod on non-dimensionless Quantity arrays'
        with pytest.raises(ValueError) as exc:
            np.cumprod(q1)
        assert exc.value.args[0] == 'cannot use cumprod on non-dimensionless Quantity arrays'

        q2 = np.array([3, 4, 5]) * u.Unit(1)
        print q2.cumprod()
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
        assert exc.value.args[0] == ("Can only apply trigonometric functions "
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
        assert exc.value.args[0] == ("Can only apply trigonometric functions "
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
        assert exc.value.args[0] == ("Can only apply trigonometric functions "
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

    def test_arctan2_invalid_units(self):
        with pytest.raises(TypeError) as exc:
            np.arctan2(np.array([1,2,3]) * u.N, 1. * u.s)
        assert "compatible dimensions" in exc.value.args[0]
        with pytest.raises(TypeError) as exc:
            np.arctan2(np.array([1,2,3]) * u.N, 1.)
        assert "dimensionless quantities if other arg" in exc.value.args[0]


class TestQuantityMathFuncs(object):
    """
    Test other mathematical functions
    """

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
        assert "raise something to a dimensionless" in exc.value.args[0]
        with pytest.raises(ValueError) as exc:
            np.power(2.*u.m, 4.6)
        assert "must be integers" in exc.value.args[0]
        with pytest.raises(TypeError) as exc:
            np.power(2.*u.m, 400.*u.cm/u.m)
        assert "dimensionless scalar quantities" in exc.value.args[0]

    @pytest.mark.parametrize('function', (np.exp, np.log, np.log2, np.log10, np.log1p))
    def test_exp_scalar(self, function):
        q = function(3. * u.m / (6. * u.m))
        assert q.unit == u.dimensionless_unscaled
        assert q.value == function(0.5)

    @pytest.mark.parametrize('function', (np.exp, np.log, np.log2, np.log10, np.log1p))
    def test_exp_array(self, function):
        q = function(np.array([2., 3., 6.]) * u.m / (6. * u.m))
        assert q.unit == u.dimensionless_unscaled
        assert np.all(q.value
                      == function(np.array([1. / 3., 1. / 2., 1.])))

    # should also work on quantities that can be made dimensionless
    @pytest.mark.parametrize('function', (np.exp, np.log, np.log2, np.log10, np.log1p))
    def test_exp_array(self, function):
        q = function(np.array([2., 3., 6.]) * u.m / (6. * u.cm))
        assert q.unit == u.dimensionless_unscaled
        assert_allclose(q.value,
                        function(np.array([100. / 3., 100. / 2., 100.])))

    @pytest.mark.parametrize('function', (np.exp, np.log, np.log2, np.log10, np.log1p))
    def test_exp_invalid_units(self, function):

        # Can't use exp() with non-dimensionless quantities
        with pytest.raises(TypeError) as exc:
            function(3. * u.m / u.s)
        assert exc.value.args[0] == ("Can only apply '{0}' function to "
                                     "dimensionless quantities"
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

    @pytest.mark.parametrize(('ufunc'), [np.absolute, np.conjugate,
                                         np.negative, np.ones_like, np.rint,
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
                                         np.maximum, np.minimum])
    def test_invariant_twoarg_scalar(self, ufunc):

        q_i1 = 4.7 * u.m
        q_i2 = 9.4 * u.km
        q_o = ufunc(q_i1, q_i2)
        assert isinstance(q_o, u.Quantity)
        assert q_o.unit == q_i1.unit
        assert np.all(q_o.value == ufunc(q_i1.value, q_i2.to(q_i1.unit).value))

    @pytest.mark.parametrize(('ufunc'), [np.add, np.subtract, np.hypot,
                                         np.maximum, np.minimum])
    def test_invariant_twoarg_array(self, ufunc):

        q_i1 = np.array([-3.3, 2.1, 10.2]) * u.kg / u.s
        q_i2 = np.array([10., -5., 1.e6]) * u.g / u.Ms
        q_o = ufunc(q_i1, q_i2)
        assert isinstance(q_o, u.Quantity)
        assert q_o.unit == q_i1.unit
        assert np.all(q_o.value == ufunc(q_i1.value, q_i2.to(q_i1.unit).value))

    @pytest.mark.parametrize(('ufunc'), [np.add, np.subtract, np.hypot,
                                         np.maximum, np.minimum])
    def test_invariant_twoarg_invalid_units(self, ufunc):

        q_i1 = 4.7 * u.m
        q_i2 = 9.4 * u.s
        with pytest.raises(TypeError) as exc:
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
    def test_invariant_twoarg_invalid_units(self, ufunc):

        q_i1 = 4.7 * u.m
        q_i2 = 9.4 * u.s
        with pytest.raises(TypeError) as exc:
            ufunc(q_i1, q_i2)
        assert "compatible dimensions" in exc.value.args[0]
