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

        q2 = (np.array([1.,2.5,4.]) * u.km / u.s) / (4000. * u.m / u.s) + 1.
        assert np.mean(q2) == 0.001625 * u.km / u.m
        assert np.mean(q2).unit == u.km / u.m
        assert np.mean(q2).value == 0.001625

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
        assert exc.value.args[0] == ("Can only apply inverse trigonometric "
                                     "functions to dimensionless and unscaled "
                                     "quantities")

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
        assert exc.value.args[0] == ("Can only apply inverse trigonometric "
                                     "functions to dimensionless and unscaled "
                                     "quantities")

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
            np.sin(np.array([1,2,3]) * u.N)
        assert exc.value.args[0] == ("Can only apply trigonometric functions "
                                     "to quantities with angle units")

    def test_arctan_invalid_units(self):
        with pytest.raises(TypeError) as exc:
            np.arctan(np.array([1,2,3]) * u.N)
        assert exc.value.args[0] == ("Can only apply inverse trigonometric "
                                     "functions to dimensionless and unscaled "
                                     "quantities")


class TestQuantityMathFuncs(object):
    """
    Test other mathematical functions
    """

    def test_sqrt_scalar(self):
        assert np.sqrt(4. * u.m) == 2. * u.m ** 0.5

    def test_sqrt_array(self):
        assert np.all(np.sqrt(np.array([1., 4., 9.]) * u.m)
                      == np.array([1., 2., 3.]) * u.m ** 0.5)

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

    @pytest.mark.parametrize('function', (np.exp, np.log, np.log2, np.log10, np.log1p))
    def test_exp_invalid_units(self, function):

        # Can't use exp() with non-dimensionless quantities
        with pytest.raises(TypeError) as exc:
            function(3. * u.m / u.s)
        assert exc.value.args[0] == ("Can only apply {0} function to dimensionless "
                                     "quantities".format(function.__name__))

        # Also can't use exp() with dimensionless scaled quantities
        with pytest.raises(TypeError) as exc:
            function(3. * u.m / u.km)
        assert exc.value.args[0] == ("Can only apply {0} function to dimensionless "
                                     "quantities".format(function.__name__))
