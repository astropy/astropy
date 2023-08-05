# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests that relate to using quantities/units on parameters of models.
"""
import numpy as np
import pytest

from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.core import Fittable1DModel, InputParameterError
from astropy.modeling.models import (
    BlackBody,
    Const1D,
    Gaussian1D,
    Pix2Sky_TAN,
    RotateNative2Celestial,
    Rotation2D,
)
from astropy.modeling.parameters import Parameter, ParameterDefinitionError
from astropy.tests.helper import assert_quantity_allclose
from astropy.units import UnitsError


class BaseTestModel(Fittable1DModel):
    @staticmethod
    def evaluate(x, a):
        return x


def test_parameter_quantity():
    """
    Basic tests for initializing general models (that do not require units)
    with parameters that have units attached.
    """
    g = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)
    assert g.amplitude.value == 1.0
    assert g.amplitude.unit is u.J
    assert g.mean.value == 1.0
    assert g.mean.unit is u.m
    assert g.stddev.value == 0.1
    assert g.stddev.unit is u.m


def test_parameter_set_quantity():
    """
    Make sure that parameters that start off as quantities can be set to any
    other quantity, regardless of whether the units of the new quantity are
    compatible with the original ones.

    We basically leave it up to the evaluate method to raise errors if there
    are issues with incompatible units, and we don't check for consistency
    at the parameter level.
    """

    g = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)

    # Try equivalent units

    g.amplitude = 4 * u.kJ
    assert_quantity_allclose(g.amplitude, 4 * u.kJ)

    g.mean = 3 * u.km
    assert_quantity_allclose(g.mean, 3 * u.km)

    g.stddev = 2 * u.mm
    assert_quantity_allclose(g.stddev, 2 * u.mm)

    # Try different units

    g.amplitude = 2 * u.s
    assert_quantity_allclose(g.amplitude, 2 * u.s)

    g.mean = 2 * u.Jy
    assert_quantity_allclose(g.mean, 2 * u.Jy)


def test_parameter_lose_units():
    """
    Check that parameters that have been set to a quantity that are then set to
    a value with no units raise an exception. We do this because setting a
    parameter to a value with no units is ambiguous if units were set before:
    if a parameter is 1 * u.Jy and the parameter is then set to 4, does this mean
    2 without units, or 2 * u.Jy?
    """

    g = Gaussian1D(1 * u.Jy, 3, 0.1)

    MESSAGE = (
        r"The .* parameter should be given as a .* because it was originally"
        r" initialized as a .*"
    )
    with pytest.raises(UnitsError, match=MESSAGE):
        g.amplitude = 2


def test_parameter_add_units():
    """
    On the other hand, if starting from a parameter with no units, we should be
    able to add units since this is unambiguous.
    """

    g = Gaussian1D(1, 3, 0.1)

    g.amplitude = 2 * u.Jy
    assert_quantity_allclose(g.amplitude, 2 * u.Jy)


def test_parameter_change_unit():
    """
    Test that changing the unit on a parameter does not work. This is an
    ambiguous operation because it's not clear if it means that the value should
    be converted or if the unit should be changed without conversion.
    """

    g = Gaussian1D(1, 1 * u.m, 0.1 * u.m)

    # Setting a unit on a unitless parameter should not work
    MESSAGE = (
        r"Cannot attach units to parameters that were not initially specified with"
        r" units"
    )
    with pytest.raises(ValueError, match=MESSAGE):
        g.amplitude.unit = u.Jy

    # But changing to another unit should not, even if it is an equivalent unit
    MESSAGE = (
        r"Cannot change the unit attribute directly, instead change the parameter to a"
        r" new quantity"
    )
    with pytest.raises(ValueError, match=MESSAGE):
        g.mean.unit = u.cm


def test_parameter_set_value():
    """
    Test that changing the value on a parameter works as expected.
    """

    g = Gaussian1D(1 * u.Jy, 1 * u.m, 0.1 * u.m)

    # To set a parameter to a quantity, we simply do
    g.amplitude = 2 * u.Jy

    # If we try setting the value, we need to pass a non-quantity value
    # TODO: determine whether this is the desired behavior?
    g.amplitude.value = 4
    assert_quantity_allclose(g.amplitude, 4 * u.Jy)
    assert g.amplitude.value == 4
    assert g.amplitude.unit is u.Jy

    # If we try setting it to a Quantity, we raise an error
    MESSAGE = (
        r"The .value property on parameters should be set to unitless values, not"
        r" Quantity objects.*"
    )
    with pytest.raises(TypeError, match=MESSAGE):
        g.amplitude.value = 3 * u.Jy


def test_parameter_quantity_property():
    """
    Test that the quantity property of Parameters behaves as expected
    """

    # Since parameters have a .value and .unit parameter that return just the
    # value and unit respectively, we also have a .quantity parameter that
    # returns a Quantity instance.

    g = Gaussian1D(1 * u.Jy, 1 * u.m, 0.1 * u.m)

    assert_quantity_allclose(g.amplitude.quantity, 1 * u.Jy)

    # Setting a parameter to a quantity changes the value and the default unit
    g.amplitude.quantity = 5 * u.mJy
    assert g.amplitude.value == 5
    assert g.amplitude.unit is u.mJy

    # And we can also set the parameter to a value with different units
    g.amplitude.quantity = 4 * u.s
    assert g.amplitude.value == 4
    assert g.amplitude.unit is u.s

    # But not to a value without units
    MESSAGE = r"The .quantity attribute should be set to a Quantity object"
    with pytest.raises(TypeError, match=MESSAGE):
        g.amplitude.quantity = 3


def test_parameter_default_units_match():
    # If the unit and default quantity units are different, raise an error
    MESSAGE = (
        r"parameter default 1.0 m does not have units equivalent to the required"
        r" unit Jy"
    )
    with pytest.raises(ParameterDefinitionError, match=MESSAGE):

        class TestC(Fittable1DModel):
            a = Parameter(default=1.0 * u.m, unit=u.Jy)


@pytest.mark.parametrize(("unit", "default"), ((u.m, 1.0), (None, 1 * u.m)))
def test_parameter_defaults(unit, default):
    """
    Test that default quantities are correctly taken into account
    """

    class TestModel(BaseTestModel):
        a = Parameter(default=default, unit=unit)

    # TODO: decide whether the default property should return a value or
    # a quantity?

    # The default unit and value should be set on the class
    assert TestModel.a.unit == u.m
    assert TestModel.a.default == 1.0

    # Check that the default unit and value are also set on a class instance
    m = TestModel()
    assert m.a.unit == u.m
    assert m.a.default == m.a.value == 1.0

    # If the parameter is set to a different value, the default is still the
    # internal default
    m = TestModel(2.0 * u.m)
    assert m.a.unit == u.m
    assert m.a.value == 2.0
    assert m.a.default == 1.0

    # Instantiate with a different, but compatible unit
    m = TestModel(2.0 * u.pc)
    assert m.a.unit == u.pc
    assert m.a.value == 2.0
    # The default is still in the original units
    # TODO: but how do we know what those units are if we don't return a
    # quantity?
    assert m.a.default == 1.0

    # Initialize with a completely different unit
    m = TestModel(2.0 * u.Jy)
    assert m.a.unit == u.Jy
    assert m.a.value == 2.0
    # TODO: this illustrates why the default doesn't make sense anymore
    assert m.a.default == 1.0

    # Instantiating with different units works, and just replaces the original unit
    MESSAGE = r".* requires a Quantity for parameter .*"
    with pytest.raises(InputParameterError, match=MESSAGE):
        TestModel(1.0)


def test_parameter_quantity_arithmetic():
    """
    Test that arithmetic operations with properties that have units return the
    appropriate Quantities.
    """

    g = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)

    # Addition should work if units are compatible
    assert g.mean + (1 * u.m) == 2 * u.m
    assert (1 * u.m) + g.mean == 2 * u.m

    # Multiplication by a scalar should also preserve the quantity-ness
    assert g.mean * 2 == (2 * u.m)
    assert 2 * g.mean == (2 * u.m)

    # Multiplication by a quantity should result in units being multiplied
    assert g.mean * (2 * u.m) == (2 * (u.m**2))
    assert (2 * u.m) * g.mean == (2 * (u.m**2))

    # Negation should work properly too
    assert -g.mean == (-1 * u.m)
    assert abs(-g.mean) == g.mean

    # However, addition of a quantity + scalar should not work
    MESSAGE = (
        r"Can only apply 'add' function to dimensionless quantities when other"
        r" argument .*"
    )
    with pytest.raises(UnitsError, match=MESSAGE):
        g.mean + 1
    with pytest.raises(UnitsError, match=MESSAGE):
        1 + g.mean


def test_parameter_quantity_comparison():
    """
    Basic test of comparison operations on properties with units.
    """

    g = Gaussian1D(1 * u.J, 1 * u.m, 0.1 * u.m)

    # Essentially here we are checking that parameters behave like Quantity

    assert g.mean == 1 * u.m
    assert 1 * u.m == g.mean
    assert g.mean != 1
    assert 1 != g.mean

    assert g.mean < 2 * u.m
    assert 2 * u.m > g.mean

    MESSAGE = (
        r"Can only apply 'less' function to dimensionless quantities when other"
        r" argument .*"
    )
    with pytest.raises(UnitsError, match=MESSAGE):
        g.mean < 2  # noqa: B015

    with pytest.raises(UnitsError, match=MESSAGE):
        2 > g.mean  # noqa: B015

    g = Gaussian1D([1, 2] * u.J, [1, 2] * u.m, [0.1, 0.2] * u.m)

    assert np.all(g.mean == [1, 2] * u.m)
    assert np.all([1, 2] * u.m == g.mean)
    assert np.all(g.mean != [1, 2])
    assert np.all([1, 2] != g.mean)

    with pytest.raises(UnitsError, match=MESSAGE):
        g.mean < [3, 4]  # noqa: B015

    with pytest.raises(UnitsError, match=MESSAGE):
        [3, 4] > g.mean  # noqa: B015


def test_parameters_compound_models():
    Pix2Sky_TAN()
    sky_coords = coord.SkyCoord(ra=5.6, dec=-72, unit=u.deg)
    lon_pole = 180 * u.deg
    n2c = RotateNative2Celestial(sky_coords.ra, sky_coords.dec, lon_pole)
    rot = Rotation2D(23)
    rot | n2c


def test_magunit_parameter():
    """Regression test for bug reproducer in issue #13133"""

    unit = u.ABmag
    c = -20.0 * unit
    model = Const1D(c)

    assert model(-23.0 * unit) == c


def test_log_getter():
    """Regression test for issue #14511"""

    x = 6000 * u.AA
    mdl_base = BlackBody(temperature=5000 * u.K, scale=u.Quantity(1))

    class CustomBlackBody(BlackBody):
        scale = Parameter(
            "scale",
            default=1,
            bounds=(0, None),
            getter=np.log,
            setter=np.exp,
            unit=u.dimensionless_unscaled,
        )

    mdl = CustomBlackBody(temperature=5000 * u.K, scale=u.Quantity(np.log(1)))
    assert mdl.scale == np.log(1)
    assert_quantity_allclose(mdl(x), mdl_base(x))


def test_sqrt_getter():
    """Regression test for issue #14511"""

    x = 1 * u.m
    mdl_base = Gaussian1D(mean=32 * u.m, stddev=3 * u.m)

    class CustomGaussian1D(Gaussian1D):
        mean = Parameter(
            "mean",
            default=1 * u.m,
            bounds=(0, None),
            getter=np.sqrt,
            setter=np.square,
            unit=u.m,
        )
        stddev = Parameter(
            "stddev",
            default=1 * u.m,
            bounds=(0, None),
            getter=np.sqrt,
            setter=np.square,
            unit=u.m,
        )

    mdl = CustomGaussian1D(mean=np.sqrt(32 * u.m), stddev=np.sqrt(3 * u.m))
    assert mdl.mean == np.sqrt(32 * u.m)
    assert (
        mdl.mean._internal_value == np.sqrt(32) ** 2
    )  # numerical inaccuracy results in 32.00000000000001
    assert mdl.mean._internal_unit == u.m
    assert mdl.stddev == np.sqrt(3 * u.m)
    assert (
        mdl.stddev._internal_value == np.sqrt(3) ** 2
    )  # numerical inaccuracy results in 3.0000000000000004
    assert mdl.stddev._internal_unit == u.m
    assert_quantity_allclose(mdl(x), mdl_base(x))
