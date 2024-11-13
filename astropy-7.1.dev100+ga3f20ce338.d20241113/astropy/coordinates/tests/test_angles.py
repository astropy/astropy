# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test initialization and other aspects of Angle and subclasses"""

import pickle
import threading

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

import astropy.units as u
from astropy.coordinates import (
    Angle,
    IllegalHourError,
    IllegalMinuteError,
    IllegalMinuteWarning,
    IllegalSecondError,
    IllegalSecondWarning,
    Latitude,
    Longitude,
)
from astropy.utils.compat.numpycompat import NUMPY_LT_2_0


def test_create_angles():
    """
    Tests creating and accessing Angle objects
    """

    """ The "angle" is a fundamental object. The internal
    representation is stored in radians, but this is transparent to the user.
    Units *must* be specified rather than a default value be assumed. This is
    as much for self-documenting code as anything else.

    Angle objects simply represent a single angular coordinate. More specific
    angular coordinates (e.g. Longitude, Latitude) are subclasses of Angle."""

    a1 = Angle(54.12412, unit=u.degree)
    a2 = Angle("54.12412", unit=u.degree)
    a3 = Angle("54:07:26.832", unit=u.degree)
    a4 = Angle("54.12412 deg")
    a5 = Angle("54.12412 degrees")
    a6 = Angle("54.12412°")  # because we like Unicode
    a8 = Angle("54°07'26.832\"")
    a9 = Angle([54, 7, 26.832], unit=u.degree)
    assert_allclose(a9.value, [54, 7, 26.832])
    assert a9.unit is u.degree

    a10 = Angle(3.60827466667, unit=u.hour)
    a11 = Angle("3:36:29.7888000120", unit=u.hour)

    Angle(0.944644098745, unit=u.radian)

    with pytest.raises(u.UnitsError):
        Angle(54.12412)
        # raises an exception because this is ambiguous

    with pytest.raises(u.UnitsError):
        Angle(54.12412, unit=u.m)

    with pytest.raises(ValueError):
        Angle(12.34, unit="not a unit")

    a14 = Angle("03h36m29.7888000120")  # no trailing 's', but unambiguous

    a15 = Angle("5h4m3s")  # single digits, no decimal
    assert a15.unit == u.hourangle

    a16 = Angle("1 d")
    a17 = Angle("1 degree")

    assert a16.degree == 1
    assert a17.degree == 1

    a18 = Angle("54 07.4472", unit=u.degree)
    a19 = Angle("54:07.4472", unit=u.degree)
    a20 = Angle("54d07.4472m", unit=u.degree)
    a21 = Angle("3h36m", unit=u.hour)
    a22 = Angle("3.6h", unit=u.hour)
    a23 = Angle("- 3h", unit=u.hour)
    a24 = Angle("+ 3h", unit=u.hour)
    a25 = Angle(3.0, unit=u.hour**1)

    # ensure the above angles that should match do
    assert a1 == a2 == a3 == a4 == a5 == a6 == a8 == a18 == a19 == a20
    assert_allclose(a1.radian, a2.radian)
    assert_allclose(a2.degree, a3.degree)
    assert_allclose(a3.radian, a4.radian)
    assert_allclose(a4.radian, a5.radian)
    assert_allclose(a5.radian, a6.radian)

    assert_allclose(a10.degree, a11.degree)
    assert a11 == a14
    assert a21 == a22
    assert a23 == -a24
    assert a24 == a25

    # check for illegal ranges / values
    with pytest.raises(IllegalSecondError):
        Angle("12 32 99", unit=u.degree)

    with pytest.raises(IllegalMinuteError):
        Angle("12 99 23", unit=u.degree)

    with pytest.raises(IllegalSecondError):
        Angle("12 32 99", unit=u.hour)

    with pytest.raises(IllegalMinuteError):
        Angle("12 99 23", unit=u.hour)

    with pytest.raises(IllegalHourError):
        Angle("99 25 51.0", unit=u.hour)

    with pytest.raises(ValueError):
        Angle("12 25 51.0xxx", unit=u.hour)

    with pytest.raises(ValueError):
        Angle("12h34321m32.2s")

    assert a1 is not None


def test_angle_from_view():
    q = np.arange(3.0) * u.deg
    a = q.view(Angle)
    assert type(a) is Angle
    assert a.unit is q.unit
    assert np.all(a == q)

    q2 = np.arange(4) * u.m
    with pytest.raises(u.UnitTypeError):
        q2.view(Angle)


def test_angle_ops():
    """
    Tests operations on Angle objects
    """

    # Angles can be added and subtracted. Multiplication and division by a
    # scalar is also permitted. A negative operator is also valid.  All of
    # these operate in a single dimension. Attempting to multiply or divide two
    # Angle objects will return a quantity.  An exception will be raised if it
    # is attempted to store output with a non-angular unit in an Angle [#2718].

    a1 = Angle(3.60827466667, unit=u.hour)
    a2 = Angle("54:07:26.832", unit=u.degree)
    a1 + a2  # creates new Angle object
    a1 - a2
    -a1

    assert_allclose((a1 * 2).hour, 2 * 3.6082746666700003)
    assert abs((a1 / 3.123456).hour - 3.60827466667 / 3.123456) < 1e-10

    # commutativity
    assert (2 * a1).hour == (a1 * 2).hour

    a3 = Angle(a1)  # makes a *copy* of the object, but identical content as a1
    assert_allclose(a1.radian, a3.radian)
    assert a1 is not a3

    a4 = abs(-a1)
    assert a4.radian == a1.radian

    a5 = Angle(5.0, unit=u.hour)
    assert a5 > a1
    assert a5 >= a1
    assert a1 < a5
    assert a1 <= a5

    # check operations with non-angular result give Quantity.
    a6 = Angle(45.0, u.degree)
    a7 = a6 * a5
    assert type(a7) is u.Quantity

    # but those with angular result yield Angle.
    # (a9 is regression test for #5327)
    a8 = a1 + 1.0 * u.deg
    assert type(a8) is Angle
    a9 = 1.0 * u.deg + a1
    assert type(a9) is Angle

    with pytest.raises(TypeError):
        a6 *= a5

    with pytest.raises(TypeError):
        a6 *= u.m

    with pytest.raises(TypeError):
        np.sin(a6, out=a6)


def test_angle_methods():
    # Most methods tested as part of the Quantity tests.
    # A few tests here which caused problems before: #8368
    a = Angle([0.0, 2.0], "deg")
    a_mean = a.mean()
    assert type(a_mean) is Angle
    assert a_mean == 1.0 * u.degree
    a_std = a.std()
    assert type(a_std) is Angle
    assert a_std == 1.0 * u.degree
    a_var = a.var()
    assert type(a_var) is u.Quantity
    assert a_var == 1.0 * u.degree**2
    if NUMPY_LT_2_0:
        # np.ndarray.ptp() method removed in numpy 2.0.
        a_ptp = a.ptp()
        assert type(a_ptp) is Angle
        assert a_ptp == 2.0 * u.degree
    a_max = a.max()
    assert type(a_max) is Angle
    assert a_max == 2.0 * u.degree
    a_min = a.min()
    assert type(a_min) is Angle
    assert a_min == 0.0 * u.degree


def test_angle_nan_functions():
    # Most numpy functions tested as part of the Quantity tests.
    # But check that we drop to Quantity when appropriate; see
    # https://github.com/astropy/astropy/pull/17221#discussion_r1813060768
    a = Angle([0.0, 2.0, np.nan], "deg")
    a_mean = np.nanmean(a)
    assert type(a_mean) is Angle
    assert a_mean == 1.0 * u.degree
    a_std = np.nanstd(a)
    assert type(a_std) is Angle
    assert a_std == 1.0 * u.degree
    a_var = np.nanvar(a)
    assert type(a_var) is u.Quantity
    assert a_var == 1.0 * u.degree**2
    a_max = np.nanmax(a)
    assert type(a_max) is Angle
    assert a_max == 2.0 * u.degree
    a_min = np.nanmin(a)
    assert type(a_min) is Angle
    assert a_min == 0.0 * u.degree


def test_angle_convert():
    """
    Test unit conversion of Angle objects
    """
    angle = Angle("54.12412", unit=u.degree)

    assert_allclose(angle.hour, 3.60827466667)
    assert_allclose(angle.radian, 0.944644098745)
    assert_allclose(angle.degree, 54.12412)

    assert len(angle.hms) == 3
    assert isinstance(angle.hms, tuple)
    assert angle.hms[0] == 3
    assert angle.hms[1] == 36
    assert_allclose(angle.hms[2], 29.78879999999947)
    # also check that the namedtuple attribute-style access works:
    assert angle.hms.h == 3
    assert angle.hms.m == 36
    assert_allclose(angle.hms.s, 29.78879999999947)

    assert len(angle.dms) == 3
    assert isinstance(angle.dms, tuple)
    assert angle.dms[0] == 54
    assert angle.dms[1] == 7
    assert_allclose(angle.dms[2], 26.831999999992036)
    # also check that the namedtuple attribute-style access works:
    assert angle.dms.d == 54
    assert angle.dms.m == 7
    assert_allclose(angle.dms.s, 26.831999999992036)

    assert isinstance(angle.dms[0], float)
    assert isinstance(angle.hms[0], float)

    # now make sure dms and signed_dms work right for negative angles
    negangle = Angle("-54.12412", unit=u.degree)

    assert negangle.dms.d == -54
    assert negangle.dms.m == -7
    assert_allclose(negangle.dms.s, -26.831999999992036)
    assert negangle.signed_dms.sign == -1
    assert negangle.signed_dms.d == 54
    assert negangle.signed_dms.m == 7
    assert_allclose(negangle.signed_dms.s, 26.831999999992036)


def test_angle_formatting():
    """
    Tests string formatting for Angle objects
    """

    """
    The string method of Angle has this signature:
    def string(self, unit=DEGREE, decimal=False, sep=" ", precision=5,
               pad=False):

    The "decimal" parameter defaults to False since if you need to print the
    Angle as a decimal, there's no need to use the "format" method (see
    above).
    """

    angle = Angle("54.12412", unit=u.degree)

    # __str__ is the default `format`
    assert str(angle) == angle.to_string()

    res = "Angle as HMS: 3h36m29.7888s"
    assert f"Angle as HMS: {angle.to_string(unit=u.hour)}" == res

    res = "Angle as HMS: 3:36:29.7888"
    assert f"Angle as HMS: {angle.to_string(unit=u.hour, sep=':')}" == res

    res = "Angle as HMS: 3:36:29.79"
    assert f"Angle as HMS: {angle.to_string(unit=u.hour, sep=':', precision=2)}" == res

    # Note that you can provide one, two, or three separators passed as a
    # tuple or list

    res = "Angle as HMS: 3h36m29.7888s"
    assert (
        "Angle as HMS:"
        f" {angle.to_string(unit=u.hour, sep=('h', 'm', 's'), precision=4)}" == res
    )

    res = "Angle as HMS: 3-36|29.7888"
    assert (
        f"Angle as HMS: {angle.to_string(unit=u.hour, sep=['-', '|'], precision=4)}"
        == res
    )

    res = "Angle as HMS: 3-36-29.7888"
    assert f"Angle as HMS: {angle.to_string(unit=u.hour, sep='-', precision=4)}" == res

    res = "Angle as HMS: 03h36m29.7888s"
    assert f"Angle as HMS: {angle.to_string(unit=u.hour, precision=4, pad=True)}" == res

    # Same as above, in degrees

    angle = Angle("3 36 29.78880", unit=u.degree)

    res = "Angle as DMS: 3d36m29.7888s"
    assert f"Angle as DMS: {angle.to_string(unit=u.degree)}" == res

    res = "Angle as DMS: 3:36:29.7888"
    assert f"Angle as DMS: {angle.to_string(unit=u.degree, sep=':')}" == res

    res = "Angle as DMS: 3:36:29.79"
    assert (
        f"Angle as DMS: {angle.to_string(unit=u.degree, sep=':', precision=2)}" == res
    )

    # Note that you can provide one, two, or three separators passed as a
    # tuple or list

    res = "Angle as DMS: 3d36m29.7888s"
    assert (
        f"Angle as DMS: {angle.to_string(unit=u.deg, sep=('d', 'm', 's'), precision=4)}"
        == res
    )

    res = "Angle as DMS: 3-36|29.7888"
    assert (
        f"Angle as DMS: {angle.to_string(unit=u.degree, sep=['-', '|'], precision=4)}"
        == res
    )

    res = "Angle as DMS: 3-36-29.7888"
    assert (
        f"Angle as DMS: {angle.to_string(unit=u.degree, sep='-', precision=4)}" == res
    )

    res = "Angle as DMS: 03d36m29.7888s"
    assert (
        f"Angle as DMS: {angle.to_string(unit=u.degree, precision=4, pad=True)}" == res
    )

    res = "Angle as rad: 0.0629763 rad"
    assert f"Angle as rad: {angle.to_string(unit=u.radian)}" == res

    res = "Angle as rad decimal: 0.0629763"
    assert (
        f"Angle as rad decimal: {angle.to_string(unit=u.radian, decimal=True)}" == res
    )

    # check negative angles

    angle = Angle(-1.23456789, unit=u.degree)
    angle2 = Angle(-1.23456789, unit=u.hour)

    assert angle.to_string() == "-1d14m04.444404s"
    assert angle.to_string(pad=True) == "-01d14m04.444404s"
    assert angle.to_string(unit=u.hour) == "-0h04m56.2962936s"
    assert angle2.to_string(unit=u.hour, pad=True) == "-01h14m04.444404s"
    assert angle.to_string(unit=u.radian, decimal=True) == "-0.0215473"

    # We should recognize units that are equal but not identical
    assert angle.to_string(unit=u.hour**1) == "-0h04m56.2962936s"


def test_to_string_vector():
    # Regression test for the fact that vectorize doesn't work with Numpy 1.6
    assert (
        Angle([1.0 / 7.0, 1.0 / 7.0], unit="deg").to_string()[0] == "0d08m34.28571429s"
    )
    assert Angle([1.0 / 7.0], unit="deg").to_string()[0] == "0d08m34.28571429s"
    assert Angle(1.0 / 7.0, unit="deg").to_string() == "0d08m34.28571429s"


@pytest.mark.parametrize(
    "unit, sep, expected_string",
    [
        ("deg", "fromunit", "15d00m00s"),
        ("deg", "dms", "15d00m00s"),
        ("deg", "hms", "1h00m00s"),
        ("hourangle", "fromunit", "15h00m00s"),
        ("hourangle", "dms", "225d00m00s"),
        ("hourangle", "hms", "15h00m00s"),
    ],
)
def test_angle_to_string_seps(unit, sep, expected_string):
    # see https://github.com/astropy/astropy/issues/11280
    a = Angle(15, unit)
    assert a.to_string(sep=sep) == expected_string


def test_angle_format_roundtripping():
    """
    Ensures that the string representation of an angle can be used to create a
    new valid Angle.
    """

    a1 = Angle(0, unit=u.radian)
    a2 = Angle(10, unit=u.degree)
    a3 = Angle(0.543, unit=u.degree)
    a4 = Angle("1d2m3.4s")

    assert Angle(str(a1)).degree == a1.degree
    assert Angle(str(a2)).degree == a2.degree
    assert Angle(str(a3)).degree == a3.degree
    assert Angle(str(a4)).degree == a4.degree

    # also check Longitude/Latitude
    ra = Longitude("1h2m3.4s")
    dec = Latitude("1d2m3.4s")

    assert_allclose(Angle(str(ra)).degree, ra.degree)
    assert_allclose(Angle(str(dec)).degree, dec.degree)


def test_radec():
    """
    Tests creation/operations of Longitude and Latitude objects
    """

    """
    Longitude and Latitude are objects that are subclassed from Angle. As with Angle, Longitude
    and Latitude can parse any unambiguous format (tuples, formatted strings, etc.).

    The intention is not to create an Angle subclass for every possible
    coordinate object (e.g. galactic l, galactic b). However, equatorial Longitude/Latitude
    are so prevalent in astronomy that it's worth creating ones for these
    units. They will be noted as "special" in the docs and use of the just the
    Angle class is to be used for other coordinate systems.
    """

    with pytest.raises(u.UnitsError):
        Longitude("4:08:15.162342")  # error - hours or degrees?
    with pytest.raises(u.UnitsError):
        Longitude("-4:08:15.162342")

    # the "smart" initializer allows >24 to automatically do degrees, but the
    # Angle-based one does not
    # TODO: adjust in 0.3 for whatever behavior is decided on

    # ra = Longitude("26:34:15.345634")  # unambiguous b/c hours don't go past 24
    # assert_allclose(ra.degree, 26.570929342)
    with pytest.raises(u.UnitsError):
        Longitude("26:34:15.345634")

    # ra = Longitude(68)
    with pytest.raises(u.UnitsError):
        Longitude(68)

    with pytest.raises(u.UnitsError):
        Longitude(12)

    with pytest.raises(ValueError):
        Longitude("garbage containing a d and no units")

    ra = Longitude("12h43m23s")
    assert_allclose(ra.hour, 12.7230555556)

    # Units can be specified
    ra = Longitude("4:08:15.162342", unit=u.hour)

    # TODO: this was the "smart" initializer behavior - adjust in 0.3 appropriately
    # Where Longitude values are commonly found in hours or degrees, declination is
    # nearly always specified in degrees, so this is the default.
    # dec = Latitude("-41:08:15.162342")
    with pytest.raises(u.UnitsError):
        Latitude("-41:08:15.162342")
    dec = Latitude("-41:08:15.162342", unit=u.degree)  # same as above


def test_negative_zero_dms():
    # Test for DMS parser
    a = Angle("-00:00:10", u.deg)
    assert_allclose(a.degree, -10.0 / 3600.0)

    # Unicode minus
    a = Angle("−00:00:10", u.deg)
    assert_allclose(a.degree, -10.0 / 3600.0)


def test_negative_zero_dm():
    # Test for DM parser
    a = Angle("-00:10", u.deg)
    assert_allclose(a.degree, -10.0 / 60.0)


def test_negative_zero_hms():
    # Test for HMS parser
    a = Angle("-00:00:10", u.hour)
    assert_allclose(a.hour, -10.0 / 3600.0)


def test_negative_zero_hm():
    # Test for HM parser
    a = Angle("-00:10", u.hour)
    assert_allclose(a.hour, -10.0 / 60.0)


def test_negative_sixty_hm():
    # Test for HM parser
    with pytest.warns(IllegalMinuteWarning):
        a = Angle("-00:60", u.hour)
    assert_allclose(a.hour, -1.0)


def test_plus_sixty_hm():
    # Test for HM parser
    with pytest.warns(IllegalMinuteWarning):
        a = Angle("00:60", u.hour)
    assert_allclose(a.hour, 1.0)


def test_negative_fifty_nine_sixty_dms():
    # Test for DMS parser
    with pytest.warns(IllegalSecondWarning):
        a = Angle("-00:59:60", u.deg)
    assert_allclose(a.degree, -1.0)


def test_plus_fifty_nine_sixty_dms():
    # Test for DMS parser
    with pytest.warns(IllegalSecondWarning):
        a = Angle("+00:59:60", u.deg)
    assert_allclose(a.degree, 1.0)


def test_negative_sixty_dms():
    # Test for DMS parser
    with pytest.warns(IllegalSecondWarning):
        a = Angle("-00:00:60", u.deg)
    assert_allclose(a.degree, -1.0 / 60.0)


def test_plus_sixty_dms():
    # Test for DMS parser
    with pytest.warns(IllegalSecondWarning):
        a = Angle("+00:00:60", u.deg)
    assert_allclose(a.degree, 1.0 / 60.0)


def test_angle_to_is_angle():
    with pytest.warns(IllegalSecondWarning):
        a = Angle("00:00:60", u.deg)
    assert isinstance(a, Angle)
    assert isinstance(a.to(u.rad), Angle)


def test_angle_to_quantity():
    with pytest.warns(IllegalSecondWarning):
        a = Angle("00:00:60", u.deg)
    q = u.Quantity(a)
    assert isinstance(q, u.Quantity)
    assert q.unit is u.deg


def test_quantity_to_angle():
    a = Angle(1.0 * u.deg)
    assert isinstance(a, Angle)
    with pytest.raises(u.UnitsError):
        Angle(1.0 * u.meter)
    a = Angle(1.0 * u.hour)
    assert isinstance(a, Angle)
    assert a.unit is u.hourangle
    with pytest.raises(u.UnitsError):
        Angle(1.0 * u.min)


def test_angle_string():
    with pytest.warns(IllegalSecondWarning):
        a = Angle("00:00:60", u.deg)
    assert str(a) == "0d01m00s"
    a = Angle("00:00:59S", u.deg)
    assert str(a) == "-0d00m59s"
    a = Angle("00:00:59N", u.deg)
    assert str(a) == "0d00m59s"
    a = Angle("00:00:59E", u.deg)
    assert str(a) == "0d00m59s"
    a = Angle("00:00:59W", u.deg)
    assert str(a) == "-0d00m59s"
    a = Angle("-00:00:10", u.hour)
    assert str(a) == "-0h00m10s"
    a = Angle("00:00:59E", u.hour)
    assert str(a) == "0h00m59s"
    a = Angle("00:00:59W", u.hour)
    assert str(a) == "-0h00m59s"
    a = Angle(3.2, u.radian)
    assert str(a) == "3.2 rad"
    a = Angle(4.2, u.microarcsecond)
    assert str(a) == "4.2 uarcsec"
    a = Angle("1.0uarcsec")
    assert a.value == 1.0
    assert a.unit == u.microarcsecond
    a = Angle("1.0uarcsecN")
    assert a.value == 1.0
    assert a.unit == u.microarcsecond
    a = Angle("1.0uarcsecS")
    assert a.value == -1.0
    assert a.unit == u.microarcsecond
    a = Angle("1.0uarcsecE")
    assert a.value == 1.0
    assert a.unit == u.microarcsecond
    a = Angle("1.0uarcsecW")
    assert a.value == -1.0
    assert a.unit == u.microarcsecond
    a = Angle("3d")
    assert_allclose(a.value, 3.0)
    assert a.unit == u.degree
    a = Angle("3dN")
    assert str(a) == "3d00m00s"
    assert a.unit == u.degree
    a = Angle("3dS")
    assert str(a) == "-3d00m00s"
    assert a.unit == u.degree
    a = Angle("3dE")
    assert str(a) == "3d00m00s"
    assert a.unit == u.degree
    a = Angle("3dW")
    assert str(a) == "-3d00m00s"
    assert a.unit == u.degree
    a = Angle('10"')
    assert_allclose(a.value, 10.0)
    assert a.unit == u.arcsecond
    a = Angle("10'N")
    assert_allclose(a.value, 10.0)
    assert a.unit == u.arcminute
    a = Angle("10'S")
    assert_allclose(a.value, -10.0)
    assert a.unit == u.arcminute
    a = Angle("10'E")
    assert_allclose(a.value, 10.0)
    assert a.unit == u.arcminute
    a = Angle("10'W")
    assert_allclose(a.value, -10.0)
    assert a.unit == u.arcminute
    a = Angle("45°55′12″N")
    assert str(a) == "45d55m12s"
    assert_allclose(a.value, 45.92)
    assert a.unit == u.deg
    a = Angle("45°55′12″S")
    assert str(a) == "-45d55m12s"
    assert_allclose(a.value, -45.92)
    assert a.unit == u.deg
    a = Angle("45°55′12″E")
    assert str(a) == "45d55m12s"
    assert_allclose(a.value, 45.92)
    assert a.unit == u.deg
    a = Angle("45°55′12″W")
    assert str(a) == "-45d55m12s"
    assert_allclose(a.value, -45.92)
    assert a.unit == u.deg
    with pytest.raises(ValueError):
        Angle("00h00m10sN")
    with pytest.raises(ValueError):
        Angle("45°55′12″NS")


def test_angle_repr():
    assert "Angle" in repr(Angle(0, u.deg))
    assert "Longitude" in repr(Longitude(0, u.deg))
    assert "Latitude" in repr(Latitude(0, u.deg))

    a = Angle(0, u.deg)
    repr(a)


def test_large_angle_representation():
    """Test that angles above 360 degrees can be output as strings,
    in repr, str, and to_string.  (regression test for #1413)"""
    a = Angle(350, u.deg) + Angle(350, u.deg)
    a.to_string()
    a.to_string(u.hourangle)
    repr(a)
    repr(a.to(u.hourangle))
    str(a)
    str(a.to(u.hourangle))


def test_wrap_at_inplace():
    a = Angle([-20, 150, 350, 360] * u.deg)
    out = a.wrap_at("180d", inplace=True)
    assert out is None
    assert np.all(a.degree == np.array([-20.0, 150.0, -10.0, 0.0]))


def test_latitude():
    """Test input validation for setting Latitude angles."""

    lim_exc = r"Latitude angle\(s\) must be within -90 deg <= angle <= 90 deg, got"
    with pytest.raises(ValueError, match=rf"{lim_exc} \[91. 89.\] deg"):
        Latitude([91, 89] * u.deg)
    with pytest.raises(ValueError, match=f"{lim_exc} -91.0 deg"):
        Latitude("-91d")

    lat = Latitude(["90d", "89d"])
    # check that one can get items
    assert lat[0] == 90 * u.deg
    assert lat[1] == 89 * u.deg
    # and that comparison with angles works
    assert np.all(lat == Angle(["90d", "89d"]))
    # check setitem works
    lat[1] = 45.0 * u.deg
    assert np.all(lat == Angle(["90d", "45d"]))
    # but not with values out of range
    with pytest.raises(ValueError, match=f"{lim_exc} 90.001 deg"):
        lat[0] = 90.001 * u.deg
    with pytest.raises(ValueError, match=f"{lim_exc} -90.001 deg"):
        lat[0] = -90.001 * u.deg
    # these should also not destroy input (#1851)
    assert np.all(lat == Angle(["90d", "45d"]))

    # Check error message for long-ish input arrays (#13994).
    with pytest.raises(
        ValueError, match=rf"{lim_exc} -143.239\d* deg <= angle <= 171.887\d* deg"
    ):
        lat = Latitude([0, 1, 2, 3, -2.5, -1, -0.5] * u.radian)


def test_latitude_manipulation():
    """Test value manipulation operations on Latitude angles."""

    lat = Latitude(["90d", "-29d"])
    # conserve type on unit change (closes #1423)
    angle = lat.to("radian")
    assert type(angle) is Latitude
    # but not on calculations
    angle = lat - 190 * u.deg
    assert type(angle) is Angle
    assert_allclose(angle, [-100, -219] * u.deg)

    lat = Latitude("80d")
    angle = lat / 2.0
    assert type(angle) is Angle
    assert angle == 40 * u.deg

    angle = lat * 2.0
    assert type(angle) is Angle
    assert angle == 160 * u.deg

    angle = -lat
    assert type(angle) is Angle
    assert angle == -80 * u.deg


def test_lon_as_lat():
    """Test validation when trying to interoperate with longitudes."""

    lon = Longitude(10, "deg")
    with pytest.raises(
        TypeError, match="A Latitude angle cannot be created from a Longitude angle"
    ):
        Latitude(lon)

    with pytest.raises(
        TypeError, match="A Longitude angle cannot be assigned to a Latitude angle"
    ):
        lat = Latitude([20], "deg")
        lat[0] = lon

    # Check we can work around the Lat vs Long checks by casting explicitly to Angle.
    lat = Latitude(Angle(lon))
    assert lat.value == 10.0

    # Check setitem.
    lat = Latitude([20], "deg")
    lat[0] = Angle(lon)
    assert lat.value[0] == 10.0


@pytest.mark.parametrize("lon", ["12.3dW", "12h13m12sE", ["1d", "1dW"]])
def test_lon_as_lat_str(lon):
    with pytest.raises(TypeError, match="Latitude.*cannot be created from a Longitude"):
        Latitude(lon)


@pytest.mark.parametrize("lat", ["12.3dN", "12d13m12sS", ["1d", "1dS"]])
def test_lat_as_lon_str(lat):
    with pytest.raises(TypeError, match="Longitude.*cannot be created from a Latitude"):
        Longitude(lat)


def test_longitude():
    """Test setting and manipulation operations on Longitude angles."""

    # Default wrapping at 360d with an array input
    lon = Longitude(["370d", "88d"])
    assert np.all(lon == Longitude(["10d", "88d"]))
    assert np.all(lon == Angle(["10d", "88d"]))

    # conserve type on unit change and keep wrap_angle (closes #1423)
    angle = lon.to("hourangle")
    assert type(angle) is Longitude
    assert angle.wrap_angle == lon.wrap_angle
    angle = lon[0]
    assert type(angle) is Longitude
    assert angle.wrap_angle == lon.wrap_angle
    angle = lon[1:]
    assert type(angle) is Longitude
    assert angle.wrap_angle == lon.wrap_angle

    # but not on calculations
    angle = lon / 2.0
    assert np.all(angle == Angle(["5d", "44d"]))
    assert type(angle) is Angle
    assert not hasattr(angle, "wrap_angle")

    angle = lon * 2.0 + 400 * u.deg
    assert np.all(angle == Angle(["420d", "576d"]))
    assert type(angle) is Angle

    # Test setting a mutable value and having it wrap
    lon[1] = -10 * u.deg
    assert np.all(lon == Angle(["10d", "350d"]))

    # Test wrapping and try hitting some edge cases
    lon = Longitude(np.array([0, 0.5, 1.0, 1.5, 2.0]) * np.pi, unit=u.radian)
    assert np.all(lon.degree == np.array([0.0, 90, 180, 270, 0]))

    lon = Longitude(
        np.array([0, 0.5, 1.0, 1.5, 2.0]) * np.pi, unit=u.radian, wrap_angle="180d"
    )
    assert np.all(lon.degree == np.array([0.0, 90, -180, -90, 0]))

    # Wrap on setting wrap_angle property (also test auto-conversion of wrap_angle to an Angle)
    lon = Longitude(np.array([0, 0.5, 1.0, 1.5, 2.0]) * np.pi, unit=u.radian)
    lon.wrap_angle = "180d"
    assert np.all(lon.degree == np.array([0.0, 90, -180, -90, 0]))

    lon = Longitude("460d")
    assert lon == Angle("100d")
    lon.wrap_angle = "90d"
    assert lon == Angle("-260d")

    # check that if we initialize a longitude with another longitude,
    # wrap_angle is kept by default
    lon2 = Longitude(lon)
    assert lon2.wrap_angle == lon.wrap_angle
    # but not if we explicitly set it
    lon3 = Longitude(lon, wrap_angle="180d")
    assert lon3.wrap_angle == 180 * u.deg

    # check that wrap_angle is always an Angle
    lon = Longitude(lon, wrap_angle=Longitude(180 * u.deg))
    assert lon.wrap_angle == 180 * u.deg
    assert lon.wrap_angle.__class__ is Angle

    # check that wrap_angle is not copied
    wrap_angle = 180 * u.deg
    lon = Longitude(lon, wrap_angle=wrap_angle)
    assert lon.wrap_angle == 180 * u.deg
    assert np.may_share_memory(lon.wrap_angle, wrap_angle)

    # check for problem reported in #2037 about Longitude initializing to -0
    lon = Longitude(0, u.deg)
    lonstr = lon.to_string()
    assert not lonstr.startswith("-")

    # also make sure dtype is correctly conserved
    assert Longitude(0, u.deg, dtype=float).dtype == np.dtype(float)
    assert Longitude(0, u.deg, dtype=int).dtype == np.dtype(int)


def test_lat_as_lon():
    """Test validation when trying to interoperate with latitudes."""

    lat = Latitude(10, "deg")
    # Test errors when trying to interoperate with latitudes.
    with pytest.raises(
        TypeError, match="A Longitude angle cannot be created from a Latitude angle"
    ):
        Longitude(lat)

    with pytest.raises(
        TypeError, match="A Latitude angle cannot be assigned to a Longitude angle"
    ):
        lon = Longitude([20], "deg")
        lon[0] = lat

    # Check we can work around the Lat vs Long checks by casting explicitly to Angle.
    lon = Longitude(Angle(lat))
    assert lon.value == 10.0

    # Check setitem.
    lon = Longitude([20], "deg")
    lon[0] = Angle(lat)
    assert lon.value[0] == 10.0


def test_wrap_at():
    a = Angle([-20, 150, 350, 360] * u.deg)
    assert np.all(a.wrap_at(360 * u.deg).degree == np.array([340.0, 150.0, 350.0, 0.0]))
    assert np.all(
        a.wrap_at(Angle(360, unit=u.deg)).degree == np.array([340.0, 150.0, 350.0, 0.0])
    )
    assert np.all(a.wrap_at("360d").degree == np.array([340.0, 150.0, 350.0, 0.0]))
    assert np.all(a.wrap_at("180d").degree == np.array([-20.0, 150.0, -10.0, 0.0]))
    assert np.all(
        a.wrap_at(np.pi * u.rad).degree == np.array([-20.0, 150.0, -10.0, 0.0])
    )

    # Test wrapping a scalar Angle
    a = Angle("190d")
    assert a.wrap_at("180d") == Angle("-170d")

    a = Angle(np.arange(-1000.0, 1000.0, 0.125), unit=u.deg)
    for wrap_angle in (270, 0.2, 0.0, 360.0, 500, -2000.125):
        aw = a.wrap_at(wrap_angle * u.deg)
        assert np.all(aw.degree >= wrap_angle - 360.0)
        assert np.all(aw.degree < wrap_angle)

        aw = a.to(u.rad).wrap_at(wrap_angle * u.deg)
        assert np.all(aw.degree >= wrap_angle - 360.0)
        assert np.all(aw.degree < wrap_angle)


def test_is_within_bounds():
    a = Angle([-20, 150, 350] * u.deg)
    assert a.is_within_bounds("0d", "360d") is False
    assert a.is_within_bounds(None, "360d") is True
    assert a.is_within_bounds(-30 * u.deg, None) is True

    a = Angle("-20d")
    assert a.is_within_bounds("0d", "360d") is False
    assert a.is_within_bounds(None, "360d") is True
    assert a.is_within_bounds(-30 * u.deg, None) is True


def test_angle_mismatched_unit():
    a = Angle("+6h7m8s", unit=u.degree)
    assert_allclose(a.value, 91.78333333333332)


def test_regression_formatting_negative():
    # Regression test for a bug that caused:
    #
    # >>> Angle(-1., unit='deg').to_string()
    # '-1d00m-0s'
    assert Angle(-0.0, unit="deg").to_string() == "-0d00m00s"
    assert Angle(-1.0, unit="deg").to_string() == "-1d00m00s"
    assert Angle(-0.0, unit="hour").to_string() == "-0h00m00s"
    assert Angle(-1.0, unit="hour").to_string() == "-1h00m00s"


def test_regression_formatting_default_precision():
    # Regression test for issue #11140
    assert Angle("10:20:30.12345678d").to_string() == "10d20m30.12345678s"
    assert Angle("10d20m30.123456784564s").to_string() == "10d20m30.12345678s"
    assert Angle("10d20m30.123s").to_string() == "10d20m30.123s"


def test_empty_sep():
    a = Angle("05h04m31.93830s")

    assert a.to_string(sep="", precision=2, pad=True) == "050431.94"


@pytest.mark.parametrize("angle_class", [Angle, Longitude])
@pytest.mark.parametrize("unit", [u.hourangle, u.hour, None])
def test_create_tuple_fail(angle_class, unit):
    """Creating an angle from an (h,m,s) tuple should fail."""
    with pytest.raises(TypeError, match="no longer supported"):
        angle_class((12, 14, 52), unit=unit)


def test_list_of_quantities():
    a1 = Angle([1 * u.deg, 1 * u.hourangle])
    assert a1.unit == u.deg
    assert_allclose(a1.value, [1, 15])

    a2 = Angle([1 * u.hourangle, 1 * u.deg], u.deg)
    assert a2.unit == u.deg
    assert_allclose(a2.value, [15, 1])


def test_multiply_divide():
    # Issue #2273
    a1 = Angle([1, 2, 3], u.deg)
    a2 = Angle([4, 5, 6], u.deg)
    a3 = a1 * a2
    assert_allclose(a3.value, [4, 10, 18])
    assert a3.unit == (u.deg * u.deg)

    a3 = a1 / a2
    assert_allclose(a3.value, [0.25, 0.4, 0.5])
    assert a3.unit == u.dimensionless_unscaled


def test_mixed_string_and_quantity():
    a1 = Angle(["1d", 1.0 * u.deg])
    assert_array_equal(a1.value, [1.0, 1.0])
    assert a1.unit == u.deg

    a2 = Angle(["1d", 1 * u.rad * np.pi, "3d"])
    assert_array_equal(a2.value, [1.0, 180.0, 3.0])
    assert a2.unit == u.deg


def test_array_angle_tostring():
    aobj = Angle([1, 2], u.deg)
    assert aobj.to_string().dtype.kind == "U"
    assert np.all(aobj.to_string() == ["1d00m00s", "2d00m00s"])


def test_wrap_at_without_new():
    """
    Regression test for subtle bugs from situations where an Angle is
    created via numpy channels that don't do the standard __new__ but instead
    depend on array_finalize to set state.  Longitude is used because the
    bug was in its _wrap_angle not getting initialized correctly
    """
    l1 = Longitude([1] * u.deg)
    l2 = Longitude([2] * u.deg)

    l = np.concatenate([l1, l2])
    assert l._wrap_angle is not None


def test__str__():
    """
    Check the __str__ method used in printing the Angle
    """

    # scalar angle
    scangle = Angle("10.2345d")
    strscangle = scangle.__str__()
    assert strscangle == "10d14m04.2s"

    # non-scalar array angles
    arrangle = Angle(["10.2345d", "-20d"])
    strarrangle = arrangle.__str__()

    assert strarrangle == "[10d14m04.2s -20d00m00s]"

    # summarizing for large arrays, ... should appear
    bigarrangle = Angle(np.ones(10000), u.deg)
    assert "..." in bigarrangle.__str__()


def test_repr_latex():
    """
    Check the _repr_latex_ method, used primarily by IPython notebooks
    """

    # try with both scalar
    scangle = Angle(2.1, u.deg)
    rlscangle = scangle._repr_latex_()

    # and array angles
    arrangle = Angle([1, 2.1], u.deg)
    rlarrangle = arrangle._repr_latex_()

    assert rlscangle == r"$2^\circ06{}^\prime00{}^{\prime\prime}$"
    assert rlscangle.split("$")[1] in rlarrangle

    # make sure the ... appears for large arrays
    bigarrangle = Angle(np.ones(50000) / 50000.0, u.deg)
    assert "..." in bigarrangle._repr_latex_()


def test_angle_with_cds_units_enabled():
    """Regression test for #5350

    Especially the example in
    https://github.com/astropy/astropy/issues/5350#issuecomment-248770151
    """
    # the problem is with the parser, so remove it temporarily
    from astropy.coordinates.angles.formats import _AngleParser
    from astropy.units import cds

    del _AngleParser._thread_local._parser
    with cds.enable():
        Angle("5d")
    del _AngleParser._thread_local._parser
    Angle("5d")


def test_longitude_nan():
    # Check that passing a NaN to Longitude doesn't raise a warning
    Longitude([0, np.nan, 1] * u.deg)


def test_latitude_nan():
    # Check that passing a NaN to Latitude doesn't raise a warning
    Latitude([0, np.nan, 1] * u.deg)


def test_angle_wrap_at_nan():
    # Check that no attempt is made to wrap a NaN angle
    angle = Angle([0, np.nan, 1] * u.deg)
    angle.flags.writeable = False  # to force an error if a write is attempted
    angle.wrap_at(180 * u.deg, inplace=True)


def test_angle_multithreading():
    """
    Regression test for issue #7168
    """
    angles = ["00:00:00"] * 10000

    def parse_test(i=0):
        Angle(angles, unit="hour")

    for i in range(10):
        threading.Thread(target=parse_test, args=(i,)).start()


@pytest.mark.parametrize("cls", [Angle, Longitude, Latitude])
@pytest.mark.parametrize(
    "input, expstr, exprepr",
    [
        (np.nan * u.deg, "nan", "nan deg"),
        ([np.nan, 5, 0] * u.deg, "[nan 5d00m00s 0d00m00s]", "[nan, 5., 0.] deg"),
        ([6, np.nan, 0] * u.deg, "[6d00m00s nan 0d00m00s]", "[6., nan, 0.] deg"),
        ([np.nan, np.nan, np.nan] * u.deg, "[nan nan nan]", "[nan, nan, nan] deg"),
        (np.nan * u.hour, "nan", "nan hourangle"),
        ([np.nan, 5, 0] * u.hour, "[nan 5h00m00s 0h00m00s]", "[nan, 5., 0.] hourangle"),
        ([6, np.nan, 0] * u.hour, "[6h00m00s nan 0h00m00s]", "[6., nan, 0.] hourangle"),
        (
            [np.nan, np.nan, np.nan] * u.hour,
            "[nan nan nan]",
            "[nan, nan, nan] hourangle",
        ),
        (np.nan * u.rad, "nan", "nan rad"),
        ([np.nan, 1, 0] * u.rad, "[nan 1 rad 0 rad]", "[nan, 1., 0.] rad"),
        ([1.50, np.nan, 0] * u.rad, "[1.5 rad nan 0 rad]", "[1.5, nan, 0.] rad"),
        ([np.nan, np.nan, np.nan] * u.rad, "[nan nan nan]", "[nan, nan, nan] rad"),
    ],
)
def test_str_repr_angles_nan(cls, input, expstr, exprepr):
    """
    Regression test for issue #11473
    """
    q = cls(input)
    assert str(q) == expstr
    # Deleting whitespaces since repr appears to be adding them for some values
    # making the test fail.
    assert repr(q).replace(" ", "") == f"<{cls.__name__}{exprepr}>".replace(" ", "")


@pytest.mark.parametrize("sign", (-1, 1))
@pytest.mark.parametrize(
    "value,expected_value,dtype,expected_dtype",
    [
        (np.pi * 2, 0.0, None, np.float64),
        (np.pi * 2, 0.0, np.float64, np.float64),
        (np.float32(2 * np.pi), np.float32(0.0), None, np.float32),
        (np.float32(2 * np.pi), np.float32(0.0), np.float32, np.float32),
    ],
)
def test_longitude_wrap(value, expected_value, dtype, expected_dtype, sign):
    """
    Test that the wrapping of the Longitude value range in radians works
    in both float32 and float64.
    """
    # This prevents upcasting to float64 as sign * value would do.
    if sign < 0:
        value = -value
        expected_value = -expected_value

    result = Longitude(value, u.rad, dtype=dtype)
    assert result.value == expected_value
    assert result.dtype == expected_dtype
    assert result.unit == u.rad


@pytest.mark.parametrize("sign", (-1, 1))
@pytest.mark.parametrize(
    "value,expected_value,dtype,expected_dtype",
    [
        (np.pi / 2, np.pi / 2, None, np.float64),
        (np.pi / 2, np.pi / 2, np.float64, np.float64),
        (np.float32(np.pi / 2), np.float32(np.pi / 2), None, np.float32),
        (np.float32(np.pi / 2), np.float32(np.pi / 2), np.float32, np.float32),
        # these cases would require coercing the float32 value to the float64 value
        # making validate have side effects, so it's not implemented for now
        # (np.float32(np.pi / 2), np.pi / 2, np.float64, np.float64),
        # (np.float32(-np.pi / 2), -np.pi / 2, np.float64, np.float64),
    ],
)
def test_latitude_limits(value, expected_value, dtype, expected_dtype, sign):
    """
    Test that the validation of the Latitude value range in radians works
    in both float32 and float64.

    As discussed in issue #13708, before, the float32 representation of pi/2
    was rejected as invalid because the comparison always used the float64
    representation.
    """
    # this prevents upcasting to float64 as sign * value would do
    if sign < 0:
        value = -value
        expected_value = -expected_value

    result = Latitude(value, u.rad, dtype=dtype)
    assert result.value == expected_value
    assert result.dtype == expected_dtype
    assert result.unit == u.rad


@pytest.mark.parametrize(
    "value,dtype",
    [
        (0.50001 * np.pi, np.float32),
        (np.float32(0.50001 * np.pi), np.float32),
        (np.float32(-0.50001 * np.pi), np.float32),
        (0.50001 * np.pi, np.float64),
    ],
)
def test_latitude_out_of_limits(value, dtype):
    """
    Test that values slightly larger than pi/2 are rejected for different dtypes.
    Test cases for issue #13708
    """
    with pytest.raises(
        ValueError,
        match=r"Latitude angle\(s\) must be within -90 deg "
        r"<= angle <= 90 deg, got -?90.001\d* deg",
    ):
        Latitude(value, u.rad, dtype=dtype)


def test_angle_pickle_to_string():
    """
    Ensure that after pickling we can still do to_string on hourangle.

    Regression test for gh-13923.
    """
    angle = Angle(0.25 * u.hourangle)
    expected = angle.to_string()
    via_pickle = pickle.loads(pickle.dumps(angle))
    via_pickle_string = via_pickle.to_string()  # This used to fail.
    assert via_pickle_string == expected
