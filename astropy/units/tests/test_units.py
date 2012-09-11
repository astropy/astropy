"""
Regression tests for the units package
"""

from __future__ import absolute_import, unicode_literals, division, print_function
import pytest
from numpy.testing import assert_allclose
from astropy.tests.helper import raises

from ... import units as u


def test_getting_started():
    """
    Corresponds to "Getting Started" section in the docs.
    """
    speed_unit = u.cm / u.s
    x = speed_unit.to(u.mile / u.hour, 1)
    assert_allclose(x, 0.02236936292054402)
    speed_converter = speed_unit.get_converter("mile hour^-1")
    x = speed_converter([1., 1000., 5000.])
    assert_allclose(x, [2.23693629e-02, 2.23693629e+01, 1.11846815e+02])


@raises(ValueError)
def test_invalid_power():
    u.m ** 0.3


def test_invalid_compare():
    assert not (u.m == u.s)


def test_convert():
    assert u.h.get_converter(u.s)(1) == 3600


def test_convert_fail():
    with pytest.raises(u.UnitsException):
        u.cm.to(u.s, 1)
    with pytest.raises(u.UnitsException):
        (u.cm / u.s).to(u.m, 1)


def test_is_equivalent():
    assert u.m.is_equivalent(u.inch)
    assert not (u.Hz.is_equivalent(u.J))
    assert u.Hz.is_equivalent(u.J, u.sp())
    assert u.J.is_equivalent(u.Hz, u.sp())


def test_composite():
    assert (u.cm / u.s * u.h).get_converter(u.m)(1) == 36
    assert u.cm * u.cm == u.cm ** 2

    assert u.cm * u.cm * u.cm == u.cm ** 3

    assert (1 / (u.cm * u.cm)) == u.cm ** -2

    assert u.m / 4.0 == 0.25 * u.m

    assert u.Hz.to(1000 * u.Hz, 1) == 0.001


def test_str():
    assert str(u.cm) == "cm"


def test_repr():
    assert repr(u.cm) == 'Unit("cm")'


def test_spectral():
    a = u.AA.to(u.Hz, 1, u.sp())
    assert_allclose(a, 2.9979245799999995e+18)
    b = u.Hz.to(u.AA, a, u.sp())
    assert_allclose(b, 1)

    a = u.AA.to(u.MHz, 1, u.sp())
    assert_allclose(a, 2.9979245799999995e+12)
    b = u.MHz.to(u.AA, a, u.sp())
    assert_allclose(b, 1)

    a = u.m.to(u.Hz, 1, u.sp())
    assert_allclose(a, 2.9979245799999995e+8)
    b = u.Hz.to(u.m, a, u.sp())
    assert_allclose(b, 1)


def test_spectral2():
    a = u.nm.to(u.J, 500, u.sp())
    assert_allclose(a, 3.972891366538605e-19)
    b = u.J.to(u.nm, a, u.sp())
    assert_allclose(b, 500)

    a = u.AA.to(u.Hz, 1, u.sp())
    b = u.Hz.to(u.J, a, u.sp())
    c = u.AA.to(u.J, 1, u.sp())
    assert_allclose(b, c)


def test_spectral3():
    a = u.nm.to(u.Hz, [1000, 2000], u.sp())
    assert_allclose(a, [2.99792458e+14, 1.49896229e+14])


def test_spectraldensity():
    a = u.AA.to(u.Jy, 1, u.sd(u.eV, 2.2))
    assert_allclose(a, 1059416252057.8357, rtol=1e-4)
    b = u.Jy.to(u.AA, a, u.sd(u.eV, 2.2))
    assert_allclose(b, 1)


def test_spectraldensity2():
    flambda = u.erg / u.angstrom / u.cm ** 2 / u.s
    fnu = u.erg / u.Hz / u.cm ** 2 / u.s

    a = flambda.to(fnu, 1, u.sd(u.AA, 3500))
    assert_allclose(a, 4.086160166177361e-12)


def test_units_conversion():
    assert_allclose(u.kpc.to(u.Mpc), 0.001)
    assert_allclose(u.Mpc.to(u.kpc), 1000)
    assert_allclose(u.yr.to(u.Myr), 1.e-6)
    assert_allclose(u.AU.to(u.pc), 4.84813681e-6)


def test_units_manipulation():
    # Just do some manipulation and check it's happy
    (u.kpc * u.yr) ** (1, 3) / u.Myr
    (u.AA * u.erg) ** 9


def test_decompose():
    assert u.Ry == u.Ry.decompose()


def test_get_equivalent_units():
    assert u.pound in u.get_equivalent_units(u.g)
    assert u.J in u.get_equivalent_units(u.Hz, u.sp())


@raises(ValueError)
def test_unknown_unit():
    u.Unit("foo")


def test_register():
    try:
        u.def_unit("foo", u.m ** 3, register=True)
        assert hasattr(u, 'foo')
    finally:
        if hasattr(u, 'foo'):
            del u.foo


def test_in_units():
    speed_unit = u.cm / u.s
    x = speed_unit.in_units(u.mile / u.hour, 1)

