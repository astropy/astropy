# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Regression tests for the units package
"""

from __future__ import absolute_import, unicode_literals, division, print_function

import warnings

from ...tests.helper import pytest, raises
from ...tests.compat import assert_allclose
from ...utils.compat.fractions import Fraction

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
    assert u.Hz.is_equivalent(u.J, u.spectral())
    assert u.J.is_equivalent(u.Hz, u.spectral())


def test_composite():
    assert (u.cm / u.s * u.h).get_converter(u.m)(1) == 36
    assert u.cm * u.cm == u.cm ** 2

    assert u.cm * u.cm * u.cm == u.cm ** 3

    assert u.Hz.to(1000 * u.Hz, 1) == 0.001


def test_str():
    assert str(u.cm) == "cm"


def test_repr():
    assert repr(u.cm) == 'Unit("cm")'


def test_spectral():
    a = u.AA.to(u.Hz, 1, u.spectral())
    assert_allclose(a, 2.9979245799999995e+18)
    b = u.Hz.to(u.AA, a, u.spectral())
    assert_allclose(b, 1)

    a = u.AA.to(u.MHz, 1, u.spectral())
    assert_allclose(a, 2.9979245799999995e+12)
    b = u.MHz.to(u.AA, a, u.spectral())
    assert_allclose(b, 1)

    a = u.m.to(u.Hz, 1, u.spectral())
    assert_allclose(a, 2.9979245799999995e+8)
    b = u.Hz.to(u.m, a, u.spectral())
    assert_allclose(b, 1)


def test_spectral2():
    a = u.nm.to(u.J, 500, u.spectral())
    assert_allclose(a, 3.972891366538605e-19)
    b = u.J.to(u.nm, a, u.spectral())
    assert_allclose(b, 500)

    a = u.AA.to(u.Hz, 1, u.spectral())
    b = u.Hz.to(u.J, a, u.spectral())
    c = u.AA.to(u.J, 1, u.spectral())
    assert_allclose(b, c)


def test_spectral3():
    a = u.nm.to(u.Hz, [1000, 2000], u.spectral())
    assert_allclose(a, [2.99792458e+14, 1.49896229e+14])


def test_spectraldensity():

    a = u.AA.to(u.Jy, 1, u.spectral_density(u.eV, 2.2))
    assert_allclose(a, 1059416252057.8357, rtol=1e-4)

    b = u.Jy.to(u.AA, a, u.spectral_density(u.eV, 2.2))
    assert_allclose(b, 1)


def test_spectraldensity2():
    flambda = u.erg / u.angstrom / u.cm ** 2 / u.s
    fnu = u.erg / u.Hz / u.cm ** 2 / u.s

    a = flambda.to(fnu, 1, u.spectral_density(u.AA, 3500))
    assert_allclose(a, 4.086160166177361e-12)


def test_spectraldensity3():

    # Define F_nu in Jy
    f_nu = u.Jy

    # Convert to ergs / cm^2 / s / Hz
    assert_allclose(f_nu.to(u.erg / u.cm ** 2 / u.s / u.Hz, 1.), 1.e-23, 10)

    # Convert to ergs / cm^2 / s at 10 Ghz
    assert_allclose(f_nu.to(u.erg / u.cm ** 2 / u.s, 1., equivalencies=u.spectral_density(u.GHz, 10)), 1.e-13, 10)

    # Convert to ergs / cm^2 / s / micron at 1 Ghz
    assert_allclose(f_nu.to(u.erg / u.cm ** 2 / u.s / u.micron, 1., equivalencies=u.spectral_density(u.Hz, 1.e9)), 3.335640951981521e-20, 10)

    # Define F_lambda in ergs / cm^2 / s / micron
    f_lambda = u.erg / u.cm ** 2 / u.s / u.micron

    # Convert to Jy at 1 Ghz
    assert_allclose(f_lambda.to(u.Jy, 1., equivalencies=u.spectral_density(u.Hz, 1.e9)), 1. / 3.335640951981521e-20, 10)

    # Convert to ergs / cm^2 / s at 10 microns
    assert_allclose(f_lambda.to(u.erg / u.cm ** 2 / u.s, 1., equivalencies=u.spectral_density(u.micron, 10.)), 10., 10)


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


def test_equivalent_units():
    assert u.pound in u.g.find_equivalent_units()
    assert u.J in u.Hz.find_equivalent_units(u.spectral())


def test_unknown_unit():
    with warnings.catch_warnings(record=True) as warning_lines:
        warnings.resetwarnings()
        warnings.simplefilter("always", u.UnitsWarning, append=True)
        u.Unit("FOO", parse_strict='warn')

    print(dir(warning_lines[0]))
    assert warning_lines[0].category == u.UnitsWarning
    assert 'FOO' in str(warning_lines[0].message)


def test_unknown_unit2():
    with warnings.catch_warnings(record=True) as warning_lines:
        warnings.resetwarnings()
        warnings.simplefilter("always", u.UnitsWarning, append=True)
        assert u.Unit("m/s/kg", parse_strict='warn').to_string() == 'm/s/kg'

    print(dir(warning_lines[0]))
    assert warning_lines[0].category == u.UnitsWarning
    assert 'm/s/kg' in str(warning_lines[0].message)


def test_unknown_unit3():
    unit = u.Unit("FOO", parse_strict='silent')
    assert isinstance(unit, u.UnrecognizedUnit)
    assert unit.name == "FOO"


def test_register():
    try:
        u.def_unit("foo", u.m ** 3, register=True)
        assert hasattr(u, 'foo')
    finally:
        if hasattr(u, 'foo'):
            del u.UnitBase._registry[u.UnitBase._registry.index(u.foo)]
            del u.foo


def test_in_units():
    speed_unit = u.cm / u.s
    x = speed_unit.in_units(u.mile / u.hour, 1)


def test_null_unit():
    assert (u.m / u.m) == u.Unit(1)


def test_unrecognized_equivalency():
    assert u.m.is_equivalent('foo') == False
    assert u.m.is_equivalent('foot') == True


@raises(TypeError)
def test_unit_noarg():
    u.Unit()


def test_convertible_exception():
    try:
        u.AA.to(u.h * u.s ** 2)
    except u.UnitsException as e:
        assert "length" in str(e)


def test_convertible_exception2():
    try:
        u.m.to(u.s)
    except u.UnitsException as e:
        assert 'time' not in str(e)


@raises(TypeError)
def test_invalid_type():
    class A(object):
        pass

    u.Unit(A())


def test_steradian():
    """
    Issue #599
    """
    assert u.sr.is_equivalent(u.rad * u.rad)
