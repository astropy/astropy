# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Regression tests for the units package
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import warnings

import numpy as np

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
    assert u.pc.is_equivalent(u.arcsecond, u.parallax())
    assert u.arcminute.is_equivalent(u.au, u.parallax())


def test_composite():
    assert (u.cm / u.s * u.h).get_converter(u.m)(1) == 36
    assert u.cm * u.cm == u.cm ** 2

    assert u.cm * u.cm * u.cm == u.cm ** 3

    assert u.Hz.to(1000 * u.Hz, 1) == 0.001


def test_str():
    assert str(u.cm) == "cm"


def test_repr():
    assert repr(u.cm) == 'Unit("cm")'


def test_parallax():
    a = u.arcsecond.to(u.pc, 10, u.parallax())
    assert_allclose(a, 0.10)
    b = u.pc.to(u.arcsecond, a, u.parallax())
    assert_allclose(b, 10)
    
    a = u.arcminute.to(u.au, 1, u.parallax())
    assert_allclose(a, 3437.7467916)
    b = u.au.to(u.arcminute, a, u.parallax())
    assert_allclose(b, 1)


def test_parallax2():
    a = u.arcsecond.to(u.pc, [0.1, 2.5], u.parallax())
    assert_allclose(a, [10, 0.4])
    

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
    assert_allclose(f_nu.to(u.erg / u.cm ** 2 / u.s, 1.,
                    equivalencies=u.spectral_density(u.GHz, 10)), 1.e-13, 10)

    # Convert to ergs / cm^2 / s / micron at 1 Ghz
    assert_allclose(f_nu.to(u.erg / u.cm ** 2 / u.s / u.micron, 1.,
                    equivalencies=u.spectral_density(u.Hz, 1.e9)),
                    3.335640951981521e-20, 10)

    # Define F_lambda in ergs / cm^2 / s / micron
    f_lambda = u.erg / u.cm ** 2 / u.s / u.micron

    # Convert to Jy at 1 Ghz
    assert_allclose(f_lambda.to(u.Jy, 1.,
                    equivalencies=u.spectral_density(u.Hz, 1.e9)),
                    1. / 3.335640951981521e-20, 10)

    # Convert to ergs / cm^2 / s at 10 microns
    assert_allclose(f_lambda.to(u.erg / u.cm ** 2 / u.s, 1.,
                    equivalencies=u.spectral_density(u.micron, 10.)), 10., 10)


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

    assert warning_lines[0].category == u.UnitsWarning
    assert 'FOO' in str(warning_lines[0].message)


def test_unknown_unit2():
    with warnings.catch_warnings(record=True) as warning_lines:
        warnings.resetwarnings()
        warnings.simplefilter("always", u.UnitsWarning, append=True)
        assert u.Unit("m/s/kg", parse_strict='warn').to_string() == 'm/s/kg'

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
            del u.UnitBase._registry['foo']
            del u.foo


def test_in_units():
    speed_unit = u.cm / u.s
    x = speed_unit.in_units(u.mile / u.hour, 1)


def test_null_unit():
    assert (u.m / u.m) == u.Unit(1)


def test_unrecognized_equivalency():
    assert u.m.is_equivalent('foo') is False
    assert u.m.is_equivalent('foot') is True


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
        assert "length" in str(e)


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

    results = u.sr.compose(units=u.cgs.bases)
    assert results[0].bases[0] is u.rad


def test_decompose_bases():
    """
    From issue #576
    """

    from .. import cgs
    from ...constants import e

    d = e.esu.unit.decompose(bases=cgs.bases)
    assert d._bases == [u.cm, u.g, u.s]
    assert d._powers == [Fraction(3, 2), Fraction(1, 2), -1]
    assert d._scale == 1.0


def test_complex_compose():
    complex = u.cd * u.sr * u.Wb
    composed = complex.compose()

    assert set(composed[0]._bases) == set([u.lm, u.Wb])


def test_equiv_compose():
    composed = u.m.compose(equivalencies=u.spectral())
    assert u.Hz in composed


def test_compose_roundtrip():
    def _test_compose_roundtrip(unit):
        composed_list = unit.decompose().compose()
        found = False
        for composed in composed_list:
            if len(composed._bases):
                if composed._bases[0] is unit:
                    found = True
                    break
            elif len(unit._bases) == 0:
                found = True
                break
        assert found

    from ... import units as u

    for val in u.__dict__.values():
        if (isinstance(val, u.UnitBase) and
                not isinstance(val, u.PrefixUnit)):
            yield _test_compose_roundtrip, val


def test_compose_cgs_to_si():
    def _test_compose_cgs_to_si(unit):
        unit.to_system(u.si)

    for val in u.cgs.__dict__.values():
        if (isinstance(val, u.UnitBase) and
                not isinstance(val, u.PrefixUnit)):
            yield _test_compose_cgs_to_si, val


def test_compose_si_to_cgs():
    def _test_compose_si_to_cgs(unit):
        # Can't convert things with Ampere to CGS without more context
        try:
            unit.to_system(u.cgs)
        except u.UnitsException:
            if u.A in unit.decompose().bases:
                pass
            else:
                raise

    for val in u.si.__dict__.values():
        if (isinstance(val, u.UnitBase) and
                not isinstance(val, u.PrefixUnit)):
            yield _test_compose_si_to_cgs, val


def test_to_cgs():
    assert u.Pa.to_system(u.cgs)[0]._bases[0] is u.Ba
    assert u.Pa.to_system(u.cgs)[0]._scale == 10.0


def test_decompose_to_cgs():
    from .. import cgs
    assert u.m.decompose(bases=cgs.bases)._bases[0] is cgs.cm


def test_compose_issue_579():
    unit = u.kg * u.s ** 2 / u.m

    result = unit.compose(units=[u.N, u.s, u.m])

    assert len(result) == 1
    assert result[0]._bases == [u.s, u.N, u.m]
    assert result[0]._powers == [4, 1, -2]


def test_self_compose():
    unit = u.kg * u.s

    assert len(unit.compose(units=[u.g, u.s])) == 1


@raises(u.UnitsException)
def test_compose_failed():
    unit = u.kg

    result = unit.compose(units=[u.N])


def test_compose_fractional_powers():
    x = (u.kg / u.s ** 3 * u.au ** 2.5 / u.yr ** 0.5 / u.sr ** 2)

    factored = x.compose()

    for unit in factored:
        assert x.decompose() == unit.decompose()

    factored = x.compose(units=u.cgs)

    for unit in factored:
        assert x.decompose() == unit.decompose()

    factored = x.compose(units=u.si)

    for unit in factored:
        assert x.decompose() == unit.decompose()


def test_compose_best_unit_first():
    results = u.l.compose()
    assert len(results[0].bases) == 1
    assert results[0].bases[0] is u.l

    results = (u.s ** -1).compose()
    assert results[0].bases[0] is u.Hz

    results = (u.Ry.decompose()).compose()
    assert results[0].bases[0] is u.Ry


def test_compose_no_duplicates():
    new = u.kg / u.s ** 3 * u.au ** 2.5 / u.yr ** 0.5 / u.sr ** 2
    composed = new.compose(units=u.cgs.bases)
    assert len(composed) == 1


def test_long_int():
    """
    Issue #672
    """
    sigma = 10 ** 21 * u.M_p / u.cm ** 2
    sigma.to(u.M_sun / u.pc ** 2)


def test_endian_independence():
    """
    Regression test for #744

    A logic issue in the units code meant that big endian arrays could not be
    converted because the dtype is '>f4', not 'float32', and the code was
    looking for the strings 'float' or 'int'.
    """
    for endian in ['<', '>']:
        for ntype in ['i', 'f']:
            for byte in ['4', '8']:
                # Note, we have to use encode() because we've imported
                # unicode_literals from __future__, and Numpy 1.4.1 crashes if
                # a unicode dtype is passed.
                x = np.array([1,2,3], dtype=(endian + ntype + byte).encode('ascii'))
                u.m.to(u.cm, x)


def test_radian_base():
    """
    Issue #863
    """
    assert (1 * u.degree).si.unit == u.rad


def test_no_as():
    # We don't define 'as', since it is a keyword, but we
    # do want to define the long form (`attosecond`).
    assert not hasattr(u, 'as')
    assert hasattr(u, 'attosecond')


def test_pickling():
    import cPickle

    p = cPickle.dumps(u.m)
    other = cPickle.loads(p)

    assert other is u.m

    new_unit = u.IrreducibleUnit(['foo'], register=False, format={'baz': 'bar'})
    p = cPickle.dumps(new_unit)
    new_unit_copy = cPickle.loads(p)
    assert new_unit is not new_unit_copy
    assert new_unit_copy.names == ['foo']
    assert new_unit_copy.get_format_name('baz') == 'bar'
