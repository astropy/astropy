# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Regression tests for the units package
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import warnings

import numpy as np

from ...tests.helper import pytest, raises, catch_warnings
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


def test_invalid_power():
    x = u.m ** (1, 3)
    assert isinstance(x.powers[0], Fraction)

    x = u.m ** (1, 2)
    assert isinstance(x.powers[0], float)


def test_invalid_compare():
    assert not (u.m == u.s)


def test_convert():
    assert u.h.get_converter(u.s)(1) == 3600


def test_convert_fail():
    with pytest.raises(u.UnitsError):
        u.cm.to(u.s, 1)
    with pytest.raises(u.UnitsError):
        (u.cm / u.s).to(u.m, 1)


def test_composite():
    assert (u.cm / u.s * u.h).get_converter(u.m)(1) == 36
    assert u.cm * u.cm == u.cm ** 2

    assert u.cm * u.cm * u.cm == u.cm ** 3

    assert u.Hz.to(1000 * u.Hz, 1) == 0.001


def test_str():
    assert str(u.cm) == "cm"


def test_repr():
    assert repr(u.cm) == 'Unit("cm")'


def test_unicode():
    assert unicode(u.m / u.s) == u' m\n \u2500\n s'


def test_units_conversion():
    assert_allclose(u.kpc.to(u.Mpc), 0.001)
    assert_allclose(u.Mpc.to(u.kpc), 1000)
    assert_allclose(u.yr.to(u.Myr), 1.e-6)
    assert_allclose(u.AU.to(u.pc), 4.84813681e-6)
    assert_allclose(u.cycle.to(u.rad), 6.283185307179586)


def test_units_manipulation():
    # Just do some manipulation and check it's happy
    (u.kpc * u.yr) ** (1, 3) / u.Myr
    (u.AA * u.erg) ** 9


def test_decompose():
    assert u.Ry == u.Ry.decompose()


def test_dimensionless_to_si():
    """
    Issue #1150: Test for conversion of dimensionless quantities
                 to the SI system
    """

    testunit = ((1.0 * u.kpc) / (1.0 * u.Mpc))

    assert testunit.unit.physical_type == 'dimensionless'
    assert_allclose(testunit.si, 0.001)


def test_dimensionless_to_cgs():
    """
    Issue #1150: Test for conversion of dimensionless quantities
                 to the CGS system
    """

    testunit = ((1.0 * u.m) / (1.0 * u.km))

    assert testunit.unit.physical_type == 'dimensionless'
    assert_allclose(testunit.cgs, 0.001)


def test_unknown_unit():
    with catch_warnings(u.UnitsWarning) as warning_lines:
        u.Unit("FOO", parse_strict='warn')

    assert 'FOO' in str(warning_lines[0].message)


def test_unknown_unit2():
    with catch_warnings(u.UnitsWarning) as warning_lines:
        assert u.Unit("m/s/kg", parse_strict='warn').to_string() == 'm/s/kg'

    assert 'm/s/kg' in str(warning_lines[0].message)


def test_unknown_unit3():
    unit = u.Unit("FOO", parse_strict='silent')
    assert isinstance(unit, u.UnrecognizedUnit)
    assert unit.name == "FOO"

    unit2 = u.Unit("FOO", parse_strict='silent')
    assert unit == unit2
    assert unit.is_equivalent(unit2)

    unit3 = u.Unit("BAR", parse_strict='silent')
    assert unit != unit3
    assert not unit.is_equivalent(unit3)

    with pytest.raises(ValueError):
        unit.get_converter(unit3)

    x = unit.to_string('latex')
    y = unit2.to_string('cgs')

    with pytest.raises(ValueError):
        unit4 = u.Unit("BAR", parse_strict='strict')

    with pytest.raises(TypeError):
        unit5 = u.Unit(None)


@raises(TypeError)
def test_invalid_scale():
    x = ['a', 'b', 'c'] * u.m


def test_cds_power():
    unit = u.Unit("10+22/cm2", format="cds", parse_strict='silent')
    assert unit.scale == 1e22


def test_register():
    from .. import core
    try:
        u.def_unit("foo", u.m ** 3, register=True)
        assert hasattr(u, 'foo')
    finally:
        u.foo.deregister(remove_from_namespace=True)


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
    except u.UnitsError as e:
        assert "length" in str(e)


def test_convertible_exception2():
    try:
        u.m.to(u.s)
    except u.UnitsError as e:
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

    results = u.sr.compose(units=u.cgs.__dict__)
    assert results[0].bases[0] is u.sr


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


def test_empty_compose():
    with pytest.raises(u.UnitsException):
        composed = u.m.compose(units=[])


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
        except u.UnitsError:
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


@raises(u.UnitsError)
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
    assert results[0].bases[0] in (u.Hz, u.Bq)

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

    new_unit = u.IrreducibleUnit(['foo'], register=True, format={'baz': 'bar'})
    new_unit.deregister()
    p = cPickle.dumps(new_unit)
    new_unit_copy = cPickle.loads(p)
    assert new_unit_copy.names == ['foo']
    assert new_unit_copy.get_format_name('baz') == 'bar'


@raises(ValueError)
def test_duplicate_define():
    u.def_unit('m')


def test_all_units():
    from ...units.core import _UnitRegistry
    assert len(_UnitRegistry().all_units) > len(_UnitRegistry().non_prefix_units)


def test_repr_latex():
    assert u.m._repr_latex_() == u.m.to_string('latex')


def test_operations_with_strings():
    assert u.m / '5s' == (u.m / (5.0 * u.s))

    assert u.m * '5s' == (5.0 * u.m * u.s)


def test_comparison():
    assert u.m > u.cm
    assert u.m >= u.cm
    assert u.cm < u.m
    assert u.cm <= u.m

    with pytest.raises(u.UnitsException):
        u.m > u.kg
