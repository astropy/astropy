"""
Regression tests for the units.format package
"""

from __future__ import absolute_import, unicode_literals, division, print_function

from astropy.tests.helper import raises, assert_allclose

from ... import units as u
from .. import core
from .. import format
from ... import wcs


def test_unit_grammar():
    def _test_unit_grammar(s, unit):
        print(s)
        unit2 = format.Generic().parse(s)
        assert unit2 == unit

    data = [
        (["m s", "m*s", "m.s"], u.m * u.s),
        (["m/s", "m*s**-1"], u.m / u.s),
        (["m**2", "m2", "m**(2)", "m**+2", "m+2", "m^(+2)"], u.m ** 2),
        (["m**-3", "m-3", "m^(-3)", "/m3"], u.m ** -3),
        (["m**(1.5)", "m(3/2)", "m**(3/2)", "m^(3/2)"], u.m ** 1.5),
        (["2.54cm"], u.cm * 2.54),
        (["10+8m"], u.m * 1e8),
        # # This is the VOUnits documentation, but doesn't seem to follow the unity grammar
        # (["3.45 10**(-4)Jy"], 3.45 * 1e-4 * u.Jy) #
        (["sqrt(m)"], u.m ** -2)
        ]

    for strings, unit in data:
        for s in strings:
            yield _test_unit_grammar, s, unit


def test_roundtrip():
    def _test_roundtrip(unit):
        a = core.Unit(unit.to_string('generic'), format='generic')
        b = core.Unit(unit.decompose().to_string('generic'), format='generic')
        assert_allclose(a.decompose().scale, unit.decompose().scale, rtol=1e-2)
        assert_allclose(b.decompose().scale, unit.decompose().scale, rtol=1e-2)

    for key, val in u.__dict__.items():
        if isinstance(val, core.Unit) and not isinstance(val, core.PrefixUnit):
            yield _test_roundtrip, val


def test_roundtrip_vo_unit():
    def _test_roundtrip_vo_unit(unit):
        a = core.Unit(unit.to_string('vounit'), format='vounit')
        b = core.Unit(unit.decompose().to_string('vounit'), format='vounit')
        assert_allclose(a.decompose().scale, unit.decompose().scale, rtol=1e-2)
        assert_allclose(b.decompose().scale, unit.decompose().scale, rtol=1e-2)

    x = format.VOUnit()
    for key, val in x._units.items():
        if isinstance(val, core.Unit) and not isinstance(val, core.PrefixUnit):
            yield _test_roundtrip_vo_unit, val


def test_roundtrip_fits():
    def _test_roundtrip_fits(unit):
        try:
            s = unit.to_string('fits')
        except ValueError:
            return
        a = core.Unit(s, format='fits')
        assert_allclose(a.decompose().scale, unit.decompose().scale, rtol=1e-2)

    for key, val in u.__dict__.items():
        if isinstance(val, core.Unit) and not isinstance(val, core.PrefixUnit):
            yield _test_roundtrip_fits, val


def test_fits_units_available():
    format.Fits()


def test_vo_units_available():
    format.VOUnit()


def test_latex():
    fluxunit = u.erg / (u.cm ** 2 * u.s)
    assert fluxunit.to_string('latex') == r'$\mathrm{\frac{erg}{s\ cm^{2}}}$'


def test_wcs_parse():
    """
    Tests that the output of to_string('fits') is also parsed by
    wcslib.  Even if we deprecated access to wcslib's unit parser, we
    may want to keep it around and hidden for this test.
    """
    def _test_wcs_parse(unit):
        try:
            fits_string = unit.decompose().to_string('fits')
        except ValueError:
            return
        wcs.UnitConverter(fits_string, fits_string)

    for key, val in u.__dict__.items():
        if isinstance(val, core.Unit) and not isinstance(val, core.PrefixUnit):
            yield _test_wcs_parse, val


def test_flatten_to_known():
    myunit = u.def_unit("FOOBAR", u.erg / u.Hz)
    assert myunit.to_string('fits') == 'erg Hz-1'
    myunit2 = myunit * u.bit ** 3
    assert myunit2.to_string('fits') == 'bit3 erg Hz-1'


@raises(ValueError)
def test_flatten_impossible():
    myunit = u.def_unit("FOOBAR")
    myunit.to_string('fits')
