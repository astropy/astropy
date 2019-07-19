# -*- coding: utf-8 -*-

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Regression tests for the units.format package
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from astropy.tests.helper import catch_warnings
from astropy import units as u
from astropy.constants import si
from astropy.units import core
from astropy.units import format as u_format
from astropy.units.utils import is_effectively_unity


@pytest.mark.parametrize('strings, unit', [
    (["m s", "m*s", "m.s"], u.m * u.s),
    (["m/s", "m*s**-1", "m /s", "m / s", "m/ s"], u.m / u.s),
    (["m**2", "m2", "m**(2)", "m**+2", "m+2", "m^(+2)"], u.m ** 2),
    (["m**-3", "m-3", "m^(-3)", "/m3"], u.m ** -3),
    (["m**(1.5)", "m(3/2)", "m**(3/2)", "m^(3/2)"], u.m ** 1.5),
    (["2.54 cm"], u.Unit(u.cm * 2.54)),
    (["10+8m"], u.Unit(u.m * 1e8)),
    # This is the VOUnits documentation, but doesn't seem to follow the
    # unity grammar (["3.45 10**(-4)Jy"], 3.45 * 1e-4 * u.Jy)
    (["sqrt(m)"], u.m ** 0.5),
    (["dB(mW)", "dB (mW)"], u.DecibelUnit(u.mW)),
    (["mag"], u.mag),
    (["mag(ct/s)"], u.MagUnit(u.ct / u.s)),
    (["dex"], u.dex),
    (["dex(cm s**-2)", "dex(cm/s2)"], u.DexUnit(u.cm / u.s**2))])
def test_unit_grammar(strings, unit):
    for s in strings:
        print(s)
        unit2 = u_format.Generic.parse(s)
        assert unit2 == unit


@pytest.mark.parametrize('string', ['sin( /pixel /s)', 'mag(mag)',
                                    'dB(dB(mW))', 'dex()'])
def test_unit_grammar_fail(string):
    with pytest.raises(ValueError):
        print(string)
        u_format.Generic.parse(string)


@pytest.mark.parametrize('strings, unit', [
    (["0.1nm"], u.AA),
    (["mW/m2"], u.Unit(u.erg / u.cm ** 2 / u.s)),
    (["mW/(m2)"], u.Unit(u.erg / u.cm ** 2 / u.s)),
    (["km/s", "km.s-1"], u.km / u.s),
    (["10pix/nm"], u.Unit(10 * u.pix / u.nm)),
    (["1.5x10+11m"], u.Unit(1.5e11 * u.m)),
    (["1.5×10+11m"], u.Unit(1.5e11 * u.m)),
    (["m2"], u.m ** 2),
    (["10+21m"], u.Unit(u.m * 1e21)),
    (["2.54cm"], u.Unit(u.cm * 2.54)),
    (["20%"], 0.20 * u.dimensionless_unscaled),
    (["10+9"], 1.e9 * u.dimensionless_unscaled),
    (["2x10-9"], 2.e-9 * u.dimensionless_unscaled),
    (["---"], u.dimensionless_unscaled),
    (["ma"], u.ma),
    (["mAU"], u.mAU),
    (["uarcmin"], u.uarcmin),
    (["uarcsec"], u.uarcsec),
    (["kbarn"], u.kbarn),
    (["Gbit"], u.Gbit),
    (["Gibit"], 2 ** 30 * u.bit),
    (["kbyte"], u.kbyte),
    (["mRy"], 0.001 * u.Ry),
    (["mmag"], u.mmag),
    (["Mpc"], u.Mpc),
    (["Gyr"], u.Gyr),
    (["°"], u.degree),
    (["°/s"], u.degree / u.s),
    (["Å"], u.AA),
    (["Å/s"], u.AA / u.s),
    (["\\h"], si.h)])
def test_cds_grammar(strings, unit):
    for s in strings:
        print(s)
        unit2 = u_format.CDS.parse(s)
        assert unit2 == unit


@pytest.mark.parametrize('string', [
    '0.1 nm',
    'solMass(3/2)',
    'km / s',
    'km s-1',
    'pix0.1nm',
    'pix/(0.1nm)',
    'km*s',
    'km**2',
    '5x8+3m',
    '0.1---',
    '---m',
    'm---',
    'mag(s-1)',
    'dB(mW)',
    'dex(cm s-2)'])
def test_cds_grammar_fail(string):
    with pytest.raises(ValueError):
        print(string)
        u_format.CDS.parse(string)


# These examples are taken from the EXAMPLES section of
# https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/
@pytest.mark.parametrize('strings, unit', [
    (["count /s", "count/s", "count s**(-1)", "count / s", "count /s "],
     u.count / u.s),
    (["/pixel /s", "/(pixel * s)"], (u.pixel * u.s) ** -1),
    (["count /m**2 /s /eV", "count m**(-2) * s**(-1) * eV**(-1)",
      "count /(m**2 * s * eV)"],
     u.count * u.m ** -2 * u.s ** -1 * u.eV ** -1),
    (["erg /pixel /s /GHz", "erg /s /GHz /pixel", "erg /pixel /(s * GHz)"],
     u.erg / (u.s * u.GHz * u.pixel)),
    (["keV**2 /yr /angstrom", "10**(10) keV**2 /yr /m"],
      # Though this is given as an example, it seems to violate the rules
      # of not raising scales to powers, so I'm just excluding it
      # "(10**2 MeV)**2 /yr /m"
     u.keV**2 / (u.yr * u.angstrom)),
    (["10**(46) erg /s", "10**46 erg /s", "10**(39) J /s", "10**(39) W",
      "10**(15) YW", "YJ /fs"],
     10**46 * u.erg / u.s),
    (["10**(-7) J /cm**2 /MeV", "10**(-9) J m**(-2) eV**(-1)",
      "nJ m**(-2) eV**(-1)", "nJ /m**2 /eV"],
     10 ** -7 * u.J * u.cm ** -2 * u.MeV ** -1),
    (["sqrt(erg /pixel /s /GHz)", "(erg /pixel /s /GHz)**(0.5)",
      "(erg /pixel /s /GHz)**(1/2)",
      "erg**(0.5) pixel**(-0.5) s**(-0.5) GHz**(-0.5)"],
     (u.erg * u.pixel ** -1 * u.s ** -1 * u.GHz ** -1) ** 0.5),
    (["(count /s) (/pixel /s)", "(count /s) * (/pixel /s)",
      "count /pixel /s**2"],
     (u.count / u.s) * (1.0 / (u.pixel * u.s)))])
def test_ogip_grammar(strings, unit):
    for s in strings:
        print(s)
        unit2 = u_format.OGIP.parse(s)
        assert unit2 == unit


@pytest.mark.parametrize('string', [
    'log(photon /m**2 /s /Hz)',
    'sin( /pixel /s)',
    'log(photon /cm**2 /s /Hz) /(sin( /pixel /s))',
    'log(photon /cm**2 /s /Hz) (sin( /pixel /s))**(-1)',
    'dB(mW)', 'dex(cm/s**2)'])
def test_ogip_grammar_fail(string):
    with pytest.raises(ValueError):
        print(string)
        u_format.OGIP.parse(string)


class RoundtripBase:
    deprecated_units = set()

    def check_roundtrip(self, unit):
        with catch_warnings() as w:
            s = unit.to_string(self.format_)
            a = core.Unit(s, format=self.format_)

        if s in self.deprecated_units:
            assert w
            assert 'deprecated' in str(w[0])
        else:
            assert not w

        assert_allclose(a.decompose().scale, unit.decompose().scale, rtol=1e-9)

    def check_roundtrip_decompose(self, unit):
        ud = unit.decompose()
        s = ud.to_string(self.format_)
        assert '  ' not in s
        a = core.Unit(s, format=self.format_)
        assert_allclose(a.decompose().scale, ud.scale, rtol=1e-5)


class TestRoundtripGeneric(RoundtripBase):
    format_ = 'generic'

    @pytest.mark.parametrize('unit', [
        unit for unit in u.__dict__.values()
        if (isinstance(unit, core.UnitBase) and
            not isinstance(unit, core.PrefixUnit))])
    def test_roundtrip(self, unit):
        self.check_roundtrip(unit)
        self.check_roundtrip_decompose(unit)


class TestRoundtripVOUnit(RoundtripBase):
    format_ = 'vounit'
    deprecated_units = u_format.VOUnit._deprecated_units

    @pytest.mark.parametrize('unit', [
        unit for unit in u_format.VOUnit._units.values()
        if (isinstance(unit, core.UnitBase) and
            not isinstance(unit, core.PrefixUnit))])
    def test_roundtrip(self, unit):
        self.check_roundtrip(unit)
        if unit not in (u.mag, u.dB):
            self.check_roundtrip_decompose(unit)


class TestRoundtripFITS(RoundtripBase):
    format_ = 'fits'
    deprecated_units = u_format.Fits._deprecated_units

    @pytest.mark.parametrize('unit', [
        unit for unit in u_format.Fits._units.values()
        if (isinstance(unit, core.UnitBase) and
            not isinstance(unit, core.PrefixUnit))])
    def test_roundtrip(self, unit):
        self.check_roundtrip(unit)


class TestRoundtripCDS(RoundtripBase):
    format_ = 'cds'

    @pytest.mark.parametrize('unit', [
        unit for unit in u_format.CDS._units.values()
        if (isinstance(unit, core.UnitBase) and
            not isinstance(unit, core.PrefixUnit))])
    def test_roundtrip(self, unit):
        self.check_roundtrip(unit)
        if unit == u.mag:
            # Skip mag: decomposes into dex, which is unknown to CDS.
            return

        self.check_roundtrip_decompose(unit)


class TestRoundtripOGIP(RoundtripBase):
    format_ = 'ogip'
    deprecated_units = u_format.OGIP._deprecated_units | {'d'}

    @pytest.mark.parametrize('unit', [
        unit for unit in u_format.OGIP._units.values()
        if (isinstance(unit, core.UnitBase) and
            not isinstance(unit, core.PrefixUnit))])
    def test_roundtrip(self, unit):
        if str(unit) in ('d', '0.001 Crab'):
            # Special-case day, which gets auto-converted to hours, and mCrab,
            # which the default check does not recognize as a deprecated unit.
            with catch_warnings() as w:
                s = unit.to_string(self.format_)
                a = core.Unit(s, format=self.format_)
            assert w
            assert 'UnitsWarning' in str(w[0])
            assert_allclose(a.decompose().scale, unit.decompose().scale, rtol=1e-9)
        else:
            self.check_roundtrip(unit)
        if str(unit) in ('mag', 'byte', 'Crab'):
            # Skip mag and byte, which decompose into dex and bit, resp.,
            # both of which are unknown to OGIP, as well as Crab, which does
            # not decompose, and thus gives a depecated unit warning.
            return

        with catch_warnings() as w:
            self.check_roundtrip_decompose(unit)

        power_of_ten = np.log10(unit.decompose().scale)
        if abs(power_of_ten - round(power_of_ten)) > 1e-3:
            assert w
            assert 'power of 10' in str(w[0])
        elif str(unit) == '0.001 Crab':
            assert w
            assert 'deprecated' in str(w[0])
        else:
            assert not w


def test_fits_units_available():
    u_format.Fits._units


def test_vo_units_available():
    u_format.VOUnit._units


def test_cds_units_available():
    u_format.CDS._units


def test_cds_non_ascii_unit():
    """Regression test for #5350.  This failed with a decoding error as
    μas could not be represented in ascii."""
    from astropy.units import cds
    with cds.enable():
        u.radian.find_equivalent_units(include_prefix_units=True)


def test_latex():
    fluxunit = u.erg / (u.cm ** 2 * u.s)
    assert fluxunit.to_string('latex') == r'$\mathrm{\frac{erg}{s\,cm^{2}}}$'


def test_new_style_latex():
    fluxunit = u.erg / (u.cm ** 2 * u.s)
    assert f"{fluxunit:latex}" == r'$\mathrm{\frac{erg}{s\,cm^{2}}}$'


def test_latex_scale():
    fluxunit = u.Unit(1.e-24 * u.erg / (u.cm ** 2 * u.s * u.Hz))
    latex = r'$\mathrm{1 \times 10^{-24}\,\frac{erg}{Hz\,s\,cm^{2}}}$'
    assert fluxunit.to_string('latex') == latex


def test_latex_inline_scale():
    fluxunit = u.Unit(1.e-24 * u.erg / (u.cm ** 2 * u.s * u.Hz))
    latex_inline = (r'$\mathrm{1 \times 10^{-24}\,erg'
                    r'\,Hz^{-1}\,s^{-1}\,cm^{-2}}$')
    assert fluxunit.to_string('latex_inline') == latex_inline


@pytest.mark.parametrize('format_spec, string', [
    ('generic', 'erg / (cm2 s)'),
    ('s', 'erg / (cm2 s)'),
    ('console', '  erg  \n ------\n s cm^2'),
    ('latex', '$\\mathrm{\\frac{erg}{s\\,cm^{2}}}$'),
    ('latex_inline', '$\\mathrm{erg\\,s^{-1}\\,cm^{-2}}$'),
    ('>20s', '       erg / (cm2 s)')])
def test_format_styles(format_spec, string):
    fluxunit = u.erg / (u.cm ** 2 * u.s)
    assert format(fluxunit, format_spec) == string


def test_flatten_to_known():
    myunit = u.def_unit("FOOBAR_One", u.erg / u.Hz)
    assert myunit.to_string('fits') == 'erg Hz-1'
    myunit2 = myunit * u.bit ** 3
    assert myunit2.to_string('fits') == 'bit3 erg Hz-1'


def test_flatten_impossible():
    myunit = u.def_unit("FOOBAR_Two")
    with u.add_enabled_units(myunit), pytest.raises(ValueError):
        myunit.to_string('fits')


def test_console_out():
    """
    Issue #436.
    """
    u.Jy.decompose().to_string('console')


def test_flexible_float():
    assert u.min._represents.to_string('latex') == r'$\mathrm{60\,s}$'


def test_fraction_repr():
    area = u.cm ** 2.0
    assert '.' not in area.to_string('latex')

    fractional = u.cm ** 2.5
    assert '5/2' in fractional.to_string('latex')

    assert fractional.to_string('unicode') == 'cm⁵⸍²'


def test_scale_effectively_unity():
    """Scale just off unity at machine precision level is OK.
    Ensures #748 does not recur
    """
    a = (3. * u.N).cgs
    assert is_effectively_unity(a.unit.scale)
    assert len(a.__repr__().split()) == 3


def test_percent():
    """Test that the % unit is properly recognized.  Since % is a special
    symbol, this goes slightly beyond the round-tripping tested above."""
    assert u.Unit('%') == u.percent == u.Unit(0.01)

    assert u.Unit('%', format='cds') == u.Unit(0.01)
    assert u.Unit(0.01).to_string('cds') == '%'

    with pytest.raises(ValueError):
        u.Unit('%', format='fits')

    with pytest.raises(ValueError):
        u.Unit('%', format='vounit')


def test_scaled_dimensionless():
    """Test that scaled dimensionless units are properly recognized in generic
    and CDS, but not in fits and vounit."""
    assert u.Unit('0.1') == u.Unit(0.1) == 0.1 * u.dimensionless_unscaled
    assert u.Unit('1.e-4') == u.Unit(1.e-4)

    assert u.Unit('10-4', format='cds') == u.Unit(1.e-4)
    assert u.Unit('10+8').to_string('cds') == '10+8'

    with pytest.raises(ValueError):
        u.Unit(0.15).to_string('fits')

    assert u.Unit(0.1).to_string('fits') == '10**-1'

    with pytest.raises(ValueError):
        u.Unit(0.1).to_string('vounit')


def test_deprecated_did_you_mean_units():
    with pytest.raises(ValueError) as exc_info:
        u.Unit('ANGSTROM', format='fits')
    assert 'Did you mean Angstrom or angstrom?' in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        u.Unit('crab', format='ogip')
    assert 'Crab (deprecated)' in str(exc_info.value)
    assert 'mCrab (deprecated)' in str(exc_info.value)

    with catch_warnings() as w:
        u.Unit('ANGSTROM', format='vounit')

    assert len(w) == 1
    assert 'angstrom (deprecated)' in str(w[0].message)
    assert '0.1nm' in str(w[0].message)
    assert str(w[0].message).count('0.1nm') == 1

    with catch_warnings() as w:
        u.Unit('angstrom', format='vounit')
    assert len(w) == 1
    assert '0.1nm' in str(w[0].message)


@pytest.mark.parametrize('string', ['mag(ct/s)', 'dB(mW)', 'dex(cm s**-2)'])
def test_fits_function(string):
    # Function units cannot be written, so ensure they're not parsed either.
    with pytest.raises(ValueError):
        print(string)
        u_format.Fits().parse(string)


@pytest.mark.parametrize('string', ['mag(ct/s)', 'dB(mW)', 'dex(cm s**-2)'])
def test_vounit_function(string):
    # Function units cannot be written, so ensure they're not parsed either.
    with pytest.raises(ValueError):
        print(string)
        with catch_warnings():  # ct, dex also raise warnings - irrelevant here.
            u_format.VOUnit().parse(string)


def test_vounit_binary_prefix():
    u.Unit('KiB', format='vounit') == u.Unit('1024 B')
    u.Unit('Kibyte', format='vounit') == u.Unit('1024 B')
    u.Unit('Kibit', format='vounit') == u.Unit('1024 B')
    with catch_warnings() as w:
        u.Unit('kibibyte', format='vounit')
    assert len(w) == 1


def test_vounit_unknown():
    assert u.Unit('unknown', format='vounit') is None
    assert u.Unit('UNKNOWN', format='vounit') is None
    assert u.Unit('', format='vounit') is u.dimensionless_unscaled


def test_vounit_details():
    with catch_warnings() as w:
        assert u.Unit('Pa', format='vounit') is u.Pascal
    assert len(w) == 1
    assert 'deprecated' in str(w[0].message)

    # The da- prefix is not allowed, and the d- prefix is discouraged
    assert u.dam.to_string('vounit') == '10m'
    assert u.Unit('dam dag').to_string('vounit') == '100g m'


def test_vounit_custom():
    x = u.Unit("'foo' m", format='vounit')
    x_vounit = x.to_string('vounit')
    assert x_vounit == "'foo' m"
    x_string = x.to_string()
    assert x_string == "foo m"

    x = u.Unit("m'foo' m", format='vounit')
    assert x.bases[1]._represents.scale == 0.001
    x_vounit = x.to_string('vounit')
    assert x_vounit == "m m'foo'"
    x_string = x.to_string()
    assert x_string == 'm mfoo'


def test_vounit_implicit_custom():
    # Yikes, this becomes "femto-urlong"...  But at least there's a warning.
    with catch_warnings() as w:
        x = u.Unit("furlong/week", format="vounit")
    assert x.bases[0]._represents.scale == 1e-15
    assert x.bases[0]._represents.bases[0].name == 'urlong'
    assert len(w) == 2
    assert 'furlong' in str(w[0].message)
    assert 'week' in str(w[1].message)


@pytest.mark.parametrize('scale, number, string',
                         [('10+2', 100, '10**2'),
                          ('10(+2)', 100, '10**2'),
                          ('10**+2', 100, '10**2'),
                          ('10**(+2)', 100, '10**2'),
                          ('10^+2', 100, '10**2'),
                          ('10^(+2)', 100, '10**2'),
                          ('10**2', 100, '10**2'),
                          ('10**(2)', 100, '10**2'),
                          ('10^2', 100, '10**2'),
                          ('10^(2)', 100, '10**2'),
                          ('10-20', 10**(-20), '10**-20'),
                          ('10(-20)', 10**(-20), '10**-20'),
                          ('10**-20', 10**(-20), '10**-20'),
                          ('10**(-20)', 10**(-20), '10**-20'),
                          ('10^-20', 10**(-20), '10**-20'),
                          ('10^(-20)', 10**(-20), '10**-20'),
                          ])
def test_fits_scale_factor(scale, number, string):

    x = u.Unit(scale + ' erg/(s cm**2 Angstrom)', format='fits')
    assert x == number * (u.erg / u.s / u.cm ** 2 / u.Angstrom)
    assert x.to_string(format='fits') == string + ' Angstrom-1 cm-2 erg s-1'

    x = u.Unit(scale + '*erg/(s cm**2 Angstrom)', format='fits')
    assert x == number * (u.erg / u.s / u.cm ** 2 / u.Angstrom)
    assert x.to_string(format='fits') == string + ' Angstrom-1 cm-2 erg s-1'


def test_fits_scale_factor_errors():
    with pytest.raises(ValueError):
        x = u.Unit('1000 erg/(s cm**2 Angstrom)', format='fits')

    with pytest.raises(ValueError):
        x = u.Unit('12 erg/(s cm**2 Angstrom)', format='fits')

    x = u.Unit(1.2 * u.erg)
    with pytest.raises(ValueError):
        x.to_string(format='fits')

    x = u.Unit(100.0 * u.erg)
    assert x.to_string(format='fits') == '10**2 erg'
