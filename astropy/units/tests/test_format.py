# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Regression tests for the units.format package
"""

from __future__ import annotations

import re
import warnings
from contextlib import nullcontext
from fractions import Fraction
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy import units as u
from astropy.constants import si
from astropy.units import (
    PrefixUnit,
    Unit,
    UnitBase,
    UnitParserWarning,
    UnitsWarning,
    cds,
    dex,
)
from astropy.units import format as u_format
from astropy.units.utils import is_effectively_unity
from astropy.utils.exceptions import AstropyDeprecationWarning

if TYPE_CHECKING:
    from collections.abc import Iterable


class FormatStringPair(NamedTuple):
    format: str
    string: str


class StringUnitPair(NamedTuple):
    string: str
    unit: UnitBase


def list_format_string_pairs(*test_cases: tuple[str, str]) -> list[FormatStringPair]:
    return [FormatStringPair(format, string) for format, string in test_cases]


def list_string_unit_pairs(
    *test_cases: tuple[Iterable[str], UnitBase],
) -> list[StringUnitPair]:
    return [
        StringUnitPair(string, unit)
        for strings, unit in test_cases
        for string in strings
    ]


@pytest.mark.parametrize(
    "test_pair",
    list_string_unit_pairs(
        (["m s", "m*s", "m.s"], u.m * u.s),
        (["m/s", "m*s**-1", "m /s", "m / s", "m/ s"], u.m / u.s),
        (["m**2", "m2", "m**(2)", "m**+2", "m+2", "m^(+2)"], u.m**2),
        (["m**-3", "m-3", "m^(-3)", "/m3"], u.m**-3),
        (["m**(1.5)", "m(3/2)", "m**(3/2)", "m^(3/2)"], u.m**1.5),
        (["2.54 cm"], u.Unit(u.cm * 2.54)),
        (["10+8m"], u.Unit(u.m * 1e8)),
        # This is the VOUnits documentation, but doesn't seem to follow the
        # unity grammar (["3.45 10**(-4)Jy"], 3.45 * 1e-4 * u.Jy)
        (["sqrt(m)"], u.m**0.5),
        (["dB(mW)", "dB (mW)"], u.DecibelUnit(u.mW)),
        (["mag"], u.mag),
        (["mag(ct/s)"], u.MagUnit(u.ct / u.s)),
        (["dex"], u.dex),
        (["dex(cm s**-2)", "dex(cm/s2)"], u.DexUnit(u.cm / u.s**2)),
    ),
    ids=lambda x: x.string,
)
def test_unit_grammar(test_pair: StringUnitPair):
    assert u_format.Generic.parse(test_pair.string) == test_pair.unit


@pytest.mark.parametrize(
    "string", ["sin( /pixel /s)", "mag(mag)", "dB(dB(mW))", "dex()"]
)
def test_unit_grammar_fail(string):
    with pytest.raises(ValueError):
        u_format.Generic.parse(string)


@pytest.mark.parametrize(
    "test_pair",
    list_string_unit_pairs(
        (["0.1nm"], u.AA),
        (["mW/m2"], u.Unit(u.erg / u.cm**2 / u.s)),
        (["mW/(m2)"], u.Unit(u.erg / u.cm**2 / u.s)),
        (["km/s", "km.s-1"], u.km / u.s),
        (["km/s/Mpc"], u.km / u.s / u.Mpc),
        (["km/(s.Mpc)"], u.km / u.s / u.Mpc),
        (["10+3J/m/s/kpc2"], u.Unit(1e3 * u.W / (u.m * u.kpc**2))),
        (["10pix/nm"], u.Unit(10 * u.pix / u.nm)),
        (["1.5x10+11m"], u.Unit(1.5e11 * u.m)),
        (["1.5×10+11/m"], u.Unit(1.5e11 / u.m)),
        (["/s"], u.s**-1),
        (["m2"], u.m**2),
        (["10+21m"], u.Unit(u.m * 1e21)),
        (["2.54cm"], u.Unit(u.cm * 2.54)),
        (["20%"], 0.20 * u.dimensionless_unscaled),
        (["10+9"], 1.0e9 * u.dimensionless_unscaled),
        (["2x10-9"], 2.0e-9 * u.dimensionless_unscaled),
        (["---"], u.dimensionless_unscaled),
        (["ma"], u.ma),
        (["mAU"], u.mAU),
        (["uarcmin"], u.uarcmin),
        (["uarcsec"], u.uarcsec),
        (["kbarn"], u.kbarn),
        (["Gbit"], u.Gbit),
        (["Gibit"], 2**30 * u.bit),
        (["kbyte"], u.kbyte),
        (["mRy"], 0.001 * u.Ry),
        (["mmag"], u.mmag),
        (["Mpc"], u.Mpc),
        (["Gyr"], u.Gyr),
        (["°"], u.degree),
        (["°/s"], u.degree / u.s),
        (["Å"], u.AA),
        (["Å/s"], u.AA / u.s),
        (["\\h"], si.h),
        (["[cm/s2]"], dex(u.cm / u.s**2)),
        (["[K]"], dex(u.K)),
        (["[-]"], dex(u.dimensionless_unscaled)),
        (["eps0/mu0"], cds.eps0 / cds.mu0),
        (["a0.s"], cds.a0 * u.s),
    ),
    ids=lambda x: x.string,
)
def test_cds_grammar(test_pair: StringUnitPair):
    assert u_format.CDS.parse(test_pair.string) == test_pair.unit


@pytest.mark.parametrize(
    "string",
    [
        "0.1 nm",
        "solMass(3/2)",
        "km / s",
        "km s-1",
        "km/s.Mpc-1",
        "/s.Mpc",
        "pix0.1nm",
        "pix/(0.1nm)",
        "km*s",
        "km**2",
        "5x8+3m",
        "0.1---",
        "---m",
        "m---",
        "--",
        "0.1-",
        "-m",
        "m-",
        "mag(s-1)",
        "dB(mW)",
        "dex(cm s-2)",
        "[--]",
    ],
)
def test_cds_grammar_fail(string):
    with pytest.raises(ValueError):
        u_format.CDS.parse(string)


def test_cds_dimensionless():
    assert u.Unit("---", format="cds") == u.dimensionless_unscaled
    assert u.dimensionless_unscaled.to_string(format="cds") == "---"


def test_cds_log10_dimensionless():
    assert u.Unit("[-]", format="cds") == u.dex(u.dimensionless_unscaled)
    assert u.dex(u.dimensionless_unscaled).to_string(format="cds") == "[-]"


def test_cds_angstrom_str():
    # Regression test for a problem noticed in
    # https://github.com/astropy/astropy/pull/17527#discussion_r1880555481
    # that the string representation of the cds version of Angstrom was "AA".
    assert str(u.cds.Angstrom) == str(u.Angstrom) == "Angstrom"
    # Since this is a NamedUnit, let's check the name for completeness.
    assert u.cds.Angstrom.name == "Angstrom"


def test_cds_solMass_str():
    # CDS allows writing solar mass as Msun or solMass,
    # but cds.solMass and u.solMass should be consistent.
    assert u.solMass.to_string("cds") == "solMass"
    assert u.cds.solMass.to_string("cds") == "solMass"


# These examples are taken from the EXAMPLES section of
# https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/
@pytest.mark.parametrize(
    "test_pair",
    list_string_unit_pairs(
        (
            ["count /s", "count/s", "count s**(-1)", "count / s", "count /s "],
            u.count / u.s,
        ),
        (
            ["/pixel /s", "/(pixel * s)"],
            (u.pixel * u.s) ** -1,
        ),
        (
            [
                "count /m**2 /s /eV",
                "count m**(-2) * s**(-1) * eV**(-1)",
                "count /(m**2 * s * eV)",
            ],
            u.count * u.m**-2 * u.s**-1 * u.eV**-1,
        ),
        (
            ["erg /pixel /s /GHz", "erg /s /GHz /pixel", "erg /pixel /(s * GHz)"],
            u.erg / (u.s * u.GHz * u.pixel),
        ),
        (
            ["keV**2 /yr /angstrom", "10**(10) keV**2 /yr /m"],
            # Though this is given as an example, it seems to violate the rules
            # of not raising scales to powers, so I'm just excluding it
            # "(10**2 MeV)**2 /yr /m"
            u.keV**2 / (u.yr * u.angstrom),
        ),
        (
            [
                "10**(46) erg /s",
                "10**46 erg /s",
                "10**(39) J /s",
                "10**(39) W",
                "10**(15) YW",
                "YJ /fs",
            ],
            10**46 * u.erg / u.s,
        ),
        (
            [
                "10**(-7) J /cm**2 /MeV",
                "10**(-9) J m**(-2) eV**(-1)",
                "nJ m**(-2) eV**(-1)",
                "nJ /m**2 /eV",
            ],
            10**-7 * u.J * u.cm**-2 * u.MeV**-1,
        ),
        (
            [
                "sqrt(erg /pixel /s /GHz)",
                "(erg /pixel /s /GHz)**(0.5)",
                "(erg /pixel /s /GHz)**(1/2)",
                "erg**(0.5) pixel**(-0.5) s**(-0.5) GHz**(-0.5)",
            ],
            (u.erg * u.pixel**-1 * u.s**-1 * u.GHz**-1) ** 0.5,
        ),
        (
            [
                "(count /s) (/pixel /s)",
                "(count /s) * (/pixel /s)",
                "count /pixel /s**2",
            ],
            (u.count / u.s) * (1.0 / (u.pixel * u.s)),
        ),
    ),
    ids=lambda x: x.string,
)
def test_ogip_grammar(test_pair: StringUnitPair):
    assert u_format.OGIP.parse(test_pair.string) == test_pair.unit


@pytest.mark.parametrize(
    "string",
    [
        "log(photon /m**2 /s /Hz)",
        "sin( /pixel /s)",
        "log(photon /cm**2 /s /Hz) /(sin( /pixel /s))",
        "log(photon /cm**2 /s /Hz) (sin( /pixel /s))**(-1)",
        "dB(mW)",
        "dex(cm/s**2)",
    ],
)
def test_ogip_grammar_fail(string):
    with pytest.raises(ValueError):
        u_format.OGIP.parse(string)


@pytest.mark.xfail(reason="'acos' is not understood by astropy", raises=ValueError)
def test_ogip_unusable_function():
    u_format.OGIP.parse("acos(m)**2")


@pytest.mark.parametrize("string", ["sqrt(m)**3", "sqrt(m**3)", "(sqrt(m))**3"])
def test_ogip_sqrt(string):
    # Regression test for #16743 - sqrt(m)**3 caused a ValueError
    assert u_format.OGIP.parse(string) == u.m ** Fraction(3, 2)


@pytest.mark.parametrize(
    "string,message,unit",
    [
        pytest.param(
            "m(s)**2",
            (
                r"^if 'm\(s\)\*\*2' was meant to be a multiplication, "
                r"it should have been written as 'm \(s\)\*\*2'.$"
            ),
            u.m * u.s**2,
            id="m(s)**2",
        ),
        pytest.param(
            "m(s)",
            (
                r"^if 'm\(s\)' was meant to be a multiplication, "
                r"it should have been written as 'm \(s\)'.$"
            ),
            u.m * u.s,
            id="m(s)",
        ),
    ],
)
def test_ogip_invalid_multiplication(string, message, unit):
    # Regression test for #16749
    with pytest.warns(UnitParserWarning, match=message):
        assert u_format.OGIP.parse(string) == unit


@pytest.mark.parametrize(
    "string,unit,power",
    [
        pytest.param("s**-1", u.s**-1, "-1", id="int_unit_power"),
        pytest.param("m**-2.0", u.m**-2, "-2.0", id="float_unit_power"),
        pytest.param("10**-3 kg", u.g, "-3", id="int_scale_power"),
    ],
)
def test_ogip_negative_exponent_parenthesis(string, unit, power):
    # Regression test for #16788 - negative powers require parenthesis
    with pytest.warns(
        UnitParserWarning,
        match=(
            r"^negative exponents must be enclosed in parenthesis\. "
            rf"Expected '\*\*\({power}\)' instead of '\*\*{power}'\.$"
        ),
    ):
        assert u_format.OGIP.parse(string) == unit


def test_ogip_ohm():
    # Regression test for #17200 - OGIP converted u.ohm to 'V / A'
    assert u_format.OGIP.to_string(u.ohm) == "ohm"


class RoundtripBase:
    def check_roundtrip(self, unit, output_format=None):
        if output_format is None:
            output_format = self.format_.name
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # Same warning shows up multiple times
            s = unit.to_string(output_format)

        if s in self.format_._deprecated_units:
            with pytest.warns(UnitsWarning, match="deprecated") as w:
                a = Unit(s, format=self.format_)
            assert len(w) == 1
        else:
            a = Unit(s, format=self.format_)  # No warning

        assert_allclose(a.decompose().scale, unit.decompose().scale, rtol=1e-9)

    def check_roundtrip_decompose(self, unit):
        ud = unit.decompose()
        s = ud.to_string(self.format_)
        assert "  " not in s
        a = Unit(s, format=self.format_)
        assert_allclose(a.decompose().scale, ud.scale, rtol=1e-5)


class TestRoundtripGeneric(RoundtripBase):
    format_ = u_format.Generic

    @pytest.mark.parametrize(
        "unit",
        [
            unit
            for unit in u.__dict__.values()
            if (isinstance(unit, UnitBase) and not isinstance(unit, PrefixUnit))
        ],
        ids=str,
    )
    def test_roundtrip(self, unit):
        self.check_roundtrip(unit)
        self.check_roundtrip(unit, output_format="unicode")
        self.check_roundtrip_decompose(unit)


class TestRoundtripVOUnit(RoundtripBase):
    format_ = u_format.VOUnit

    @pytest.mark.parametrize(
        "unit",
        [u for u in u_format.VOUnit._units.values() if not isinstance(u, PrefixUnit)],
        ids=str,
    )
    def test_roundtrip(self, unit):
        self.check_roundtrip(unit)
        if unit not in (u.mag, u.dB):
            self.check_roundtrip_decompose(unit)


class TestRoundtripFITS(RoundtripBase):
    format_ = u_format.FITS

    @pytest.mark.parametrize(
        "unit",
        [u for u in u_format.FITS._units.values() if not isinstance(u, PrefixUnit)],
        ids=str,
    )
    def test_roundtrip(self, unit):
        self.check_roundtrip(unit)


class TestRoundtripCDS(RoundtripBase):
    format_ = u_format.CDS

    @pytest.mark.parametrize(
        "unit",
        [u for u in u_format.CDS._units.values() if not isinstance(u, PrefixUnit)],
        ids=str,
    )
    def test_roundtrip(self, unit):
        self.check_roundtrip(unit)
        if unit == u.mag:
            # Skip mag: decomposes into dex, which is unknown to CDS.
            return

        self.check_roundtrip_decompose(unit)

    @pytest.mark.parametrize(
        "unit", [u.dex(unit) for unit in (u.cm / u.s**2, u.K, u.Lsun)], ids=str
    )
    def test_roundtrip_dex(self, unit):
        string = unit.to_string(format="cds")
        recovered = u.Unit(string, format="cds")
        assert recovered == unit


class TestRoundtripOGIP(RoundtripBase):
    format_ = u_format.OGIP

    @pytest.mark.parametrize(
        "unit",
        [
            unit
            for unit in u_format.OGIP._units.values()
            if (isinstance(unit, UnitBase) and not isinstance(unit, PrefixUnit))
        ],
        ids=str,
    )
    def test_roundtrip(self, unit):
        if str(unit) == "0.001 Crab":
            # Special-case mCrab, which the default check does not recognize
            # as a deprecated unit.
            with pytest.warns(UnitsWarning):
                s = unit.to_string(self.format_)
                a = Unit(s, format=self.format_)
            assert_allclose(a.decompose().scale, unit.decompose().scale, rtol=1e-9)
        else:
            self.check_roundtrip(unit)
        if str(unit) in ("mag", "byte", "Crab"):
            # Skip mag and byte, which decompose into dex and bit, resp.,
            # both of which are unknown to OGIP, as well as Crab, which does
            # not decompose, and thus gives a deprecated unit warning.
            return

        power_of_ten = np.log10(unit.decompose().scale)
        if abs(power_of_ten - round(power_of_ten)) > 1e-3:
            ctx = pytest.warns(UnitsWarning, match="power of 10")
        elif str(unit) == "0.001 Crab":
            ctx = pytest.warns(UnitsWarning, match="deprecated")
        else:
            ctx = nullcontext()
        with ctx:
            self.check_roundtrip_decompose(unit)


@pytest.mark.parametrize(
    "unit_formatter_class,n_units",
    [(u_format.FITS, 765), (u_format.VOUnit, 1303), (u_format.CDS, 3326)],
)
def test_units_available(unit_formatter_class, n_units):
    assert len(unit_formatter_class._units) == n_units


def test_cds_non_ascii_unit():
    """Regression test for #5350.  This failed with a decoding error as
    μas could not be represented in ascii."""
    with cds.enable():
        u.radian.find_equivalent_units(include_prefix_units=True)


def test_latex():
    fluxunit = u.erg / (u.cm**2 * u.s)
    assert fluxunit.to_string("latex") == r"$\mathrm{\frac{erg}{s\,cm^{2}}}$"


def test_new_style_latex():
    fluxunit = u.erg / (u.cm**2 * u.s)
    assert f"{fluxunit:latex}" == r"$\mathrm{\frac{erg}{s\,cm^{2}}}$"


def test_latex_scale():
    fluxunit = u.Unit(1.0e-24 * u.erg / (u.cm**2 * u.s * u.Hz))
    latex = r"$\mathrm{1 \times 10^{-24}\,\frac{erg}{Hz\,s\,cm^{2}}}$"
    assert fluxunit.to_string("latex") == latex


def test_latex_inline_scale():
    fluxunit = u.Unit(1.0e-24 * u.erg / (u.cm**2 * u.s * u.Hz))
    latex_inline = r"$\mathrm{1 \times 10^{-24}\,erg\,Hz^{-1}\,s^{-1}\,cm^{-2}}$"
    assert fluxunit.to_string("latex_inline") == latex_inline


@pytest.mark.parametrize(
    "format_spec, string, decomposed",
    [
        ("generic", "erg / (Angstrom s cm2)", "1e+07 kg / (m s3)"),
        ("s", "erg / (Angstrom s cm2)", "1e+07 kg / (m s3)"),
        ("console", "erg Angstrom^-1 s^-1 cm^-2", "10000000 kg m^-1 s^-3"),
        (
            "latex",
            r"$\mathrm{\frac{erg}{\mathring{A}\,s\,cm^{2}}}$",
            r"$\mathrm{10000000\,\frac{kg}{m\,s^{3}}}$",
        ),
        (
            "latex_inline",
            r"$\mathrm{erg\,\mathring{A}^{-1}\,s^{-1}\,cm^{-2}}$",
            r"$\mathrm{10000000\,kg\,m^{-1}\,s^{-3}}$",
        ),
        ("unicode", "erg Å⁻¹ s⁻¹ cm⁻²", "10000000 kg m⁻¹ s⁻³"),
        (">25s", "   erg / (Angstrom s cm2)", "        1e+07 kg / (m s3)"),
        ("cds", "erg.Angstrom-1.s-1.cm-2", "10000000kg.m-1.s-3"),
        ("ogip", "erg / (angstrom s cm**2)", "1e+07 kg / (m s**3)"),
        ("fits", "erg Angstrom-1 s-1 cm-2", "10**7 kg m-1 s-3"),
        ("vounit", "erg.Angstrom**-1.s**-1.cm**-2", "10000000kg.m**-1.s**-3"),
        # TODO: make fits and vounit less awful!
    ],
)
def test_format_styles(format_spec, string, decomposed):
    fluxunit = u.erg / (u.cm**2 * u.s * u.Angstrom)
    if format_spec == "vounit":
        # erg and Angstrom are deprecated in vounit.
        with pytest.warns(UnitsWarning, match="deprecated"):
            formatted = format(fluxunit, format_spec)
    else:
        formatted = format(fluxunit, format_spec)
    assert formatted == string
    # Decomposed mostly to test that scale factors are dealt with properly
    # in the various formats.
    assert format(fluxunit.decompose(), format_spec) == decomposed


@pytest.mark.parametrize(
    "format_spec, fraction, string, decomposed",
    [
        ("generic", False, "erg s-1 cm-2", "0.001 kg s-3"),
        (
            "console",
            "multiline",
            " erg  \n------\ns cm^2",
            "      kg \n0.001 ---\n      s^3",
        ),
        ("console", "inline", "erg / (s cm^2)", "0.001 kg / s^3"),
        ("unicode", "multiline", " erg \n─────\ns cm²", "      kg\n0.001 ──\n      s³"),
        ("unicode", "inline", "erg / (s cm²)", "0.001 kg / s³"),
        (
            "latex",
            False,
            r"$\mathrm{erg\,s^{-1}\,cm^{-2}}$",
            r"$\mathrm{0.001\,kg\,s^{-3}}$",
        ),
        (
            "latex",
            "inline",
            r"$\mathrm{erg / (s\,cm^{2})}$",
            r"$\mathrm{0.001\,kg / s^{3}}$",
        ),
        # TODO: make generic with fraction=False less awful!
    ],
)
def test_format_styles_non_default_fraction(format_spec, fraction, string, decomposed):
    fluxunit = u.erg / (u.cm**2 * u.s)
    assert fluxunit.to_string(format_spec, fraction=fraction) == string
    assert fluxunit.decompose().to_string(format_spec, fraction=fraction) == decomposed


@pytest.mark.parametrize("format_spec", u_format.Base.registry)
def test_multiline_fraction_different_if_available(format_spec):
    fluxunit = u.W / u.m**2
    inline_format = fluxunit.to_string(format_spec, fraction="inline")
    if format_spec in ["generic", "cds", "fits", "ogip", "vounit"]:
        with pytest.warns(UnitsWarning, match="does not support multiline"):
            multiline_format = fluxunit.to_string(format_spec, fraction="multiline")
        assert multiline_format == inline_format
    else:
        multiline_format = fluxunit.to_string(format_spec, fraction="multiline")
        assert multiline_format != inline_format


@pytest.mark.parametrize("format_spec", u_format.Base.registry)
def test_unknown_fraction_style(format_spec):
    fluxunit = u.W / u.m**2
    msg = "fraction can only be False, 'inline', or 'multiline', not 'parrot'"
    with pytest.raises(ValueError, match=msg):
        fluxunit.to_string(format_spec, fraction="parrot")


def test_flatten_to_known():
    myunit = u.def_unit("FOOBAR_One", u.erg / u.Hz)
    assert myunit.to_string("fits") == "erg Hz-1"
    myunit2 = myunit * u.bit**3
    assert myunit2.to_string("fits") == "bit3 erg Hz-1"


def test_flatten_impossible():
    myunit = u.def_unit("FOOBAR_Two")
    with u.add_enabled_units(myunit), pytest.raises(ValueError):
        myunit.to_string("fits")


def test_console_out():
    """
    Issue #436.
    """
    u.Jy.decompose().to_string("console")


@pytest.mark.parametrize(
    "test_pair",
    list_format_string_pairs(
        ("generic", "10"),
        ("console", "10"),
        ("unicode", "10"),
        ("cds", "10"),
        ("latex", r"$\mathrm{10}$"),
    ),
    ids=lambda x: x.format,
)
def test_scale_only(test_pair: FormatStringPair):
    assert u.Unit(10).to_string(test_pair.format) == test_pair.string


def test_flexible_float():
    assert u.min._represents.to_string("latex") == r"$\mathrm{60\,s}$"


def test_fits_to_string_function_error():
    """Test function raises TypeError on bad input.

    This instead of returning None, see gh-11825.
    """

    with pytest.raises(TypeError, match="unit argument must be"):
        u_format.FITS.to_string(None)


def test_fraction_repr():
    area = u.cm**2.0
    assert "." not in area.to_string("latex")

    fractional = u.cm**2.5
    assert "5/2" in fractional.to_string("latex")

    assert fractional.to_string("unicode") == "cm⁵⸍²"


def test_scale_effectively_unity():
    """Scale just off unity at machine precision level is OK.
    Ensures #748 does not recur
    """
    a = (3.0 * u.N).cgs
    assert is_effectively_unity(a.unit.scale)
    assert len(a.__repr__().split()) == 3


def test_percent():
    """Test that the % unit is properly recognized.  Since % is a special
    symbol, this goes slightly beyond the round-tripping tested above."""
    assert u.Unit("%") == u.percent == u.Unit(0.01)

    assert u.Unit("%", format="cds") == u.Unit(0.01)
    assert u.Unit(0.01).to_string("cds") == "%"

    with pytest.raises(ValueError):
        u.Unit("%", format="fits")

    with pytest.raises(ValueError):
        u.Unit("%", format="vounit")


def test_scaled_dimensionless():
    """Test that scaled dimensionless units are properly recognized in generic
    and CDS, but not in fits and vounit."""
    assert u.Unit("0.1") == u.Unit(0.1) == 0.1 * u.dimensionless_unscaled
    assert u.Unit("1.e-4") == u.Unit(1.0e-4)

    assert u.Unit("10-4", format="cds") == u.Unit(1.0e-4)
    assert u.Unit("10+8").to_string("cds") == "10+8"

    with pytest.raises(ValueError):
        u.Unit(0.15).to_string("fits")

    assert u.Unit(0.1).to_string("fits") == "10**-1"

    with pytest.raises(ValueError):
        u.Unit(0.1).to_string("vounit")


def test_deprecated_did_you_mean_units():
    with pytest.raises(ValueError) as exc_info:
        u.Unit("ANGSTROM", format="fits")
    assert "Did you mean Angstrom or angstrom?" in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        u.Unit("crab", format="ogip")
    assert "Crab (deprecated)" in str(exc_info.value)
    assert "mCrab (deprecated)" in str(exc_info.value)

    with pytest.raises(
        ValueError,
        match=(
            r"Did you mean 0\.1nm, Angstrom \(deprecated\) or angstrom \(deprecated\)\?"
        ),
    ):
        u.Unit("ANGSTROM", format="vounit")

    with pytest.warns(UnitsWarning, match=r".* 0\.1nm\.") as w:
        u.Unit("angstrom", format="vounit")
    assert len(w) == 1


@pytest.mark.parametrize("string", ["mag(ct/s)", "dB(mW)", "dex(cm s**-2)"])
def test_fits_function(string):
    # Function units cannot be written, so ensure they're not parsed either.
    with pytest.raises(ValueError):
        u_format.FITS().parse(string)


@pytest.mark.parametrize("string", ["mag(ct/s)", "dB(mW)", "dex(cm s**-2)"])
def test_vounit_function(string):
    # Function units cannot be written, so ensure they're not parsed either.
    with pytest.raises(ValueError), warnings.catch_warnings():
        # ct, dex also raise warnings - irrelevant here.
        warnings.simplefilter("ignore")
        u_format.VOUnit().parse(string)


def test_vounit_binary_prefix():
    assert u.Unit("KiB", format="vounit") == u.Unit("1024 B")
    assert u.Unit("Kibyte", format="vounit") == u.Unit("1024 B")
    assert u.Unit("Kibit", format="vounit") == u.Unit("128 B")
    with pytest.raises(ValueError, match="not supported by the VOUnit standard"):
        u.Unit("kibibyte", format="vounit")


def test_vounit_unknown():
    assert u.Unit("unknown", format="vounit") is None
    assert u.Unit("UNKNOWN", format="vounit") is None
    assert u.Unit("", format="vounit") is u.dimensionless_unscaled


def test_vounit_details():
    assert u.Unit("Pa", format="vounit") is u.Pascal
    assert u.Unit("ka", format="vounit") == u.Unit("1000 yr")
    assert u.Unit("pix", format="vounit") == u.Unit("pixel", format="vounit")

    # Regression test for astropy/astroquery#2480
    assert u.Unit("Sun", format="vounit") is u.Sun

    # Test that adding a prefix to a simple units raises a warning
    with pytest.warns(
        UnitsWarning, match="Unit 'kdB' not supported by the VOUnit standard.*"
    ):
        u.Unit("kdB", format="vounit", parse_strict="warn")

    # The da- prefix is not allowed, and the d- prefix is discouraged
    assert u.dam.to_string("vounit") == "10m"
    assert u.Unit("dam dag").to_string("vounit") == "100g.m"

    # Parse round-trip
    with pytest.warns(UnitsWarning, match="deprecated"):
        flam = u.erg / u.cm / u.cm / u.s / u.AA
        x = u.format.VOUnit.to_string(flam)
        assert x == "erg.Angstrom**-1.s**-1.cm**-2"
        new_flam = u.format.VOUnit.parse(x)
        assert new_flam == flam


@pytest.mark.parametrize(
    "unit, vounit, number, scale, voscale",
    [
        ("nm", "nm", 0.1, "10^-1", "0.1"),
        ("fm", "fm", 100.0, "10+2", "100"),
        ("m^2", "m**2", 100.0, "100.0", "100"),
        ("cm", "cm", 2.54, "2.54", "2.54"),
        ("kg", "kg", 1.898124597e27, "1.898124597E27", "1.8981246e+27"),
        ("m/s", "m.s**-1", 299792458.0, "299792458", "2.9979246e+08"),
        ("cm2", "cm**2", 1.0e-20, "10^(-20)", "1e-20"),
    ],
)
def test_vounit_scale_factor(unit, vounit, number, scale, voscale):
    x = u.Unit(f"{scale} {unit}")
    assert x == number * u.Unit(unit)
    assert x.to_string(format="vounit") == voscale + vounit


@pytest.mark.parametrize(
    "unit, vounit",
    [
        ("m s^-1", "m/s"),
        ("s^-1", "1/s"),
        ("100 s^-2", "100/s**2"),
        ("kg m-1 s-2", "kg/(m.s**2)"),
    ],
)
@pytest.mark.parametrize("fraction", [True, "inline"])
def test_vounit_fraction(unit, vounit, fraction):
    x = u.Unit(unit)
    assert x.to_string(format="vounit", fraction=fraction) == vounit


@pytest.mark.parametrize(
    "unit, vounit",
    [
        ("m^2", "m**2"),
        ("s^-1", "s**-1"),
        ("s(0.333)", "s**(0.333)"),
        ("s(-0.333)", "s**(-0.333)"),
        ("s(1/3)", "s**(1/3)"),
        ("s(-1/3)", "s**(-1/3)"),
    ],
)
def test_vounit_power(unit, vounit):
    x = u.Unit(unit)
    assert x.to_string(format="vounit") == vounit


def test_vounit_custom():
    x = u.Unit("'foo' m", format="vounit")
    x_vounit = x.to_string("vounit")
    assert x_vounit == "'foo'.m"
    x_string = x.to_string()
    assert x_string == "foo m"

    x = u.Unit("m'foo' m", format="vounit")
    assert x.bases[1]._represents.scale == 0.001
    x_vounit = x.to_string("vounit")
    assert x_vounit == "m.m'foo'"
    x_string = x.to_string()
    assert x_string == "m mfoo"


def test_vounit_implicit_custom():
    # Yikes, this becomes "femto-urlong"...  But at least there's a warning.
    with pytest.warns(UnitsWarning) as w:
        x = u.Unit("furlong/week", format="vounit", parse_strict="warn")
    assert x.bases[0]._represents.scale == 1e-15
    assert x.bases[0]._represents.bases[0].name == "urlong"
    assert len(w) == 2
    assert "furlong" in str(w[0].message)
    assert "week" in str(w[1].message)


@pytest.mark.parametrize(
    "scale, number, string",
    [
        ("10+2", 100, "10**2"),
        ("10(+2)", 100, "10**2"),
        ("10**+2", 100, "10**2"),
        ("10**(+2)", 100, "10**2"),
        ("10^+2", 100, "10**2"),
        ("10^(+2)", 100, "10**2"),
        ("10**2", 100, "10**2"),
        ("10**(2)", 100, "10**2"),
        ("10^2", 100, "10**2"),
        ("10^(2)", 100, "10**2"),
        ("10-20", 10 ** (-20), "10**-20"),
        ("10(-20)", 10 ** (-20), "10**-20"),
        ("10**-20", 10 ** (-20), "10**-20"),
        ("10**(-20)", 10 ** (-20), "10**-20"),
        ("10^-20", 10 ** (-20), "10**-20"),
        ("10^(-20)", 10 ** (-20), "10**-20"),
    ],
)
def test_fits_scale_factor(scale, number, string):
    x = u.Unit(scale + " erg/(s cm**2 Angstrom)", format="fits")
    assert x == number * (u.erg / u.s / u.cm**2 / u.Angstrom)
    assert x.to_string(format="fits") == string + " erg Angstrom-1 s-1 cm-2"

    x = u.Unit(scale + "*erg/(s cm**2 Angstrom)", format="fits")
    assert x == number * (u.erg / u.s / u.cm**2 / u.Angstrom)
    assert x.to_string(format="fits") == string + " erg Angstrom-1 s-1 cm-2"


def test_fits_scale_factor_errors():
    with pytest.raises(ValueError):
        x = u.Unit("1000 erg/(s cm**2 Angstrom)", format="fits")

    with pytest.raises(ValueError):
        x = u.Unit("12 erg/(s cm**2 Angstrom)", format="fits")

    x = u.Unit(1.2 * u.erg)
    with pytest.raises(ValueError):
        x.to_string(format="fits")

    x = u.Unit(100.0 * u.erg)
    assert x.to_string(format="fits") == "10**2 erg"


@pytest.mark.parametrize(
    "unit, latex, unicode",
    [
        (u.deg, r"$\mathrm{{}^{\circ}}$", "°"),
        (u.deg**2, r"$\mathrm{deg^{2}}$", "deg²"),
        (u.arcmin, r"$\mathrm{{}^{\prime}}$", "′"),
        (u.arcmin**2, r"$\mathrm{arcmin^{2}}$", "arcmin²"),
        (u.arcsec, r"$\mathrm{{}^{\prime\prime}}$", "″"),
        (u.arcsec**2, r"$\mathrm{arcsec^{2}}$", "arcsec²"),
        (u.hourangle, r"$\mathrm{{}^{h}}$", "ʰ"),
        (u.hourangle**2, r"$\mathrm{hourangle^{2}}$", "hourangle²"),
        (u.electron, r"$\mathrm{e^{-}}$", "e⁻"),
        (u.electron**2, r"$\mathrm{electron^{2}}$", "electron²"),
    ],
)
def test_double_superscript(unit, latex, unicode):
    """Regression test for #5870, #8699, #9218, #14403; avoid double superscripts."""
    assert unit.to_string("latex") == latex
    assert unit.to_string("unicode") == unicode


def test_no_prefix_superscript():
    """Regression test for gh-911 and #14419."""
    assert u.mdeg.to_string("latex") == r"$\mathrm{mdeg}$"
    assert u.narcmin.to_string("latex") == r"$\mathrm{narcmin}$"
    assert u.parcsec.to_string("latex") == r"$\mathrm{parcsec}$"
    assert u.mdeg.to_string("unicode") == "mdeg"
    assert u.narcmin.to_string("unicode") == "narcmin"
    assert u.parcsec.to_string("unicode") == "parcsec"


@pytest.mark.parametrize(
    "power,expected",
    (
        (1.0, "m"),
        (2.0, "m2"),
        (-10, "1 / m10"),
        (1.5, "m(3/2)"),
        (2 / 3, "m(2/3)"),
        (7 / 11, "m(7/11)"),
        (-1 / 64, "1 / m(1/64)"),
        (1 / 100, "m(1/100)"),
        (2 / 101, "m(0.019801980198019802)"),
        (Fraction(2, 101), "m(2/101)"),
    ),
)
def test_powers(power, expected):
    """Regression test for #9279 - powers should not be oversimplified."""
    unit = u.m**power
    s = unit.to_string()
    assert s == expected
    assert unit == s


@pytest.mark.parametrize(
    "string,unit",
    [
        ("\N{MICRO SIGN}g", u.microgram),
        ("\N{GREEK SMALL LETTER MU}g", u.microgram),
        ("g\N{MINUS SIGN}1", u.g ** (-1)),
        ("m\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}", u.m**-1),
        ("m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}", u.m / u.s),
        ("m\N{SUPERSCRIPT TWO}", u.m**2),
        ("m\N{SUPERSCRIPT PLUS SIGN}\N{SUPERSCRIPT TWO}", u.m**2),
        ("m\N{SUPERSCRIPT THREE}", u.m**3),
        ("m\N{SUPERSCRIPT ONE}\N{SUPERSCRIPT ZERO}", u.m**10),
        ("\N{GREEK CAPITAL LETTER OMEGA}", u.ohm),
        ("\N{OHM SIGN}", u.ohm),  # deprecated but for compatibility
        ("\N{MICRO SIGN}\N{GREEK CAPITAL LETTER OMEGA}", u.microOhm),
        ("\N{ANGSTROM SIGN}", u.Angstrom),
        ("\N{ANGSTROM SIGN} \N{OHM SIGN}", u.Angstrom * u.Ohm),
        ("\N{LATIN CAPITAL LETTER A WITH RING ABOVE}", u.Angstrom),
        ("\N{LATIN CAPITAL LETTER A}\N{COMBINING RING ABOVE}", u.Angstrom),
        ("m\N{ANGSTROM SIGN}", u.milliAngstrom),
        ("°C", u.deg_C),
        ("°", u.deg),
        ("M⊙", u.Msun),  # \N{CIRCLED DOT OPERATOR}
        ("L☉", u.Lsun),  # \N{SUN}
        ("M⊕", u.Mearth),  # normal earth symbol = \N{CIRCLED PLUS}
        ("M♁", u.Mearth),  # be generous with \N{EARTH}
        ("R♃", u.Rjup),  # \N{JUPITER}
        ("′", u.arcmin),  # \N{PRIME}
        ("R∞", u.Ry),
        ("Mₚ", u.M_p),
    ],
)
def test_unicode(string, unit):
    assert u_format.Generic.parse(string) == unit
    assert u.Unit(string) == unit
    # Should work in composites too.
    assert u.Unit(f"{string}/s") == unit / u.s
    assert u.Unit(f"m {string}") == u.m * unit
    assert u.Unit(f"{string} {string}") == unit**2
    # Not obvious that "°2" should be "deg**2", but not easy to reject,
    # and "R♃²" should work.  But don't run on examples with a space or that
    # already end in a number.
    if re.match(r"^\S*[^\d⁰¹²³⁴⁵⁶⁷⁸⁹]$", string):
        assert u.Unit(f"{string}2") == unit**2
    assert u.Unit(f"{string}/{string}") == u.dimensionless_unscaled
    # Finally, check round-trip
    assert u.Unit(unit.to_string("unicode")) == unit


@pytest.mark.parametrize(
    "string",
    [
        "g\N{MICRO SIGN}",
        "g\N{MINUS SIGN}",
        "m\N{SUPERSCRIPT MINUS}1",
        "m+\N{SUPERSCRIPT ONE}",
        "m\N{MINUS SIGN}\N{SUPERSCRIPT ONE}",
        "k\N{ANGSTROM SIGN}",
    ],
)
def test_unicode_failures(string):
    with pytest.raises(ValueError):
        u.Unit(string)


@pytest.mark.parametrize("format_", ("unicode", "latex", "latex_inline"))
def test_parse_error_message_for_output_only_format(format_):
    with pytest.raises(NotImplementedError, match="not parse"):
        u.Unit("m", format=format_)


@pytest.mark.parametrize(
    "parser,error_type,err_msg_start",
    [
        pytest.param("foo", ValueError, "Unknown format 'foo'", id="ValueError"),
        pytest.param(
            {}, TypeError, "Expected a formatter name, not {}", id="TypeError"
        ),
    ],
)
def test_unknown_parser(parser, error_type, err_msg_start):
    with pytest.raises(
        error_type,
        match=(
            f"^{err_msg_start}\\.\nValid parser names are: "
            "'cds', 'generic', 'fits', 'ogip', 'vounit'$"
        ),
    ):
        u.Unit("m", format=parser)


@pytest.mark.parametrize(
    "formatter,error_type,err_msg_start",
    [
        pytest.param("abc", ValueError, "Unknown format 'abc'", id="ValueError"),
        pytest.param(
            float,
            TypeError,
            "Expected a formatter name, not <class 'float'>",
            id="TypeError",
        ),
    ],
)
def test_unknown_output_format(formatter, error_type, err_msg_start):
    with pytest.raises(
        error_type,
        match=(
            f"^{err_msg_start}\\.\nValid formatter names are: "
            "'cds', 'console', 'generic', 'fits', 'latex', 'latex_inline', 'ogip', "
            "'unicode', 'vounit'$"
        ),
    ):
        u.m.to_string(formatter)


def test_celsius_fits():
    assert u.Unit("Celsius", format="fits") == u.deg_C
    assert u.Unit("deg C", format="fits") == u.deg_C

    # check that compounds do what we expect: what do we expect?
    assert u.Unit("deg C kg-1", format="fits") == u.C * u.deg / u.kg
    assert u.Unit("Celsius kg-1", format="fits") == u.deg_C / u.kg
    assert u.deg_C.to_string("fits") == "Celsius"


@pytest.mark.parametrize(
    "test_pair",
    list_format_string_pairs(
        ("generic", "dB(1 / m)"),
        ("latex", r"$\mathrm{dB}$$\mathrm{\left( \mathrm{\frac{1}{m}} \right)}$"),
        ("latex_inline", r"$\mathrm{dB}$$\mathrm{\left( \mathrm{m^{-1}} \right)}$"),
        ("console", "dB(m^-1)"),
        ("unicode", "dB(m⁻¹)"),
    ),
    ids=lambda x: x.format,
)
def test_function_format_styles(test_pair: FormatStringPair):
    dbunit = u.decibel(u.m**-1)
    assert dbunit.to_string(test_pair.format) == test_pair.string
    assert f"{dbunit:{test_pair.format}}" == test_pair.string


@pytest.mark.parametrize(
    "format_spec, fraction, string",
    [
        ("console", "multiline", "   1\ndB(-)\n   m"),
        ("console", "inline", "dB(1 / m)"),
        ("unicode", "multiline", "   1\ndB(─)\n   m"),
        ("unicode", "inline", "dB(1 / m)"),
        ("latex", False, r"$\mathrm{dB}$$\mathrm{\left( \mathrm{m^{-1}} \right)}$"),
        ("latex", "inline", r"$\mathrm{dB}$$\mathrm{\left( \mathrm{1 / m} \right)}$"),
    ],
)
def test_function_format_styles_non_default_fraction(format_spec, fraction, string):
    dbunit = u.decibel(u.m**-1)
    assert dbunit.to_string(format_spec, fraction=fraction) == string


@pytest.mark.parametrize(
    "test_pair",
    list_format_string_pairs(
        ("", "1"),
        (".1g", "1"),
        (".3g", "1"),
        (".1e", "1.0"),
        (".1f", "1.0"),
        (".3e", "1.000"),
    ),
    ids=lambda x: repr(x.format),
)
def test_format_latex_one(test_pair: FormatStringPair):
    # see https://github.com/astropy/astropy/issues/12571
    assert (
        u_format.Latex.format_exponential_notation(1, test_pair.format)
        == test_pair.string
    )


def test_Fits_name_deprecation():
    with pytest.warns(
        AstropyDeprecationWarning,
        match=(
            r'^The class "Fits" has been renamed to "FITS" in version 7\.0\. '
            r"The old name is deprecated and may be removed in a future version\.\n"
            r"        Use FITS instead\.$"
        ),
    ):
        from astropy.units.format import Fits
    assert Fits is u.format.FITS
