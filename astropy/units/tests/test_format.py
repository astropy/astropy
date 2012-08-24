from astropy.units import format
from astropy import units as u


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


def test_fits_units_available():
    format.Fits()


def test_vo_units_available():
    format.VOUnit()


def test_latex():
    fluxunit = u.erg / (u.cm ** 2 * u.s)
    assert fluxunit.to_string('latex') == r'$\mathrm{\frac{erg}{s\ cm^{2}}}$'
