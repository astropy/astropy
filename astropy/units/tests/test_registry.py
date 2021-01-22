import pytest

import astropy.io.fits as fits
import astropy.units as u
from astropy.utils.data import get_pkg_data_filename


@pytest.mark.parametrize(
    "aliases,bad_unit_string,good_unit_string",
    [
        ({"counts": u.count}, "10**(-6) counts/s", "10**(-6) count/s"),
        ({"Angstroms": u.AA}, "Angstroms", "Angstrom"),
        (
            {"ergs": u.erg, "Angstroms": u.AA},
            "10^-17 ergs/(s.cm^2.Angstroms)",
            "10^-17 erg/(s.cm^2.Angstrom)"
        ),
    ]
)
def test_unit_aliases_helper(aliases, bad_unit_string, good_unit_string):
    with u.add_unit_aliases(aliases):
        assert u.Unit(bad_unit_string) == u.Unit(good_unit_string)
