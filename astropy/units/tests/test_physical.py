# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Regression tests for the physical_type support in the units package
"""

import pytest

from astropy import units as u
from astropy.units import physical
from astropy.constants import hbar
from astropy.tests.helper import raises


unit_physical_type_pairs = [
    (u.m, "length"),
    (u.cm ** 3, "volume"),
    (u.km / u.h, "speed"),
    (u.barn * u.Mpc, "volume"),
    (u.m * u.s, "unknown"),
    (u.m / u.m, "dimensionless"),
    (hbar.unit, "angular momentum"),
    (u.erg / (u.cm**2 * u.s * u.AA), "spectral flux density wav"),  # flam
    (u.photon / (u.cm ** 2 * u.s * u.AA), "photon flux density wav"),  # photlam
    (u.photon / (u.cm ** 2 * u.s * u.Hz), "photon flux density"),  # photnu
    (u.byte, "data quantity"),
    (u.bit, "data quantity"),
    (u.J * u.m ** -2 * u.s ** -1, "energy flux"),
    (u.cm ** -3 * u.hr ** -1, "volumetric rate"),
    (u.imperial.mi / u.week, "speed"),
    (u.erg / u.s, "power"),
    (u.C / u.s, "electrical current"),
    (u.C / u.s / u.cm ** 2, "electrical current density"),
    (u.T * u.m ** 2, "magnetic flux"),
    (u.N * u.m, "energy"),
    (u.rad / u.ms, "angular speed"),
]


@pytest.mark.parametrize("unit, physical_type", unit_physical_type_pairs)
def test_physical_types(unit, physical_type):
    """
    Test that the `physical_type` attribute of `Unit` objects provides
    the expected physical type for various units.
    """
    if unit.physical_type != physical_type:
        pytest.fail(
            f"{repr(unit)}.physical_type was expected to return "
            f"{repr(physical_type)}, but instead returned "
            f"{unit.physical_type}."
        )


@raises(ValueError)
def test_redundant_physical_type():
    physical.def_physical_type(u.m, 'kessel run')
