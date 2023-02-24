# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

import pytest

from astropy import units as u
from astropy.cosmology import Planck15
from astropy.cosmology.units import with_H0
from astropy.units import equivalencies as eq

asdf = pytest.importorskip("asdf", minversion="2.3.0.dev0")

from asdf.exceptions import AsdfDeprecationWarning

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        category=AsdfDeprecationWarning,
        message=r"asdf.tests.helpers is deprecated.*",
    )
    from asdf.tests.helpers import assert_roundtrip_tree


def get_equivalencies():
    """
    Return a list of example equivalencies for testing serialization.
    """
    return [
        eq.plate_scale(0.3 * u.deg / u.mm),
        eq.pixel_scale(0.5 * u.deg / u.pix),
        eq.pixel_scale(100.0 * u.pix / u.cm),
        eq.spectral_density(350 * u.nm, factor=2),
        eq.spectral_density(350 * u.nm),
        eq.spectral(),
        eq.brightness_temperature(500 * u.GHz),
        eq.brightness_temperature(500 * u.GHz, beam_area=23 * u.sr),
        eq.temperature_energy(),
        eq.temperature(),
        eq.thermodynamic_temperature(300 * u.Hz),
        eq.thermodynamic_temperature(140 * u.GHz, Planck15.Tcmb0),
        eq.beam_angular_area(3 * u.sr),
        eq.mass_energy(),
        eq.molar_mass_amu(),
        eq.doppler_relativistic(2 * u.m),
        eq.doppler_optical(2 * u.nm),
        eq.doppler_radio(2 * u.Hz),
        eq.parallax(),
        eq.logarithmic(),
        eq.dimensionless_angles(),
        eq.spectral() + eq.temperature(),
        (
            eq.spectral_density(35 * u.nm)
            + eq.brightness_temperature(5 * u.Hz, beam_area=2 * u.sr)
        ),
        (
            eq.spectral()
            + eq.spectral_density(35 * u.nm)
            + eq.brightness_temperature(5 * u.Hz, beam_area=2 * u.sr)
        ),
        with_H0(),
    ]


@pytest.mark.parametrize("equiv", get_equivalencies())
@pytest.mark.filterwarnings(
    "ignore:`with_H0` is deprecated from `astropy.units.equivalencies` "
    "since astropy 5.0 and may be removed in a future version. "
    "Use `astropy.cosmology.units.with_H0` instead."
)
def test_equivalencies(tmpdir, equiv):
    tree = {"equiv": equiv}
    assert_roundtrip_tree(tree, tmpdir)
