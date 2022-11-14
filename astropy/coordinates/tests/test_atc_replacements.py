# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Test replacements for ERFA functions atciqz and aticq."""

import erfa
import pytest

import astropy.units as u
from astropy.coordinates import SphericalRepresentation
from astropy.coordinates.builtin_frames.utils import atciqz, aticq, get_jd12
from astropy.tests.helper import assert_quantity_allclose as assert_allclose
from astropy.time import Time

# Hard-coded random values
sph = SphericalRepresentation(
    lon=[15.0, 214.0] * u.deg, lat=[-12.0, 64.0] * u.deg, distance=[1, 1.0]
)


@pytest.mark.parametrize(
    "t", [Time("2014-06-25T00:00"), Time(["2014-06-25T00:00", "2014-09-24"])]
)
@pytest.mark.parametrize("pos", [sph[0], sph])
def test_atciqz_aticq(t, pos):
    """Check replacements against erfa versions for consistency."""
    jd1, jd2 = get_jd12(t, "tdb")
    astrom, _ = erfa.apci13(jd1, jd2)

    ra = pos.lon.to_value(u.rad)
    dec = pos.lat.to_value(u.rad)
    assert_allclose(erfa.atciqz(ra, dec, astrom), atciqz(pos, astrom))
    assert_allclose(erfa.aticq(ra, dec, astrom), aticq(pos, astrom))
