"""
Tests for the new `_OriginMismatch` StrEnum used in BaseCoordinateFrame.
"""

from __future__ import annotations

import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates.baseframe import _OriginMismatch


@pytest.mark.parametrize("value", [_OriginMismatch.IGNORE, "ignore"])
def test_origin_mismatch_valid(value):
    """
    The helper should accept both the enum and the legacy strings
    and return identical spherical components for the same coordinate.
    """
    sc = SkyCoord("00h00m00s", "00d00m00s")
    fr = sc.frame

    lon1, lat1, lon2, lat2 = fr._prepare_unit_sphere_coords(sc, value)

    # Same coordinate → components should match
    assert lon1 == lon2 and lat1 == lat2
    # Units must be degrees
    assert lon1.unit is u.deg and lat1.unit is u.deg


def test_origin_mismatch_invalid():
    """
    Passing an invalid value must raise the original ValueError
    with the exact message expected by existing tests.
    """
    sc = SkyCoord("00h00m00s", "00d00m00s")
    fr = sc.frame

    msg = (
        r"^origin_mismatch='bogus' is invalid\. "
        r"Allowed values are 'ignore', 'warn' or 'error'\.$"
    )
    with pytest.raises(ValueError, match=msg):
        fr._prepare_unit_sphere_coords(sc, "bogus")
