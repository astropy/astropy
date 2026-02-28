# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for projected separation / elliptical distance feature.
"""

import numpy as np

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord


def test_projected_components_simple():
    c1 = SkyCoord(ra=0 * u.deg, dec=0 * u.deg)
    c2 = SkyCoord(ra=1 * u.deg, dec=0 * u.deg)

    # project onto PA=90 deg (eastwards major axis)
    dmaj, dmin = c1.frame.separation_projected(
        c2, pa=Angle(90, "deg"), return_components=True
    )
    assert np.allclose(dmaj.to(u.deg).value, 1.0, atol=1e-8)
    assert np.allclose(dmin.to(u.deg).value, 0.0, atol=1e-8)


def test_elliptical_distance_basic():
    c1 = SkyCoord(ra=0 * u.deg, dec=0 * u.deg)
    c2 = SkyCoord(ra=1 * u.deg, dec=1 * u.deg)

    # get components
    dmaj, dmin = c1.frame.separation_projected(
        c2, pa=Angle(0, "deg"), return_components=True
    )
    assert hasattr(dmaj, "to") and hasattr(dmin, "to")

    # elliptical distance with b_over_a = 0.5 should return an Angle
    ed = c1.frame.separation_projected(c2, pa=Angle(0, "deg"), b_over_a=0.5)
    assert hasattr(ed, "to")
    assert ed.unit == u.deg
