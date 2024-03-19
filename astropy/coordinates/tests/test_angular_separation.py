# Licensed under a 3-clause BSD style license - see LICENSE.rst


"""
Tests for the projected separation stuff
"""

import numpy as np
import pytest

from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates.builtin_frames import ICRS

# lon1, lat1, lon2, lat2 in degrees
coords = [
    (1, 0, 0, 0),
    (0, 1, 0, 0),
    (0, 0, 1, 0),
    (0, 0, 0, 1),
    (0, 0, 10, 0),
    (0, 0, 90, 0),
    (0, 0, 180, 0),
    (0, 45, 0, -45),
    (0, 60, 0, -30),
    (-135, -15, 45, 15),
    (100, -89, -80, 89),
    (0, 0, 0, 0),
    (0, 0, 1.0 / 60.0, 1.0 / 60.0),
]
correct_seps = [1, 1, 1, 1, 10, 90, 180, 90, 90, 180, 180, 0, 0.023570225877234643]
correctness_margin = 2e-10


def test_angsep():
    """
    Tests that the angular separation object also behaves correctly.
    """
    from astropy.coordinates import angular_separation

    # check it both works with floats in radians, Quantities, or Angles
    for conv in (np.deg2rad, lambda x: u.Quantity(x, "deg"), lambda x: Angle(x, "deg")):
        for (lon1, lat1, lon2, lat2), corrsep in zip(coords, correct_seps):
            angsep = angular_separation(conv(lon1), conv(lat1), conv(lon2), conv(lat2))
            assert np.fabs(angsep - conv(corrsep)) < conv(correctness_margin)


def test_proj_separations():
    """
    Test angular separation functionality
    """
    c1 = ICRS(ra=0 * u.deg, dec=0 * u.deg)
    c2 = ICRS(ra=0 * u.deg, dec=1 * u.deg)

    # these operations have ambiguous interpretations for points on a sphere
    with pytest.raises(TypeError):
        c1 + c2
    with pytest.raises(TypeError):
        c1 - c2
