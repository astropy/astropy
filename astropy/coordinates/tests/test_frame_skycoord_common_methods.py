# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tests for methods that are implemented in `BaseCoordinateFrame`,
but are also exposed by `SkyCoord` instances.
"""

import pytest

from astropy import units as u
from astropy.coordinates import FK5, ICRS, SkyCoord
from astropy.tests.helper import assert_quantity_allclose


@pytest.mark.parametrize("self_class", [SkyCoord, ICRS])
@pytest.mark.parametrize("other_class", [SkyCoord, ICRS])
def test_position_angle_array(self_class, other_class):
    c1 = self_class(0 * u.deg, [0, 1, 2] * u.deg)
    c2 = other_class([-1, -2, -3] * u.deg, [0.1, 1.1, 2.1] * u.deg)
    assert_quantity_allclose(
        c1.position_angle(c2), [275.710887, 272.880921, 271.963607] * u.deg
    )


@pytest.mark.parametrize("self_class", [SkyCoord, ICRS])
@pytest.mark.parametrize(
    "other_coord,value",
    [
        (SkyCoord(1 * u.deg, 0 * u.deg), 90 * u.deg),
        (SkyCoord(1 * u.deg, 0.1 * u.deg), 84.289113 * u.deg),
        (SkyCoord(1 * u.deg, 1 * u.deg), 44.995636455344844 * u.deg),
        (ICRS(0 * u.deg, 1 * u.deg), 0 * u.deg),
        (FK5(1 * u.deg, 0 * u.deg), 90.000139 * u.deg),
    ],
)
def test_position_angle_scalar(self_class, other_coord, value):
    assert_quantity_allclose(
        self_class(0 * u.deg, 0 * u.deg).position_angle(other_coord), value
    )
