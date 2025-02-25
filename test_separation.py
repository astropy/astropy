import numpy as np
import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord


def projected_separation(coord1, coord2, axis_angle):
    sep = coord1.separation(coord2)
    return sep * np.cos(axis_angle.to(u.rad))


def elliptical_separation(coord1, coord2, ellipticity):
    sep = coord1.separation(coord2)
    return sep * (1 + ellipticity)


def test_standard_separation():
    coord1 = SkyCoord(ra=10 * u.deg, dec=20 * u.deg)
    coord2 = SkyCoord(ra=15 * u.deg, dec=25 * u.deg)
    sep = coord1.separation(coord2)

    assert np.isclose(sep.deg, 6.8056, atol=1e-3)


def test_projected_separation():
    coord1 = SkyCoord(ra=10 * u.deg, dec=20 * u.deg)
    coord2 = SkyCoord(ra=15 * u.deg, dec=25 * u.deg)
    axis_angle = 45 * u.deg
    sep = projected_separation(coord1, coord2, axis_angle)
    assert sep.unit == u.deg
    assert sep.value > 0


# Test for elliptical separation
def test_elliptical_separation():
    coord1 = SkyCoord(ra=10 * u.deg, dec=20 * u.deg)
    coord2 = SkyCoord(ra=15 * u.deg, dec=25 * u.deg)
    ellipticity = 0.5
    sep = elliptical_separation(coord1, coord2, ellipticity)
    assert sep.unit == u.deg
    assert sep.value > 0


if __name__ == "__main__":
    pytest.main()
