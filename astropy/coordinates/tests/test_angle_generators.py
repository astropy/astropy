"""Unit tests for the :mod:`astropy.coordinates.angles.utils` module."""

import pytest
from numpy.testing import assert_allclose

import astropy.units as u
from astropy.coordinates import (
    golden_spiral_grid,
    uniform_spherical_random_surface,
    uniform_spherical_random_volume,
)
from astropy.utils import NumpyRNGContext


@pytest.mark.parametrize(
    "size,lon,lat",
    [
        (3, [1.94161104, 5.82483312, 3.42486989], [0.72972766, 0.0, -0.72972766]),
        (
            6,
            [1.94161104, 5.82483312, 3.4248699, 1.0249067, 4.908129, 2.5081655],
            [0.9851108, 0.5235988, 0.16744808, -0.16744808, -0.5235988, -0.9851108],
        ),
    ],
)
def test_golden_spiral_grid_input(size, lon, lat):
    grid = golden_spiral_grid(size)
    assert_allclose(grid.lon.to_value(u.rad), lon)
    assert_allclose(grid.lat.to_value(u.rad), lat)


@pytest.mark.parametrize(
    "func", [uniform_spherical_random_surface, uniform_spherical_random_volume]
)
def test_uniform_spherical_random_input(func):
    with NumpyRNGContext(42):
        sph = func(size=100)
        assert len(sph) == 100


def test_uniform_spherical_random_volume_input():
    with NumpyRNGContext(42):
        sph = uniform_spherical_random_volume(size=100, max_radius=1)
        assert len(sph) == 100
        assert sph.distance.unit == u.dimensionless_unscaled
        assert sph.distance.max() <= 1.0

        sph = uniform_spherical_random_volume(size=100, max_radius=4 * u.pc)
        assert len(sph) == 100
        assert sph.distance.max() <= 4 * u.pc
