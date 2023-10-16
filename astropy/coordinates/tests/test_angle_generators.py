"""Unit tests for the :mod:`astropy.coordinates.angles.utils` module."""

import pytest

import astropy.units as u
from astropy.coordinates import (
    golden_spiral_grid,
    uniform_spherical_random_surface,
    uniform_spherical_random_volume,
)
from astropy.utils import NumpyRNGContext


def test_golden_spiral_grid_input():
    usph = golden_spiral_grid(size=100)
    assert len(usph) == 100


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
