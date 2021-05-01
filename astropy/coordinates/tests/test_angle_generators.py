"""Unit tests for the astropy.coordinates.angle_utilities module"""

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates.angle_utilities import (
    golden_spiral_grid,
    uniform_spherical_random_surface,
    uniform_spherical_random_volume
)


def test_golden_spiral_grid_input():
    usph = golden_spiral_grid(size=100)
    assert len(usph) == 100


@pytest.mark.parametrize("func", [uniform_spherical_random_surface,
                                  uniform_spherical_random_volume])
@pytest.mark.parametrize("rng", [None, np.random.default_rng(42)])
def test_uniform_spherical_random_input(func, rng):
    sph = func(size=100, rng=rng)
    assert len(sph) == 100


def test_uniform_spherical_random_volume_input():
    sph = uniform_spherical_random_volume(size=100, distance_scale=1)
    assert len(sph) == 100
    assert sph.distance.unit == u.dimensionless_unscaled
    assert sph.distance.max() <= 1.

    sph = uniform_spherical_random_volume(size=100, distance_scale=4*u.pc)
    assert len(sph) == 100
    assert sph.distance.max() <= 4*u.pc
