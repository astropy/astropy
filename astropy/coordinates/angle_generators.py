# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains functions for generating angles, either randomly or gridded
"""

# Third-party
import numpy as np

# Astropy
import astropy.units as u
from astropy.coordinates.representation import (
    UnitSphericalRepresentation,
    SphericalRepresentation)


def golden_spiral_grid(size):
    """Generate a grid of points on the surface of the unit sphere using the
    Fibonacci or Golden Spiral method.

    .. seealso::

        `Evenly distributing points on a sphere <https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere>`_

    Parameters
    ----------
    size : int
        The number of points to generate.

    Returns
    -------
    rep : `~astropy.coordinates.UnitSphericalRepresentation`
        The grid of points.
    """
    golden_r = (1 + 5**0.5) / 2

    grid = np.arange(0, size, dtype=float) + 0.5
    lon = 2*np.pi / golden_r * grid * u.rad
    lat = np.arcsin(1 - 2 * grid / size) * u.rad

    return UnitSphericalRepresentation(lon, lat)


def uniform_spherical_random_surface(size=1, rng=None):
    """Generate a random sampling of points on the surface of the unit sphere.

    Parameters
    ----------
    size : int
        The number of points to generate.
    rng : `numpy.random.Generator`, optional
        A random number generator instance.

    Returns
    -------
    rep : `~astropy.coordinates.UnitSphericalRepresentation`
        The random points.
    """

    if rng is None:
        rng = np.random.default_rng()

    lon = rng.uniform(0, 2*np.pi, size) * u.rad
    lat = 90*u.deg - np.arccos(2 * rng.uniform(size=size) - 1) * u.rad

    return UnitSphericalRepresentation(lon, lat)


@u.quantity_input(distance_scale=u.pc)
def uniform_spherical_random_volume(size=1, distance_scale=1*u.pc, rng=None):
    """Generate a random sampling of points in a spherical volume.

    Parameters
    ----------
    size : int
        The number of points to generate.
    distance_scale : `~astropy.units.Quantity`, optional
        A unit-ful factor to scale the random distances.
    rng : `numpy.random.Generator`, optional
        A random number generator instance.

    Returns
    -------
    rep : `~astropy.coordinates.SphericalRepresentation`
        The random points.
    """
    if rng is None:
        rng = np.random.default_rng()

    usph = uniform_spherical_random_surface(size=size, rng=rng)

    r = np.cbrt(rng.uniform(size=size)) * distance_scale
    return SphericalRepresentation(
        usph.lon, usph.lat, r)
