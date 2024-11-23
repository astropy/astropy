# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains utility functions for working with angles. These are both
used internally in ``astropy.coordinates.angles``, and of possible use externally.
"""

__all__ = [
    "angular_separation",
    "golden_spiral_grid",
    "offset_by",
    "position_angle",
    "uniform_spherical_random_surface",
    "uniform_spherical_random_volume",
]

# Third-party
import numpy as np

# Astropy
import astropy.units as u
from astropy.coordinates.representation import (
    SphericalRepresentation,
    UnitSphericalRepresentation,
)
from astropy.utils.compat import COPY_IF_NEEDED

_TWOPI = 2 * np.pi


def angular_separation(lon1, lat1, lon2, lat2):
    """
    Angular separation between two points on a sphere.

    Parameters
    ----------
    lon1, lat1, lon2, lat2 : `~astropy.coordinates.Angle`, `~astropy.units.Quantity` or float
        Longitude and latitude of the two points. Quantities should be in
        angular units; floats in radians.

    Returns
    -------
    angular separation : `~astropy.units.Quantity` ['angle'] or float
        Type depends on input; ``Quantity`` in angular units, or float in
        radians.

    Notes
    -----
    The angular separation is calculated using the Vincenty formula [1]_,
    which is slightly more complex and computationally expensive than
    some alternatives, but is stable at at all distances, including the
    poles and antipodes.

    .. [1] https://en.wikipedia.org/wiki/Great-circle_distance
    """
    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.arctan2(np.hypot(num1, num2), denominator)


def position_angle(lon1, lat1, lon2, lat2):
    """
    Position Angle (East of North) between two points on a sphere.

    Parameters
    ----------
    lon1, lat1, lon2, lat2 : `~astropy.coordinates.Angle`, `~astropy.units.Quantity` or float
        Longitude and latitude of the two points. Quantities should be in
        angular units; floats in radians.

    Returns
    -------
    pa : `~astropy.coordinates.Angle`
        The (positive) position angle of the vector pointing from position 1 to
        position 2.  If any of the angles are arrays, this will contain an array
        following the appropriate `numpy` broadcasting rules.

    """
    from .core import Angle

    deltalon = lon2 - lon1
    colat = np.cos(lat2)

    x = np.sin(lat2) * np.cos(lat1) - colat * np.sin(lat1) * np.cos(deltalon)
    y = np.sin(deltalon) * colat

    return Angle(np.arctan2(y, x), u.radian).wrap_at(360 * u.deg)


def offset_by(lon, lat, posang, distance):
    """
    Point with the given offset from the given point.

    Parameters
    ----------
    lon, lat, posang, distance : `~astropy.coordinates.Angle`, `~astropy.units.Quantity` or float
        Longitude and latitude of the starting point,
        position angle and distance to the final point.
        Quantities should be in angular units; floats in radians.
        Polar points at lat= +/-90 are treated as limit of +/-(90-epsilon) and same lon.

    Returns
    -------
    lon, lat : `~astropy.coordinates.Angle`
        The position of the final point.  If any of the angles are arrays,
        these will contain arrays following the appropriate `numpy` broadcasting rules.
        0 <= lon < 2pi.
    """
    from .core import Angle

    # Calculations are done using the spherical trigonometry sine and cosine rules
    # of the triangle A at North Pole,   B at starting point,   C at final point
    # with angles     A (change in lon), B (posang),            C (not used, but negative reciprocal posang)
    # with sides      a (distance),      b (final co-latitude), c (starting colatitude)
    # B, a, c are knowns; A and b are unknowns
    # https://en.wikipedia.org/wiki/Spherical_trigonometry

    cos_a = np.cos(distance)
    sin_a = np.sin(distance)
    cos_c = np.sin(lat)
    sin_c = np.cos(lat)
    cos_B = np.cos(posang)
    sin_B = np.sin(posang)

    # cosine rule: Know two sides: a,c and included angle: B; get unknown side b
    cos_b = cos_c * cos_a + sin_c * sin_a * cos_B
    # sin_b = np.sqrt(1 - cos_b**2)
    # sine rule and cosine rule for A (using both lets arctan2 pick quadrant).
    # multiplying both sin_A and cos_A by x=sin_b * sin_c prevents /0 errors
    # at poles.  Correct for the x=0 multiplication a few lines down.
    # sin_A/sin_a == sin_B/sin_b    # Sine rule
    xsin_A = sin_a * sin_B * sin_c
    # cos_a == cos_b * cos_c + sin_b * sin_c * cos_A  # cosine rule
    xcos_A = cos_a - cos_b * cos_c

    A = Angle(np.arctan2(xsin_A, xcos_A), u.radian)
    # Treat the poles as if they are infinitesimally far from pole but at given lon
    small_sin_c = sin_c < 1e-12
    if small_sin_c.any():
        # For south pole (cos_c = -1), A = posang; for North pole, A=180 deg - posang
        A_pole = (90 * u.deg + cos_c * (90 * u.deg - Angle(posang, u.radian))).to(u.rad)
        if A.shape:
            # broadcast to ensure the shape is like that of A, which is also
            # affected by the (possible) shapes of lat, posang, and distance.
            small_sin_c = np.broadcast_to(small_sin_c, A.shape)
            A[small_sin_c] = A_pole[small_sin_c]
        else:
            A = A_pole

    outlon = (Angle(lon, u.radian) + A).wrap_at(360.0 * u.deg).to(u.deg)
    outlat = Angle(np.arcsin(cos_b), u.radian).to(u.deg)

    return outlon, outlat


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
    lon = _TWOPI / golden_r * grid * u.rad
    lat = np.arcsin(1 - 2 * grid / size) * u.rad

    return UnitSphericalRepresentation(lon, lat)


def uniform_spherical_random_surface(size=1):
    """Generate a random sampling of points on the surface of the unit sphere.

    Parameters
    ----------
    size : int
        The number of points to generate.

    Returns
    -------
    rep : `~astropy.coordinates.UnitSphericalRepresentation`
        The random points.
    """
    rng = np.random  # can maybe switch to this being an input later - see #11628

    lon = rng.uniform(0, _TWOPI, size) * u.rad
    lat = np.arcsin(rng.uniform(-1, 1, size=size)) * u.rad

    return UnitSphericalRepresentation(lon, lat)


def uniform_spherical_random_volume(size=1, max_radius=1):
    """Generate a random sampling of points that follow a uniform volume
    density distribution within a sphere.

    Parameters
    ----------
    size : int
        The number of points to generate.
    max_radius : number, quantity-like, optional
        A dimensionless or unit-ful factor to scale the random distances.

    Returns
    -------
    rep : `~astropy.coordinates.SphericalRepresentation`
        The random points.
    """
    rng = np.random  # can maybe switch to this being an input later - see #11628

    usph = uniform_spherical_random_surface(size=size)

    r = np.cbrt(rng.uniform(size=size)) * u.Quantity(max_radius, copy=COPY_IF_NEEDED)
    return SphericalRepresentation(usph.lon, usph.lat, r)
