# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to "observed" systems from
ITRS.  These are direct, topocentric transforms: the target ITRS position is
kept fixed and the observer's ITRS position is subtracted.  The resulting
horizontal/equatorial coordinates are purely geometric and do **not** include
refraction, aberration, or precession-nutation.  They are intended for nearby or
Earth-fixed objects (satellites, aircraft, ground features) rather than for
distant celestial sources.

.. todo::
    Refraction parity with CIRS/ICRS observed paths.  The CIRS and ICRS
    ``observed_transforms`` modules apply refraction when the ``pressure``
    attribute on the observed frame is non-zero.  This direct ITRS path
    currently ignores ``pressure``, ``temperature``, ``relative_humidity``,
    and ``obswl`` and always returns topocentric coordinates.
"""

import numpy as np

from astropy import units as u
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.matrix_utilities import matrix_product, matrix_transpose, rotation_matrix
from astropy.coordinates.representation import (
    CartesianRepresentation,
    SphericalRepresentation,
    UnitSphericalRepresentation,
)
from astropy.coordinates.transformations import FunctionTransformWithFiniteDifference

from .altaz import AltAz
from .hadec import HADec
from .itrs import ITRS


def _itrs_to_enu_mat(lon, lat):
    """Rotation from ITRS to local East-North-Up (ENU) basis.

    The returned matrix maps an ITRS column vector into the observed-frame
    Cartesian basis in which the x-axis points north (azimuth 0), the y-axis
    points east (azimuth 90 degrees), and the z-axis points up (altitude 90
    degrees).  The matrix rows are the geodetic north, east, and up unit
    vectors expressed in ITRS.
    """
    slon = np.sin(lon.to_value(u.radian))
    clon = np.cos(lon.to_value(u.radian))
    slat = np.sin(lat.to_value(u.radian))
    clat = np.cos(lat.to_value(u.radian))
    zero = np.zeros_like(slon)

    # Geodetic north, east, and up unit vectors in ITRS.
    mat = np.stack([
        np.stack([-clon*slat, -slon*slat, clat], axis=-1),
        np.stack([-slon, clon, zero], axis=-1),
        np.stack([clon*clat, slon*clat, slat], axis=-1),
    ], axis=-2)
    return mat


def _itrs_to_hadec_mat(lon, lat):
    """Rotation from ITRS to local Hour-Angle-Declination basis.

    The returned matrix maps an ITRS column vector into the observed-frame
    Cartesian basis used by `~astropy.coordinates.SphericalRepresentation`
    with longitude = hour angle and latitude = declination:

    * x-axis: hour angle 0, declination 0 (the south point on the local
      horizon, i.e. along the local meridian towards decreasing geodetic
      latitude).
    * y-axis: hour angle +6 h, declination 0 (west on the horizon).
    * z-axis: declination +90 deg (the north celestial pole / ITRS +z axis).

    This is a right-handed basis with hour angle increasing westwards.
    """
    slon = np.sin(lon.to_value(u.radian))
    clon = np.cos(lon.to_value(u.radian))
    zero = np.zeros_like(slon)
    one = np.ones_like(slon)

    # HA=0 points to the local meridian at declination 0 (south on horizon).
    # HA=+6 h points west.  Dec=+90 deg points to the ITRS +z axis.
    mat = np.stack([
        np.stack([-clon, -slon, zero], axis=-1),
        np.stack([slon, -clon, zero], axis=-1),
        np.stack([zero, zero, one], axis=-1),
    ], axis=-2)
    return mat


def itrs_to_observed_mat(observed_frame):
    """Build the ITRS -> observed-frame rotation matrix.

    Parameters
    ----------
    observed_frame : `~astropy.coordinates.AltAz` or `~astropy.coordinates.HADec`
        The target observed frame, carrying ``location`` (the observer) and
        ``obstime``.

    Returns
    -------
    mat : `~numpy.ndarray`
        A ``(3, 3)`` rotation matrix, or a broadcasted stack of shape
        ``(..., 3, 3)`` for array-valued frame attributes.
    """
    geodetic = observed_frame.location.to_geodetic('WGS84')
    lon = geodetic.lon
    lat = geodetic.lat

    if isinstance(observed_frame, AltAz):
        return _itrs_to_enu_mat(lon, lat)
    else:
        return _itrs_to_hadec_mat(lon, lat)


def _is_unitspherical(coo):
    """Return whether the coordinate data is unit-spherical (dimensionless)."""
    return (isinstance(coo.data, UnitSphericalRepresentation) or
            coo.cartesian.x.unit == u.one)


def itrs_to_observed(itrs_coo, observed_frame):
    """Transform an ITRS coordinate to an observed frame.

    The ITRS Cartesian vector is treated as time-invariant: it is **not**
    rotated to the observed frame's ``obstime``.  Only the observer's ITRS
    position is evaluated at ``observed_frame.obstime``.
    """
    is_unitspherical = _is_unitspherical(itrs_coo)

    # Rotate the ITRS vector into the local observed basis.  For unit-
    # spherical inputs the coordinate is a direction only; we follow the
    # convention used by icrs_observed_transforms and cirs_observed_transforms
    # and do not subtract the observer position (which would change the
    # direction).  For spherical/Cartesian inputs we form the true topocentric
    # vector (target minus observer).
    mat = itrs_to_observed_mat(observed_frame)
    if is_unitspherical:
        observed_cart = itrs_coo.cartesian.transform(mat)
        srepr = observed_cart.represent_as(UnitSphericalRepresentation)
    else:
        # Observer's ITRS position at the output obstime.
        obs_itrs = observed_frame.location.get_itrs(observed_frame.obstime).cartesian
        topocentric = itrs_coo.cartesian - obs_itrs
        observed_cart = topocentric.transform(mat)
        srepr = observed_cart.represent_as(SphericalRepresentation)

    return observed_frame.realize_frame(srepr)


def observed_to_itrs(observed_coo, itrs_frame):
    """Transform an observed coordinate back to ITRS.

    The inverse rotation is applied and the observer's ITRS position at
    ``observed_coo.obstime`` is added back.  The resulting Cartesian vector is
    realized in ``itrs_frame``; its numerical value does not depend on
    ``itrs_frame.obstime``.
    """
    is_unitspherical = _is_unitspherical(observed_coo)

    # Inverse rotation from observed basis to ITRS.
    mat = itrs_to_observed_mat(observed_coo)
    itrs_cart = observed_coo.cartesian.transform(matrix_transpose(mat))

    if is_unitspherical:
        # The output direction vector is dimensionless; no observer offset is
        # added, so the result remains a unit direction.
        rep = itrs_cart.represent_as(UnitSphericalRepresentation)
    else:
        # Add back the observer's ITRS position.
        obs_itrs = observed_coo.location.get_itrs(observed_coo.obstime).cartesian
        itrs_cart = itrs_cart + obs_itrs
        rep = itrs_cart

    return itrs_frame.realize_frame(rep)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, AltAz)
@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, HADec)
def _itrs_to_observed(itrs_coo, observed_frame):
    return itrs_to_observed(itrs_coo, observed_frame)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, AltAz, ITRS)
@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, HADec, ITRS)
def _observed_to_itrs(observed_coo, itrs_frame):
    return observed_to_itrs(observed_coo, itrs_frame)


# Note: loopback transformations through ITRS are intentionally omitted here.
# Self-transforms for AltAz and HADec are already registered through ICRS in
# icrs_observed_transforms.py, and the transform graph only permits one direct
# self-transform per frame pair.  Earth-fixed loopback semantics for
# AltAz/HADec therefore cannot be provided without first removing the ICRS
# loopbacks.
