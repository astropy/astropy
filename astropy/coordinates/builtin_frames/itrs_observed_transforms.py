# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to/from observed frames
(AltAz, HADec) directly from ITRS without going through CIRS. This avoids
applying stellar aberration corrections that are inappropriate for nearby
(terrestrial/near-Earth) objects in ITRS.

The transformation is purely geometric: it converts between Earth-Centered
Earth-Fixed (ECEF) Cartesian coordinates and local horizontal / equatorial
coordinates using only rotation matrices based on the observer's geodetic
longitude and latitude.  No ERA, polar motion, or sidereal time is involved,
so the result is independent of ``obstime``.
"""
import numpy as np

from astropy import units as u
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import FunctionTransformWithFiniteDifference
from astropy.coordinates.representation import CartesianRepresentation
from astropy.coordinates.matrix_utilities import matrix_transpose, rotation_matrix

from .altaz import AltAz
from .hadec import HADec
from .itrs import ITRS

# Latitude of the north pole (used for rotation_matrix calls).
NORTH_POLE = 90.0 * u.deg


def itrs_to_altaz_mat(lon, lat):
    """Compute the rotation matrix from ITRS Cartesian to AltAz Cartesian.

    The AltAz Cartesian frame is left-handed (azimuth measured clockwise).
    Its axes correspond to (North, East, Up) in the SphericalRepresentation
    convention where lon=az (0=North, 90°=East) and lat=alt.

    Parameters
    ----------
    lon : `~astropy.units.Quantity`
        Observer's geodetic longitude.
    lat : `~astropy.units.Quantity`
        Observer's geodetic latitude.

    Returns
    -------
    mat : `numpy.ndarray`, shape (3, 3)
        Rotation matrix from ITRS to AltAz.
    """
    minus_x = np.eye(3)
    minus_x[0][0] = -1.0
    return minus_x @ rotation_matrix(NORTH_POLE - lat, 'y') @ rotation_matrix(lon, 'z')


def itrs_to_hadec_mat(lon):
    """Compute the rotation matrix from ITRS Cartesian to HADec Cartesian.

    The HADec frame is left-handed (hour angle increases Westward).

    Parameters
    ----------
    lon : `~astropy.units.Quantity`
        Observer's geodetic longitude.

    Returns
    -------
    mat : `numpy.ndarray`, shape (3, 3)
        Rotation matrix from ITRS to HADec.
    """
    minus_y = np.eye(3)
    minus_y[1][1] = -1.0
    return minus_y @ rotation_matrix(lon, 'z')


def altaz_to_hadec_mat(lat):
    """Compute the rotation matrix from AltAz Cartesian to HADec Cartesian.

    Parameters
    ----------
    lat : `~astropy.units.Quantity`
        Observer's geodetic latitude.

    Returns
    -------
    mat : `numpy.ndarray`, shape (3, 3)
        Rotation matrix from AltAz to HADec.
    """
    z180 = np.eye(3)
    z180[0][0] = -1.0
    z180[1][1] = -1.0
    return z180 @ rotation_matrix(NORTH_POLE - lat, 'y')


def itrs_to_observed_mat(observed_frame):
    """Return the ITRS→observed rotation matrix for the given observed frame.

    Parameters
    ----------
    observed_frame : `~astropy.coordinates.AltAz` or `~astropy.coordinates.HADec`
        The target observed frame (must have a ``location`` attribute).

    Returns
    -------
    mat : `numpy.ndarray`, shape (3, 3)
        Rotation matrix from ITRS Cartesian to the observed frame's Cartesian.
    """
    lon, lat, height = observed_frame.location.to_geodetic('WGS84')
    if isinstance(observed_frame, AltAz):
        return itrs_to_altaz_mat(lon, lat)
    else:
        return itrs_to_hadec_mat(lon)


def _get_cart(location):
    """Return the ECEF CartesianRepresentation for an EarthLocation."""
    return CartesianRepresentation(location.x, location.y, location.z)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, AltAz)
@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, HADec)
def itrs_to_observed(itrs_coo, observed_frame):
    """Transform from ITRS to AltAz or HADec.

    This is a purely geometric, time-invariant transform.  No CIRS rotation,
    no ERA, and no stellar aberration is applied.  This makes it correct for
    near-Earth objects and terrestrial objects whose positions are naturally
    described in Earth-fixed (ITRS) coordinates.

    The ITRS ``location`` attribute determines whether the coordinate data
    is geocentric or topocentric:

    * **Geocentric ITRS** (``location`` == geocenter, the default):
      the absolute ECEF position is used; the observer's ECEF position is
      subtracted to form the topocentric vector before rotating.
    * **Topocentric ITRS** (``location`` == observer location):
      the data already represents the displacement from the observer;
      only the rotation is applied.
    """
    pmat = itrs_to_observed_mat(observed_frame)

    obs_loc = observed_frame.location
    itrs_loc = itrs_coo.location

    if np.any(itrs_loc != obs_loc):
        # Convert from whatever reference location the ITRS uses to the
        # observer's topocentric frame.  This is purely geometric — no CIRS.
        # absolute ECEF = data + ITRS_location_offset
        # topocentric   = absolute ECEF - observer_ECEF
        loc_cart = _get_cart(itrs_loc)
        obs_cart = _get_cart(obs_loc)
        topo_cart = itrs_coo.cartesian + loc_cart - obs_cart
        crepr = topo_cart.transform(pmat)
    else:
        # Locations match: data is already the topocentric vector.
        # Pure rotation — time-invariant.
        crepr = itrs_coo.cartesian.transform(pmat)

    return observed_frame.realize_frame(crepr)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, AltAz, ITRS)
@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, HADec, ITRS)
def observed_to_itrs(observed_coo, itrs_frame):
    """Transform from AltAz or HADec to ITRS.

    This is the inverse of ``itrs_to_observed``.  The result is expressed as
    a displacement relative to the target ITRS frame's ``location``.
    """
    pmat = matrix_transpose(itrs_to_observed_mat(observed_coo))
    # Topocentric vector relative to observer
    topo_cart = observed_coo.cartesian.transform(pmat)

    obs_loc = observed_coo.location
    target_loc = itrs_frame.location

    if np.any(target_loc != obs_loc):
        # Re-express the topocentric vector in the target ITRS location frame.
        # displacement relative to target = (absolute ECEF) - target_ECEF
        # = (topo_cart + obs_ECEF) - target_ECEF
        obs_cart = _get_cart(obs_loc)
        target_cart = _get_cart(target_loc)
        crepr = topo_cart + obs_cart - target_cart
    else:
        crepr = topo_cart

    return itrs_frame.realize_frame(crepr)
