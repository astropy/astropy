# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to "observed" systems from ITRS.

This module provides a direct approach to ITRS<->Observed transformations that stays
entirely within the ITRS frame. Unlike the CIRS-based transformations, these treat
ITRS coordinates as time-invariant and perform simple rotations based on the
observer's location.

This approach is particularly useful for observing nearby objects (satellites,
airplanes, etc.) where staying within ITRS avoids issues with geocentric vs
topocentric aberration.
"""

import numpy as np

from astropy import units as u
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import FunctionTransformWithFiniteDifference
from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_transpose

from .altaz import AltAz
from .hadec import HADec
from .itrs import ITRS
from .utils import PIOVER2


def itrs_to_observed_mat(observed_frame):
    """
    Compute the transformation matrix from ITRS to observed frame (AltAz or HADec).

    Parameters
    ----------
    observed_frame : AltAz or HADec
        The observed frame with location information

    Returns
    -------
    mat : ndarray
        The 3x3 rotation matrix for transforming from ITRS to the observed frame
    """
    lon, lat, height = observed_frame.location.to_geodetic('WGS84')
    elong = lon.to_value(u.radian)

    if isinstance(observed_frame, AltAz):
        # Form ITRS to AltAz matrix
        elat = lat.to_value(u.radian)
        # AltAz frame is left-handed
        minus_x = np.eye(3)
        minus_x[0][0] = -1.0
        mat = (minus_x
               @ rotation_matrix(PIOVER2 - elat, 'y', unit=u.radian)
               @ rotation_matrix(elong, 'z', unit=u.radian))
    else:
        # Form ITRS to HADec matrix
        # HADec frame is left-handed
        minus_y = np.eye(3)
        minus_y[1][1] = -1.0
        mat = (minus_y
               @ rotation_matrix(elong, 'z', unit=u.radian))

    return mat


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, AltAz)
@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, HADec)
def itrs_to_observed(itrs_coo, observed_frame):
    """
    Transform from ITRS to observed frame (AltAz or HADec).

    This transformation treats ITRS coordinates as time-invariant. The coordinates
    are simply rotated based on the observer's location to produce topocentric
    coordinates. This is appropriate for nearby objects where the ITRS coordinates
    are tied to Earth's surface.

    Parameters
    ----------
    itrs_coo : ITRS
        The ITRS coordinate to transform
    observed_frame : AltAz or HADec
        The target observed frame

    Returns
    -------
    observed : AltAz or HADec
        The coordinate in the observed frame
    """
    # Trying to synchronize the obstimes here makes no sense. In fact,
    # it's a real gotcha as doing an ITRS->ITRS transform references
    # ITRS coordinates, which should be tied to the Earth, to the SSB.
    # Instead, we treat ITRS coordinates as time-invariant here.

    # Form the topocentric ITRS position
    topocentric_itrs_repr = (itrs_coo.cartesian
                             - observed_frame.location.get_itrs().cartesian)

    # Transform to observed frame
    rep = topocentric_itrs_repr.transform(itrs_to_observed_mat(observed_frame))

    return observed_frame.realize_frame(rep)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, AltAz, ITRS)
@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, HADec, ITRS)
def observed_to_itrs(observed_coo, itrs_frame):
    """
    Transform from observed frame (AltAz or HADec) to ITRS.

    This transformation treats ITRS coordinates as time-invariant. The observed
    coordinates are rotated back to ITRS and then the observer's location is added
    to produce geocentric ITRS coordinates.

    Parameters
    ----------
    observed_coo : AltAz or HADec
        The observed coordinate to transform
    itrs_frame : ITRS
        The target ITRS frame

    Returns
    -------
    itrs : ITRS
        The coordinate in ITRS frame
    """
    # Form the topocentric ITRS position
    topocentric_itrs_repr = observed_coo.cartesian.transform(
        matrix_transpose(itrs_to_observed_mat(observed_coo)))

    # Form the geocentric ITRS position
    rep = topocentric_itrs_repr + observed_coo.location.get_itrs().cartesian

    return itrs_frame.realize_frame(rep)


# Create loopback transformations for AltAz and HADec through ITRS
# This provides an alternative path for transformations between these frames
# that avoids the CIRS-based path and is more appropriate for nearby objects
frame_transform_graph._add_merged_transform(AltAz, ITRS, AltAz)
frame_transform_graph._add_merged_transform(HADec, ITRS, HADec)
