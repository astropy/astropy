# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains the transformation functions for getting to "observed" systems from ITRS.

This module provides direct transformations between ITRS and observed frames
(AltAz, HADec) that stay entirely within the ITRS, avoiding the geocentric vs
topocentric aberration issues that arise when transforming through ICRS/GCRS/CIRS.

These transformations treat ITRS coordinates as time-invariant, which is appropriate
for nearby objects like satellites, airplanes, mountains, and buildings. The obstime
of the output frame is simply adopted without attempting to transform the ITRS
coordinates themselves across different times (which would incorrectly reference
them to the solar system barycenter rather than keeping them tied to the Earth).
"""

import numpy as np

from astropy import units as u
from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_transpose
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import FunctionTransformWithFiniteDifference

from .altaz import AltAz
from .hadec import HADec
from .itrs import ITRS
from .utils import PIOVER2


def itrs_to_observed_mat(observed_frame):
    """
    Create the rotation matrix for transforming from ITRS to an observed frame.

    This is a purely geometric transformation that converts from geocentric
    Cartesian ITRS coordinates to either topocentric AltAz or HADec coordinates.

    Parameters
    ----------
    observed_frame : AltAz or HADec
        The observed frame to transform to, which provides the observer location.

    Returns
    -------
    mat : numpy.ndarray
        A 3x3 rotation matrix for the transformation.
    """
    lon, lat, height = observed_frame.location.to_geodetic('WGS84')
    elong = lon.to_value(u.radian)

    if isinstance(observed_frame, AltAz):
        # Form ITRS to AltAz matrix
        # AltAz frame is left-handed (azimuth increases eastward)
        elat = lat.to_value(u.radian)
        # Create left-handed coordinate system by negating x
        minus_x = np.eye(3)
        minus_x[0][0] = -1.0
        mat = (minus_x
               @ rotation_matrix(PIOVER2 - elat, 'y', unit=u.radian)
               @ rotation_matrix(elong, 'z', unit=u.radian))
    else:
        # Form ITRS to HADec matrix
        # HADec frame is left-handed (hour angle increases westward)
        # Create left-handed coordinate system by negating y
        minus_y = np.eye(3)
        minus_y[1][1] = -1.0
        mat = (minus_y
               @ rotation_matrix(elong, 'z', unit=u.radian))

    return mat


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, AltAz)
@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, HADec)
def itrs_to_observed(itrs_coo, observed_frame):
    """
    Transform ITRS coordinates to observed coordinates (AltAz or HADec).

    This transformation stays entirely within the ITRS and treats ITRS positions
    as time-invariant. This is appropriate for nearby objects that are fixed to
    or moving near the Earth's surface (satellites, aircraft, ground features).

    Trying to synchronize obstimes here makes no sense - doing an ITRS->ITRS
    transform for differing obstimes would incorrectly reference ITRS coordinates
    (which should be tied to the Earth) to the solar system barycenter. Instead,
    we treat ITRS coordinates as time-invariant and simply adopt the obstime of
    the output frame.

    Parameters
    ----------
    itrs_coo : ITRS
        The input ITRS coordinates.
    observed_frame : AltAz or HADec
        The observed frame to transform to.

    Returns
    -------
    observed : AltAz or HADec
        The transformed coordinates in the observed frame.
    """
    # Form the topocentric ITRS position
    # Subtract the observer's position to get topocentric coordinates
    topocentric_itrs_repr = (itrs_coo.cartesian
                             - observed_frame.location.get_itrs().cartesian)

    # Apply the rotation matrix to convert to observed frame
    rep = topocentric_itrs_repr.transform(itrs_to_observed_mat(observed_frame))

    return observed_frame.realize_frame(rep)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, AltAz, ITRS)
@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, HADec, ITRS)
def observed_to_itrs(observed_coo, itrs_frame):
    """
    Transform observed coordinates (AltAz or HADec) to ITRS coordinates.

    This transformation stays entirely within the ITRS and treats ITRS positions
    as time-invariant. The reverse of itrs_to_observed.

    Parameters
    ----------
    observed_coo : AltAz or HADec
        The input observed coordinates.
    itrs_frame : ITRS
        The ITRS frame to transform to.

    Returns
    -------
    itrs : ITRS
        The transformed coordinates in ITRS.
    """
    # Form the topocentric ITRS position
    # Apply inverse rotation (transpose) to convert from observed to ITRS
    topocentric_itrs_repr = observed_coo.cartesian.transform(
        matrix_transpose(itrs_to_observed_mat(observed_coo)))

    # Form the geocentric ITRS position
    # Add the observer's position to convert from topocentric to geocentric
    rep = topocentric_itrs_repr + observed_coo.location.get_itrs().cartesian

    return itrs_frame.realize_frame(rep)
