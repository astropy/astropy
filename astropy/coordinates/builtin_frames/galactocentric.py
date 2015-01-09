# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import numpy as np

from ... import units as u
from ..angles import Angle
from ..representation import CartesianRepresentation, UnitSphericalRepresentation
from ..baseframe import (BaseCoordinateFrame, FrameAttribute,
                         frame_transform_graph)
from ..transformations import FunctionTransform
from ..errors import ConvertError

from .icrs import ICRS

# Measured by minimizing the difference between a plane of coordinates along
#   l=0, b=[-90,90] and the Galactocentric x-z plane
ROLL0 = Angle(58.5986320306*u.degree)

class Galactocentric(BaseCoordinateFrame):
    """
    A coordinate or frame in the Galactocentric system. This frame
    requires specifying the Sun-Galactic center distance, and optionally
    the height of the Sun above the Galactic midplane.

    The position of the Sun is assumed to be on the x axis of the final,
    right-handed system. That is, the x axis points from the position of
    the Sun projected to the Galactic midplane to the Galactic center --
    roughly towards ``l,b = (0º,0º)``. For the default transformation
    (``roll=0º``), the y axis points roughly towards Galactic longitude
    ``l=90º``, and the z axis points roughly towards the North Galactic
    Pole (``b=90º``).

    The default position of the Galactic Center in ICRS coordinates is
    taken from Reid et al. 2004,
    http://adsabs.harvard.edu/abs/2004ApJ...616..872R.

        RA = 17:45:37.224 (hours)
        Dec = -28:56:10.23 (degrees)

    The default distance to the Galactic Center is 8.3 kpc, e.g.,
    Gillessen et al. 2009,
    http://adsabs.harvard.edu/abs/2009ApJ...692.1075G.

    The default height of the Sun above the Galactic midplane is taken to
    be 27 pc, as measured by
    http://adsabs.harvard.edu/abs/2001ApJ...553..184C.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other
        keywords)
    galcen_distance : `~astropy.units.Quantity`, optional, must be keyword
        The distance from the Sun to the Galactic center.
    galcen_ra : `Angle`, optional, must be keyword
        The Right Ascension (RA) of the Galactic center in the ICRS frame.
    galcen_dec : `Angle`, optional, must be keyword
        The Declination (Dec) of the Galactic center in the ICRS frame.
    z_sun : `~astropy.units.Quantity`, optional, must be keyword
        The distance from the Sun to the Galactic midplane.
    roll : `Angle`, optional, must be keyword
        The angle to rotate about the final x-axis, relative to the
        orientation for Galactic. For example, if this roll angle is 0,
        the final x-z plane will align with the Galactic coordinates x-z
        plane. Unless you really know what this means, you probably should
        not change this!

    """
    default_representation = CartesianRepresentation

    # TODO: these can all become QuantityFrameAttribute's once #3217 is merged
    galcen_distance = FrameAttribute(default=8.3*u.kpc)
    galcen_ra = FrameAttribute(default=Angle(266.4051*u.degree))
    galcen_dec = FrameAttribute(default=Angle(-28.936175*u.degree))
    z_sun = FrameAttribute(default=27.*u.pc)
    roll = FrameAttribute(default=0.*u.deg)

# ICRS to/from Galactocentric ----------------------->
@frame_transform_graph.transform(FunctionTransform, ICRS, Galactocentric)
def icrs_to_galactocentric(icrs_coord, galactocentric_frame):
    from ..representation import CartesianRepresentation
    from ..angles import rotation_matrix

    if isinstance(icrs_coord.data, UnitSphericalRepresentation):
        raise ConvertError("Transforming to a Galactocentric frame requires "
                           "a 3D coordinates, e.g. (angle, angle, distance) or"
                           " (x, y, z).")

    xyz = icrs_coord.cartesian.xyz

    # define rotation matrix to align x(ICRS) with the vector to the Galactic center
    mat1 = rotation_matrix(-galactocentric_frame.galcen_dec, 'y')
    mat2 = rotation_matrix(galactocentric_frame.galcen_ra, 'z')
    R1 = mat1 * mat2

    # extra roll away from the Galactic x-z plane
    R2 = rotation_matrix(ROLL0 - galactocentric_frame.roll, 'x')

    # construct transformation matrix
    R = R2*R1

    # some reshape hacks to handle ND arrays
    orig_shape = xyz.shape
    xyz = R.dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)

    # translate by Sun-Galactic center distance along new x axis
    xyz[0] = xyz[0] - galactocentric_frame.galcen_distance

    # rotate about y' to account for tilt due to Sun's height above the plane
    z_d = (galactocentric_frame.z_sun / galactocentric_frame.galcen_distance).decompose()
    R = rotation_matrix(-np.arcsin(z_d), 'y')
    xyz = R.dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)

    representation = CartesianRepresentation(xyz)
    return galactocentric_frame.realize_frame(representation)

@frame_transform_graph.transform(FunctionTransform, Galactocentric, ICRS)
def galactocentric_to_icrs(galactocentric_coord, icrs_frame):
    from ..representation import CartesianRepresentation
    from ..angles import rotation_matrix

    if isinstance(galactocentric_coord.data, UnitSphericalRepresentation):
        raise ConvertError("Transforming from a Galactocentric frame requires "
                           "a 3D coordinate, e.g. (angle, angle, distance) or"
                           " (x, y, z).")

    xyz = galactocentric_coord.cartesian.xyz

    # rotate about y' to account for tilt due to Sun's height above the plane
    z_d = (galactocentric_coord.z_sun / galactocentric_coord.galcen_distance).decompose()
    R = rotation_matrix(np.arcsin(z_d), 'y')

    # some reshape hacks to handle ND arrays
    orig_shape = xyz.shape
    xyz = R.dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)

    # translate by Sun-Galactic center distance along x axis
    xyz[0] = xyz[0] + galactocentric_coord.galcen_distance

    # define inverse rotation matrix that aligns x(ICRS) with the vector to the Galactic center
    mat1 = rotation_matrix(-galactocentric_coord.galcen_dec, 'y')
    mat2 = rotation_matrix(galactocentric_coord.galcen_ra, 'z')
    R1 = mat1 * mat2

    # extra roll away from the Galactic x-z plane
    R2 = rotation_matrix(ROLL0 - galactocentric_coord.roll, 'x')

    # construct transformation matrix
    R = R2*R1

    # rotate into ICRS frame
    xyz = np.linalg.inv(R).dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)

    representation = CartesianRepresentation(xyz)
    return icrs_frame.realize_frame(representation)
