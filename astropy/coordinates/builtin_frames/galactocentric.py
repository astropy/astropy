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
# This is not used directly, but accessed via `get_roll0`.  We define it here to
# prevent having to create new Angle objects every time `get_roll0` is called.
_ROLL0 = Angle(58.5986320306*u.degree)

class Galactocentric(BaseCoordinateFrame):
    r"""
    A coordinate or frame in the Galactocentric system. This frame
    requires specifying the Sun-Galactic center distance, and optionally
    the height of the Sun above the Galactic midplane.

    The position of the Sun is assumed to be on the x axis of the final,
    right-handed system. That is, the x axis points from the position of
    the Sun projected to the Galactic midplane to the Galactic center --
    roughly towards :math:`(l,b) = (0^\circ,0^\circ)`. For the default
    transformation (:math:`{\rm roll}=0^\circ`), the y axis points roughly
    towards Galactic longitude :math:`l=90^\circ`, and the z axis points
    roughly towards the North Galactic Pole (:math:`b=90^\circ`).

    The default position of the Galactic Center in ICRS coordinates is
    taken from Reid et al. 2004,
    http://adsabs.harvard.edu/abs/2004ApJ...616..872R.

    .. math::

        {\rm RA} = 17:45:37.224~{\rm hr}\\
        {\rm Dec} = -28:56:10.23~{\rm deg}

    The default distance to the Galactic Center is 8.3 kpc, e.g.,
    Gillessen et al. 2009,
    http://adsabs.harvard.edu/abs/2009ApJ...692.1075G.

    The default height of the Sun above the Galactic midplane is taken to
    be 27 pc, as measured by
    http://adsabs.harvard.edu/abs/2001ApJ...553..184C.

    For a more detailed look at the math behind this transformation, see
    the document :ref:`coordinates-galactocentric`.

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

    Examples
    --------
    To transform to the Galactocentric frame with the default
    frame attributes, pass the uninstantiated class name to the
    ``transform_to()`` method of a coordinate frame or
    `~astropy.coordinates.SkyCoord` object::

        >>> import astropy.units as u
        >>> import astropy.coordinates as coord
        >>> c = coord.ICRS(ra=[158.3122, 24.5] * u.degree,
        ...                dec=[-17.3, 81.52] * u.degree,
        ...                distance=[11.5, 24.12] * u.kpc)
        >>> c.transform_to(coord.Galactocentric) # doctest: +FLOAT_CMP
        <Galactocentric Coordinate (galcen_distance=8.3 kpc, galcen_ra=266d24m18.36s, galcen_dec=-28d56m10.23s, z_sun=27.0 pc, roll=0.0 deg): (x, y, z) in kpc
            [(-9.6083818980977, -9.400621883358546, 6.520560663896347),
             (-21.283023068029138, 18.763340128812384, 7.846938548636718)]>

    To specify a custom set of parameters, you have to include extra keyword
    arguments when initializing the Galactocentric frame object::

        >>> c.transform_to(coord.Galactocentric(galcen_distance=8.1*u.kpc)) # doctest: +FLOAT_CMP
        <Galactocentric Coordinate (galcen_distance=8.1 kpc, galcen_ra=266d24m18.36s, galcen_dec=-28d56m10.23s, z_sun=27.0 pc, roll=0.0 deg): (x, y, z) in kpc
            [(-9.407859235565343, -9.400621883358546, 6.520665737962164),
             (-21.08239383088295, 18.763340128812384, 7.84798134569032)]>

    Similarly, transforming from the Galactocentric frame to another coordinate frame::

        >>> c = coord.Galactocentric(x=[-8.3, 4.5] * u.kpc,
        ...                          y=[0., 81.52] * u.kpc,
        ...                          z=[0.027, 24.12] * u.kpc)
        >>> c.transform_to(coord.ICRS) # doctest: +FLOAT_CMP
        <ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)
            [(86.22349058727241, 28.8389413808627, 4.391577882957292e-05),
             (289.6680265194508, 49.88763881149547, 85.96407345372828)]>

    Or, with custom specification of the Galactic center::

        >>> c = coord.Galactocentric(x=[-8.0, 4.5] * u.kpc,
        ...                          y=[0., 81.52] * u.kpc,
        ...                          z=[21.0, 24120.0] * u.pc,
        ...                          z_sun=21 * u.pc, galcen_distance=8. * u.kpc)
        >>> c.transform_to(coord.ICRS) # doctest: +FLOAT_CMP
        <ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)
            [(86.25852490164378, 28.85773187391088, 2.7562547481200286e-05),
             (289.77285254989323, 50.062904565432014, 85.92160096237191)]>

    """
    default_representation = CartesianRepresentation

    # TODO: these can all become QuantityFrameAttribute's once #3217 is merged
    galcen_distance = FrameAttribute(default=8.3*u.kpc)
    galcen_ra = FrameAttribute(default=Angle(266.4051*u.degree))
    galcen_dec = FrameAttribute(default=Angle(-28.936175*u.degree))
    z_sun = FrameAttribute(default=27.*u.pc)
    roll = FrameAttribute(default=0.*u.deg)

    @classmethod
    def get_roll0(cls):
        """
        The additional roll angle (about the final x axis) necessary to align
        the final z axis to match the Galactic yz-plane.  Setting the ``roll``
        frame attribute to  -this method's return value removes this rotation,
        allowing the use of the `Galactocentric` frame in more general contexts.
        """
        # not that the actual value is defined at the module level.  We make at
        # a property here because this module isn't actually part of the public
        # API, so it's better for it to be accessable from Galactocentric
        return _ROLL0

# ICRS to/from Galactocentric ----------------------->
@frame_transform_graph.transform(FunctionTransform, ICRS, Galactocentric)
def icrs_to_galactocentric(icrs_coord, galactocentric_frame):
    from ..representation import CartesianRepresentation
    from ..angles import rotation_matrix

    if isinstance(icrs_coord.data, UnitSphericalRepresentation):
        raise ConvertError("Transforming to a Galactocentric frame requires "
                           "a 3D coordinate, e.g. (angle, angle, distance) or"
                           " (x, y, z).")

    xyz = icrs_coord.cartesian.xyz

    # define rotation matrix to align x(ICRS) with the vector to the Galactic center
    mat1 = rotation_matrix(-galactocentric_frame.galcen_dec, 'y')
    mat2 = rotation_matrix(galactocentric_frame.galcen_ra, 'z')
    R1 = mat1 * mat2

    # extra roll away from the Galactic x-z plane
    R2 = rotation_matrix(galactocentric_frame.get_roll0() - galactocentric_frame.roll, 'x')

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
    R2 = rotation_matrix(galactocentric_coord.get_roll0() - galactocentric_coord.roll, 'x')

    # construct transformation matrix
    R = R2*R1

    # rotate into ICRS frame
    xyz = np.linalg.inv(R).dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)

    representation = CartesianRepresentation(xyz)
    return icrs_frame.realize_frame(representation)
