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

class Astrometric(BaseCoordinateFrame):
    r"""
    A coordinate or frame in the Astrometric system. This frame
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
    the document :ref:`coordinates-astrometric`.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other
        keywords)
    origin_distance : `~astropy.units.Quantity`, optional, must be keyword
        The distance from the Sun to the Galactic center.
    origin_ra : `Angle`, optional, must be keyword
        The Right Ascension (RA) of the Galactic center in the ICRS frame.
    origin_dec : `Angle`, optional, must be keyword
        The Declination (Dec) of the Galactic center in the ICRS frame.
    z_sun : `~astropy.units.Quantity`, optional, must be keyword
        The distance from the Sun to the Galactic midplane.

    Examples
    --------
    To transform to the Astrometric frame with the default
    frame attributes, pass the uninstantiated class name to the
    ``transform_to()`` method of a coordinate frame or
    `~astropy.coordinates.SkyCoord` object::

        >>> import astropy.units as u
        >>> import astropy.coordinates as coord
        >>> c = coord.ICRS(ra=[158.3122, 24.5] * u.degree,
        ...                dec=[-17.3, 81.52] * u.degree,
        ...                distance=[11.5, 24.12] * u.kpc)
        >>> c.transform_to(coord.Astrometric) # doctest: +FLOAT_CMP
        <Astrometric Coordinate (origin_distance=8.3 kpc, origin_ra=266d24m18.36s, origin_dec=-28d56m10.23s, z_sun=27.0 pc, roll=0.0 deg): (x, y, z) in kpc
            [(-9.6083818980977, -9.400621883358546, 6.520560663896347),
             (-21.283023068029138, 18.763340128812384, 7.846938548636718)]>

    To specify a custom set of parameters, you have to include extra keyword
    arguments when initializing the Astrometric frame object::

        >>> c.transform_to(coord.Astrometric(origin_distance=8.1*u.kpc)) # doctest: +FLOAT_CMP
        <Astrometric Coordinate (origin_distance=8.1 kpc, origin_ra=266d24m18.36s, origin_dec=-28d56m10.23s, z_sun=27.0 pc, roll=0.0 deg): (x, y, z) in kpc
            [(-9.407859235565343, -9.400621883358546, 6.520665737962164),
             (-21.08239383088295, 18.763340128812384, 7.84798134569032)]>

    Similarly, transforming from the Astrometric frame to another coordinate frame::

        >>> c = coord.Astrometric(x=[-8.3, 4.5] * u.kpc,
        ...                          y=[0., 81.52] * u.kpc,
        ...                          z=[0.027, 24.12] * u.kpc)
        >>> c.transform_to(coord.ICRS) # doctest: +FLOAT_CMP
        <ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)
            [(86.22349058727241, 28.8389413808627, 4.391577882957292e-05),
             (289.6680265194508, 49.88763881149547, 85.96407345372828)]>

    Or, with custom specification of the Galactic center::

        >>> c = coord.Astrometric(x=[-8.0, 4.5] * u.kpc,
        ...                          y=[0., 81.52] * u.kpc,
        ...                          z=[21.0, 24120.0] * u.pc,
        ...                          z_sun=21 * u.pc, origin_distance=8. * u.kpc)
        >>> c.transform_to(coord.ICRS) # doctest: +FLOAT_CMP
        <ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)
            [(86.25852490164378, 28.85773187391088, 2.7562547481200286e-05),
             (289.77285254989323, 50.062904565432014, 85.92160096237191)]>

    """
    default_representation = CartesianRepresentation

    # TODO: these can all become QuantityFrameAttribute's once #3217 is merged
    
    origin_distance = QuantityFrameAttribute(default=0, units=u.kpc)
    origin_ra = QuantityFrameAttribute(default=0, units=u.degree)
    origin_dec = QuantitFrameAttribute(default=0, units=u.degree)


# ICRS to/from Astrometric ----------------------->
@frame_transform_graph.transform(FunctionTransform, ICRS, Astrometric)
def icrs_to_astrometric(icrs_coord, astrometric_frame):
    from ..representation import CartesianRepresentation
    from ..angles import rotation_matrix

    if isinstance(icrs_coord.data, UnitSphericalRepresentation):
        raise ConvertError("Transforming to a Astrometric frame requires "
                           "a 3D coordinate, e.g. (angle, angle, distance) or"
                           " (x, y, z).")

    xyz = icrs_coord.cartesian.xyz

    # define rotation matrix to align x(ICRS) with the vector to the Galactic center
    mat1 = rotation_matrix(-astrometric_frame.origin_dec, 'y')
    mat2 = rotation_matrix(astrometric_frame.origin_ra, 'z')
    R = mat1 * mat2

    # some reshape hacks to handle ND arrays
    orig_shape = xyz.shape
    xyz = R.dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)

    # translate by Sun-Galactic center distance along new x axis
    xyz[0] = xyz[0] - astrometric_frame.origin_distance

    # rotate about y' to account for tilt due to Sun's height above the plane
    z_d = (astrometric_frame.z_sun / astrometric_frame.origin_distance).decompose()
    R = rotation_matrix(-np.arcsin(z_d), 'y')
    xyz = R.dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)

    representation = CartesianRepresentation(xyz)
    return astrometric_frame.realize_frame(representation)

@frame_transform_graph.transform(FunctionTransform, Astrometric, ICRS)
def astrometric_to_icrs(astrometric_coord, icrs_frame):
    from ..representation import CartesianRepresentation
    from ..angles import rotation_matrix

    if isinstance(astrometric_coord.data, UnitSphericalRepresentation):
        raise ConvertError("Transforming from a Astrometric frame requires "
                           "a 3D coordinate, e.g. (angle, angle, distance) or"
                           " (x, y, z).")

    xyz = astrometric_coord.cartesian.xyz

    # rotate about y' to account for tilt due to Sun's height above the plane
    z_d = (astrometric_coord.z_sun / astrometric_coord.origin_distance).decompose()
    R = rotation_matrix(np.arcsin(z_d), 'y')

    # some reshape hacks to handle ND arrays
    orig_shape = xyz.shape
    xyz = R.dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)

    # translate by Sun-Galactic center distance along x axis
    xyz[0] = xyz[0] + astrometric_coord.origin_distance

    # define inverse rotation matrix that aligns x(ICRS) with the vector to the Galactic center
    mat1 = rotation_matrix(-astrometric_coord.origin_dec, 'y')
    mat2 = rotation_matrix(astrometric_coord.origin_ra, 'z')
    R = mat1 * mat2

    # rotate into ICRS frame
    xyz = np.linalg.inv(R).dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)

    representation = CartesianRepresentation(xyz)
    return icrs_frame.realize_frame(representation)
