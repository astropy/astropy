# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import warnings

import numpy as np

from ... import units as u
from ...utils.exceptions import AstropyDeprecationWarning
from ..angles import Angle
from ..matrix_utilities import rotation_matrix, matrix_product, matrix_transpose
from .. import representation as r
from ..baseframe import (BaseCoordinateFrame, frame_transform_graph,
                         RepresentationMapping)
from ..attributes import (Attribute, CoordinateAttribute,
                          QuantityAttribute,
                          DifferentialAttribute)
from ..transformations import AffineTransform
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
    Gillessen et al. (2009),
    https://ui.adsabs.harvard.edu/#abs/2009ApJ...692.1075G/abstract

    The default height of the Sun above the Galactic midplane is taken to
    be 27 pc, as measured by Chen et al. (2001),
    https://ui.adsabs.harvard.edu/#abs/2001ApJ...553..184C/abstract

    The default solar motion relative to the Galactic center is taken from a
    combination of Schönrich et al. (2010) [for the peculiar velocity] and
    Bovy (2015) [for the circular velocity at the solar radius],
    https://ui.adsabs.harvard.edu/#abs/2010MNRAS.403.1829S/abstract
    https://ui.adsabs.harvard.edu/#abs/2015ApJS..216...29B/abstract

    For a more detailed look at the math behind this transformation, see
    the document :ref:`coordinates-galactocentric`.

    The frame attributes are listed under **Other Parameters**.

    Parameters
    ----------
    representation : `~astropy.coordinates.representation.BaseRepresentation` or None
        A representation object or None to have no data (or use the other
        keywords)

    x : `~astropy.units.Quantity`, optional
        Cartesian, Galactocentric :math:`x` position component.
    y : `~astropy.units.Quantity`, optional
        Cartesian, Galactocentric :math:`y` position component.
    z : `~astropy.units.Quantity`, optional
        Cartesian, Galactocentric :math:`z` position component.

    v_x : `~astropy.units.Quantity`, optional
        Cartesian, Galactocentric :math:`v_x` velocity component.
    v_y : `~astropy.units.Quantity`, optional
        Cartesian, Galactocentric :math:`v_y` velocity component.
    v_z : `~astropy.units.Quantity`, optional
        Cartesian, Galactocentric :math:`v_z` velocity component.

    copy : bool, optional
        If `True` (default), make copies of the input coordinate arrays.
        Can only be passed in as a keyword argument.

    differential_cls : `BaseDifferential`, dict, optional
        A differential class or dictionary of differential classes (currently
        only a velocity differential with key 's' is supported). This sets
        the expected input differential class, thereby changing the expected
        keyword arguments of the data passed in. For example, passing
        ``differential_cls=CartesianDifferential`` will make the classes
        expect velocity data with the argument names ``v_x, v_y, v_z``.

    Other parameters
    ----------------
    galcen_coord : `ICRS`, optional, must be keyword
        The ICRS coordinates of the Galactic center.
    galcen_distance : `~astropy.units.Quantity`, optional, must be keyword
        The distance from the sun to the Galactic center.
    galcen_v_sun : `~astropy.coordinates.representation.CartesianDifferential`, optional, must be keyword
        The velocity of the sun *in the Galactocentric frame* as Cartesian
        velocity components.
    z_sun : `~astropy.units.Quantity`, optional, must be keyword
        The distance from the sun to the Galactic midplane.
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
        <Galactocentric Coordinate (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
            ( 266.4051, -28.936175)>, galcen_distance=8.3 kpc, galcen_v_sun=( 11.1,  232.24,  7.25) km / s, z_sun=27.0 pc, roll=0.0 deg): (x, y, z) in kpc
            [( -9.6083819 ,  -9.40062188,  6.52056066),
             (-21.28302307,  18.76334013,  7.84693855)]>

    To specify a custom set of parameters, you have to include extra keyword
    arguments when initializing the Galactocentric frame object::

        >>> c.transform_to(coord.Galactocentric(galcen_distance=8.1*u.kpc)) # doctest: +FLOAT_CMP
        <Galactocentric Coordinate (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
            ( 266.4051, -28.936175)>, galcen_distance=8.1 kpc, galcen_v_sun=( 11.1,  232.24,  7.25) km / s, z_sun=27.0 pc, roll=0.0 deg): (x, y, z) in kpc
            [( -9.40785924,  -9.40062188,  6.52066574),
             (-21.08239383,  18.76334013,  7.84798135)]>

    Similarly, transforming from the Galactocentric frame to another coordinate frame::

        >>> c = coord.Galactocentric(x=[-8.3, 4.5] * u.kpc,
        ...                          y=[0., 81.52] * u.kpc,
        ...                          z=[0.027, 24.12] * u.kpc)
        >>> c.transform_to(coord.ICRS) # doctest: +FLOAT_CMP
        <ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)
            [(  86.22349059, 28.83894138,  4.39157788e-05),
             ( 289.66802652, 49.88763881,  8.59640735e+01)]>

    Or, with custom specification of the Galactic center::

        >>> c = coord.Galactocentric(x=[-8.0, 4.5] * u.kpc,
        ...                          y=[0., 81.52] * u.kpc,
        ...                          z=[21.0, 24120.0] * u.pc,
        ...                          z_sun=21 * u.pc, galcen_distance=8. * u.kpc)
        >>> c.transform_to(coord.ICRS) # doctest: +FLOAT_CMP
        <ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)
            [(  86.2585249 ,  28.85773187,  2.75625475e-05),
             ( 289.77285255,  50.06290457,  8.59216010e+01)]>

    """
    frame_specific_representation_info = {
        r.CartesianDifferential: [
            RepresentationMapping('d_x', 'v_x', u.km/u.s),
            RepresentationMapping('d_y', 'v_y', u.km/u.s),
            RepresentationMapping('d_z', 'v_z', u.km/u.s),
        ],
    }

    default_representation = r.CartesianRepresentation
    default_differential = r.CartesianDifferential

    # frame attributes
    galcen_coord = CoordinateAttribute(default=ICRS(ra=266.4051*u.degree,
                                                    dec=-28.936175*u.degree),
                                       frame=ICRS)
    galcen_distance = QuantityAttribute(default=8.3*u.kpc)

    galcen_v_sun = DifferentialAttribute(
        default=r.CartesianDifferential([11.1, 220+12.24, 7.25] * u.km/u.s),
        allowed_classes=[r.CartesianDifferential])

    z_sun = QuantityAttribute(default=27.*u.pc)
    roll = QuantityAttribute(default=0.*u.deg)

    def __init__(self, *args, **kwargs):

        # backwards-compatibility
        if ('galcen_ra' in kwargs or 'galcen_dec' in kwargs):
            warnings.warn("The arguments 'galcen_ra', and 'galcen_dec' are "
                          "deprecated in favor of specifying the sky coordinate"
                          " as a CoordinateAttribute using the 'galcen_coord' "
                          "argument", AstropyDeprecationWarning)

            galcen_kw = dict()
            galcen_kw['ra'] = kwargs.pop('galcen_ra', self.galcen_coord.ra)
            galcen_kw['dec'] = kwargs.pop('galcen_dec', self.galcen_coord.dec)
            kwargs['galcen_coord'] = ICRS(**galcen_kw)

        super(Galactocentric, self).__init__(*args, **kwargs)

    @property
    def galcen_ra(self):
        warnings.warn("The attribute 'galcen_ra' is deprecated. Use "
                      "'.galcen_coord.ra' instead.", AstropyDeprecationWarning)
        return self.galcen_coord.ra

    @property
    def galcen_dec(self):
        warnings.warn("The attribute 'galcen_dec' is deprecated. Use "
                      "'.galcen_coord.dec' instead.", AstropyDeprecationWarning)
        return self.galcen_coord.dec

    @classmethod
    def get_roll0(cls):
        """
        The additional roll angle (about the final x axis) necessary to align
        the final z axis to match the Galactic yz-plane.  Setting the ``roll``
        frame attribute to  -this method's return value removes this rotation,
        allowing the use of the `Galactocentric` frame in more general contexts.
        """
        # note that the actual value is defined at the module level.  We make at
        # a property here because this module isn't actually part of the public
        # API, so it's better for it to be accessable from Galactocentric
        return _ROLL0

# ICRS to/from Galactocentric ----------------------->


def get_matrix_vectors(galactocentric_frame, inverse=False):
    """
    Use the ``inverse`` argument to get the inverse transformation, matrix and
    offsets to go from Galactocentric to ICRS.
    """
    # shorthand
    gcf = galactocentric_frame

    # rotation matrix to align x(ICRS) with the vector to the Galactic center
    mat1 = rotation_matrix(-gcf.galcen_coord.dec, 'y')
    mat2 = rotation_matrix(gcf.galcen_coord.ra, 'z')
    # extra roll away from the Galactic x-z plane
    mat0 = rotation_matrix(gcf.get_roll0() - gcf.roll, 'x')

    # construct transformation matrix and use it
    R = matrix_product(mat0, mat1, mat2)

    # Now need to translate by Sun-Galactic center distance around x' and
    # rotate about y' to account for tilt due to Sun's height above the plane
    translation = r.CartesianRepresentation(gcf.galcen_distance * [1., 0., 0.])
    z_d = gcf.z_sun / gcf.galcen_distance
    H = rotation_matrix(-np.arcsin(z_d), 'y')

    # compute total matrices
    A = matrix_product(H, R)

    # Now we re-align the translation vector to account for the Sun's height
    # above the midplane
    offset = -translation.transform(H)

    if inverse:
        # the inverse of a rotation matrix is a transpose, which is much faster
        #   and more stable to compute
        A = matrix_transpose(A)
        offset = (-offset).transform(A)
        offset_v = r.CartesianDifferential.from_cartesian(
            (-gcf.galcen_v_sun).to_cartesian().transform(A))
        offset = offset.with_differentials(offset_v)

    else:
        offset = offset.with_differentials(gcf.galcen_v_sun)

    return A, offset


def _check_coord_repr_diff_types(c):
    if isinstance(c.data, r.UnitSphericalRepresentation):
        raise ConvertError("Transforming to/from a Galactocentric frame "
                           "requires a 3D coordinate, e.g. (angle, angle, "
                           "distance) or (x, y, z).")

    if ('s' in c.data.differentials and
            isinstance(c.data.differentials['s'],
                       (r.UnitSphericalDifferential,
                        r.UnitSphericalCosLatDifferential,
                        r.RadialDifferential))):
        raise ConvertError("Transforming to/from a Galactocentric frame "
                           "requires a 3D velocity, e.g., proper motion "
                           "components and radial velocity.")


@frame_transform_graph.transform(AffineTransform, ICRS, Galactocentric)
def icrs_to_galactocentric(icrs_coord, galactocentric_frame):
    _check_coord_repr_diff_types(icrs_coord)
    return get_matrix_vectors(galactocentric_frame)


@frame_transform_graph.transform(AffineTransform, Galactocentric, ICRS)
def galactocentric_to_icrs(galactocentric_coord, icrs_frame):
    _check_coord_repr_diff_types(galactocentric_coord)
    return get_matrix_vectors(galactocentric_coord, inverse=True)
