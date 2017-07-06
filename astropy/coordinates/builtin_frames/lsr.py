# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u
from ...time import Time
from .. import representation as r
from ..baseframe import (BaseCoordinateFrame, RepresentationMapping,
                         frame_transform_graph)
from ..transformations import AffineTransform
from ..attributes import DifferentialAttribute

from .baseradec import _base_radec_docstring, BaseRADecFrame
from .icrs import ICRS
from .galactic import Galactic

# For speed
J2000 = Time('J2000')

v_bary_Schoenrich2010 = r.CartesianDifferential([11.1, 12.24, 7.25]*u.km/u.s)

__all__ = ['LSR', 'GalacticLSR']


class LSR(BaseRADecFrame):
    r"""A coordinate or frame in the Local Standard of Rest (LSR).

    This coordinate frame is axis-aligned and co-spatial with `ICRS`, but has
    a velocity offset relative to the solar system barycenter to remove the
    peculiar motion of the sun relative to the LSR. Roughly, the LSR is the mean
    velocity of the stars in the solar neighborhood, but the precise definition
    of which depends on the study. As defined in Schönrich et al. (2010):
    "The LSR is the rest frame at the location of the Sun of a star that would
    be on a circular orbit in the gravitational potential one would obtain by
    azimuthally averaging away non-axisymmetric features in the actual Galactic
    potential." No such orbit truly exists, but it is still a commonly used
    velocity frame.

    We use default values from Schönrich et al. (2010) for the barycentric
    velocity relative to the LSR, which is defined in Galactic (right-handed)
    cartesian velocity components
    :math:`(U, V, W) = (11.1, 12.24, 7.25)~{{\rm km}}~{{\rm s}}^{{-1}}`. These
    values are customizable via the ``v_bary`` argument which specifies the
    velocity of the solar system barycenter with respect to the LSR.

    The frame attributes are listed under **Other Parameters**.

    {params}

    Other parameters
    ----------------
    v_bary : `~astropy.coordinates.representation.CartesianDifferential`
        The velocity of the solar system barycenter with respect to the LSR, in
        Galactic cartesian velocity components.

    """

    # frame attributes:
    v_bary = DifferentialAttribute(default=v_bary_Schoenrich2010,
                                   allowed_classes=[r.CartesianDifferential])


LSR.__doc__ = LSR.__doc__.format(params=_base_radec_docstring)


@frame_transform_graph.transform(AffineTransform, ICRS, LSR)
def icrs_to_lsr(icrs_coord, lsr_frame):
    v_bary_gal = Galactic(lsr_frame.v_bary.to_cartesian())
    v_bary_icrs = v_bary_gal.transform_to(icrs_coord)
    v_offset = v_bary_icrs.data.represent_as(r.CartesianDifferential)
    offset = r.CartesianRepresentation([0, 0, 0]*u.au, differentials=v_offset)
    return None, offset


@frame_transform_graph.transform(AffineTransform, LSR, ICRS)
def lsr_to_icrs(lsr_coord, icrs_frame):
    v_bary_gal = Galactic(lsr_coord.v_bary.to_cartesian())
    v_bary_icrs = v_bary_gal.transform_to(icrs_frame)
    v_offset = v_bary_icrs.data.represent_as(r.CartesianDifferential)
    offset = r.CartesianRepresentation([0, 0, 0]*u.au, differentials=-v_offset)
    return None, offset

# ------------------------------------------------------------------------------


class GalacticLSR(BaseCoordinateFrame):
    r"""A coordinate or frame in the Local Standard of Rest (LSR), axis-aligned
    to the `Galactic` frame.

    This coordinate frame is axis-aligned and co-spatial with `ICRS`, but has
    a velocity offset relative to the solar system barycenter to remove the
    peculiar motion of the sun relative to the LSR. Roughly, the LSR is the mean
    velocity of the stars in the solar neighborhood, but the precise definition
    of which depends on the study. As defined in Schönrich et al. (2010):
    "The LSR is the rest frame at the location of the Sun of a star that would
    be on a circular orbit in the gravitational potential one would obtain by
    azimuthally averaging away non-axisymmetric features in the actual Galactic
    potential." No such orbit truly exists, but it is still a commonly used
    velocity frame.

    We use default values from Schönrich et al. (2010) for the barycentric
    velocity relative to the LSR, which is defined in Galactic (right-handed)
    cartesian velocity components
    :math:`(U, V, W) = (11.1, 12.24, 7.25)~{{\rm km}}~{{\rm s}}^{{-1}}`. These
    values are customizable via the ``v_bary`` argument which specifies the
    velocity of the solar system barycenter with respect to the LSR.

    The frame attributes are listed under **Other Parameters**.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)

    l : `Angle`, optional, must be keyword
        The Galactic longitude for this object (``b`` must also be given and
        ``representation`` must be None).
    b : `Angle`, optional, must be keyword
        The Galactic latitude for this object (``l`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).

    pm_l_cosb : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in Galactic longitude (including the ``cos(b)`` term)
        for this object (``pm_b`` must also be given).
    pm_b : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in Galactic latitude for this object (``pm_l_cosb``
        must also be given).
    radial_velocity : :class:`~astropy.units.Quantity`, optional, must be keyword
        The radial velocity of this object.

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
    v_bary : `~astropy.coordinates.representation.CartesianDifferential`
        The velocity of the solar system barycenter with respect to the LSR, in
        Galactic cartesian velocity components.
    """

    frame_specific_representation_info = {
        r.SphericalRepresentation: [
            RepresentationMapping('lon', 'l'),
            RepresentationMapping('lat', 'b')
        ],
        r.SphericalCosLatDifferential: [
            RepresentationMapping('d_lon_coslat', 'pm_l_cosb', u.mas/u.yr),
            RepresentationMapping('d_lat', 'pm_b', u.mas/u.yr),
            RepresentationMapping('d_distance', 'radial_velocity', u.km/u.s)
        ],
        r.SphericalDifferential: [
            RepresentationMapping('d_lon', 'pm_l', u.mas/u.yr),
            RepresentationMapping('d_lat', 'pm_b', u.mas/u.yr),
            RepresentationMapping('d_distance', 'radial_velocity', u.km/u.s)
        ],
        r.CartesianDifferential: [
            RepresentationMapping('d_x', 'v_x', u.km/u.s),
            RepresentationMapping('d_y', 'v_y', u.km/u.s),
            RepresentationMapping('d_z', 'v_z', u.km/u.s)
        ],
    }
    frame_specific_representation_info[r.UnitSphericalRepresentation] = \
        frame_specific_representation_info[r.SphericalRepresentation]
    frame_specific_representation_info[r.UnitSphericalCosLatDifferential] = \
        frame_specific_representation_info[r.SphericalCosLatDifferential]
    frame_specific_representation_info[r.UnitSphericalDifferential] = \
        frame_specific_representation_info[r.SphericalDifferential]

    default_representation = r.SphericalRepresentation
    default_differential = r.SphericalCosLatDifferential

    # frame attributes:
    v_bary = DifferentialAttribute(default=v_bary_Schoenrich2010)


@frame_transform_graph.transform(AffineTransform, Galactic, GalacticLSR)
def galactic_to_galacticlsr(galactic_coord, lsr_frame):
    v_bary_gal = Galactic(lsr_frame.v_bary.to_cartesian())
    v_offset = v_bary_gal.data.represent_as(r.CartesianDifferential)
    offset = r.CartesianRepresentation([0, 0, 0]*u.au, differentials=v_offset)
    return None, offset


@frame_transform_graph.transform(AffineTransform, GalacticLSR, Galactic)
def galacticlsr_to_galactic(lsr_coord, galactic_frame):
    v_bary_gal = Galactic(lsr_coord.v_bary.to_cartesian())
    v_offset = v_bary_gal.data.represent_as(r.CartesianDifferential)
    offset = r.CartesianRepresentation([0, 0, 0]*u.au, differentials=-v_offset)
    return None, offset
