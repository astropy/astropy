# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy import units as u
from astropy.utils.decorators import format_doc
from astropy.time import Time
from astropy.coordinates import representation as r
from astropy.coordinates.baseframe import (BaseCoordinateFrame,
                                           RepresentationMapping,
                                           frame_transform_graph, base_doc)
from astropy.coordinates.transformations import AffineTransform
from astropy.coordinates.attributes import DifferentialAttribute

from .baseradec import BaseRADecFrame, doc_components as doc_components_radec
from .icrs import ICRS
from .galactic import Galactic

# For speed
J2000 = Time('J2000')

v_bary_Schoenrich2010 = r.CartesianDifferential([11.1, 12.24, 7.25]*u.km/u.s)

__all__ = ['LSR', 'GalacticLSR']


doc_footer_lsr = """
    Other parameters
    ----------------
    v_bary : `~astropy.coordinates.representation.CartesianDifferential`
        The velocity of the solar system barycenter with respect to the LSR, in
        Galactic cartesian velocity components.
"""


@format_doc(base_doc, components=doc_components_radec, footer=doc_footer_lsr)
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
    """

    # frame attributes:
    v_bary = DifferentialAttribute(default=v_bary_Schoenrich2010,
                                   allowed_classes=[r.CartesianDifferential])


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


doc_components_gal = """
    l : `~astropy.coordinates.Angle`, optional, must be keyword
        The Galactic longitude for this object (``b`` must also be given and
        ``representation`` must be None).
    b : `~astropy.coordinates.Angle`, optional, must be keyword
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
"""


@format_doc(base_doc, components=doc_components_gal, footer=doc_footer_lsr)
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
    """

    frame_specific_representation_info = {
        r.SphericalRepresentation: [
            RepresentationMapping('lon', 'l'),
            RepresentationMapping('lat', 'b')
        ]
    }

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
