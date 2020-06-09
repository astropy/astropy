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
from .fk4 import FK4

# For speed
J2000 = Time('J2000')

v_bary_Schoenrich2010 = r.CartesianDifferential([11.1, 12.24, 7.25]*u.km/u.s)

__all__ = ['LSR', 'GalacticLSR', 'LSRK', 'LSRD']


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


# ------------------------------------------------------------------------------

# The LSRK velocity frame, defined as having a velocity of 20 km/s towards
# RA=270 Dec=30 (B1900) relative to the solar system Barycenter. This is defined
# in:
#
#   Gordon 1975, Methods of Experimental Physics: Volume 12:
#   Astrophysics, Part C: Radio Observations - Section 6.1.5.


class LSRK(BaseRADecFrame):
    r"""
    A coordinate or frame in the Kinematic Local Standard of Rest (LSR).

    This frame is defined as having a velocity of 20 km/s towards RA=270 Dec=30
    (B1900) relative to the solar system Barycenter. This is defined in:

        Gordon 1975, Methods of Experimental Physics: Volume 12:
        Astrophysics, Part C: Radio Observations - Section 6.1.5.

    This coordinate frame is axis-aligned and co-spatial with `ICRS`, but has
    a velocity offset relative to the solar system barycenter to remove the
    peculiar motion of the sun relative to the LSRK.
    """


# NOTE: To avoid a performance penalty at import time, we hard-code the ICRS
# offsets here. The code to generate the offsets is provided for reproducibility.
# GORDON1975_V_BARY = 20*u.km/u.s
# GORDON1975_DIRECTION = FK4(ra=270*u.deg, dec=30*u.deg, equinox='B1900')
# V_OFFSET_LSRK = ((GORDON1975_V_BARY * GORDON1975_DIRECTION.transform_to(ICRS()).data)
#                  .represent_as(r.CartesianDifferential))

V_OFFSET_LSRK = r.CartesianDifferential([0.28999706839034606,
                                         -17.317264789717928,
                                         10.00141199546947]*u.km/u.s)

ICRS_LSRK_OFFSET = r.CartesianRepresentation([0, 0, 0]*u.au, differentials=V_OFFSET_LSRK)
LSRK_ICRS_OFFSET = r.CartesianRepresentation([0, 0, 0]*u.au, differentials=-V_OFFSET_LSRK)


@frame_transform_graph.transform(AffineTransform, ICRS, LSRK)
def icrs_to_lsrk(icrs_coord, lsr_frame):
    return None, ICRS_LSRK_OFFSET


@frame_transform_graph.transform(AffineTransform, LSRK, ICRS)
def lsrk_to_icrs(lsr_coord, icrs_frame):
    return None, LSRK_ICRS_OFFSET


# ------------------------------------------------------------------------------

# The LSRD velocity frame, defined as a velocity of U=9 km/s, V=12 km/s,
# and W=7 km/s in Galactic coordinates or 16.552945 km/s
# towards l=53.13 b=25.02. This is defined in:
#
#   Delhaye 1965, Solar Motion and Velocity Distribution of
#   Common Stars.


class LSRD(BaseRADecFrame):
    r"""
    A coordinate or frame in the Dynamical Local Standard of Rest (LSRD)

    This frame is defined as a velocity of U=9 km/s, V=12 km/s,
    and W=7 km/s in Galactic coordinates or 16.552945 km/s
    towards l=53.13 b=25.02. This is defined in:

       Delhaye 1965, Solar Motion and Velocity Distribution of
       Common Stars.

    This coordinate frame is axis-aligned and co-spatial with `ICRS`, but has
    a velocity offset relative to the solar system barycenter to remove the
    peculiar motion of the sun relative to the LSRD.
    """


# NOTE: To avoid a performance penalty at import time, we hard-code the ICRS
# offsets here. The code to generate the offsets is provided for reproducibility.
# V_BARY_DELHAYE1965 = r.CartesianDifferential([9, 12, 7] * u.km/u.s)
# V_OFFSET_LSRD = (Galactic(V_BARY_DELHAYE1965.to_cartesian()).transform_to(ICRS()).data
#                  .represent_as(r.CartesianDifferential))

V_OFFSET_LSRD = r.CartesianDifferential([-0.6382306360182073,
                                         -14.585424483191094,
                                         7.8011572411006815]*u.km/u.s)

ICRS_LSRD_OFFSET = r.CartesianRepresentation([0, 0, 0]*u.au, differentials=V_OFFSET_LSRD)
LSRD_ICRS_OFFSET = r.CartesianRepresentation([0, 0, 0]*u.au, differentials=-V_OFFSET_LSRD)


@frame_transform_graph.transform(AffineTransform, ICRS, LSRD)
def icrs_to_lsrd(icrs_coord, lsr_frame):
    return None, ICRS_LSRD_OFFSET


@frame_transform_graph.transform(AffineTransform, LSRD, ICRS)
def lsrd_to_icrs(lsr_coord, icrs_frame):
    return None, LSRD_ICRS_OFFSET
