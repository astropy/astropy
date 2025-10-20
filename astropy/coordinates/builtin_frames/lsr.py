# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy import units as u
from astropy.coordinates import representation as r
from astropy.coordinates.attributes import DifferentialAttribute
from astropy.coordinates.baseframe import (
    BaseCoordinateFrame,
    RepresentationMapping,
    base_doc,
    frame_transform_graph,
)
from astropy.coordinates.transformations import AffineTransform
from astropy.time import Time
from astropy.utils.decorators import format_doc

from .baseradec import BaseRADecFrame
from .baseradec import doc_components as doc_components_radec
from .galactic import Galactic
from .icrs import ICRS

# For speed
J2000 = Time("J2000")

v_bary_Schoenrich2010 = r.CartesianDifferential([11.1, 12.24, 7.25] * u.km / u.s)

__all__ = ["LSR", "LSRD", "LSRK", "GalacticLSR"]


doc_footer_lsr = """
    Other parameters
    ----------------
    v_bary : `~astropy.coordinates.CartesianDifferential`
        The velocity of the solar system barycenter with respect to the LSR, in
        Galactic cartesian velocity components.
"""


@format_doc(base_doc, components=doc_components_radec, footer=doc_footer_lsr)
class LSR(BaseRADecFrame):
    r"""A frame in the Local Standard of Rest (LSR).

    For Earth-bound observers it is often convenient to use a reference
    frame that is tied to the Solar System barycenter, but such frames
    are not very useful for describing galactic dynamics. The dynamical
    LSR is instead tied to the circular velocity at the Sun's location,
    but defining a circular velocity in a non-axisymmetric galaxy
    requires non-trivial averaging. The kinematic LSR is understood as a
    frame in which the average motion of the stars in the solar
    neighborhood is zero, but in practice that is not straightforward
    either because the average motion is different for different
    spectral types.

    The default parameters of this frame are those of the dynamical LSR
    of Schönrich et al. (2010), meaning the Galactic (right-handed)
    Cartesian velocity components of the solar motion are
    :math:`(U, V, W) = (11.1, 12.24, 7.25)~{{\rm km}}~{{\rm s}}^{{-1}}`,
    but a different solar motion can be specified with the ``v_bary``
    argument. The frame is axis-aligned and co-spatial with
    `~astropy.coordinates.ICRS`.

    The frame attributes are listed under **Other Parameters**.

    """

    # frame attributes:
    v_bary = DifferentialAttribute(
        default=v_bary_Schoenrich2010,
        allowed_classes=[r.CartesianDifferential],
        doc="The relative velocity of the solar-system barycenter",
    )


@frame_transform_graph.transform(AffineTransform, ICRS, LSR)
def icrs_to_lsr(icrs_coord, lsr_frame):
    v_bary_gal = Galactic(lsr_frame.v_bary.to_cartesian())
    v_bary_icrs = v_bary_gal.transform_to(icrs_coord)
    v_offset = v_bary_icrs.data.represent_as(r.CartesianDifferential)
    offset = r.CartesianRepresentation([0, 0, 0] * u.au, differentials=v_offset)
    return None, offset


@frame_transform_graph.transform(AffineTransform, LSR, ICRS)
def lsr_to_icrs(lsr_coord, icrs_frame):
    v_bary_gal = Galactic(lsr_coord.v_bary.to_cartesian())
    v_bary_icrs = v_bary_gal.transform_to(icrs_frame)
    v_offset = v_bary_icrs.data.represent_as(r.CartesianDifferential)
    offset = r.CartesianRepresentation([0, 0, 0] * u.au, differentials=-v_offset)
    return None, offset


# ------------------------------------------------------------------------------


doc_components_gal = """
    l : `~astropy.coordinates.Angle`, optional, keyword-only
        The Galactic longitude for this object (``b`` must also be given and
        ``representation`` must be None).
    b : `~astropy.coordinates.Angle`, optional, keyword-only
        The Galactic latitude for this object (``l`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity` ['length'], optional, keyword-only
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).

    pm_l_cosb : `~astropy.units.Quantity` ['angular speed'], optional, keyword-only
        The proper motion in Galactic longitude (including the ``cos(b)`` term)
        for this object (``pm_b`` must also be given).
    pm_b : `~astropy.units.Quantity` ['angular speed'], optional, keyword-only
        The proper motion in Galactic latitude for this object (``pm_l_cosb``
        must also be given).
    radial_velocity : `~astropy.units.Quantity` ['speed'], optional, keyword-only
        The radial velocity of this object.
"""


@format_doc(base_doc, components=doc_components_gal, footer=doc_footer_lsr)
class GalacticLSR(BaseCoordinateFrame):
    r"""A frame in the Local Standard of Rest (LSR), aligned to the Galactic frame.

    For Earth-bound observers it is often convenient to use a reference
    frame that is tied to the Solar System barycenter, but such frames
    are not very useful for describing galactic dynamics. The dynamical
    LSR is instead tied to the circular velocity at the Sun's location,
    but defining a circular velocity in a non-axisymmetric galaxy
    requires non-trivial averaging. The kinematic LSR is understood as a
    frame in which the average motion of the stars in the solar
    neighborhood is zero, but in practice that is not straightforward
    either because the average motion is different for different
    spectral types.

    The default parameters of this frame are those of the dynamical LSR
    of Schönrich et al. (2010), meaning the Galactic (right-handed)
    Cartesian velocity components of the solar motion are
    :math:`(U, V, W) = (11.1, 12.24, 7.25)~{{\rm km}}~{{\rm s}}^{{-1}}`,
    but a different solar motion can be specified with the ``v_bary``
    argument. The frame is rotated relative to the
    `~astropy.coordinates.ICRS` so that it is axis-aligned and
    co-spatial with the `~astropy.coordinates.Galactic` frame.

    The frame attributes are listed under **Other Parameters**.

    """

    frame_specific_representation_info = {
        r.SphericalRepresentation: [
            RepresentationMapping("lon", "l"),
            RepresentationMapping("lat", "b"),
        ]
    }

    default_representation = r.SphericalRepresentation
    default_differential = r.SphericalCosLatDifferential

    # frame attributes:
    v_bary = DifferentialAttribute(
        default=v_bary_Schoenrich2010,
        doc="The relative velocity of the solar-system barycenter",
    )


@frame_transform_graph.transform(AffineTransform, Galactic, GalacticLSR)
def galactic_to_galacticlsr(galactic_coord, lsr_frame):
    v_bary_gal = Galactic(lsr_frame.v_bary.to_cartesian())
    v_offset = v_bary_gal.data.represent_as(r.CartesianDifferential)
    offset = r.CartesianRepresentation([0, 0, 0] * u.au, differentials=v_offset)
    return None, offset


@frame_transform_graph.transform(AffineTransform, GalacticLSR, Galactic)
def galacticlsr_to_galactic(lsr_coord, galactic_frame):
    v_bary_gal = Galactic(lsr_coord.v_bary.to_cartesian())
    v_offset = v_bary_gal.data.represent_as(r.CartesianDifferential)
    offset = r.CartesianRepresentation([0, 0, 0] * u.au, differentials=-v_offset)
    return None, offset


# ------------------------------------------------------------------------------


class LSRK(BaseRADecFrame):
    """A frame in the Kinematic Local Standard of Rest (LSR).

    Conceptually the kinematic LSR is a frame where the average motion
    of the stars in the solar neighborhood is zero. In practice, the
    observed average motion is different for different spectral types,
    which has historically justified using convenient rounded values for
    the solar motion relative to the LSR. This LSRK frame uses the
    definition from

        Gordon 1975, Methods of Experimental Physics: Volume 12:
        Astrophysics, Part C: Radio Observations - Section 6.1.5.

    meaning the solar motion is 20 km/s towards RA=270 Dec=30 (B1900).
    The frame is axis-aligned and co-spatial with `~astropy.coordinates.ICRS`.

    """


# NOTE: To avoid a performance penalty at import time, we hard-code the ICRS
# offsets here. The code to generate the offsets is provided for reproducibility.
# GORDON1975_V_BARY = 20*u.km/u.s
# GORDON1975_DIRECTION = FK4(ra=270*u.deg, dec=30*u.deg, equinox='B1900')
# V_OFFSET_LSRK = ((GORDON1975_V_BARY * GORDON1975_DIRECTION.transform_to(ICRS()).data)
#                  .represent_as(r.CartesianDifferential))

V_OFFSET_LSRK = r.CartesianDifferential(
    [0.28999706839034606, -17.317264789717928, 10.00141199546947] * u.km / u.s
)

ICRS_LSRK_OFFSET = r.CartesianRepresentation(
    [0, 0, 0] * u.au, differentials=V_OFFSET_LSRK
)
LSRK_ICRS_OFFSET = r.CartesianRepresentation(
    [0, 0, 0] * u.au, differentials=-V_OFFSET_LSRK
)


@frame_transform_graph.transform(AffineTransform, ICRS, LSRK)
def icrs_to_lsrk(icrs_coord, lsr_frame):
    return None, ICRS_LSRK_OFFSET


@frame_transform_graph.transform(AffineTransform, LSRK, ICRS)
def lsrk_to_icrs(lsr_coord, icrs_frame):
    return None, LSRK_ICRS_OFFSET


# ------------------------------------------------------------------------------


class LSRD(BaseRADecFrame):
    r"""A frame in the Dynamical Local Standard of Rest (LSR).

    Conceptually the dynamical LSR is a frame moving at the circular
    velocity at the Sun's location. In practice, the concept of a
    circular velocity in a non-axisymmetric galaxy is not trivial.
    This LSRD frame uses the historical definition from

       Delhaye 1965, Solar Motion and Velocity Distribution of
       Common Stars - Section 2.1.

    meaning the solar motion is
    :math:`(U, V, W) = (9, 12, 7)~{{\rm km}}~{{\rm s}}^{{-1}}`,
    or 16.5 km/s towards l=53 b=25. The frame is axis-aligned and
    co-spatial with `~astropy.coordinates.ICRS`.

    """


# NOTE: To avoid a performance penalty at import time, we hard-code the ICRS
# offsets here. The code to generate the offsets is provided for reproducibility.
# V_BARY_DELHAYE1965 = r.CartesianDifferential([9, 12, 7] * u.km/u.s)
# V_OFFSET_LSRD = (Galactic(V_BARY_DELHAYE1965.to_cartesian()).transform_to(ICRS()).data
#                  .represent_as(r.CartesianDifferential))

V_OFFSET_LSRD = r.CartesianDifferential(
    [-0.6382306360182073, -14.585424483191094, 7.8011572411006815] * u.km / u.s
)

ICRS_LSRD_OFFSET = r.CartesianRepresentation(
    [0, 0, 0] * u.au, differentials=V_OFFSET_LSRD
)
LSRD_ICRS_OFFSET = r.CartesianRepresentation(
    [0, 0, 0] * u.au, differentials=-V_OFFSET_LSRD
)


@frame_transform_graph.transform(AffineTransform, ICRS, LSRD)
def icrs_to_lsrd(icrs_coord, lsr_frame):
    return None, ICRS_LSRD_OFFSET


@frame_transform_graph.transform(AffineTransform, LSRD, ICRS)
def lsrd_to_icrs(lsr_coord, icrs_frame):
    return None, LSRD_ICRS_OFFSET


# ------------------------------------------------------------------------------

# Create loopback transformations
frame_transform_graph._add_merged_transform(LSR, ICRS, LSR)
frame_transform_graph._add_merged_transform(GalacticLSR, Galactic, GalacticLSR)
