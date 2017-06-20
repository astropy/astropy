# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u
from ...time import Time
from ..representation import (CartesianRepresentation,
                              CartesianDifferential,
                              SphericalRepresentation,
                              UnitSphericalRepresentation,
                              SphericalCosLatDifferential,
                              UnitSphericalCosLatDifferential)
from ..baseframe import (BaseCoordinateFrame, RepresentationMapping,
                         frame_transform_graph, BaseRADecFrame)
from ..transformations import AffineTransform
from ..frame_attributes import DifferentialFrameAttribute

from .icrs import ICRS
from .galactic import Galactic

# For speed
J2000 = Time('J2000')

v_bary_Schoenrich2010 = CartesianDifferential([-11.1, 12.24, 7.25]*u.km/u.s)

__all__ = ['LSR', 'GalacticLSR']

class LSR(BaseRADecFrame):
    """
    A coordinate or frame in the Local Standard of Rest (LSR).

    TODO: more words
    - axis-aligned with ICRS
    - co-spatial (SS barycenter)
    - v_bary assumed to be Galactic
    """

    # frame attributes:
    v_bary = DifferentialFrameAttribute(default=v_bary_Schoenrich2010)

LSR.__doc__ += BaseRADecFrame.__doc__

@frame_transform_graph.transform(AffineTransform, ICRS, LSR)
def icrs_to_lsr(icrs_coord, lsr_frame):
    v_bary_gal = Galactic(lsr_frame.v_bary.to_cartesian())
    v_bary_icrs = v_bary_gal.transform_to(icrs_coord)
    v_offset = v_bary_icrs.data.represent_as(CartesianDifferential)
    offset = CartesianRepresentation([0,0,0]*u.au, differentials=v_offset)
    return None, offset

@frame_transform_graph.transform(AffineTransform, LSR, ICRS)
def lsr_to_icrs(lsr_coord, icrs_frame):
    v_bary_gal = Galactic(lsr_coord.v_bary.to_cartesian())
    v_bary_icrs = v_bary_gal.transform_to(icrs_frame)
    v_offset = v_bary_icrs.data.represent_as(CartesianDifferential)
    offset = CartesianRepresentation([0,0,0]*u.au, differentials=-v_offset)
    return None, offset


class GalacticLSR(BaseCoordinateFrame):
    """
    A coordinate or frame in the Local Standard of Rest (LSR).

    TODO: more words
    - axis-aligned with Galactic
    - co-spatial (SS barycenter)

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
    copy : bool, optional
        If `True` (default), make copies of the input coordinate arrays.
        Can only be passed in as a keyword argument.
    """

    frame_specific_representation_info = {
        SphericalRepresentation: [
            RepresentationMapping('lon', 'l'),
            RepresentationMapping('lat', 'b')
        ],
        SphericalCosLatDifferential: [
            RepresentationMapping('d_lon_coslat', 'pm_l', u.mas/u.yr),
            RepresentationMapping('d_lat', 'pm_b', u.mas/u.yr),
            RepresentationMapping('d_distance', 'radial_velocity', u.km/u.s)
        ],
        CartesianDifferential: [
            RepresentationMapping('d_x', 'v_x', u.km/u.s),
            RepresentationMapping('d_y', 'v_y', u.km/u.s),
            RepresentationMapping('d_z', 'v_z', u.km/u.s)
        ],
    }
    frame_specific_representation_info[UnitSphericalRepresentation] = \
        frame_specific_representation_info[SphericalRepresentation]
    frame_specific_representation_info[UnitSphericalCosLatDifferential] = \
        frame_specific_representation_info[SphericalCosLatDifferential]

    default_representation = SphericalRepresentation
    default_differential = SphericalCosLatDifferential

    # frame attributes:
    v_bary = DifferentialFrameAttribute(default=v_bary_Schoenrich2010)


@frame_transform_graph.transform(AffineTransform, Galactic, GalacticLSR)
def galactic_to_galacticlsr(galactic_coord, lsr_frame):
    v_bary_gal = Galactic(lsr_frame.v_bary.to_cartesian())
    v_offset = v_bary_gal.data.represent_as(CartesianDifferential)
    offset = CartesianRepresentation([0,0,0]*u.au, differentials=v_offset)
    return None, offset

@frame_transform_graph.transform(AffineTransform, GalacticLSR, Galactic)
def galacticlsr_to_galactic(lsr_coord, galactic_frame):
    v_bary_gal = Galactic(lsr_coord.v_bary.to_cartesian())
    v_offset = v_bary_gal.data.represent_as(CartesianDifferential)
    offset = CartesianRepresentation([0,0,0]*u.au, differentials=-v_offset)
    return None, offset
