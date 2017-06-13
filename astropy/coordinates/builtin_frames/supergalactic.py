# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u
from ..representation import (CartesianRepresentation,
                              CartesianDifferential,
                              SphericalRepresentation,
                              UnitSphericalRepresentation,
                              SphericalCosLatDifferential,
                              UnitSphericalCosLatDifferential)
from ..baseframe import BaseCoordinateFrame, RepresentationMapping
from .galactic import Galactic


class Supergalactic(BaseCoordinateFrame):
    """
    Supergalactic Coordinates
    (see Lahav et al. 2000, <http://adsabs.harvard.edu/abs/2000MNRAS.312..166L>,
    and references therein).

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    sgl : `Angle`, optional, must be keyword
        The supergalactic longitude for this object (``sgb`` must also be given and
        ``representation`` must be None).
    sgb : `Angle`, optional, must be keyword
        The supergalactic latitude for this object (``sgl`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
    copy : bool, optional
        If `True` (default), make copies of the input coordinate arrays.
        Can only be passed in as a keyword argument.
    """

    frame_specific_representation_info = {
        SphericalRepresentation: [
            RepresentationMapping('lon', 'sgl'),
            RepresentationMapping('lat', 'sgb')
        ],
        CartesianRepresentation: [
            RepresentationMapping('x', 'sgx'),
            RepresentationMapping('y', 'sgy'),
            RepresentationMapping('z', 'sgz')
        ],
        SphericalCosLatDifferential: [
            RepresentationMapping('d_lon_coslat', 'pm_sgl'),
            RepresentationMapping('d_lat', 'pm_sgb'),
            RepresentationMapping('d_distance', 'radial_velocity'),
        ],
        CartesianDifferential: [
            RepresentationMapping('d_x', 'v_x'),
            RepresentationMapping('d_y', 'v_y'),
            RepresentationMapping('d_z', 'v_z')
        ],
    }
    frame_specific_representation_info[UnitSphericalRepresentation] = \
        frame_specific_representation_info[SphericalRepresentation]
    frame_specific_representation_info[UnitSphericalCosLatDifferential] = \
        frame_specific_representation_info[SphericalCosLatDifferential]

    default_representation = SphericalRepresentation
    default_differential = SphericalCosLatDifferential

    # North supergalactic pole in Galactic coordinates.
    # Needed for transformations to/from Galactic coordinates.
    _nsgp_gal = Galactic(l=47.37*u.degree, b=+6.32*u.degree)
