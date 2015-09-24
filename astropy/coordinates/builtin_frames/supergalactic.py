# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u
from ..angles import Angle
from ..representation import SphericalRepresentation
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
    """

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'sgl'),
                      RepresentationMapping('lat', 'sgb')],
        'cartesian': [RepresentationMapping('x', 'sgx'),
                      RepresentationMapping('y', 'sgy'),
                      RepresentationMapping('z', 'sgz')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation

    # North supergalactic pole in Galactic coordinates.
    # Needed for transformations to/from Galactic coordinates.
    _nsgp_gal = Galactic(l=47.37*u.degree, b=+6.32*u.degree)
