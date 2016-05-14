# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u
from ..angles import Angle
from ..representation import SphericalRepresentation
from ..baseframe import BaseCoordinateFrame, RepresentationMapping

# these are needed for defining the NGP
from .fk5 import FK5
from .fk4 import FK4NoETerms


class Galactic(BaseCoordinateFrame):
    """
    Galactic Coordinates.

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
    copy : bool, optional
        If `True` (default), make copies of the input coordinate arrays.
        Can only be passed in as a keyword argument.
    """

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'l'),
                      RepresentationMapping('lat', 'b')],
        'cartesian': [RepresentationMapping('x', 'w'),
                      RepresentationMapping('y', 'u'),
                      RepresentationMapping('z', 'v')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation

    # North galactic pole and zeropoint of l in FK4/FK5 coordinates. Needed for
    # transformations to/from FK4/5

    # These are from the IAU's definition of galactic coordinates
    _ngp_B1950 = FK4NoETerms(ra=192.25*u.degree, dec=27.4*u.degree)
    _lon0_B1950 = Angle(123, u.degree)

    # These are *not* from Reid & Brunthaler 2004 - instead, they were
    # derived by doing:
    #
    # >>> FK4NoETerms(ra=192.25*u.degree, dec=27.4*u.degree).transform_to(FK5)
    #
    # This gives better consistency with other codes than using the values
    # from Reid & Brunthaler 2004 and the best self-consistency between FK5
    # -> Galactic and FK5 -> FK4 -> Galactic. The lon0 angle was found by
    # optimizing the self-consistency.
    _ngp_J2000 = FK5(ra=192.8594812065348*u.degree, dec=27.12825118085622*u.degree)
    _lon0_J2000 = Angle(122.9319185680026, u.degree)
