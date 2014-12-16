# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..representation import SphericalRepresentation
from ..baseframe import BaseCoordinateFrame, RepresentationMapping, TimeFrameAttribute
from .utils import DEFAULT_OBSTIME

class CIRS(BaseCoordinateFrame):
    """
    A coordinate or frame in the Celestial Intermediate Reference System (CIRS).

    This frame has one frame attribute:

    * ``obstime``
        The time at which the observation is taken.  Used for determining the
        position of the Earth and its precession.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    ra : `Angle`, optional, must be keyword
        The RA for this object (``dec`` must also be given and ``representation``
        must be None).
    dec : `Angle`, optional, must be keyword
        The Declination for this object (``ra`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).
    """

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'ra'),
                      RepresentationMapping('lat', 'dec')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation

    obstime = TimeFrameAttribute(default=DEFAULT_OBSTIME)

#The "self-transform" is defined in icrs_cirs_transformations.py, because in
#the current implementation it goes through ICRS (like GCRS)
