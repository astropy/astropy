# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..representation import SphericalRepresentation
from ..baseframe import BaseCoordinateFrame, RepresentationMapping


class ICRS(BaseCoordinateFrame):
    """
    A coordinate or frame in the ICRS system.

    If you're looking for "J2000" coordinates, and aren't sure if you want to
    use this or `~astropy.coordinates.FK5`, you probably want to use ICRS. It's
    more well-defined as a catalog coordinate and is an inertial system, and is
    very close (within tens of milliarcseconds) to J2000 equatorial.

    For more background on the ICRS and related coordinate transformations, see the
    references provided in the  :ref:`astropy-coordinates-seealso` section of the
    documentation.

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
    copy : bool, optional
        If `True` (default), make copies of the input coordinate arrays.
        Can only be passed in as a keyword argument.
    """

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'ra'),
                      RepresentationMapping('lat', 'dec')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation
