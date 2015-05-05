# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..representation import SphericalRepresentation
from ..baseframe import BaseCoordinateFrame, RepresentationMapping, TimeFrameAttribute
from .utils import EQUINOX_J2000


class GeocentricEcliptic(BaseCoordinateFrame):
    """
    Geocentric ecliptic coordinates.  These origin of the coordinates are the
    geocenter (Earth), with the x axis pointing to the *true* (not mean) equinox
    at the time specified by the ``equinox`` attribute, and the xy-plane in the
    plane of the ecliptic for that date.

    Be aware that the definition of "geocentric" here means that this frame
    *includes* light deflection from the sun, abberation, etc when transfoming
    to/from e.g. ICRS.

    This frame has one frame attribute:

    * ``equinox``
        The date to assume for the .  Determines the location of
        the x-axis and the location of the Earth (necessary for transformation
        to non-geocentric systems)

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    lambda : `Angle`, optional, must be keyword
        The ecliptic longitude for this object (``beta`` must also be given and
        ``representation`` must be None).
    beta : `Angle`, optional, must be keyword
        The ecliptic latitde for this object (``lambda`` must also be given and
        ``representation`` must be None).
    delta : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object from the geocenter.
        (``representation`` must be None).
    """

    frame_specific_representation_info = {
        'unitspherical': [RepresentationMapping('lon', 'lambda'),
                          RepresentationMapping('lat', 'beta')]
    }
    frame_specific_representation_info['spherical'] = \
        frame_specific_representation_info['unitspherical'] + \
        [RepresentationMapping('distance', 'delta')]

    default_representation = SphericalRepresentation

    equinox = TimeFrameAttribute(default=EQUINOX_J2000)


class BarycentricEcliptic(BaseCoordinateFrame):
    """
    Barycentric ecliptic coordinates.  These origin of the coordinates are the
    barycenter of the solar system, with the x axis pointing in the direction of
    the *true* (not mean) equinox as at the time specified by the ``equinox``
    attribute (as seen from Earth), and the xy-plane in the plane of the
    ecliptic for that date.

    This frame has one frame attribute:

    * ``equinox``
        The date to assume for the .  Determines the location of
        the x-axis and the location of the Earth (necessary for transformation
        to non-geocentric systems)

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    l : `Angle`, optional, must be keyword
        The ecliptic longitude for this object (``beta`` must also be given and
        ``representation`` must be None).
    b : `Angle`, optional, must be keyword
        The ecliptic latitde for this object (``lambda`` must also be given and
        ``representation`` must be None).
    r : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object from the sun's center.
        (``representation`` must be None).
    """

    frame_specific_representation_info = {
        'unitspherical': [RepresentationMapping('lon', 'l'),
                          RepresentationMapping('lat', 'b')]
    }
    frame_specific_representation_info['spherical'] = \
        frame_specific_representation_info['unitspherical'] + \
        [RepresentationMapping('distance', 'r')]

    default_representation = SphericalRepresentation

    equinox = TimeFrameAttribute(default=EQUINOX_J2000)
