# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..representation import SphericalRepresentation
from ..baseframe import BaseCoordinateFrame, RepresentationMapping, TimeFrameAttribute
from .utils import EQUINOX_J2000


__all__ = ['GeocentricTrueEcliptic', 'BarycentricTrueEcliptic', 'HeliocentricTrueEcliptic']


class GeocentricTrueEcliptic(BaseCoordinateFrame):
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
        The date to assume for this frame.  Determines the location of the
        x-axis and the location of the Earth (necessary for transformation to
        non-geocentric systems).

    .. warning::
        In the current version of astropy, the ecliptic frames do not yet have
        stringent accuracy tests.  We recommend you test to "known-good" cases
        to ensure this frames are what you are looking for. (and then ideally
        you would contribute these tests to Astropy!)

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    lon : `Angle`, optional, must be keyword
        The ecliptic longitude for this object (``lat`` must also be given and
        ``representation`` must be None).
    lat : `Angle`, optional, must be keyword
        The ecliptic latitde for this object (``lon`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object from the geocenter.
        (``representation`` must be None).
    """
    default_representation = SphericalRepresentation

    equinox = TimeFrameAttribute(default=EQUINOX_J2000)


class BarycentricTrueEcliptic(BaseCoordinateFrame):
    """
    Barycentric ecliptic coordinates.  These origin of the coordinates are the
    barycenter of the solar system, with the x axis pointing in the direction of
    the *true* (not mean) equinox as at the time specified by the ``equinox``
    attribute (as seen from Earth), and the xy-plane in the plane of the
    ecliptic for that date.

    This frame has one frame attribute:

    * ``equinox``
        The date to assume for this frame.  Determines the location of the
        x-axis and the location of the Earth and Sun.

    .. warning::
        In the current version of astropy, the ecliptic frames do not yet have
        stringent accuracy tests.  We recommend you test to "known-good" cases
        to ensure this frames are what you are looking for. (and then ideally
        you would contribute these tests to Astropy!)

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    l : `Angle`, optional, must be keyword
        The ecliptic longitude for this object (``b`` must also be given and
        ``representation`` must be None).
    b : `Angle`, optional, must be keyword
        The ecliptic latitde for this object (``l`` must also be given and
        ``representation`` must be None).
    r : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object from the sun's center.
        (``representation`` must be None).
    """
    default_representation = SphericalRepresentation

    equinox = TimeFrameAttribute(default=EQUINOX_J2000)


class HeliocentricTrueEcliptic(BaseCoordinateFrame):
    """
    Heliocentric ecliptic coordinates.  These origin of the coordinates are the
    center of the sun, with the x axis pointing in the direction of
    the *true* (not mean) equinox as at the time specified by the ``equinox``
    attribute (as seen from Earth), and the xy-plane in the plane of the
    ecliptic for that date.

    This frame has one frame attribute:

    * ``equinox``
        The date to assume for this frame.  Determines the location of the
        x-axis and the location of the Earth and Sun.

    .. warning::
        In the current version of astropy, the ecliptic frames do not yet have
        stringent accuracy tests.  We recommend you test to "known-good" cases
        to ensure this frames are what you are looking for. (and then ideally
        you would contribute these tests to Astropy!)

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    l : `Angle`, optional, must be keyword
        The ecliptic longitude for this object (``b`` must also be given and
        ``representation`` must be None).
    b : `Angle`, optional, must be keyword
        The ecliptic latitde for this object (``l`` must also be given and
        ``representation`` must be None).
    r : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object from the sun's center.
        (``representation`` must be None).
    """
    default_representation = SphericalRepresentation

    equinox = TimeFrameAttribute(default=EQUINOX_J2000)
