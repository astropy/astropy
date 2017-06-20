# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u
from ..representation import (CartesianDifferential,
                              SphericalRepresentation,
                              UnitSphericalRepresentation,
                              SphericalCosLatDifferential,
                              UnitSphericalCosLatDifferential)
from ..baseframe import BaseCoordinateFrame, RepresentationMapping
from ..frame_attributes import TimeFrameAttribute
from .utils import EQUINOX_J2000

__all__ = ['GeocentricTrueEcliptic', 'BarycentricTrueEcliptic',
           'HeliocentricTrueEcliptic']


class BaseEclipticFrame(BaseCoordinateFrame):
    """
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
        The ecliptic latitude for this object (``lon`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The distance for this object from the {0}.
        (``representation`` must be None).

    pm_lon : `Angle`, optional, must be keyword
        The proper motion in the ecliptic longitude for this object (``pm_lat``
        must also be given).
    pm_lat : `Angle`, optional, must be keyword
        The proper motion in the ecliptic latitude for this object (``pm_lon``
        must also be given).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The distance for this object from the {0}.
        (``representation`` must be None).

    copy : bool, optional
        If `True` (default), make copies of the input coordinate arrays.
        Can only be passed in as a keyword argument.
    """

    frame_specific_representation_info = {
        SphericalCosLatDifferential: [
            RepresentationMapping('d_lon_coslat', 'pm_lon', u.mas/u.yr),
            RepresentationMapping('d_lat', 'pm_alt', u.mas/u.yr),
            RepresentationMapping('d_distance', 'radial_velocity', u.km/u.s),
        ]
    }

    frame_specific_representation_info[UnitSphericalCosLatDifferential] = \
        frame_specific_representation_info[SphericalCosLatDifferential]

    default_representation = SphericalRepresentation
    default_differential = SphericalCosLatDifferential


class GeocentricTrueEcliptic(BaseEclipticFrame):
    """
    Geocentric ecliptic coordinates.  These origin of the coordinates are the
    geocenter (Earth), with the x axis pointing to the *true* (not mean) equinox
    at the time specified by the ``equinox`` attribute, and the xy-plane in the
    plane of the ecliptic for that date.

    Be aware that the definition of "geocentric" here means that this frame
    *includes* light deflection from the sun, aberration, etc when transforming
    to/from e.g. ICRS.

    This frame has one frame attribute:

    * ``equinox``
        The date to assume for this frame.  Determines the location of the
        x-axis and the location of the Earth (necessary for transformation to
        non-geocentric systems).
    """

    equinox = TimeFrameAttribute(default=EQUINOX_J2000)


GeocentricTrueEcliptic.__doc__ += BaseEclipticFrame.__doc__.format("geocenter")


class BarycentricTrueEcliptic(BaseEclipticFrame):
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
    """

    equinox = TimeFrameAttribute(default=EQUINOX_J2000)


BarycentricTrueEcliptic.__doc__ += BaseEclipticFrame.__doc__.format("sun's center")


class HeliocentricTrueEcliptic(BaseEclipticFrame):
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
    """
    default_representation = SphericalRepresentation
    default_differential = SphericalCosLatDifferential

    equinox = TimeFrameAttribute(default=EQUINOX_J2000)


HeliocentricTrueEcliptic.__doc__ += BaseEclipticFrame.__doc__.format("sun's center")
