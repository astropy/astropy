# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u
from .. import representation as r
from ..baseframe import BaseCoordinateFrame, RepresentationMapping
from ..attributes import TimeAttribute
from .utils import EQUINOX_J2000, DEFAULT_OBSTIME

__all__ = ['GeocentricTrueEcliptic', 'BarycentricTrueEcliptic',
           'HeliocentricTrueEcliptic', 'BaseEclipticFrame']

_base_ecliptic_docstring = """.. warning::
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

    pm_lon_coslat : `Angle`, optional, must be keyword
        The proper motion in the ecliptic longitude (including the ``cos(lat)``
        factor) for this object (``pm_lat`` must also be given).
    pm_lat : `Angle`, optional, must be keyword
        The proper motion in the ecliptic latitude for this object
        (``pm_lon_coslat`` must also be given).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The distance for this object from the {0}.
        (``representation`` must be None).

    copy : bool, optional
        If `True` (default), make copies of the input coordinate arrays.
        Can only be passed in as a keyword argument.

    differential_cls : `BaseDifferential`, dict, optional
        A differential class or dictionary of differential classes (currently
        only a velocity differential with key 's' is supported). This sets
        the expected input differential class, thereby changing the expected
        keyword arguments of the data passed in. For example, passing
        ``differential_cls=CartesianDifferential`` will make the classes
        expect velocity data with the argument names ``v_x, v_y, v_z``.
"""


class BaseEclipticFrame(BaseCoordinateFrame):
    """
    A base class for frames that have names and conventions like that of
    ecliptic frames.

    {params}
    """

    frame_specific_representation_info = {
        r.SphericalCosLatDifferential: [
            RepresentationMapping('d_lon_coslat', 'pm_lon_coslat', u.mas/u.yr),
            RepresentationMapping('d_lat', 'pm_lat', u.mas/u.yr),
            RepresentationMapping('d_distance', 'radial_velocity', u.km/u.s),
        ],
        r.SphericalDifferential: [
            RepresentationMapping('d_lon', 'pm_lon', u.mas/u.yr),
            RepresentationMapping('d_lat', 'pm_lat', u.mas/u.yr),
            RepresentationMapping('d_distance', 'radial_velocity', u.km/u.s),
        ],
        r.CartesianDifferential: [
            RepresentationMapping('d_x', 'v_x', u.km/u.s),
            RepresentationMapping('d_y', 'v_y', u.km/u.s),
            RepresentationMapping('d_z', 'v_z', u.km/u.s),
        ],
    }

    frame_specific_representation_info[r.UnitSphericalCosLatDifferential] = \
        frame_specific_representation_info[r.SphericalCosLatDifferential]
    frame_specific_representation_info[r.UnitSphericalDifferential] = \
        frame_specific_representation_info[r.SphericalDifferential]

    default_representation = r.SphericalRepresentation
    default_differential = r.SphericalCosLatDifferential


BaseEclipticFrame.__doc__ = BaseEclipticFrame.__doc__.format(
    params=_base_ecliptic_docstring)


class GeocentricTrueEcliptic(BaseEclipticFrame):
    """
    Geocentric ecliptic coordinates.  These origin of the coordinates are the
    geocenter (Earth), with the x axis pointing to the *true* (not mean) equinox
    at the time specified by the ``equinox`` attribute, and the xy-plane in the
    plane of the ecliptic for that date.

    Be aware that the definition of "geocentric" here means that this frame
    *includes* light deflection from the sun, aberration, etc when transforming
    to/from e.g. ICRS.

    The frame attributes are listed under **Other Parameters**.

    {params}

    Other parameters
    ----------------
    equinox : `~astropy.time.Time`, optional
        The date to assume for this frame.  Determines the location of the
        x-axis and the location of the Earth (necessary for transformation to
        non-geocentric systems). Defaults to the 'J2000' equinox.
    """

    equinox = TimeAttribute(default=EQUINOX_J2000)


GeocentricTrueEcliptic.__doc__ = GeocentricTrueEcliptic.__doc__.format(
    params=_base_ecliptic_docstring.format("geocenter"))


class BarycentricTrueEcliptic(BaseEclipticFrame):
    """
    Barycentric ecliptic coordinates.  These origin of the coordinates are the
    barycenter of the solar system, with the x axis pointing in the direction of
    the *true* (not mean) equinox as at the time specified by the ``equinox``
    attribute (as seen from Earth), and the xy-plane in the plane of the
    ecliptic for that date.

    The frame attributes are listed under **Other Parameters**.

    {params}

    Other parameters
    ----------------
    equinox : `~astropy.time.Time`, optional
        The date to assume for this frame.  Determines the location of the
        x-axis and the location of the Earth and Sun.
        Defaults to the 'J2000' equinox.
    """

    equinox = TimeAttribute(default=EQUINOX_J2000)


BarycentricTrueEcliptic.__doc__ = BarycentricTrueEcliptic.__doc__.format(
    params=_base_ecliptic_docstring.format("sun's center"))


class HeliocentricTrueEcliptic(BaseEclipticFrame):
    """
    Heliocentric ecliptic coordinates.  These origin of the coordinates are the
    center of the sun, with the x axis pointing in the direction of
    the *true* (not mean) equinox as at the time specified by the ``equinox``
    attribute (as seen from Earth), and the xy-plane in the plane of the
    ecliptic for that date.

    The frame attributes are listed under **Other Parameters**.

    {params}

    Other parameters
    ----------------
    equinox : `~astropy.time.Time`, optional
        The date to assume for this frame.  Determines the location of the
        x-axis and the location of the Earth and Sun.
        Defaults to the 'J2000' equinox.
    """

    equinox = TimeAttribute(default=EQUINOX_J2000)
    obstime = TimeAttribute(default=DEFAULT_OBSTIME)


HeliocentricTrueEcliptic.__doc__ = HeliocentricTrueEcliptic.__doc__.format(
    params=_base_ecliptic_docstring.format("sun's center"))
