# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ... import units as u
from ...utils.decorators import format_doc
from .. import representation as r
from ..baseframe import BaseCoordinateFrame, RepresentationMapping, base_doc
from ..attributes import TimeAttribute
from .utils import EQUINOX_J2000, DEFAULT_OBSTIME

__all__ = ['GeocentricTrueEcliptic', 'BarycentricTrueEcliptic',
           'HeliocentricTrueEcliptic', 'BaseEclipticFrame']


doc_components_ecl = """
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
    radial_velocity : :class:`~astropy.units.Quantity`, optional, must be keyword
        The radial velocity of this object.
"""


@format_doc(base_doc,
            components=doc_components_ecl.format('specified location'),
            footer="")
class BaseEclipticFrame(BaseCoordinateFrame):
    """
    A base class for frames that have names and conventions like that of
    ecliptic frames.

    .. warning::
            In the current version of astropy, the ecliptic frames do not yet have
            stringent accuracy tests.  We recommend you test to "known-good" cases
            to ensure this frames are what you are looking for. (and then ideally
            you would contribute these tests to Astropy!)
    """

    default_representation = r.SphericalRepresentation
    default_differential = r.SphericalCosLatDifferential


doc_footer_geo = """
    Other parameters
    ----------------
    equinox : `~astropy.time.Time`, optional
        The date to assume for this frame.  Determines the location of the
        x-axis and the location of the Earth (necessary for transformation to
        non-geocentric systems). Defaults to the 'J2000' equinox.
"""


@format_doc(base_doc, components=doc_components_ecl.format('geocenter'),
            footer=doc_footer_geo)
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
    """

    equinox = TimeAttribute(default=EQUINOX_J2000)


doc_footer_bary = """
    Other parameters
    ----------------
    equinox : `~astropy.time.Time`, optional
        The date to assume for this frame.  Determines the location of the
        x-axis and the location of the Earth and Sun.
        Defaults to the 'J2000' equinox.
"""


@format_doc(base_doc, components=doc_components_ecl.format("barycenter"),
            footer=doc_footer_bary)
class BarycentricTrueEcliptic(BaseEclipticFrame):
    """
    Barycentric ecliptic coordinates.  These origin of the coordinates are the
    barycenter of the solar system, with the x axis pointing in the direction of
    the *true* (not mean) equinox as at the time specified by the ``equinox``
    attribute (as seen from Earth), and the xy-plane in the plane of the
    ecliptic for that date.

    The frame attributes are listed under **Other Parameters**.
    """

    equinox = TimeAttribute(default=EQUINOX_J2000)


doc_footer_helio = """
    Other parameters
    ----------------
    equinox : `~astropy.time.Time`, optional
        The date to assume for this frame.  Determines the location of the
        x-axis and the location of the Earth and Sun.
        Defaults to the 'J2000' equinox.
"""


@format_doc(base_doc, components=doc_components_ecl.format("sun's center"),
            footer=doc_footer_helio)
class HeliocentricTrueEcliptic(BaseEclipticFrame):
    """
    Heliocentric ecliptic coordinates.  These origin of the coordinates are the
    center of the sun, with the x axis pointing in the direction of
    the *true* (not mean) equinox as at the time specified by the ``equinox``
    attribute (as seen from Earth), and the xy-plane in the plane of the
    ecliptic for that date.

    The frame attributes are listed under **Other Parameters**.

    {params}


    """

    equinox = TimeAttribute(default=EQUINOX_J2000)
    obstime = TimeAttribute(default=DEFAULT_OBSTIME)
