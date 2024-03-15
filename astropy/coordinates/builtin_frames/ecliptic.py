# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy import units as u
from astropy.coordinates import representation as r
from astropy.coordinates.attributes import QuantityAttribute, TimeAttribute
from astropy.coordinates.baseframe import BaseCoordinateFrame, base_doc
from astropy.utils.decorators import format_doc

from .utils import DEFAULT_OBSTIME, EQUINOX_J2000

__all__ = [
    "GeocentricMeanEcliptic",
    "BarycentricMeanEcliptic",
    "HeliocentricMeanEcliptic",
    "BaseEclipticFrame",
    "GeocentricTrueEcliptic",
    "BarycentricTrueEcliptic",
    "HeliocentricTrueEcliptic",
    "HeliocentricEclipticIAU76",
    "CustomBarycentricEcliptic",
]


doc_components_ecl = """
    lon : `~astropy.coordinates.Angle`, optional, keyword-only
        The ecliptic longitude for this object (``lat`` must also be given and
        ``representation`` must be None).
    lat : `~astropy.coordinates.Angle`, optional, keyword-only
        The ecliptic latitude for this object (``lon`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity` ['length'], optional, keyword-only
        The distance for this object from the {0}.
        (``representation`` must be None).

    pm_lon_coslat : `~astropy.units.Quantity` ['angualar speed'], optional, keyword-only
        The proper motion in the ecliptic longitude (including the ``cos(lat)``
        factor) for this object (``pm_lat`` must also be given).
    pm_lat : `~astropy.units.Quantity` ['angular speed'], optional, keyword-only
        The proper motion in the ecliptic latitude for this object
        (``pm_lon_coslat`` must also be given).
    radial_velocity : `~astropy.units.Quantity` ['speed'], optional, keyword-only
        The radial velocity of this object.
"""


@format_doc(
    base_doc, components=doc_components_ecl.format("specified location"), footer=""
)
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
    obstime : `~astropy.time.Time`, optional
        The time at which the observation is taken.  Used for determining the
        position of the Earth. Defaults to J2000.
"""


@format_doc(
    base_doc, components=doc_components_ecl.format("geocenter"), footer=doc_footer_geo
)
class GeocentricMeanEcliptic(BaseEclipticFrame):
    """
    Geocentric mean ecliptic coordinates.  These origin of the coordinates are the
    geocenter (Earth), with the x axis pointing to the *mean* (not true) equinox
    at the time specified by the ``equinox`` attribute, and the xy-plane in the
    plane of the ecliptic for that date.

    Be aware that the definition of "geocentric" here means that this frame
    *includes* light deflection from the sun, aberration, etc when transforming
    to/from e.g. ICRS.

    The frame attributes are listed under **Other Parameters**.
    """

    equinox = TimeAttribute(default=EQUINOX_J2000, doc="The equinox time")
    obstime = TimeAttribute(
        default=DEFAULT_OBSTIME, doc="The reference time (e.g., time of observation)"
    )


@format_doc(
    base_doc, components=doc_components_ecl.format("geocenter"), footer=doc_footer_geo
)
class GeocentricTrueEcliptic(BaseEclipticFrame):
    """
    Geocentric true ecliptic coordinates.  These origin of the coordinates are the
    geocenter (Earth), with the x axis pointing to the *true* (not mean) equinox
    at the time specified by the ``equinox`` attribute, and the xy-plane in the
    plane of the ecliptic for that date.

    Be aware that the definition of "geocentric" here means that this frame
    *includes* light deflection from the sun, aberration, etc when transforming
    to/from e.g. ICRS.

    The frame attributes are listed under **Other Parameters**.
    """

    equinox = TimeAttribute(default=EQUINOX_J2000, doc="The equinox time")
    obstime = TimeAttribute(
        default=DEFAULT_OBSTIME, doc="The reference time (e.g., time of observation)"
    )


doc_footer_bary = """
    Other parameters
    ----------------
    equinox : `~astropy.time.Time`, optional
        The date to assume for this frame.  Determines the location of the
        x-axis and the location of the Earth and Sun.
        Defaults to the 'J2000' equinox.
"""


@format_doc(
    base_doc, components=doc_components_ecl.format("barycenter"), footer=doc_footer_bary
)
class BarycentricMeanEcliptic(BaseEclipticFrame):
    """
    Barycentric mean ecliptic coordinates.  These origin of the coordinates are the
    barycenter of the solar system, with the x axis pointing in the direction of
    the *mean* (not true) equinox as at the time specified by the ``equinox``
    attribute (as seen from Earth), and the xy-plane in the plane of the
    ecliptic for that date.

    The frame attributes are listed under **Other Parameters**.
    """

    equinox = TimeAttribute(default=EQUINOX_J2000, doc="The equinox time")


@format_doc(
    base_doc, components=doc_components_ecl.format("barycenter"), footer=doc_footer_bary
)
class BarycentricTrueEcliptic(BaseEclipticFrame):
    """
    Barycentric true ecliptic coordinates.  These origin of the coordinates are the
    barycenter of the solar system, with the x axis pointing in the direction of
    the *true* (not mean) equinox as at the time specified by the ``equinox``
    attribute (as seen from Earth), and the xy-plane in the plane of the
    ecliptic for that date.

    The frame attributes are listed under **Other Parameters**.
    """

    equinox = TimeAttribute(default=EQUINOX_J2000, doc="The equinox time")


doc_footer_helio = """
    Other parameters
    ----------------
    equinox : `~astropy.time.Time`, optional
        The date to assume for this frame.  Determines the location of the
        x-axis and the location of the Earth and Sun.
        Defaults to the 'J2000' equinox.
    obstime : `~astropy.time.Time`, optional
        The time at which the observation is taken.  Used for determining the
        position of the Sun. Defaults to J2000.
"""


@format_doc(
    base_doc,
    components=doc_components_ecl.format("sun's center"),
    footer=doc_footer_helio,
)
class HeliocentricMeanEcliptic(BaseEclipticFrame):
    """
    Heliocentric mean ecliptic coordinates.  These origin of the coordinates are the
    center of the sun, with the x axis pointing in the direction of
    the *mean* (not true) equinox as at the time specified by the ``equinox``
    attribute (as seen from Earth), and the xy-plane in the plane of the
    ecliptic for that date.

    The frame attributes are listed under **Other Parameters**.

    {params}


    """

    equinox = TimeAttribute(default=EQUINOX_J2000, doc="The equinox time")
    obstime = TimeAttribute(
        default=DEFAULT_OBSTIME, doc="The reference time (e.g., time of observation)"
    )


@format_doc(
    base_doc,
    components=doc_components_ecl.format("sun's center"),
    footer=doc_footer_helio,
)
class HeliocentricTrueEcliptic(BaseEclipticFrame):
    """
    Heliocentric true ecliptic coordinates.  These origin of the coordinates are the
    center of the sun, with the x axis pointing in the direction of
    the *true* (not mean) equinox as at the time specified by the ``equinox``
    attribute (as seen from Earth), and the xy-plane in the plane of the
    ecliptic for that date.

    The frame attributes are listed under **Other Parameters**.

    {params}


    """

    equinox = TimeAttribute(default=EQUINOX_J2000, doc="The equinox time")
    obstime = TimeAttribute(
        default=DEFAULT_OBSTIME, doc="The reference time (e.g., time of observation)"
    )


@format_doc(base_doc, components=doc_components_ecl.format("sun's center"), footer="")
class HeliocentricEclipticIAU76(BaseEclipticFrame):
    """
    Heliocentric mean (IAU 1976) ecliptic coordinates.  These origin of the coordinates are the
    center of the sun, with the x axis pointing in the direction of
    the *mean* (not true) equinox of J2000, and the xy-plane in the plane of the
    ecliptic of J2000 (according to the IAU 1976/1980 obliquity model).
    It has, therefore, a fixed equinox and an older obliquity value
    than the rest of the frames.

    The frame attributes are listed under **Other Parameters**.

    {params}


    """

    obstime = TimeAttribute(
        default=DEFAULT_OBSTIME, doc="The reference time (e.g., time of observation)"
    )


@format_doc(base_doc, components=doc_components_ecl.format("barycenter"), footer="")
class CustomBarycentricEcliptic(BaseEclipticFrame):
    """
    Barycentric ecliptic coordinates with custom obliquity.
    These origin of the coordinates are the
    barycenter of the solar system, with the x axis pointing in the direction of
    the *mean* (not true) equinox of J2000, and the xy-plane in the plane of the
    ecliptic tilted a custom obliquity angle.

    The frame attributes are listed under **Other Parameters**.
    """

    obliquity = QuantityAttribute(
        default=84381.448 * u.arcsec, unit=u.arcsec, doc="The obliquity of the ecliptic"
    )
