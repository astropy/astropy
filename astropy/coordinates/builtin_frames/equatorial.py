# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Coordinate frames tied to the Equator and Equinox of Earth.

TEME is a True equator, Mean Equinox coordinate frame used in NORAD TLE
satellite files.

TETE is a True equator, True Equinox coordinate frame often called the
"apparent" coordinates. It is the same frame as used by JPL Horizons
and can be combined with Local Apparent Sidereal Time to calculate the
hour angle.
"""

from astropy.coordinates.attributes import EarthLocationAttribute, TimeAttribute
from astropy.coordinates.baseframe import BaseCoordinateFrame, base_doc
from astropy.coordinates.builtin_frames.baseradec import BaseRADecFrame, doc_components
from astropy.coordinates.representation import (
    CartesianDifferential,
    CartesianRepresentation,
)
from astropy.utils.decorators import format_doc

from .utils import DEFAULT_OBSTIME, EARTH_CENTER

__all__ = ["TEME", "TETE"]

doc_footer_teme = """
    Other parameters
    ----------------
    obstime : `~astropy.time.Time`
        The time at which the frame is defined.  Used for determining the
        position of the Earth.
"""

doc_footer_tete = """
    Other parameters
    ----------------
    obstime : `~astropy.time.Time`
        The time at which the observation is taken.  Used for determining the
        position of the Earth.
    location : `~astropy.coordinates.EarthLocation`
        The location on the Earth.  This can be specified either as an
        `~astropy.coordinates.EarthLocation` object or as anything that can be
        transformed to an `~astropy.coordinates.ITRS` frame. The default is the
        centre of the Earth.
"""


@format_doc(base_doc, components=doc_components, footer=doc_footer_tete)
class TETE(BaseRADecFrame):
    """
    An equatorial coordinate or frame using the True Equator and True Equinox (TETE).

    Equatorial coordinate frames measure RA with respect to the equinox and declination
    with with respect to the equator. The location of the equinox and equator vary due
    the gravitational torques on the oblate Earth. This variation is split into precession
    and nutation, although really they are two aspects of a single phenomena. The smooth,
    long term variation is known as precession, whilst smaller, periodic components are
    called nutation.

    Calculation of the true equator and equinox involves the application of both precession
    and nutation, whilst only applying precession gives a mean equator and equinox.

    TETE coordinates are often referred to as "apparent" coordinates, or
    "apparent place". TETE is the apparent coordinate system used by JPL Horizons
    and is the correct coordinate system to use when combining the right ascension
    with local apparent sidereal time to calculate the apparent (TIRS) hour angle.

    For more background on TETE, see the references provided in the
    :ref:`astropy:astropy-coordinates-seealso` section of the documentation.
    Of particular note are Sections 5 and 6 of
    `USNO Circular 179 <https://arxiv.org/abs/astro-ph/0602086>`_) and
    especially the diagram at the top of page 57.

    This frame also includes frames that are defined *relative* to the center of the Earth,
    but that are offset (in both position and velocity) from the center of the Earth. You
    may see such non-geocentric coordinates referred to as "topocentric".

    The frame attributes are listed under **Other Parameters**.
    """

    obstime = TimeAttribute(
        default=DEFAULT_OBSTIME, doc="The reference time (e.g., time of observation)"
    )
    location = EarthLocationAttribute(
        default=EARTH_CENTER, doc="The location on Earth of the observer"
    )


# Self transform goes through ICRS and is defined in icrs_cirs_transforms.py


@format_doc(base_doc, components="", footer=doc_footer_teme)
class TEME(BaseCoordinateFrame):
    """
    A coordinate or frame in the True Equator Mean Equinox frame (TEME).

    This frame is a geocentric system similar to CIRS or geocentric apparent place,
    except that the mean sidereal time is used to rotate from TIRS. TEME coordinates
    are most often used in combination with orbital data for satellites in the
    two-line-ephemeris format.

    Different implementations of the TEME frame exist. For clarity, this frame follows the
    conventions and relations to other frames that are set out in Vallado et al (2006).

    For more background on TEME, see the references provided in the
    :ref:`astropy:astropy-coordinates-seealso` section of the documentation.
    """

    default_representation = CartesianRepresentation
    default_differential = CartesianDifferential

    obstime = TimeAttribute(doc="The reference time (e.g., time of observation)")


# Transformation functions for getting to/from TEME and ITRS are in
# intermediate rotation transforms.py
