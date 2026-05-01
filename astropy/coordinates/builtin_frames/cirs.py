# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.coordinates.attributes import EarthLocationAttribute, TimeAttribute
from astropy.coordinates.baseframe import base_doc
from astropy.utils.decorators import format_doc

from .baseradec import BaseRADecFrame, doc_components
from .utils import DEFAULT_OBSTIME, EARTH_CENTER

__all__ = ["CIRS"]


doc_footer = """
    Other parameters
    ----------------
    obstime : `~astropy.time.Time`
        The time at which the observation is taken.  Used for determining the
        position of the Earth and its precession.
    location : `~astropy.coordinates.EarthLocation`
        The location on the Earth.  This can be specified either as an
        `~astropy.coordinates.EarthLocation` object or as anything that can be
        transformed to an `~astropy.coordinates.ITRS` frame. The default is the
        centre of the Earth.
"""


@format_doc(base_doc, components=doc_components, footer=doc_footer)
class CIRS(BaseRADecFrame):
    """
    A coordinate or frame in the Celestial Intermediate Reference System (CIRS).

    CIRS is a geocentric reference system that serves as the *intermediate*
    step between the quasi-inertial `~astropy.coordinates.ICRS` (fixed to
    distant quasars) and the Earth-fixed `~astropy.coordinates.ITRS`.  It was
    introduced by the International Astronomical Union (IAU) 2000 and 2006
    resolutions to replace the older equinox-based "apparent place" systems
    with a cleaner, higher-precision framework.

    **Definition**

    CIRS is defined by two fundamental concepts:

    * **Celestial Intermediate Pole (CIP)** — the pole of the system.  Its
      motion in the `~astropy.coordinates.GCRS` is described by the IAU 2006
      precession model combined with the IAU 2000A nutation model, giving
      micro-arcsecond accuracy.

    * **Celestial Intermediate Origin (CIO)** — the origin of right ascension
      in CIRS.  Unlike the classical vernal equinox (which is defined by the
      intersection of the Earth's equator and the ecliptic and therefore
      "swings" with precession and nutation), the CIO is a *non-rotating
      origin*: it is defined so that it has no rotation component around the
      CIP.  This cleanly separates the Earth's spin (described by the
      **Earth Rotation Angle**, ERA) from its wobble (precession and nutation).

    **Relation to other frames**

    * **ICRS → CIRS**: applying the IAU precession-nutation model rotates the
      ICRS axes to the intermediate (CIP/CIO) axes.  Astropy performs this
      via the ERFA ``apco``/``atciqz`` routines.

    * **CIRS → ITRS**: a single rotation by the Earth Rotation Angle (ERA)
      around the CIP axis.  ERA is a linear function of UT1 and replaces the
      more complex "equation of the equinoxes" needed in the older
      equinox-based systems.

    * **CIRS vs. GCRS**: both share the same pole (CIP) and the same
      geocentric origin, but they differ in their right-ascension origin.
      GCRS uses the equinox; CIRS uses the CIO.  The angular offset between
      the two origins along the intermediate equator is called the
      *equation of the origins* (EO).

    **Astropy implementation**

    Astropy implements CIRS using the `ERFA
    <https://github.com/liberfa/erfa>`_ library (the open-source version of
    IAU SOFA), which provides the IAU 2006/2000A precession-nutation model.
    The default ``location`` is the geocentre; supply an
    `~astropy.coordinates.EarthLocation` to get the topocentric variant
    (which additionally accounts for diurnal aberration).

    For more background see Section 2.6 of `USNO Circular 179
    <https://arxiv.org/abs/astro-ph/0602086>`_ and the references in the
    :ref:`astropy:astropy-coordinates-seealso` section of the documentation.

    The frame attributes are listed under **Other Parameters**.
    """

    obstime = TimeAttribute(
        default=DEFAULT_OBSTIME, doc="The reference time (e.g., time of observation"
    )
    location = EarthLocationAttribute(
        default=EARTH_CENTER, doc="The location on Earth of the observer"
    )


# The "self-transform" is defined in icrs_cirs_transformations.py, because in
# the current implementation it goes through ICRS (like GCRS)
