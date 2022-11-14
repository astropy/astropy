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

    The frame attributes are listed under **Other Parameters**.
    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)
    location = EarthLocationAttribute(default=EARTH_CENTER)


# The "self-transform" is defined in icrs_cirs_transformations.py, because in
# the current implementation it goes through ICRS (like GCRS)
