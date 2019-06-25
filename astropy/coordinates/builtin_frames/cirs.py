# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.decorators import format_doc
from astropy.coordinates.attributes import TimeAttribute
from astropy.coordinates.baseframe import base_doc
from .baseradec import doc_components, BaseRADecFrame
from .utils import DEFAULT_OBSTIME

__all__ = ['CIRS']


doc_footer = """
    Other parameters
    ----------------
    obstime : `~astropy.time.Time`
        The time at which the observation is taken.  Used for determining the
        position of the Earth and its precession.
"""


@format_doc(base_doc, components=doc_components, footer=doc_footer)
class CIRS(BaseRADecFrame):
    """
    A coordinate or frame in the Celestial Intermediate Reference System (CIRS).

    The frame attributes are listed under **Other Parameters**.
    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)


# The "self-transform" is defined in icrs_cirs_transformations.py, because in
# the current implementation it goes through ICRS (like GCRS)
