# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.decorators import format_doc
from astropy.coordinates.attributes import TimeAttribute
from .utils import DEFAULT_OBSTIME
from astropy.coordinates.baseframe import base_doc
from .baseradec import BaseRADecFrame, doc_components

__all__ = ['HCRS']


doc_footer = """
    Other parameters
    ----------------
    obstime : `~astropy.time.Time`
        The time at which the observation is taken.  Used for determining the
        position of the Sun.
"""


@format_doc(base_doc, components=doc_components, footer=doc_footer)
class HCRS(BaseRADecFrame):
    """
    A coordinate or frame in a Heliocentric system, with axes aligned to ICRS.

    The ICRS has an origin at the Barycenter and axes which are fixed with
    respect to space.

    This coordinate system is distinct from ICRS mainly in that it is relative
    to the Sun's center-of-mass rather than the solar system Barycenter.
    In principle, therefore, this frame should include the effects of
    aberration (unlike ICRS), but this is not done, since they are very small,
    of the order of 8 milli-arcseconds.

    For more background on the ICRS and related coordinate transformations, see
    the references provided in the :ref:`astropy-coordinates-seealso` section of
    the documentation.

    The frame attributes are listed under **Other Parameters**.
    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

# Transformations are defined in icrs_circ_transforms.py
