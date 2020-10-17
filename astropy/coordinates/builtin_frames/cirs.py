# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy import units as u
from astropy.utils.decorators import format_doc
from astropy.coordinates.attributes import (TimeAttribute,
                                            CartesianRepresentationAttribute)
from astropy.coordinates.baseframe import base_doc
from .baseradec import doc_components, BaseRADecFrame
from .utils import DEFAULT_OBSTIME

__all__ = ['CIRS']


doc_footer = """
    Other parameters
    ----------------
    obstime : `~astropy.time.Time`
        The time at which the observation is taken.  Used for determining the
        position of the Earth and its precession
    obsgeoloc : `~astropy.coordinates.CartesianRepresentation`, `~astropy.units.Quantity`
        The position of the observer relative to the center-of-mass of the
        Earth, oriented the same as GCRS. Either [0, 0, 0],
        `~astropy.coordinates.CartesianRepresentation`, or proper input for one,
        i.e., a `~astropy.units.Quantity` with shape (3, ...) and length units.
        Defaults to [0, 0, 0], meaning a geocentric observer.
    obsgeovel : `~astropy.coordinates.CartesianRepresentation`, `~astropy.units.Quantity`
        The velocity of the observer relative to the center-of-mass of the
        Earth, oriented the same as GCRS. Either [0, 0, 0],
        `~astropy.coordinates.CartesianRepresentation`, or proper input for one,
        i.e., a `~astropy.units.Quantity` with shape (3, ...) and velocity
        units.  Defaults to [0, 0, 0], meaning a geocentric observer.
"""


@format_doc(base_doc, components=doc_components, footer=doc_footer)
class CIRS(BaseRADecFrame):
    """
    A coordinate or frame in the Celestial Intermediate Reference System (CIRS).

    The frame attributes are listed under **Other Parameters**.
    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)
    obsgeoloc = CartesianRepresentationAttribute(default=[0, 0, 0],
                                                 unit=u.m)
    obsgeovel = CartesianRepresentationAttribute(default=[0, 0, 0],
                                                 unit=u.m/u.s)

# The "self-transform" is defined in icrs_cirs_transformations.py, because in
# the current implementation it goes through ICRS (like GCRS)
