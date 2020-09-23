# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy import units as u
from astropy.utils.decorators import format_doc
from astropy.coordinates.attributes import (TimeAttribute,
                                            CartesianRepresentationAttribute)
from .utils import DEFAULT_OBSTIME, EQUINOX_J2000
from astropy.coordinates.baseframe import base_doc
from .baseradec import BaseRADecFrame, doc_components

__all__ = ['GCRS', 'PrecessedGeocentric']


doc_footer_gcrs = """
    Other parameters
    ----------------
    obstime : `~astropy.time.Time`
        The time at which the observation is taken.  Used for determining the
        position of the Earth.
    obsgeoloc : `~astropy.coordinates.CartesianRepresentation`, `~astropy.units.Quantity`
        The position of the observer relative to the center-of-mass of the
        Earth, oriented the same as BCRS/ICRS. Either [0, 0, 0],
        `~astropy.coordinates.CartesianRepresentation`, or proper input for one,
        i.e., a `~astropy.units.Quantity` with shape (3, ...) and length units.
        Defaults to [0, 0, 0], meaning "true" GCRS.
    obsgeovel : `~astropy.coordinates.CartesianRepresentation`, `~astropy.units.Quantity`
        The velocity of the observer relative to the center-of-mass of the
        Earth, oriented the same as BCRS/ICRS. Either [0, 0, 0],
        `~astropy.coordinates.CartesianRepresentation`, or proper input for one,
        i.e., a `~astropy.units.Quantity` with shape (3, ...) and velocity
        units.  Defaults to [0, 0, 0], meaning "true" GCRS.
"""


@format_doc(base_doc, components=doc_components, footer=doc_footer_gcrs)
class GCRS(BaseRADecFrame):
    """
    A coordinate or frame in the Geocentric Celestial Reference System (GCRS).

    GCRS is distinct form ICRS mainly in that it is relative to the Earth's
    center-of-mass rather than the solar system Barycenter.  That means this
    frame includes the effects of aberration (unlike ICRS). For more background
    on the GCRS, see the references provided in the
    :ref:`astropy-coordinates-seealso` section of the documentation. (Of
    particular note is Section 1.2 of
    `USNO Circular 179 <https://arxiv.org/abs/astro-ph/0602086>`_)

    This frame also includes frames that are defined *relative* to the Earth,
    but that are offset (in both position and velocity) from the Earth.

    The frame attributes are listed under **Other Parameters**.
    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)
    obsgeoloc = CartesianRepresentationAttribute(default=[0, 0, 0],
                                                 unit=u.m)
    obsgeovel = CartesianRepresentationAttribute(default=[0, 0, 0],
                                                 unit=u.m/u.s)


# The "self-transform" is defined in icrs_cirs_transformations.py, because in
# the current implementation it goes through ICRS (like CIRS)


doc_footer_prec_geo = """
    Other parameters
    ----------------
    equinox : `~astropy.time.Time`
        The (mean) equinox to precess the coordinates to.
    obstime : `~astropy.time.Time`
        The time at which the observation is taken.  Used for determining the
        position of the Earth.
    obsgeoloc : `~astropy.coordinates.CartesianRepresentation`, `~astropy.units.Quantity`
        The position of the observer relative to the center-of-mass of the
        Earth, oriented the same as BCRS/ICRS. Either [0, 0, 0],
        `~astropy.coordinates.CartesianRepresentation`, or proper input for one,
        i.e., a `~astropy.units.Quantity` with shape (3, ...) and length units.
        Defaults to [0, 0, 0], meaning "true" Geocentric.
    obsgeovel : `~astropy.coordinates.CartesianRepresentation`, `~astropy.units.Quantity`
        The velocity of the observer relative to the center-of-mass of the
        Earth, oriented the same as BCRS/ICRS. Either 0,
        `~astropy.coordinates.CartesianRepresentation`, or proper input for one,
        i.e., a `~astropy.units.Quantity` with shape (3, ...) and velocity
        units. Defaults to [0, 0, 0], meaning "true" Geocentric.
"""


@format_doc(base_doc, components=doc_components, footer=doc_footer_prec_geo)
class PrecessedGeocentric(BaseRADecFrame):
    """
    A coordinate frame defined in a similar manner as GCRS, but precessed to a
    requested (mean) equinox.  Note that this does *not* end up the same as
    regular GCRS even for J2000 equinox, because the GCRS orientation is fixed
    to that of ICRS, which is not quite the same as the dynamical J2000
    orientation.

    The frame attributes are listed under **Other Parameters**
    """

    equinox = TimeAttribute(default=EQUINOX_J2000)
    obstime = TimeAttribute(default=DEFAULT_OBSTIME)
    obsgeoloc = CartesianRepresentationAttribute(default=[0, 0, 0], unit=u.m)
    obsgeovel = CartesianRepresentationAttribute(default=[0, 0, 0], unit=u.m/u.s)
