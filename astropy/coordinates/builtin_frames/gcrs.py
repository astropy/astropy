# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u
from ..representation import SphericalRepresentation
from ..baseframe import (BaseCoordinateFrame, RepresentationMapping,
                         TimeFrameAttribute, QuantityFrameAttribute)
from .utils import DEFAULT_OBSTIME


class GCRS(BaseCoordinateFrame):
    """
    A coordinate or frame in the Geocentric Celestial Reference System (GCRS).

    GCRS is distinct form ICRS mainly in that it is relative to the Earth's
    center-of-mass rather than the solar system Barycenter.  That means this
    frame includes the effects of abberation (unlike ICRS). For more background
    on the GCRS, see the references provided in the :ref:`my-reference-label`
    section of the documentation. (Of particular note is Section 1.2 of
    `USNO Circular 179 <http://aa.usno.navy.mil/publications/docs/Circular_179.php>`_)

    This frame also includes frames that are defined *relative* to the Earth,
    but that are offset (in both position and velocity) from the Earth.



    This frame has three frame attributes:

    * ``obstime``
        The time at which the observation is taken.  Used for determining the
        position of the Earth.
    * ``obsgeoloc``
        3-vector giving the position of the observer relative to the
        center-of-mass of the Earth, oriented the same as BCRS/ICRS.
        Defaults to 0,  meaning "true" GCRS.
    * ``obsgeovel``
        3-vector giving the celocity of the observer relative to the
        center-of-mass of the Earth, oriented the same as BCRS/ICRS.
        Defaults to 0, meaning "true" GCRS.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    ra : `Angle`, optional, must be keyword
        The RA for this object (``dec`` must also be given and ``representation``
        must be None).
    dec : `Angle`, optional, must be keyword
        The Declination for this object (``ra`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).
    """

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'ra'),
                      RepresentationMapping('lat', 'dec')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation

    obstime = TimeFrameAttribute(default=DEFAULT_OBSTIME)
    obsgeoloc = QuantityFrameAttribute(default=0, unit=u.m, shape=(3,))
    obsgeovel = QuantityFrameAttribute(default=0, unit=u.m/u.s, shape=(3,))

#The "self-transform" is defined in icrs_cirs_transformations.py, because in
#the current implementation it goes through ICRS (like CIRS)
