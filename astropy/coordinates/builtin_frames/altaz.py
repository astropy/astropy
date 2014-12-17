# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ... import units as u
from ..representation import SphericalRepresentation
from ..baseframe import (BaseCoordinateFrame, FrameAttribute,
                         TimeFrameAttribute, QuantityFrameAttribute,
                         RepresentationMapping, EarthLocationAttribute)


class AltAz(BaseCoordinateFrame):
    """
    A coordinate or frame in the Altitude-Azimuth system (i.e., Horizontal
    coordinates).  This frame is assumed to *include* refraction effects if
    the ``pressure`` frame attribute is non-zero.

    This frame has the following frame attributes, which are necessary for
    transforming from AltAz to some other system:

    * ``obstime``
        The time at which the observation is taken.  Used for determining the
        position and orientation of the Earth.
    * ``location``
        The location on the Earth.  This can be specified either as an
        `~astropy.coordinates.EarthLocation` object or as anything that can be
        transformed to an `~astropy.coordinates.ITRS` frame.
    * ``pressure``
        The atmospheric pressure as an `~astropy.units.Quantity` with pressure
        units.  This is necessary for performing refraction corrections.
        Setting this to 0 (the default) will disable refraction calculations
        when transforming to/from this frame.
    * ``temperature``
        The ground-level temperature as an `~astropy.units.Quantity` in
        deg C.  This is necessary for performing refraction corrections.
    * ``relative_humidity``
        The relative humidity as a number from 0 to 1.  This is necessary for
        performing refraction corrections.
    * ``obswl``
        The ground-level temperature as an `~astropy.units.Quantity` with length
        units.  This is necessary for performing refraction corrections.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    az : `Angle`, optional, must be keyword
        The Azimuth for this object (``alt`` must also be given and
        ``representation`` must be None).
    alt : `Angle`, optional, must be keyword
        The Altitude for this object (``az`` must also be given and
        ``representation`` must be None).
    distance : :class:`~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.

    Notes
    -----
    The refraction model is based on that implemented in ERFA, which is fast
    but becomes inaccurate for altitudes below about 5 degrees.  Near and below
    altitudes of 0, it can even give meaningless answers, and in this case
    transforming to AltAz and back to another frame can give highly discrepent
    results.  For much better numerical stability, leaving the ``pressure`` at
    ``0`` (the default), disabling the refraction correction (yielding
    "topocentric" horizontal coordinates).

    """

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'az'),
                      RepresentationMapping('lat', 'alt')],
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation

    obstime = TimeFrameAttribute(default=None)
    location = EarthLocationAttribute(default=None)
    pressure = QuantityFrameAttribute(default=0, unit=u.hPa)
    temperature = QuantityFrameAttribute(default=0, unit=u.deg_C)
    relative_humidity = FrameAttribute(default=0)
    obswl = QuantityFrameAttribute(default=1*u.micron, unit=u.micron)

    def __init__(self, *args, **kwargs):
        super(AltAz, self).__init__(*args, **kwargs)

#self-transform defined in cirs_observed_transforms.py
