# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from ... import units as u
from ...utils.decorators import format_doc
from .. import representation as r
from ..baseframe import BaseCoordinateFrame, RepresentationMapping, base_doc
from ..attributes import (Attribute, TimeAttribute,
                          QuantityAttribute, EarthLocationAttribute)

__all__ = ['AltAz']


_90DEG = 90*u.deg

doc_components = """
    az : `Angle`, optional, must be keyword
        The Azimuth for this object (``alt`` must also be given and
        ``representation`` must be None).
    alt : `Angle`, optional, must be keyword
        The Altitude for this object (``az`` must also be given and
        ``representation`` must be None).
    distance : :class:`~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.

    pm_az_cosalt : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in azimuth (including the ``cos(alt)`` factor) for
        this object (``pm_alt`` must also be given).
    pm_alt : :class:`~astropy.units.Quantity`, optional, must be keyword
        The proper motion in altitude for this object (``pm_az_cosalt`` must
        also be given).
    radial_velocity : :class:`~astropy.units.Quantity`, optional, must be keyword
        The radial velocity of this object."""

doc_footer = """
    Other parameters
    ----------------
    obstime : `~astropy.time.Time`
        The time at which the observation is taken.  Used for determining the
        position and orientation of the Earth.
    location : `~astropy.coordinates.EarthLocation`
        The location on the Earth.  This can be specified either as an
        `~astropy.coordinates.EarthLocation` object or as anything that can be
        transformed to an `~astropy.coordinates.ITRS` frame.
    pressure : `~astropy.units.Quantity`
        The atmospheric pressure as an `~astropy.units.Quantity` with pressure
        units.  This is necessary for performing refraction corrections.
        Setting this to 0 (the default) will disable refraction calculations
        when transforming to/from this frame.
    temperature : `~astropy.units.Quantity`
        The ground-level temperature as an `~astropy.units.Quantity` in
        deg C.  This is necessary for performing refraction corrections.
    relative_humidity`` : `~astropy.units.Quantity` or number.
        The relative humidity as a dimensionless quantity between 0 to 1.
        This is necessary for performing refraction corrections.
    obswl : `~astropy.units.Quantity`
        The average wavelength of observations as an `~astropy.units.Quantity`
         with length units.  This is necessary for performing refraction
         corrections.

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

@format_doc(base_doc, components=doc_components, footer=doc_footer)
class AltAz(BaseCoordinateFrame):
    """
    A coordinate or frame in the Altitude-Azimuth system (Horizontal
    coordinates).  Azimuth is oriented East of North (i.e., N=0, E=90 degrees).

    This frame is assumed to *include* refraction effects if the ``pressure``
    frame attribute is non-zero.

    The frame attributes are listed under **Other Parameters**, which are
    necessary for transforming from AltAz to some other system.
    """

    frame_specific_representation_info = {
        r.SphericalRepresentation: [
            RepresentationMapping('lon', 'az'),
            RepresentationMapping('lat', 'alt')
        ]
    }

    default_representation = r.SphericalRepresentation
    default_differential = r.SphericalCosLatDifferential

    obstime = TimeAttribute(default=None)
    location = EarthLocationAttribute(default=None)
    pressure = QuantityAttribute(default=0, unit=u.hPa)
    temperature = QuantityAttribute(default=0, unit=u.deg_C)
    relative_humidity = QuantityAttribute(default=0, unit=u.dimensionless_unscaled)
    obswl = QuantityAttribute(default=1*u.micron, unit=u.micron)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def secz(self):
        """
        Secant if the zenith angle for this coordinate, a common estimate of the
        airmass.
        """
        return 1/np.sin(self.alt)

    @property
    def zen(self):
        """
        The zenith angle for this coordinate
        """
        return _90DEG.to(self.alt.unit) - self.alt


# self-transform defined in cirs_observed_transforms.py
