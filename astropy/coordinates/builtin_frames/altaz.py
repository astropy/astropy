# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from astropy import units as u
from astropy.coordinates import representation as r
from astropy.coordinates.attributes import (
    EarthLocationAttribute,
    QuantityAttribute,
    TimeAttribute,
)
from astropy.coordinates.baseframe import (
    BaseCoordinateFrame,
    RepresentationMapping,
    base_doc,
)
from astropy.utils.decorators import format_doc

__all__ = ["AltAz"]


_90DEG = 90 * u.deg

doc_components = """
    az : `~astropy.coordinates.Angle`, optional, keyword-only
        The Azimuth for this object (``alt`` must also be given and
        ``representation`` must be None).
    alt : `~astropy.coordinates.Angle`, optional, keyword-only
        The Altitude for this object (``az`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity` ['length'], optional, keyword-only
        The Distance for this object along the line-of-sight.

    pm_az_cosalt : `~astropy.units.Quantity` ['angular speed'], optional, keyword-only
        The proper motion in azimuth (including the ``cos(alt)`` factor) for
        this object (``pm_alt`` must also be given).
    pm_alt : `~astropy.units.Quantity` ['angular speed'], optional, keyword-only
        The proper motion in altitude for this object (``pm_az_cosalt`` must
        also be given).
    radial_velocity : `~astropy.units.Quantity` ['speed'], optional, keyword-only
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
    pressure : `~astropy.units.Quantity` ['pressure']
        The atmospheric pressure as an `~astropy.units.Quantity` with pressure
        units.  This is necessary for performing refraction corrections.
        Setting this to 0 (the default) will disable refraction calculations
        when transforming to/from this frame.
    temperature : `~astropy.units.Quantity` ['temperature']
        The ground-level temperature as an `~astropy.units.Quantity` in
        deg C.  This is necessary for performing refraction corrections.
    relative_humidity : `~astropy.units.Quantity` ['dimensionless'] or number
        The relative humidity as a dimensionless quantity between 0 to 1.
        This is necessary for performing refraction corrections.
    obswl : `~astropy.units.Quantity` ['length']
        The average wavelength of observations as an `~astropy.units.Quantity`
         with length units.  This is necessary for performing refraction
         corrections.

    Notes
    -----
    The refraction model is based on that implemented in ERFA, which is fast
    but becomes inaccurate for altitudes below about 5 degrees.  Near and below
    altitudes of 0, it can even give meaningless answers, and in this case
    transforming to AltAz and back to another frame can give highly discrepant
    results.  For much better numerical stability, leave the ``pressure`` at
    ``0`` (the default), thereby disabling the refraction correction and
    yielding "topocentric" horizontal coordinates.
    """


@format_doc(base_doc, components=doc_components, footer=doc_footer)
class AltAz(BaseCoordinateFrame):
    """
    A coordinate or frame in the Altitude-Azimuth system (Horizontal
    coordinates) with respect to the WGS84 ellipsoid.  Azimuth is oriented
    East of North (i.e., N=0, E=90 degrees).  Altitude is also known as
    elevation angle, so this frame is also in the Azimuth-Elevation system.

    This frame is assumed to *include* refraction effects if the ``pressure``
    frame attribute is non-zero.

    The frame attributes are listed under **Other Parameters**, which are
    necessary for transforming from AltAz to some other system.
    """

    frame_specific_representation_info = {
        r.SphericalRepresentation: [
            RepresentationMapping("lon", "az"),
            RepresentationMapping("lat", "alt"),
        ]
    }

    default_representation = r.SphericalRepresentation
    default_differential = r.SphericalCosLatDifferential

    obstime = TimeAttribute(
        default=None, doc="The reference time (e.g., time of observation)"
    )
    location = EarthLocationAttribute(
        default=None, doc="The location on Earth of the observer"
    )
    pressure = QuantityAttribute(default=0, unit=u.hPa, doc="The atmospheric pressure")
    temperature = QuantityAttribute(
        default=0, unit=u.deg_C, doc="The ground-level temperature"
    )
    relative_humidity = QuantityAttribute(
        default=0, unit=u.dimensionless_unscaled, doc="The relative humidity"
    )
    obswl = QuantityAttribute(
        default=1 * u.micron,
        unit=u.micron,
        doc="The average wavelength of observations",
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def secz(self):
        """
        Secant of the zenith angle for this coordinate, a common estimate of
        the airmass.
        """
        return 1 / np.sin(self.alt)

    @property
    def zen(self):
        """
        The zenith angle (or zenith distance / co-altitude) for this coordinate.
        """
        return _90DEG.to(self.alt.unit) - self.alt


# self-transform defined in icrs_observed_transforms.py
