# Licensed under a 3-clause BSD style license - see LICENSE.rst

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

__all__ = ["HADec"]


doc_components = """
    ha : `~astropy.coordinates.Angle`, optional, keyword-only
        The Hour Angle for this object (``dec`` must also be given and
        ``representation`` must be None).
    dec : `~astropy.coordinates.Angle`, optional, keyword-only
        The Declination for this object (``ha`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity` ['length'], optional, keyword-only
        The Distance for this object along the line-of-sight.

    pm_ha_cosdec : `~astropy.units.Quantity` ['angular speed'], optional, keyword-only
        The proper motion in hour angle (including the ``cos(dec)`` factor) for
        this object (``pm_dec`` must also be given).
    pm_dec : `~astropy.units.Quantity` ['angular speed'], optional, keyword-only
        The proper motion in declination for this object (``pm_ha_cosdec`` must
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
    relative_humidity : `~astropy.units.Quantity` ['dimensionless'] or number.
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
    transforming to HADec and back to another frame can give highly discrepant
    results.  For much better numerical stability, leave the ``pressure`` at
    ``0`` (the default), thereby disabling the refraction correction and
    yielding "topocentric" equatorial coordinates.
    """


@format_doc(base_doc, components=doc_components, footer=doc_footer)
class HADec(BaseCoordinateFrame):
    """
    A coordinate or frame in the Hour Angle-Declination system (Equatorial
    coordinates) with respect to the WGS84 ellipsoid.  Hour Angle is oriented
    with respect to upper culmination such that the hour angle is negative to
    the East and positive to the West.

    This frame is assumed to *include* refraction effects if the ``pressure``
    frame attribute is non-zero.

    The frame attributes are listed under **Other Parameters**, which are
    necessary for transforming from HADec to some other system.
    """

    frame_specific_representation_info = {
        r.SphericalRepresentation: [
            RepresentationMapping("lon", "ha", u.hourangle),
            RepresentationMapping("lat", "dec"),
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
        if self.has_data:
            self._set_data_lon_wrap_angle(self.data)

    @staticmethod
    def _set_data_lon_wrap_angle(data):
        if hasattr(data, "lon"):
            data.lon.wrap_angle = 180.0 * u.deg
        return data

    def represent_as(self, base, s="base", in_frame_units=False):
        """
        Ensure the wrap angle for any spherical
        representations.
        """
        data = super().represent_as(base, s, in_frame_units=in_frame_units)
        self._set_data_lon_wrap_angle(data)
        return data


# self-transform defined in icrs_observed_transforms.py
