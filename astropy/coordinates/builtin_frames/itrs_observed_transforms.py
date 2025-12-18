# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Direct transforms between ITRS and observed frames."""

from __future__ import annotations

import numpy as np
import warnings

from astropy import units as u
from astropy.coordinates import Latitude, Longitude
from astropy.coordinates import representation as rep
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.representation import (CartesianRepresentation,
                                                SphericalRepresentation,
                                                UnitSphericalRepresentation)
from astropy.coordinates.transformations import FunctionTransformWithFiniteDifference
from astropy.utils.exceptions import AstropyUserWarning

from .altaz import AltAz
from .hadec import HADec
from .itrs import ITRS

__all__ = ()

_MISSING_LOCATION_ERROR = (
    "Observed frame requires a valid EarthLocation for direct ITRS<->Observed transforms. "
    "Set location=EarthLocation(...)."
)

_UNIT_SPHERICAL_WARNING = (
    "ITRS input has no distance: assuming infinite distance (ignoring topocentric parallax) "
    "for Direct ITRS->Observed."
)

_REFRACTION_WARNING = (
    "Direct ITRS<->Observed ignores atmospheric refraction; set pressure=0 or use "
    "ITRS->ICRS->Observed for ERFA-based refraction."
)

_OBSTIME_WARNING = (
    "Direct ITRS<->Observed transformation assumes static ITRS coordinates; obstime mismatch "
    "detected (from {t_in} to {t_out}). Geometry is computed purely in ITRS without aberration "
    "or time-dependent ITRS shifts. For aberration-aware transforms use ITRS->ICRS->Observed."
)

_OBSERVED_TO_ITRS_DISTANCE_ERROR = (
    "Observed->ITRS requires a distance to compute absolute ITRS coordinates. Provide distance "
    "or use geoid_fallback='wgs84'."
)


def _warn_if_obstime_mismatch(source_time, target_time):
    if source_time is None or target_time is None:
        return

    try:
        equal = np.all(source_time == target_time)
    except Exception:  # pragma: no cover - defensive against unexpected shapes
        equal = False

    if not equal:
        warnings.warn(
            _OBSTIME_WARNING.format(t_in=source_time.iso, t_out=target_time.iso),
            AstropyUserWarning,
            stacklevel=3,
        )


def _warn_if_pressure(frame):
    pressure = getattr(frame, "pressure", None)
    if pressure is None:
        return

    try:
        pressure_values = pressure.to_value(u.hPa)
    except Exception:
        return

    if np.any(pressure_values != 0):
        warnings.warn(_REFRACTION_WARNING, AstropyUserWarning, stacklevel=3)


def _ensure_location(location, target_shape):
    if location is None:
        raise ValueError(_MISSING_LOCATION_ERROR)

    loc_shape = location.shape
    try:
        broadcast_shape = np.broadcast_shapes(target_shape, loc_shape)
    except ValueError as exc:
        raise ValueError(
            f"EarthLocation shape {loc_shape} is not broadcastable to coordinate shape {target_shape}."
        ) from exc

    return broadcast_shape


def _location_to_itrs_and_rotation(location, obstime, target_shape):
    broadcast_shape = _ensure_location(location, target_shape)
    observer_itrs = location.get_itrs(obstime=obstime)
    if observer_itrs.shape != broadcast_shape:
        observer_itrs = observer_itrs.broadcast_to(broadcast_shape)

    lon, lat, _ = location.to_geodetic("WGS84")
    lon_val, lat_val = np.broadcast_arrays(
        lon.to_value(u.rad),
        lat.to_value(u.rad),
        subok=True,
    )
    if lon_val.shape != broadcast_shape:
        lon_val = np.broadcast_to(lon_val, broadcast_shape)
        lat_val = np.broadcast_to(lat_val, broadcast_shape)

    rotation = _rotation_itrs_to_local(lon_val, lat_val)

    return observer_itrs.cartesian, rotation, lat_val, broadcast_shape


def _rotation_itrs_to_local(lon_rad, lat_rad):
    sin_lon = np.sin(lon_rad)
    cos_lon = np.cos(lon_rad)
    sin_lat = np.sin(lat_rad)
    cos_lat = np.cos(lat_rad)

    shape = lon_rad.shape + (3, 3)
    matrix = np.empty(shape, dtype=float)

    matrix[..., 0, 0] = -sin_lon
    matrix[..., 0, 1] = cos_lon
    matrix[..., 0, 2] = 0.0

    matrix[..., 1, 0] = -sin_lat * cos_lon
    matrix[..., 1, 1] = -sin_lat * sin_lon
    matrix[..., 1, 2] = cos_lat

    matrix[..., 2, 0] = cos_lat * cos_lon
    matrix[..., 2, 1] = cos_lat * sin_lon
    matrix[..., 2, 2] = sin_lat

    return matrix


def _is_direction_only_itrs(itrs_coo):
    data = itrs_coo.data
    if isinstance(data, UnitSphericalRepresentation):
        return True
    cart = itrs_coo.cartesian
    return cart.x.unit == u.one and cart.y.unit == u.one and cart.z.unit == u.one


def _cartesian_to_local_components(cart, rotation):
    unit = cart.x.unit
    xyz = np.moveaxis(cart.xyz.to_value(unit), 0, -1)
    local = np.einsum("...ij,...j->...i", rotation, xyz)
    return local[..., 0] * unit, local[..., 1] * unit, local[..., 2] * unit


def _local_components_to_cartesian(east, north, up):
    return CartesianRepresentation(x=east, y=north, z=up)


def _az_alt_from_local(east, north, up):
    rho = np.sqrt(east**2 + north**2 + up**2)
    alt = Latitude(np.arcsin(np.clip((up / rho).to_value(u.one), -1.0, 1.0)), unit=u.rad)
    unit = east.unit
    east_val = east.to_value(unit)
    north_val = north.to_value(unit)
    az = Longitude(np.arctan2(east_val, north_val), unit=u.rad).wrap_at(2 * np.pi * u.rad)
    return az, alt, rho


def _ha_dec_from_altaz(az, alt, lat_rad):
    alt_val = alt.to_value(u.rad)
    az_val = az.to_value(u.rad)

    sin_alt = np.sin(alt_val)
    cos_alt = np.cos(alt_val)
    sin_az = np.sin(az_val)
    cos_az = np.cos(az_val)

    sin_lat = np.sin(lat_rad)
    cos_lat = np.cos(lat_rad)

    sin_dec = sin_alt * sin_lat + cos_alt * cos_lat * cos_az
    sin_dec = np.clip(sin_dec, -1.0, 1.0)
    dec = Latitude(np.arcsin(sin_dec), unit=u.rad)

    cos_dec = np.cos(dec.to_value(u.rad))
    eps = 1e-15

    sin_ha_num = -cos_alt * sin_az
    cos_ha_num = sin_alt * cos_lat - cos_alt * sin_lat * cos_az

    sin_ha = np.divide(sin_ha_num, cos_dec, out=np.zeros_like(sin_ha_num), where=np.abs(cos_dec) > eps)
    cos_ha = np.divide(cos_ha_num, cos_dec, out=np.ones_like(cos_ha_num), where=np.abs(cos_dec) > eps)

    ha = Longitude(np.arctan2(sin_ha, cos_ha), unit=u.rad).wrap_at(np.pi * u.rad)

    return ha, dec


def _enu_from_altaz(az, alt, distance):
    alt_val = alt.to_value(u.rad)
    az_val = az.to_value(u.rad)

    cos_alt = np.cos(alt_val)
    sin_alt = np.sin(alt_val)
    sin_az = np.sin(az_val)
    cos_az = np.cos(az_val)

    east = distance * cos_alt * sin_az
    north = distance * cos_alt * cos_az
    up = distance * sin_alt
    return east, north, up


def _enu_from_hadec(ha, dec, lat_rad, distance):
    ha_val = ha.to_value(u.rad)
    dec_val = dec.to_value(u.rad)

    sin_dec = np.sin(dec_val)
    cos_dec = np.cos(dec_val)
    sin_ha = np.sin(ha_val)
    cos_ha = np.cos(ha_val)

    sin_lat = np.sin(lat_rad)
    cos_lat = np.cos(lat_rad)

    east = -distance * cos_dec * sin_ha
    north = distance * (sin_dec * cos_lat - cos_dec * cos_ha * sin_lat)
    up = distance * (sin_dec * sin_lat + cos_dec * cos_ha * cos_lat)
    return east, north, up


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, AltAz)
@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, ITRS, HADec)
def itrs_to_observed(itrs_coo, observed_frame):
    _warn_if_pressure(observed_frame)
    _warn_if_obstime_mismatch(getattr(itrs_coo, "obstime", None), observed_frame.obstime)

    observer_cart, rotation, lat_rad, broadcast_shape = _location_to_itrs_and_rotation(
        observed_frame.location,
        observed_frame.obstime,
        itrs_coo.shape,
    )

    itrs_cart = itrs_coo.cartesian
    if itrs_cart.shape != broadcast_shape:
        itrs_cart = itrs_cart.broadcast_to(broadcast_shape)

    direction_only = _is_direction_only_itrs(itrs_coo)

    if not direction_only:
        if observer_cart.shape != broadcast_shape:
            observer_cart = observer_cart.broadcast_to(broadcast_shape)
        topo_cart = itrs_cart - observer_cart
    else:
        topo_cart = itrs_cart

    east, north, up = _cartesian_to_local_components(topo_cart, rotation)
    az, alt, rho = _az_alt_from_local(east, north, up)

    if direction_only:
        warnings.warn(_UNIT_SPHERICAL_WARNING, AstropyUserWarning, stacklevel=3)
        if isinstance(observed_frame, AltAz):
            rep_out = UnitSphericalRepresentation(az, alt)
        else:
            ha, dec = _ha_dec_from_altaz(az, alt, lat_rad)
            rep_out = UnitSphericalRepresentation(ha, dec)
    else:
        if isinstance(observed_frame, AltAz):
            rep_out = SphericalRepresentation(lon=az, lat=alt, distance=rho)
        else:
            ha, dec = _ha_dec_from_altaz(az, alt, lat_rad)
            rep_out = SphericalRepresentation(lon=ha, lat=dec, distance=rho)

    return observed_frame.realize_frame(rep_out)


@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, AltAz, ITRS)
@frame_transform_graph.transform(FunctionTransformWithFiniteDifference, HADec, ITRS)
def observed_to_itrs(observed_coo, itrs_frame):
    _warn_if_pressure(observed_coo)
    _warn_if_obstime_mismatch(observed_coo.obstime, itrs_frame.obstime)

    distance = observed_coo.distance
    if distance is None or not distance.unit.is_equivalent(u.m):
        raise ValueError(_OBSERVED_TO_ITRS_DISTANCE_ERROR)

    observer_cart, rotation, lat_rad, _ = _location_to_itrs_and_rotation(
        observed_coo.location,
        observed_coo.obstime,
        observed_coo.shape,
    )

    rotation_T = np.swapaxes(rotation, -1, -2)

    if isinstance(observed_coo, AltAz):
        sph = observed_coo.represent_as(rep.SphericalRepresentation)
        az = sph.lon.to(u.rad)
        alt = sph.lat.to(u.rad)
        east, north, up = _enu_from_altaz(az, alt, distance)
    else:
        sph = observed_coo.represent_as(rep.SphericalRepresentation)
        ha = sph.lon.to(u.rad)
        dec = sph.lat.to(u.rad)
        east, north, up = _enu_from_hadec(ha, dec, lat_rad, distance)

    enu_cart = _local_components_to_cartesian(east, north, up)
    topo_itrs = enu_cart.transform(rotation_T)

    if observer_cart.shape != topo_itrs.shape:
        observer_cart = observer_cart.broadcast_to(topo_itrs.shape)

    geo_cart = topo_itrs + observer_cart
    return itrs_frame.realize_frame(geo_cart)
