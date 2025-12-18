# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

import numpy as np
import pytest

from astropy import units as u
from astropy.coordinates import (AltAz, EarthLocation, HADec, ICRS, ITRS,
                                 CartesianRepresentation, UnitSphericalRepresentation)
from astropy.time import Time
from astropy.utils import iers
from astropy.utils.exceptions import AstropyUserWarning

iers.conf.auto_download = False


def _itrs_offset_from_enu(location, obstime, east, north, up):
    lon, lat, _ = location.to_geodetic("WGS84")
    lon = lon.to_value(u.rad)
    lat = lat.to_value(u.rad)

    sin_lon = np.sin(lon)
    cos_lon = np.cos(lon)
    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)

    matrix = np.array(
        [
            [-sin_lon, cos_lon, 0.0],
            [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat],
            [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat],
        ]
    )

    enu_vector = np.array([
        east.to_value(u.m),
        north.to_value(u.m),
        up.to_value(u.m),
    ])
    itrs_vector = matrix.T @ enu_vector

    offset = CartesianRepresentation(
        itrs_vector[0] * u.m,
        itrs_vector[1] * u.m,
        itrs_vector[2] * u.m,
    )

    return location.get_itrs(obstime=obstime).cartesian + offset


def _make_test_location():
    return EarthLocation(lat=35.0 * u.deg, lon=-105.0 * u.deg, height=2100 * u.m)


def test_round_trip_itrs_altaz_hadec():
    location = _make_test_location()
    obstime = Time("2024-03-21T12:00:00", scale="utc")

    offset = _itrs_offset_from_enu(
        location,
        obstime,
        3200 * u.m,
        -1800 * u.m,
        9700 * u.m,
    )
    target = ITRS(offset, obstime=obstime)

    altaz_frame = AltAz(location=location, obstime=obstime)
    altaz = target.transform_to(altaz_frame)
    back_from_altaz = altaz.transform_to(ITRS(obstime=obstime))

    had_frame = HADec(location=location, obstime=obstime)
    hadec = target.transform_to(had_frame)
    back_from_hadec = hadec.transform_to(ITRS(obstime=obstime))

    diff_altaz = (back_from_altaz.cartesian - target.cartesian).norm().to(u.mm)
    diff_hadec = (back_from_hadec.cartesian - target.cartesian).norm().to(u.mm)

    assert diff_altaz.value < 1.0
    assert diff_hadec.value < 1.0

    top_range = np.sqrt((3200 * u.m) ** 2 + (-1800 * u.m) ** 2 + (9700 * u.m) ** 2)
    assert altaz.distance.to(u.m).value == pytest.approx(top_range.to_value(u.m), rel=1e-9)
    assert hadec.distance.to(u.m).value == pytest.approx(top_range.to_value(u.m), rel=1e-9)


def test_overhead_target_expectations():
    location = _make_test_location()
    obstime = Time("2024-06-01T09:15:00", scale="utc")

    offset = _itrs_offset_from_enu(location, obstime, 0 * u.m, 0 * u.m, 1500 * u.m)
    target = ITRS(offset, obstime=obstime)

    altaz = target.transform_to(AltAz(location=location, obstime=obstime))
    hadec = target.transform_to(HADec(location=location, obstime=obstime))

    assert altaz.alt.to(u.deg).value == pytest.approx(90.0, abs=1e-8)

    latitude = location.to_geodetic("WGS84")[1]
    assert hadec.ha.wrap_at(180 * u.deg).to(u.arcsec).value == pytest.approx(0.0, abs=1e-6)
    assert hadec.dec.to(u.arcsec).value == pytest.approx(latitude.to(u.arcsec).value, abs=1e-6)


def test_location_shape_mismatch_raises():
    itrs = ITRS(
        CartesianRepresentation(
            [[6.5e6, 0, 0], [6.6e6, 1000, 1000]] * u.m,
            xyz_axis=-1,
        ),
        obstime=Time("2024-01-01T00:00:00"),
    )
    locations = EarthLocation(
        lon=[0, 90, 180] * u.deg,
        lat=[0, 10, 20] * u.deg,
        height=0 * u.m,
    )

    with pytest.raises(ValueError, match=r"EarthLocation shape \(3,\) is not broadcastable to coordinate shape \(2,\)"):
        itrs.transform_to(AltAz(location=locations))


def test_broadcast_array_obstime():
    location = _make_test_location()
    obstimes = Time(["2024-01-01T00:00:00", "2024-01-01T00:05:00"])  # noqa: E501

    offsets = [
        _itrs_offset_from_enu(location, t, 1000 * u.m, 200 * u.m, 800 * u.m)
        for t in obstimes
    ]
    stacked = CartesianRepresentation(
        np.stack([off.xyz.to_value(u.m) for off in offsets], axis=0) * u.m,
        xyz_axis=-1,
    )

    itrs = ITRS(stacked, obstime=obstimes)
    altaz = itrs.transform_to(AltAz(location=location, obstime=obstimes))
    assert altaz.shape == (2,)


def test_unitspherical_handling():
    location = _make_test_location()
    obstime = Time("2024-04-04T00:00:00")

    itrs = ITRS(
        UnitSphericalRepresentation(lon=45 * u.deg, lat=10 * u.deg),
        obstime=obstime,
    )
    with pytest.warns(AstropyUserWarning, match="ITRS input has no distance"):
        altaz = itrs.transform_to(AltAz(location=location, obstime=obstime))
    assert isinstance(altaz.data, UnitSphericalRepresentation)

    altaz = AltAz(az=10 * u.deg, alt=30 * u.deg, location=location, obstime=obstime)
    with pytest.raises(ValueError, match="Observed->ITRS requires a distance"):
        altaz.transform_to(ITRS(obstime=obstime))


def test_refraction_warning_triggered():
    location = _make_test_location()
    obstime = Time("2024-02-02T12:30:00")

    offset = _itrs_offset_from_enu(location, obstime, 2000 * u.m, 0 * u.m, 6000 * u.m)
    target = ITRS(offset, obstime=obstime)

    frame = AltAz(location=location, obstime=obstime, pressure=1013 * u.hPa)
    with pytest.warns(AstropyUserWarning, match="ignores atmospheric refraction"):
        target.transform_to(frame)


def test_refraction_warning_observed_to_itrs():
    location = _make_test_location()
    obstime = Time("2024-02-02T12:30:00")
    altaz = AltAz(
        az=40 * u.deg,
        alt=20 * u.deg,
        distance=12 * u.km,
        location=location,
        obstime=obstime,
        pressure=600 * u.hPa,
    )

    with pytest.warns(AstropyUserWarning, match="ignores atmospheric refraction"):
        altaz.transform_to(ITRS(obstime=obstime))


def test_obstime_mismatch_warning():
    location = _make_test_location()
    obstime_src = Time("2024-05-01T00:00:00")
    obstime_tgt = Time("2024-05-01T01:00:00")

    offset = _itrs_offset_from_enu(location, obstime_src, 1500 * u.m, 500 * u.m, 4000 * u.m)
    itrs = ITRS(offset, obstime=obstime_src)

    with pytest.warns(AstropyUserWarning, match="obstime mismatch detected"):
        itrs.transform_to(AltAz(location=location, obstime=obstime_tgt))

    altaz = AltAz(az=20 * u.deg, alt=50 * u.deg, distance=3 * u.km, location=location, obstime=obstime_src)
    with pytest.warns(AstropyUserWarning, match="obstime mismatch detected"):
        altaz.transform_to(ITRS(obstime=obstime_tgt))


def test_direct_vs_icrs_routes_differ():
    location = _make_test_location()
    obstime = Time("2024-07-01T22:00:00")

    offset = _itrs_offset_from_enu(location, obstime, 5000 * u.m, 2500 * u.m, 7500 * u.m)
    target = ITRS(offset, obstime=obstime)

    frame = AltAz(location=location, obstime=obstime)
    direct = target.transform_to(frame)

    via_icrs = target.transform_to(ICRS()).transform_to(frame)

    delta_alt = np.abs((direct.alt - via_icrs.alt).to(u.arcsec).value)
    delta_az = np.abs((direct.az - via_icrs.az).to(u.arcsec).value)

    assert delta_alt > 1e-3  # difference is measurable in arcsec
    assert delta_az > 1e-3


def test_missing_location_raises():
    itrs = ITRS(x=1 * u.km, y=2 * u.km, z=3 * u.km, obstime=Time("2024-01-01"))
    with pytest.raises(ValueError, match="requires a valid EarthLocation"):
        itrs.transform_to(AltAz())
