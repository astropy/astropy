# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Mixin columns for use in ascii/tests/test_ecsv.py, fits/tests/test_connect.py,
and misc/tests/test_hdf5.py.

All columns should have length 2.
"""

import numpy as np

from astropy import coordinates, table, time
from astropy import units as u

el = coordinates.EarthLocation(x=[1, 2] * u.km, y=[3, 4] * u.km, z=[5, 6] * u.km)
sr = coordinates.SphericalRepresentation([0, 1] * u.deg, [2, 3] * u.deg, 1 * u.kpc)
cr = coordinates.CartesianRepresentation([0, 1] * u.pc, [4, 5] * u.pc, [8, 6] * u.pc)
sd = coordinates.SphericalCosLatDifferential(
    [0, 1] * u.mas / u.yr, [0, 1] * u.mas / u.yr, 10 * u.km / u.s
)
srd = coordinates.SphericalRepresentation(sr, differentials=sd)
sc = coordinates.SkyCoord(
    [1, 2], [3, 4], unit="deg,deg", frame="fk4", obstime="J1990.5"
)
scd = coordinates.SkyCoord(
    [1, 2], [3, 4], [5, 6], unit="deg,deg,m", frame="fk4", obstime=["J1990.5"] * 2
)
scdc = scd.copy()
scdc.representation_type = "cartesian"
scpm = coordinates.SkyCoord(
    [1, 2],
    [3, 4],
    [5, 6],
    unit="deg,deg,pc",
    pm_ra_cosdec=[7, 8] * u.mas / u.yr,
    pm_dec=[9, 10] * u.mas / u.yr,
)
scpmrv = coordinates.SkyCoord(
    [1, 2],
    [3, 4],
    [5, 6],
    unit="deg,deg,pc",
    pm_ra_cosdec=[7, 8] * u.mas / u.yr,
    pm_dec=[9, 10] * u.mas / u.yr,
    radial_velocity=[11, 12] * u.km / u.s,
)
scrv = coordinates.SkyCoord(
    [1, 2], [3, 4], [5, 6], unit="deg,deg,pc", radial_velocity=[11, 12] * u.km / u.s
)
# Use UnitSpherical and Radial for a regression test for gh-16998.
icrs = coordinates.ICRS(
    [1, 2] * u.deg,
    [3, 4] * u.deg,
    radial_velocity=[5, 6] * u.km / u.s,
    representation_type="unitspherical",
    differential_type="radial",
)
altaz = coordinates.AltAz(
    [1, 2] * u.deg,
    [3, 4] * u.deg,
    obstime=time.Time("2010-11-12T13:14:15"),
    location=el,
)
so = coordinates.SkyOffsetFrame(
    [-1, -2] * u.deg,
    [-3, -4] * u.deg,
    origin=altaz,
    rotation=[5, 6] * u.deg,
)
sond = coordinates.SkyOffsetFrame(
    origin=icrs,
    rotation=[5, 6] * u.deg,
)

tm = time.Time(
    [51000.5, 51001.5], format="mjd", scale="tai", precision=5, location=el[0]
)
tm2 = time.Time(tm, precision=3, format="iso")
tm3 = time.Time(tm, location=el)
tm3.info.serialize_method["ecsv"] = "jd1_jd2"
obj = table.Column([{"a": 1}, {"b": [2]}], dtype="object")
su = table.Column(
    [(1, (1.5, 1.6)), (2, (2.5, 2.6))],
    name="su",
    dtype=[("i", np.int64), ("f", [("p1", np.float64), ("p0", np.float64)])],
)
su2 = table.Column(
    [(["snake", "c"], [1.6, 1.5]), (["eal", "a"], [2.5, 2.6])],
    dtype=[("name", "U5", (2,)), ("f", "f8", (2,))],
)
stokes = coordinates.StokesCoord(["RR", "LL"])

# NOTE: for testing, the name of the column "x" for the
# Quantity is important since it tests the fix for #10215
# (namespace clash, where "x" clashes with "el.x").
mixin_cols = {
    "tm": tm,
    "tm2": tm2,
    "tm3": tm3,
    "dt": time.TimeDelta([1, 2] * u.day),
    "sc": sc,
    "scd": scd,
    "scdc": scdc,
    "scpm": scpm,
    "scpmrv": scpmrv,
    "scrv": scrv,
    "icrs": icrs,
    "altaz": altaz,
    "so": so,
    "sond": sond,
    "x": [1, 2] * u.m,
    "qdb": [10, 20] * u.dB(u.mW),
    "qdex": [4.5, 5.5] * u.dex(u.cm / u.s**2),
    "qmag": [21, 22] * u.ABmag,
    "lat": coordinates.Latitude([1, 2] * u.deg),
    "lon": coordinates.Longitude([1, 2] * u.deg, wrap_angle=180.0 * u.deg),
    "ang": coordinates.Angle([1, 2] * u.deg),
    "el": el,
    "sr": sr,
    "cr": cr,
    "sd": sd,
    "srd": srd,
    "nd": table.NdarrayMixin([1, 2]),
    "obj": obj,
    "su": su,
    "su2": su2,
    "stokes": stokes,
}
time_attrs = [
    "value",
    "shape",
    "format",
    "scale",
    "precision",
    "in_subfmt",
    "out_subfmt",
    "location",
]
compare_attrs = {
    "tm": time_attrs,
    "tm2": time_attrs,
    "tm3": time_attrs,
    "dt": ["shape", "value", "format", "scale"],
    "sc": ["ra", "dec", "representation_type", "frame.name"],
    "scd": ["ra", "dec", "distance", "representation_type", "frame.name"],
    "scdc": ["x", "y", "z", "representation_type", "frame.name"],
    "scpm": [
        "ra",
        "dec",
        "distance",
        "pm_ra_cosdec",
        "pm_dec",
        "representation_type",
        "frame.name",
    ],
    "scpmrv": [
        "ra",
        "dec",
        "distance",
        "pm_ra_cosdec",
        "pm_dec",
        "radial_velocity",
        "representation_type",
        "frame.name",
    ],
    "scrv": [
        "ra",
        "dec",
        "distance",
        "radial_velocity",
        "representation_type",
        "frame.name",
    ],
    "icrs": [
        "ra",
        "dec",
        "radial_velocity",
        "representation_type",
        "differential_type",
    ],
    "altaz": ["alt", "az", "obstime", "location", "representation_type"],
    "so": ["lon", "lat", "rotation", "origin"],
    "sond": ["rotation", "origin"],
    "x": ["value", "unit"],
    "qdb": ["value", "unit"],
    "qdex": ["value", "unit"],
    "qmag": ["value", "unit"],
    "lon": ["value", "unit", "wrap_angle"],
    "lat": ["value", "unit"],
    "ang": ["value", "unit"],
    "el": ["x", "y", "z", "ellipsoid"],
    "nd": ["data"],
    "sr": ["lon", "lat", "distance"],
    "cr": ["x", "y", "z"],
    "sd": ["d_lon_coslat", "d_lat", "d_distance"],
    "srd": [
        "lon",
        "lat",
        "distance",
        "differentials.s.d_lon_coslat",
        "differentials.s.d_lat",
        "differentials.s.d_distance",
    ],
    "obj": [],
    "su": ["i", "f.p0", "f.p1"],
    "su2": ["name", "f"],
    "stokes": ["value"],
}
non_trivial_names = {
    "cr": ["cr.x", "cr.y", "cr.z"],
    "dt": ["dt.jd1", "dt.jd2"],
    "el": ["el.x", "el.y", "el.z"],
    "sc": ["sc.ra", "sc.dec"],
    "scd": ["scd.ra", "scd.dec", "scd.distance", "scd.obstime.jd1", "scd.obstime.jd2"],
    "scdc": ["scdc.x", "scdc.y", "scdc.z", "scdc.obstime.jd1", "scdc.obstime.jd2"],
    "scfc": ["scdc.x", "scdc.y", "scdc.z", "scdc.obstime.jd1", "scdc.obstime.jd2"],
    "scpm": [
        "scpm.ra",
        "scpm.dec",
        "scpm.distance",
        "scpm.pm_ra_cosdec",
        "scpm.pm_dec",
    ],
    "scpmrv": [
        "scpmrv.ra",
        "scpmrv.dec",
        "scpmrv.distance",
        "scpmrv.pm_ra_cosdec",
        "scpmrv.pm_dec",
        "scpmrv.radial_velocity",
    ],
    "scrv": ["scrv.ra", "scrv.dec", "scrv.distance", "scrv.radial_velocity"],
    "icrs": ["icrs.ra", "icrs.dec", "icrs.radial_velocity"],
    "altaz": [
        "altaz.az",
        "altaz.alt",
        "altaz.location.x",
        "altaz.location.y",
        "altaz.location.z",
    ],
    "so": [
        "so.lon",
        "so.lat",
        "so.rotation",
        "so.origin.az",
        "so.origin.alt",
        "so.origin.location.x",
        "so.origin.location.y",
        "so.origin.location.z",
    ],
    "sond": [
        "sond.rotation",
        "sond.origin.ra",
        "sond.origin.dec",
        "sond.origin.radial_velocity",
    ],
    "sd": ["sd.d_lon_coslat", "sd.d_lat", "sd.d_distance"],
    "sr": ["sr.lon", "sr.lat", "sr.distance"],
    "srd": [
        "srd.lon",
        "srd.lat",
        "srd.distance",
        "srd.differentials.s.d_lon_coslat",
        "srd.differentials.s.d_lat",
        "srd.differentials.s.d_distance",
    ],
    "su": ["su.i", "su.f.p1", "su.f.p0"],
    "su2": ["su2.name", "su2.f"],
    "tm": ["tm.jd1", "tm.jd2"],
    "tm2": ["tm2.jd1", "tm2.jd2"],
    "tm3": ["tm3.jd1", "tm3.jd2", "tm3.location.x", "tm3.location.y", "tm3.location.z"],
}
serialized_names = {
    name: non_trivial_names.get(name, [name]) for name in sorted(mixin_cols)
}
