import numpy as np
from .. import units as u


# What follows is a set of functions to convert between spherical and
# cartesian coordinate frames. The aim here is to be as performant as possible,
# so this includes routines that assume angles are in radians and ones that
# assume angles are in degrees, since using astropy.units can introduce
# overheads. For even better performance, the *_radian functions could be
# re-written in Cython, but it is important that all these functions work with
# n-d arrays.


def spherical_to_cartesian_radian(lon, lat, radius=1):
    cos_lat = np.cos(lat)
    x = radius * cos_lat * np.cos(lon)
    y = radius * cos_lat * np.sin(lon)
    z = radius * np.sin(lat)
    return x, y, z


def cartesian_to_spherical_radian(x, y, z):
    s = np.hypot(x, y)
    r = np.hypot(s, z)
    lon = np.arctan2(y, x)
    lat = np.arctan2(z, s)
    return lon, lat, r


def spherical_to_cartesian_degree(lon, lat, radius=1):
    return spherical_to_cartesian_radian(np.radians(lon), np.radians(lat), radius=radius)


def cartesian_to_spherical_degree(x, y, z):
    lon, lat, r = cartesian_to_spherical_radian(x, y, z)
    return np.degrees(lon), np.degrees(lat), r


def spherical_to_cartesian(lon, lat, radius=1 * u.one):
    lon = lon.to(u.rad).value
    lat = lat.to(u.rad).value
    x, y, z = spherical_to_cartesian_radian(lon, lat, radius.value)
    return x * radius.unit, y * radius.unit, z * radius.unit


def cartesian_to_spherical(x, y, z):
    xv = x.value
    yv = y.to(x.unit).value
    zv = z.to(x.unit).value
    lon, lat, r = cartesian_to_spherical_radian(xv, yv, zv)
    return lon * u.radian, lat * u.radian, r * x.unit
