# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from .. import units as u
from ..coordinates import CartesianPoints, Longitude, Latitude

try:
    # Not guaranteed available at setup time
    from . import erfa_time
except ImportError:
    if not _ASTROPY_SETUP_:
        raise

ELLIPSOIDS = {'WGS84': 1, 'GRS80': 2, 'WGS72': 3}


class EarthLocation(CartesianPoints):

    ellipsoid = 'WGS84'

    @classmethod
    def from_geodetic(cls, lon, lat, height, ellipsoid=None):
        if ellipsoid is None:
            ellipsoid = cls.ellipsoid
        if ellipsoid not in ELLIPSOIDS:
            raise ValueError('Ellipsoid {0} not among known ones ({1})'
                             .format(ellipsoid, ELLIPSOIDS.keys()))
        lon = Longitude(lon, u.degree, wrap_angle='180d')
        lat = Latitude(lat, u.degree)
        height = u.Quantity(height, u.meter)
        xyz = erfa_time.era_gd2gc(ELLIPSOIDS[ellipsoid],
                                  np.atleast_1d(lon.to(u.radian).value),
                                  np.atleast_1d(lat.to(u.radian).value),
                                  np.atleast_1d(height.value))
        self = cls(xyz.squeeze() * u.meter)
        self.ellipsoid = ellipsoid
        return self

    @property
    def geodetic(self):
        self_value = self.to(u.meter).value
        if self_value.ndim == 1:
            self_value = self_value.reshape(-1, 1)
        lon, lat, height = erfa_time.era_gc2gd(ELLIPSOIDS[self.ellipsoid],
                                               self_value)
        return (Longitude(lon.squeeze() * u.radian, u.degree,
                          wrap_angle='180d'),
                Latitude(lat.squeeze() * u.radian, u.degree),
                u.Quantity(height.squeeze(), u.meter))

    @property
    def longitude(self):
        return self.geodetic[0]

    @property
    def latitude(self):
        return self.geodetic[1]

    @property
    def height(self):
        return self.geodetic[2]
