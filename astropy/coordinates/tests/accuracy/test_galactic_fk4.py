# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os

import numpy as np

from .... import units as u
from ...builtin_frames import Galactic, FK4
from ....time import Time
from ....table import Table
from ...angle_utilities import angular_separation
from ....utils.data import get_pkg_data_fileobj

TOLERANCE = 0.5  # arcseconds


def test_galactic_fk4():
    with get_pkg_data_fileobj('galactic_fk4.csv') as f:
        t = Table.read(f, format='ascii')

    for i in range(len(t)):

        # Extract row
        r = t[i]

        # Galactic to FK4
        c1 = Galactic(l=r['lon_in']*u.deg, b=r['lat_in']*u.deg)
        c2 = c1.transform_to(FK4(equinox=Time(r['equinox_fk4'], scale='utc')))

        # Find difference
        diff = angular_separation(c2.ra.radian, c2.dec.radian,
                                  np.radians(r['ra_fk4']),
                                  np.radians(r['dec_fk4']))

        assert np.degrees(diff) * 3600. < TOLERANCE

        # FK4 to Galactic
        c1 = FK4(ra=r['lon_in']*u.deg, dec=r['lat_in']*u.deg,
                 obstime=Time(r['obstime'], scale='utc'),
                 equinox=Time(r['equinox_fk4'], scale='utc'))
        c2 = c1.transform_to(Galactic)

        # Find difference
        diff = angular_separation(c2.l.radian, c2.b.radian,
                                  np.radians(r['lon_gal']),
                                  np.radians(r['lat_gal']))

        assert np.degrees(diff) * 3600. < TOLERANCE
