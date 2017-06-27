# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


import numpy as np

from .... import units as u
from ...builtin_frames import Galactic, FK4
from ....time import Time
from ....table import Table
from ...angle_utilities import angular_separation
from ....utils.data import get_pkg_data_contents
from ....extern.six.moves import range

# the number of tests to run
from . import N_ACCURACY_TESTS

TOLERANCE = 0.3  # arcseconds


def test_galactic_fk4():
    lines = get_pkg_data_contents('galactic_fk4.csv').split('\n')
    t = Table.read(lines, format='ascii', delimiter=',', guess=False)

    if N_ACCURACY_TESTS >= len(t):
        idxs = range(len(t))
    else:
        idxs = np.random.randint(len(t), size=N_ACCURACY_TESTS)

    diffarcsec1 = []
    diffarcsec2 = []
    for i in idxs:
        # Extract row
        r = t[int(i)]  # int here is to get around a py 3.x astropy.table bug

        # Galactic to FK4
        c1 = Galactic(l=r['lon_in']*u.deg, b=r['lat_in']*u.deg)
        c2 = c1.transform_to(FK4(equinox=Time(r['equinox_fk4'], scale='utc')))

        # Find difference
        diff = angular_separation(c2.ra.radian, c2.dec.radian,
                                  np.radians(r['ra_fk4']),
                                  np.radians(r['dec_fk4']))

        diffarcsec1.append(np.degrees(diff) * 3600.)

        # FK4 to Galactic
        c1 = FK4(ra=r['lon_in']*u.deg, dec=r['lat_in']*u.deg,
                 obstime=Time(r['obstime'], scale='utc'),
                 equinox=Time(r['equinox_fk4'], scale='utc'))
        c2 = c1.transform_to(Galactic)

        # Find difference
        diff = angular_separation(c2.l.radian, c2.b.radian,
                                  np.radians(r['lon_gal']),
                                  np.radians(r['lat_gal']))

        diffarcsec2.append(np.degrees(diff) * 3600.)

    np.testing.assert_array_less(diffarcsec1, TOLERANCE)
    np.testing.assert_array_less(diffarcsec2, TOLERANCE)
