# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os

import numpy as np

from .... import units as u
from ... import Galactic, FK4
from ....time import Time
from ....table import Table
from ...angle_utilities import angular_separation

TOLERANCE = 0.5  # arcseconds

ROOT = os.path.dirname(os.path.abspath(__file__))


def test_galactic_fk4():

    t = Table.read(os.path.join(ROOT, 'galactic_fk4.csv'), format='ascii')

    for i in range(len(t)):

        # Extract row
        r = t[i]

        # FK4 to FK5
        c1 = Galactic(r['lon_in'], r['lat_in'],
                                 unit=(u.degree, u.degree),
                                 obstime=Time(r['obstime'], scale='utc'))
        c2 = c1.transform_to(FK4).precess_to(Time(r['equinox_fk4'], scale='utc'))

        # Find difference
        diff = angular_separation(c2.ra.radian, c2.dec.radian,
                                    np.radians(r['ra_fk4']), np.radians(r['dec_fk4']))

        assert np.degrees(diff) * 3600. < TOLERANCE

        # FK5 to FK4
        c1 = FK4(r['lon_in'], r['lat_in'],
                            unit=(u.degree, u.degree),
                            obstime=Time(r['obstime'], scale='utc'),
                            equinox=Time(r['equinox_fk4'], scale='utc'))
        c2 = c1.transform_to(Galactic)

        # Find difference
        diff = angular_separation(c2.l.radian, c2.b.radian,
                                    np.radians(r['lon_gal']), np.radians(r['lat_gal']))

        assert np.degrees(diff) * 3600. < TOLERANCE
