# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os

import numpy as np

from .... import units as u
from ...builtin_frames import FK4NoETerms, FK5
from ....time import Time
from ....table import Table
from ...angle_utilities import angular_separation
from ....utils.data import get_pkg_data_contents

#the number of tests to run
from . import N_ACCURACY_TESTS

TOLERANCE = 0.03  # arcseconds


def test_fk4_no_e_fk5():
    lines = get_pkg_data_contents('fk4_no_e_fk5.csv').split('\n')
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

        # FK4NoETerms to FK5
        c1 = FK4NoETerms(ra=r['ra_in']*u.deg, dec=r['dec_in']*u.deg,
                         obstime=Time(r['obstime'], scale='utc'),
                         equinox=Time(r['equinox_fk4'], scale='utc'))
        c2 = c1.transform_to(FK5(equinox=Time(r['equinox_fk5'], scale='utc')))

        # Find difference
        diff = angular_separation(c2.ra.radian, c2.dec.radian,
                                  np.radians(r['ra_fk5']),
                                  np.radians(r['dec_fk5']))

        diffarcsec1.append(np.degrees(diff) * 3600.)

        # FK5 to FK4NoETerms
        c1 = FK5(ra=r['ra_in']*u.deg, dec=r['dec_in']*u.deg,
                 equinox=Time(r['equinox_fk5'], scale='utc'))
        fk4neframe = FK4NoETerms(obstime=Time(r['obstime'], scale='utc'),
                                 equinox=Time(r['equinox_fk4'], scale='utc'))
        c2 = c1.transform_to(fk4neframe)

        # Find difference
        diff = angular_separation(c2.ra.radian, c2.dec.radian,
                                  np.radians(r['ra_fk4']),
                                  np.radians(r['dec_fk4']))

        diffarcsec2.append(np.degrees(diff) * 3600.)

    np.testing.assert_array_less(diffarcsec1, TOLERANCE)
    np.testing.assert_array_less(diffarcsec2, TOLERANCE)
