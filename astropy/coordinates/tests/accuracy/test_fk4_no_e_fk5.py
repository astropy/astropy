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
from ....utils.data import get_pkg_data_fileobj

TOLERANCE = 0.03  # arcseconds


def test_fk4_no_e_fk5():
    with get_pkg_data_fileobj('fk4_no_e_fk5.csv') as f:
        t = Table.read(f, format='ascii')

    for i in range(len(t)):

        # Extract row
        r = t[i]

        # FK4NoETerms to FK5
        c1 = FK4NoETerms(ra=r['ra_in']*u.deg, dec=r['dec_in']*u.deg,
                         obstime=Time(r['obstime'], scale='utc'),
                         equinox=Time(r['equinox_fk4'], scale='utc'))
        c2 = c1.transform_to(FK5(equinox=Time(r['equinox_fk5'], scale='utc')))

        # Find difference
        diff = angular_separation(c2.ra.radian, c2.dec.radian,
                                  np.radians(r['ra_fk5']),
                                  np.radians(r['dec_fk5']))

        #assert np.degrees(diff) * 3600. < TOLERANCE

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

        assert np.degrees(diff) * 3600. < TOLERANCE
