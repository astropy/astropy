# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

import numpy as np

from .... import units as u
from ... import FK4NoETerms, FK5
from ....time import Time
from ....table import Table
from ...angle_utilities import angular_separation

TOLERANCE = 0.03  # arcseconds

ROOT = os.path.dirname(os.path.abspath(__file__))


def test_fk4_no_e_fk5():

    t = Table.read(os.path.join(ROOT, 'fk4_no_e_fk5.csv'), format='ascii')

    for i in range(len(t)):

        # Extract row
        r = t[i]

        # FK4 to FK5
        c1 = FK4NoETerms(r['ra_in'], r['dec_in'],
                                   unit=(u.degree, u.degree),
                                   obstime=Time(r['obstime'], scale='utc'),
                                   equinox=Time(r['equinox_fk4'], scale='utc'))
        c2 = c1.transform_to(FK5).precess_to(Time(r['equinox_fk5'], scale='utc'))

        # Find difference
        diff = angular_separation(c2.ra.radian, c2.dec.radian,
                                    np.radians(r['ra_fk5']), np.radians(r['dec_fk5']))

        assert np.degrees(diff) * 3600. < TOLERANCE

        # FK5 to FK4
        c1 = FK5(r['ra_in'], r['dec_in'],
                            unit=(u.degree, u.degree),
                            obstime=Time(r['obstime'], scale='utc'),
                            equinox=Time(r['equinox_fk5'], scale='utc'))
        c2 = c1.transform_to(FK4NoETerms).precess_to(Time(r['equinox_fk4'], scale='utc'))

        # Find difference
        diff = angular_separation(c2.ra.radian, c2.dec.radian,
                                    np.radians(r['ra_fk4']), np.radians(r['dec_fk4']))

        assert np.degrees(diff) * 3600. < TOLERANCE
