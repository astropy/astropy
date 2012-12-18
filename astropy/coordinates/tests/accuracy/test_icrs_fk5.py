import os

import numpy as np

from .... import units as u
from ... import ICRSCoordinates, FK5Coordinates
from ....time import Time
from ....table import Table
from ...angle_utilities import vincenty_sphere_dist

TOLERANCE = 0.025  # arcseconds

ROOT = os.path.dirname(os.path.abspath(__file__))


def test_icrs_no_e_fk5():

    t = Table.read(os.path.join(ROOT, 'icrs_fk5.csv'), format='ascii')

    for i in range(len(t)):

        # Extract row
        r = t[i]

        # FK4 to FK5
        c1 = ICRSCoordinates(r['ra_in'], r['dec_in'],
                                   unit=(u.degree, u.degree),
                                   obstime=Time(r['obstime'], scale='utc'))
        c2 = c1.transform_to(FK5Coordinates).precess_to(Time(r['equinox_fk5'], scale='utc'))

        # Find difference
        diff = vincenty_sphere_dist(c2.dec.radians, c2.ra.radians,
                                    np.radians(r['dec_fk5']), np.radians(r['ra_fk5']))

        assert np.degrees(diff) * 3600. < TOLERANCE

        # FK5 to FK4
        c1 = FK5Coordinates(r['ra_in'], r['dec_in'],
                            unit=(u.degree, u.degree),
                            obstime=Time(r['obstime'], scale='utc'),
                            equinox=Time(r['equinox_fk5'], scale='utc'))
        c2 = c1.transform_to(ICRSCoordinates)

        # Find difference
        diff = vincenty_sphere_dist(c2.dec.radians, c2.ra.radians,
                                    np.radians(r['dec_icrs']), np.radians(r['ra_icrs']))

        assert np.degrees(diff) * 3600. < TOLERANCE
