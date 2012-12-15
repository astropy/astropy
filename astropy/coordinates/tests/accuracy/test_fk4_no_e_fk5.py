import os

import numpy as np

from .... import units as u
from .... import coordinates as coord
from ....time import Time
from ....table import Table
from ...angle_utilities import vincenty_sphere_dist

TOLERANCE = 0.1  # arcseconds

ROOT = os.path.dirname(os.path.abspath(__file__))

def test_fk4_no_e_fk5():

    t = Table.read(os.path.join(ROOT, 'fk4_no_e_fk5.csv'), format='ascii')

    for i in range(len(t)):

        # Extract row
        r = t[i]

        # FK4 to FK5
        c1 = coord.FK4NoETermCoordinates(r['ra_in'], r['dec_in'],
                                         unit=(u.degree, u.degree),
                                         obstime=Time(r['obstime'], scale='utc'),
                                         equinox=Time(r['equinox_in'], scale='utc'))
        c2 = c1.transform_to(coord.FK5Coordinates).precess_to(Time(r['equinox_out'], scale='utc'))

        # Find difference
        diff = vincenty_sphere_dist(c2.dec.radians, c2.ra.radians,
                                    np.radians(r['dec_fk5']), np.radians(r['ra_fk5']))

        assert np.degrees(diff) * 3600. < TOLERANCE
