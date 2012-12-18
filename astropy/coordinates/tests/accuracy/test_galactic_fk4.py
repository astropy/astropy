import os

import numpy as np

from .... import units as u
from ... import GalacticCoordinates, FK4Coordinates
from ....time import Time
from ....table import Table
from ...angle_utilities import vincenty_sphere_dist

TOLERANCE = 1.  # arcseconds

ROOT = os.path.dirname(os.path.abspath(__file__))


def test_gal_no_e_fk4():

    t = Table.read(os.path.join(ROOT, 'galactic_fk4.csv'), format='ascii')

    for i in range(len(t)):

        # Extract row
        r = t[i]

        # FK4 to FK5
        c1 = GalacticCoordinates(r['lon_in'], r['lat_in'],
                                 unit=(u.degree, u.degree),
                                 obstime=Time(r['obstime'], scale='utc'))
        c2 = c1.transform_to(FK4Coordinates).precess_to(Time(r['equinox_fk4'], scale='utc'))

        # Find difference
        diff = vincenty_sphere_dist(c2.dec.radians, c2.ra.radians,
                                    np.radians(r['dec_fk4']), np.radians(r['ra_fk4']))

        assert np.degrees(diff) * 3600. < TOLERANCE

        # FK5 to FK4
        c1 = FK4Coordinates(r['lon_in'], r['lat_in'],
                            unit=(u.degree, u.degree),
                            obstime=Time(r['obstime'], scale='utc'),
                            equinox=Time(r['equinox_fk4'], scale='utc'))
        c2 = c1.transform_to(GalacticCoordinates)

        # Find difference
        diff = vincenty_sphere_dist(c2.l.radians, c2.b.radians,
                                    np.radians(r['lon_gal']), np.radians(r['lat_gal']))

        assert np.degrees(diff) * 3600. < TOLERANCE
