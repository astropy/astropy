# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


import os

import numpy as np

from .... import units as u
from ... import FK4NoETerms, FK4
from ....time import Time
from ....table import Table
from ...angle_utilities import angular_separation

# It looks as though SLALIB, which AST relies on, assumes a simplified version
# of the e-terms corretion, so we have to up the tolerance a bit to get things
# to agree.
TOLERANCE = 1.e-5  # arcseconds

ROOT = os.path.dirname(os.path.abspath(__file__))


def test_fk4_no_e_fk5():

    t = Table.read(os.path.join(ROOT, 'fk4_no_e_fk4.csv'), format='ascii')

    for i in range(len(t)):

        # Extract row
        r = t[i]

        # FK4 to FK5
        c1 = FK4(r['ra_in'], r['dec_in'],
                            unit=(u.degree, u.degree),
                            obstime=Time(r['obstime'], scale='utc'))
        c2 = c1.transform_to(FK4NoETerms)

        # Find difference
        diff = angular_separation(c2.ra.radian, c2.dec.radian,
                                    np.radians(r['ra_fk4ne']), np.radians(r['dec_fk4ne']))

        assert np.degrees(diff) * 3600. < TOLERANCE

        # FK5 to FK4
        c1 = FK4NoETerms(r['ra_in'], r['dec_in'],
                                   unit=(u.degree, u.degree),
                                   obstime=Time(r['obstime'], scale='utc'))
        c2 = c1.transform_to(FK4)

        # Find difference
        diff = angular_separation(c2.ra.radian, c2.dec.radian,
                                    np.radians(r['ra_fk4']), np.radians(r['dec_fk4']))

        assert np.degrees(diff) * 3600. < TOLERANCE
