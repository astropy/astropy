# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools
import numpy as np

from .. import iers
from ....table import Table

allclose_sec = functools.partial(np.allclose, rtol=1e-15, atol=1e-9)
# 1 nanosec atol


class TestBasic():
    """Basic tests that IERS_B returns correct values"""

    def test_simple(self):
        iers.IERS_B.close()
        assert iers.IERS_B.iers_table is None
        iers_b = iers.IERS_B.open()
        assert iers.IERS_B.iers_table is not None
        assert isinstance(iers.IERS_B.iers_table, Table)
        jd1 = np.array([2456108.5, 2456108.5, 2456108.5,
                        2456109.5, 2456109.5])
        jd2 = np.array([0.49999421, 0.99997685, 0.99998843, 0., 0.5])
        ut1_utc, status = iers_b.ut1_utc(jd1, jd2)
        assert np.all(status == iers.FROM_IERS_B)
        assert allclose_sec(ut1_utc, np.array([-0.5868211, -0.5868184,
                                               -0.5868184, 0.4131816,
                                               0.41328895]))
        ut1_utc2, status2 = iers_b.ut1_utc(1e11, 0.)
        assert status2 == iers.TIME_BEYOND_IERS_RANGE
