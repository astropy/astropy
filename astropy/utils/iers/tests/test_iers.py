# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools
import numpy as np

from ....tests.helper import pytest
from .. import iers
from ....table import Table
from ....time import Time

allclose_sec = functools.partial(np.allclose, rtol=1e-15, atol=1e-9)
# 1 nanosec atol

try:
    iers.IERS_A.open()  # check if IERS_A is available
except IOError:
    HAS_IERS_A = False
else:
    HAS_IERS_A = True


class TestBasic():
    """Basic tests that IERS_B returns correct values"""

    def test_simple(self):
        iers.IERS.close()
        assert iers.IERS.iers_table is None
        iers_tab = iers.IERS.open()
        assert iers.IERS.iers_table is not None
        assert isinstance(iers.IERS.iers_table, Table)
        jd1 = np.array([2456108.5, 2456108.5, 2456108.5,
                        2456109.5, 2456109.5])
        jd2 = np.array([0.49999421, 0.99997685, 0.99998843, 0., 0.5])
        ut1_utc, status = iers_tab.ut1_utc(jd1, jd2)
        assert np.all(status == iers.FROM_IERS_B)
        assert allclose_sec(ut1_utc, np.array([-0.5868211, -0.5868184,
                                               -0.5868184, 0.4131816,
                                               0.41328895]))
        # should be future-proof; surely we've moved to another planet by then
        ut1_utc2, status2 = iers_tab.ut1_utc(1e11, 0.)
        assert status2 == iers.TIME_BEYOND_IERS_RANGE

        # check it works via Time too
        t = Time(jd1, jd2, format='jd', scale='utc')
        ut1_utc3, status3 = iers_tab.ut1_utc(t)
        assert np.all(status3 == iers.FROM_IERS_B)
        assert allclose_sec(ut1_utc3, np.array([-0.5868211, -0.5868184,
                                                -0.5868184, 0.4131816,
                                                0.41328895]))


@pytest.mark.skipif('not HAS_IERS_A')
class TestIERS_A():
    def test_simple(self):
        iers_tab = iers.IERS_A.open()
        jd1 = np.array([2456108.5, 2456108.5, 2456108.5,
                        2456109.5, 2456109.5])
        jd2 = np.array([0.49999421, 0.99997685, 0.99998843, 0., 0.5])
        ut1_utc, status = iers_tab.ut1_utc(jd1, jd2)
        assert np.all(status == iers.FROM_IERS_B)
        assert allclose_sec(ut1_utc, np.array([-0.5868211, -0.5868184,
                                               -0.5868184, 0.4131816,
                                               0.41328895]))
        ut1_utc2, status2 = iers_tab.ut1_utc(1e11, 0.)
        assert status2 == iers.TIME_BEYOND_IERS_RANGE

        tnow = Time.now()

        ut1_utc3, status3 = iers_tab.ut1_utc(tnow)
        assert status3 == iers.FROM_IERS_A_PREDICTION
        assert ut1_utc3 != 0.
