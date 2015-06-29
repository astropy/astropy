# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import numpy as np

from ....tests.helper import pytest, assert_quantity_allclose, remote_data
from .. import iers
from .... import units as u
from ....table import QTable
from ....time import Time


FILE_NOT_FOUND_ERROR = getattr(__builtins__, 'FileNotFoundError', IOError)

try:
    iers.IERS_A.open()  # check if IERS_A is available
except IOError:
    HAS_IERS_A = False
else:
    HAS_IERS_A = True

IERS_A_EXCERPT = os.path.join(os.path.dirname(__file__), 'iers_a_excerpt')
IERS_A_EXCERPT_NETWORK = "https://raw.githubusercontent.com/astropy/astropy/master/astropy/utils/iers/tests/iers_a_excerpt"

class TestBasic():
    """Basic tests that IERS_B returns correct values"""

    def test_simple(self):
        iers.IERS.close()
        assert iers.IERS.iers_table is None
        iers_tab = iers.IERS.open()
        assert iers.IERS.iers_table is not None
        assert isinstance(iers.IERS.iers_table, QTable)
        assert iers_tab['UT1_UTC'].unit is u.second
        assert iers_tab['PM_x'].unit is u.arcsecond
        assert iers_tab['PM_y'].unit is u.arcsecond
        jd1 = np.array([2456108.5, 2456108.5, 2456108.5, 2456109.5, 2456109.5])
        jd2 = np.array([0.49999421, 0.99997685, 0.99998843, 0., 0.5])
        ut1_utc = iers_tab.ut1_utc(jd1, jd2)
        assert isinstance(ut1_utc, u.Quantity)
        assert ut1_utc.unit is u.second
        assert_quantity_allclose(ut1_utc, [-0.5868211, -0.5868184, -0.5868184,
                                           0.4131816, 0.41328895] * u.s,
                                 atol=1.*u.ns)
        # should be future-proof; surely we've moved to another planet by then
        with pytest.raises(IndexError):
            ut1_utc2, status2 = iers_tab.ut1_utc(1e11, 0.)
        # also check it returns the right status
        ut1_utc2, status2 = iers_tab.ut1_utc(jd1, jd2, return_status=True)
        assert np.all(status2 == iers.FROM_IERS_B)
        ut1_utc4, status4 = iers_tab.ut1_utc(1e11, 0., return_status=True)
        assert status4 == iers.TIME_BEYOND_IERS_RANGE

        # check it works via Time too
        t = Time(jd1, jd2, format='jd', scale='utc')
        ut1_utc3 = iers_tab.ut1_utc(t)
        assert_quantity_allclose(ut1_utc3, [-0.5868211, -0.5868184, -0.5868184,
                                            0.4131816, 0.41328895] *u.s,
                                 atol=1.*u.ns)

    def test_open_filename(self):
        iers.IERS.close()
        iers.IERS.open(iers.IERS_B_FILE)
        assert iers.IERS.iers_table is not None
        assert isinstance(iers.IERS.iers_table, QTable)
        iers.IERS.close()
        with pytest.raises(FILE_NOT_FOUND_ERROR):
            iers.IERS.open('surely this does not exist')

    @remote_data
    def test_open_network_url(self):
        iers.IERS_A.close()
        iers.IERS_A.open(IERS_A_EXCERPT_NETWORK)
        assert iers.IERS_A.iers_table is not None
        assert isinstance(iers.IERS_A.iers_table, QTable)
        iers.IERS_A.close()

class TestIERS_AExcerpt():
    def test_simple(self):
        iers_tab = iers.IERS_A.open(IERS_A_EXCERPT)
        assert iers_tab['UT1_UTC'].unit is u.second
        assert 'P' in iers_tab['UT1Flag']
        assert 'I' in iers_tab['UT1Flag']
        assert 'B' in iers_tab['UT1Flag']
        assert np.all((iers_tab['UT1Flag'] == 'I') |
                      (iers_tab['UT1Flag'] == 'P') |
                      (iers_tab['UT1Flag'] == 'B'))
        assert iers_tab['PM_x'].unit is u.arcsecond
        assert iers_tab['PM_y'].unit is u.arcsecond
        assert 'P' in iers_tab['PolPMFlag']
        assert 'I' in iers_tab['PolPMFlag']
        assert 'B' in iers_tab['PolPMFlag']
        assert np.all((iers_tab['PolPMFlag'] == 'P') |
                      (iers_tab['PolPMFlag'] == 'I') |
                      (iers_tab['PolPMFlag'] == 'B'))
        t = Time([57053., 57054., 57055.], format='mjd')
        ut1_utc, status = iers_tab.ut1_utc(t, return_status=True)
        assert status[0] == iers.FROM_IERS_B
        assert np.all(status[1:] == iers.FROM_IERS_A)
        assert_quantity_allclose(ut1_utc,
                                 [-0.4916557, -0.4925323, -0.4934373] * u.s,
                                 atol=1.*u.ns)
        pm_x, pm_y, status = iers_tab.pm_xy(t, return_status=True)
        assert status[0] == iers.FROM_IERS_B
        assert np.all(status[1:] == iers.FROM_IERS_A)
        assert_quantity_allclose(pm_x,
                                 [0.003734, 0.004581, 0.004623] * u.arcsec,
                                 atol=1.*u.narcsec)
        assert_quantity_allclose(pm_y,
                                 [0.310824, 0.313150, 0.315517] * u.arcsec,
                                 atol=1.*u.narcsec)


@pytest.mark.skipif(str('not HAS_IERS_A'))
class TestIERS_A():
    def test_simple(self):
        iers_tab = iers.IERS_A.open()
        jd1 = np.array([2456108.5, 2456108.5, 2456108.5, 2456109.5, 2456109.5])
        jd2 = np.array([0.49999421, 0.99997685, 0.99998843, 0., 0.5])
        ut1_utc, status = iers_tab.ut1_utc(jd1, jd2, return_status=True)
        assert np.all(status == iers.FROM_IERS_B)
        assert_quantity_allclose(ut1_utc, [-0.5868211, -0.5868184, -0.5868184,
                                           0.4131816, 0.41328895] * u.s,
                                 atol=1.*u.ns)
        ut1_utc2, status2 = iers_tab.ut1_utc(1e11, 0., return_status=True)
        assert status2 == iers.TIME_BEYOND_IERS_RANGE

        tnow = Time.now()

        ut1_utc3, status3 = iers_tab.ut1_utc(tnow, return_status=True)
        assert status3 == iers.FROM_IERS_A_PREDICTION
        assert ut1_utc3 != 0.
