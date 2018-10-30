# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import urllib.request
import warnings

import pytest
import numpy as np

from ....tests.helper import assert_quantity_allclose, catch_warnings
from .. import iers
from .... import units as u
from ....table import QTable
from ....time import Time, TimeDelta

FILE_NOT_FOUND_ERROR = getattr(__builtins__, 'FileNotFoundError', OSError)

try:
    iers.IERS_A.open('finals2000A.all')  # check if IERS_A is available
except OSError:
    HAS_IERS_A = False
else:
    HAS_IERS_A = True

IERS_A_EXCERPT = os.path.join(os.path.dirname(__file__), 'iers_a_excerpt')


class TestBasic():
    """Basic tests that IERS_B returns correct values"""

    def test_simple(self):
        iers.IERS.close()
        assert iers.IERS.iers_table is None
        iers_tab = iers.IERS.open()
        assert iers.IERS.iers_table is not None
        assert isinstance(iers.IERS.iers_table, QTable)
        assert (iers_tab['UT1_UTC'].unit / u.second).is_unity()
        assert (iers_tab['PM_x'].unit / u.arcsecond).is_unity()
        assert (iers_tab['PM_y'].unit / u.arcsecond).is_unity()
        jd1 = np.array([2456108.5, 2456108.5, 2456108.5, 2456109.5, 2456109.5])
        jd2 = np.array([0.49999421, 0.99997685, 0.99998843, 0., 0.5])
        ut1_utc = iers_tab.ut1_utc(jd1, jd2)
        assert isinstance(ut1_utc, u.Quantity)
        assert (ut1_utc.unit / u.second).is_unity()
        # IERS files change at the 0.1 ms level; see gh-6981
        assert_quantity_allclose(ut1_utc, [-0.5868211, -0.5868184, -0.5868184,
                                           0.4131816, 0.41328895] * u.s,
                                 atol=0.1*u.ms)
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
                                            0.4131816, 0.41328895] * u.s,
                                 atol=0.1*u.ms)

        # Table behaves properly as a table (e.g. can be sliced)
        assert len(iers_tab[:2]) == 2

    def test_open_filename(self):
        iers.IERS.close()
        iers.IERS.open(iers.IERS_B_FILE)
        assert iers.IERS.iers_table is not None
        assert isinstance(iers.IERS.iers_table, QTable)
        iers.IERS.close()
        with pytest.raises(FILE_NOT_FOUND_ERROR):
            iers.IERS.open('surely this does not exist')

    def test_open_network_url(self):
        iers.IERS_A.close()
        iers.IERS_A.open("file:" + urllib.request.pathname2url(IERS_A_EXCERPT))
        assert iers.IERS_A.iers_table is not None
        assert isinstance(iers.IERS_A.iers_table, QTable)
        iers.IERS_A.close()


class TestIERS_AExcerpt():
    def test_simple(self):
        # Test the IERS A reader. It is also a regression tests that ensures
        # values do not get overridden by IERS B; see #4933.
        iers_tab = iers.IERS_A.open(IERS_A_EXCERPT)

        assert (iers_tab['UT1_UTC'].unit / u.second).is_unity()
        assert 'P' in iers_tab['UT1Flag']
        assert 'I' in iers_tab['UT1Flag']
        assert 'B' in iers_tab['UT1Flag']
        assert np.all((iers_tab['UT1Flag'] == 'I') |
                      (iers_tab['UT1Flag'] == 'P') |
                      (iers_tab['UT1Flag'] == 'B'))

        assert (iers_tab['dX_2000A'].unit / u.marcsec).is_unity()
        assert (iers_tab['dY_2000A'].unit / u.marcsec).is_unity()
        assert 'P' in iers_tab['NutFlag']
        assert 'I' in iers_tab['NutFlag']
        assert 'B' in iers_tab['NutFlag']
        assert np.all((iers_tab['NutFlag'] == 'P') |
                      (iers_tab['NutFlag'] == 'I') |
                      (iers_tab['NutFlag'] == 'B'))

        assert (iers_tab['PM_x'].unit / u.arcsecond).is_unity()
        assert (iers_tab['PM_y'].unit / u.arcsecond).is_unity()
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
        # These values are *exactly* as given in the table, so they should
        # match to double precision accuracy.
        assert_quantity_allclose(ut1_utc,
                                 [-0.4916557, -0.4925323, -0.4934373] * u.s,
                                 atol=0.1*u.ms)

        dcip_x, dcip_y, status = iers_tab.dcip_xy(t, return_status=True)
        assert status[0] == iers.FROM_IERS_B
        assert np.all(status[1:] == iers.FROM_IERS_A)
        # These values are *exactly* as given in the table, so they should
        # match to double precision accuracy.
        print(dcip_x)
        print(dcip_y)
        assert_quantity_allclose(dcip_x,
                                 [-0.086, -0.093, -0.087] * u.marcsec,
                                 atol=1.*u.narcsec)
        assert_quantity_allclose(dcip_y,
                                 [0.094, 0.081, 0.072] * u.marcsec,
                                 atol=1*u.narcsec)

        pm_x, pm_y, status = iers_tab.pm_xy(t, return_status=True)
        assert status[0] == iers.FROM_IERS_B
        assert np.all(status[1:] == iers.FROM_IERS_A)
        assert_quantity_allclose(pm_x,
                                 [0.003734, 0.004581, 0.004623] * u.arcsec,
                                 atol=0.1*u.marcsec)
        assert_quantity_allclose(pm_y,
                                 [0.310824, 0.313150, 0.315517] * u.arcsec,
                                 atol=0.1*u.marcsec)

        # Table behaves properly as a table (e.g. can be sliced)
        assert len(iers_tab[:2]) == 2


@pytest.mark.skipif(str('not HAS_IERS_A'))
class TestIERS_A():

    def test_simple(self):
        """Test that open() by default reads a 'finals2000A.all' file."""
        # Ensure we remove any cached table (gh-5131).
        iers.IERS_A.close()
        iers_tab = iers.IERS_A.open()
        jd1 = np.array([2456108.5, 2456108.5, 2456108.5, 2456109.5, 2456109.5])
        jd2 = np.array([0.49999421, 0.99997685, 0.99998843, 0., 0.5])
        ut1_utc, status = iers_tab.ut1_utc(jd1, jd2, return_status=True)
        assert np.all(status == iers.FROM_IERS_B)
        assert_quantity_allclose(ut1_utc, [-0.5868211, -0.5868184, -0.5868184,
                                           0.4131816, 0.41328895] * u.s,
                                 atol=0.1*u.ms)
        ut1_utc2, status2 = iers_tab.ut1_utc(1e11, 0., return_status=True)
        assert status2 == iers.TIME_BEYOND_IERS_RANGE

        tnow = Time.now()

        ut1_utc3, status3 = iers_tab.ut1_utc(tnow, return_status=True)
        assert status3 == iers.FROM_IERS_A_PREDICTION
        assert ut1_utc3 != 0.


class TestIERS_Auto():

    def setup_class(self):
        """Set up useful data for the tests.
        """
        self.N = 40
        self.ame = 30.0
        self.iers_a_file_1 = os.path.join(os.path.dirname(__file__), 'finals2000A-2016-02-30-test')
        self.iers_a_file_2 = os.path.join(os.path.dirname(__file__), 'finals2000A-2016-04-30-test')
        self.iers_a_url_1 = os.path.normpath('file://' + os.path.abspath(self.iers_a_file_1))
        self.iers_a_url_2 = os.path.normpath('file://' + os.path.abspath(self.iers_a_file_2))
        self.t = Time.now() + TimeDelta(10, format='jd') * np.arange(self.N)

    def teardown_method(self, method):
        """Run this after every test.
        """
        iers.IERS_Auto.close()

    def test_interpolate_error_formatting(self):
        """Regression test: make sure the error message in
        IERS_Auto._check_interpolate_indices() is formatted correctly.
        """
        with iers.conf.set_temp('iers_auto_url', self.iers_a_url_1):
            with iers.conf.set_temp('auto_max_age', self.ame):
                with pytest.raises(ValueError) as err:
                    iers_table = iers.IERS_Auto.open()
                    with warnings.catch_warnings():
                        # Ignoring this if it comes up -- IERS_Auto predictive
                        # values are older than 30.0 days but downloading the
                        # latest table did not find newer values
                        warnings.simplefilter('ignore', iers.IERSStaleWarning)
                        iers_table.ut1_utc(self.t.jd1, self.t.jd2)
        assert str(err.value) == iers.INTERPOLATE_ERROR.format(self.ame)

    def test_auto_max_age_none(self):
        """Make sure that iers.INTERPOLATE_ERROR's advice about setting
        auto_max_age = None actually works.
        """
        with iers.conf.set_temp('iers_auto_url', self.iers_a_url_1):
            with iers.conf.set_temp('auto_max_age', None):
                iers_table = iers.IERS_Auto.open()
                delta = iers_table.ut1_utc(self.t.jd1, self.t.jd2)
        assert isinstance(delta, np.ndarray)
        assert delta.shape == (self.N,)
        assert_quantity_allclose(delta, np.array([-0.2246227]*self.N)*u.s)

    def test_auto_max_age_minimum(self):
        """Check that the minimum auto_max_age is enforced.
        """
        with iers.conf.set_temp('iers_auto_url', self.iers_a_url_1):
            with iers.conf.set_temp('auto_max_age', 5.0):
                with pytest.raises(ValueError) as err:
                    iers_table = iers.IERS_Auto.open()
                    delta = iers_table.ut1_utc(self.t.jd1, self.t.jd2)
        assert str(err.value) == 'IERS auto_max_age configuration value must be larger than 10 days'

    @pytest.mark.remote_data
    def test_no_auto_download(self):
        with iers.conf.set_temp('auto_download', False):
            t = iers.IERS_Auto.open()
        assert type(t) is iers.IERS_B

    @pytest.mark.remote_data
    def test_simple(self):

        with iers.conf.set_temp('iers_auto_url', self.iers_a_url_1):

            dat = iers.IERS_Auto.open()
            assert dat['MJD'][0] == 57359.0 * u.d
            assert dat['MJD'][-1] == 57539.0 * u.d

            # Pretend we are accessing at a time 7 days after start of predictive data
            predictive_mjd = dat.meta['predictive_mjd']
            dat._time_now = Time(predictive_mjd, format='mjd') + 7 * u.d

            # Look at times before and after the test file begins.  0.1292905 is
            # the IERS-B value from MJD=57359.  The value in
            # finals2000A-2016-02-30-test has been replaced at this point.
            assert np.allclose(dat.ut1_utc(Time(50000, format='mjd').jd).value, 0.1292905)
            assert np.allclose(dat.ut1_utc(Time(60000, format='mjd').jd).value, -0.2246227)

            # Now pretend we are accessing at time 60 days after start of predictive data.
            # There will be a warning when downloading the file doesn't give new data
            # and an exception when extrapolating into the future with insufficient data.
            dat._time_now = Time(predictive_mjd, format='mjd') + 60 * u.d
            assert np.allclose(dat.ut1_utc(Time(50000, format='mjd').jd).value, 0.1292905)
            with catch_warnings(iers.IERSStaleWarning) as warns:
                with pytest.raises(ValueError) as err:
                    dat.ut1_utc(Time(60000, format='mjd').jd)
            assert 'interpolating from IERS_Auto using predictive values' in str(err)
            assert len(warns) == 1
            assert 'IERS_Auto predictive values are older' in str(warns[0].message)

            # Warning only if we are getting return status
            with catch_warnings(iers.IERSStaleWarning) as warns:
                dat.ut1_utc(Time(60000, format='mjd').jd, return_status=True)
            assert len(warns) == 1
            assert 'IERS_Auto predictive values are older' in str(warns[0].message)

            # Now set auto_max_age = None which says that we don't care how old the
            # available IERS-A file is.  There should be no warnings or exceptions.
            with iers.conf.set_temp('auto_max_age', None):
                with catch_warnings(iers.IERSStaleWarning) as warns:
                    dat.ut1_utc(Time(60000, format='mjd').jd)
                assert not warns

        # Now point to a later file with same values but MJD increased by
        # 60 days and see that things work.  dat._time_now is still the same value
        # as before, i.e. right around the start of predictive values for the new file.
        # (In other words this is like downloading the latest file online right now).
        with iers.conf.set_temp('iers_auto_url', self.iers_a_url_2):

            # Look at times before and after the test file begins.  This forces a new download.
            assert np.allclose(dat.ut1_utc(Time(50000, format='mjd').jd).value, 0.1292905)
            assert np.allclose(dat.ut1_utc(Time(60000, format='mjd').jd).value, -0.3)

            # Now the time range should be different.
            assert dat['MJD'][0] == 57359.0 * u.d
            assert dat['MJD'][-1] == (57539.0 + 60) * u.d
