# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for `astropy.vo.client.conesearch` and `astropy.vo.client.async`."""
from __future__ import absolute_import, division, print_function, unicode_literals

# STDLIB
import os

# THIRD-PARTY
import numpy as np

# LOCAL
from .. import conesearch, vos_catalog
from ..exceptions import VOSError, ConeSearchError
from ... import conf
from .... import units as u
from ....coordinates import ICRS, SkyCoord
from ....tests.helper import pytest, remote_data
from ....utils.data import get_pkg_data_filename
from ....utils import data


__doctest_skip__ = ['*']


# Global variables for TestConeSearch
SCS_RA = 0
SCS_DEC = 0
SCS_SR = 0.1
SCS_CENTER = ICRS(SCS_RA * u.degree, SCS_DEC * u.degree)
SCS_RADIUS = SCS_SR * u.degree


@remote_data
class TestConeSearch(object):
    """Test Cone Search on a pre-defined access URL.

    .. note::

        This test will fail if the URL becomes inaccessible,
        which is beyond AstroPy's control. When this happens,
        change the test to use a different URL.

        At the time this was written, ``pedantic=True`` will
        not yield any successful search.

    """
    def setup_class(self):
        # If this link is broken, use the next in database that works
        self.url = 'http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=USNO-A2&'
        self.catname = 'USNO-A2'

        # Avoid downloading the full database
        conf.conesearch_dbname = 'conesearch_simple'

        # Sometimes 3s is not enough
        data.conf.remote_timeout = 10

        self.verbose = False
        self.pedantic = False

    def test_cat_listing(self):
        assert (conesearch.list_catalogs() ==
                ['BROKEN', 'USNO ACT', 'USNO NOMAD', 'USNO-A2', 'USNO-B1'])
        assert (conesearch.list_catalogs(pattern='usno*a') ==
                ['USNO ACT', 'USNO NOMAD', 'USNO-A2'])

    def test_no_result(self):
        with pytest.raises(VOSError):
            tab_1 = conesearch.conesearch(
                SCS_CENTER, 0.001, catalog_db=self.url,
                pedantic=self.pedantic, verbose=self.verbose)

    @pytest.mark.parametrize(('center', 'radius'),
                             [((SCS_RA, SCS_DEC), SCS_SR),
                              (SCS_CENTER, SCS_RADIUS)])
    def test_one_search(self, center, radius):
        """This does not necessarily uses ``self.url`` because of
        unordered dict in JSON tree.

        """
        tab_1 = conesearch.conesearch(
            center, radius, pedantic=None, verbose=self.verbose)

        assert tab_1.array.size > 0

    def test_sky_coord(self):
        """
        Check that searching with a SkyCoord works too
        """
        sc_cen = SkyCoord(SCS_CENTER)
        tab = conesearch.conesearch(
            sc_cen, SCS_RADIUS, catalog_db=self.url,
            pedantic=self.pedantic, verbose=self.verbose)

        assert tab.array.size > 0

    def test_timeout(self):
        """Test time out error."""
        try:
            with data.conf.set_temp('remote_timeout', 0.001):
                tab = conesearch.conesearch(
                    SCS_CENTER, SCS_RADIUS, pedantic=self.pedantic,
                    verbose=self.verbose, catalog_db=self.url)
        except VOSError as e:
            assert 'timed out' in str(e), 'test_timeout failed'
        else:
            raise Exception('test_timeout failed')

    def test_searches(self):
        tab_2 = conesearch.conesearch(
            SCS_CENTER, SCS_RADIUS, catalog_db=self.url,
            pedantic=self.pedantic, verbose=self.verbose)

        tab_3 = conesearch.conesearch(
            SCS_CENTER, SCS_RADIUS,
            catalog_db=[self.catname, self.url],
            pedantic=self.pedantic, verbose=self.verbose)

        tab_4 = conesearch.conesearch(
            SCS_CENTER, SCS_RADIUS,
            catalog_db=vos_catalog.get_remote_catalog_db(
                conf.conesearch_dbname),
            pedantic=self.pedantic, verbose=self.verbose)

        assert tab_2.url == tab_3.url
        np.testing.assert_array_equal(tab_2.array, tab_3.array)

        # If this fails, it is because of dict hashing, no big deal.
        if tab_2.url == tab_4.url:
            np.testing.assert_array_equal(tab_2.array, tab_4.array)
        else:
            pytest.xfail('conesearch_simple.json used a different URL')

    @pytest.mark.parametrize(('center', 'radius'),
                             [((SCS_RA, SCS_DEC), SCS_SR),
                              (SCS_CENTER, SCS_RADIUS)])
    def test_search_all(self, center, radius):
        all_results = conesearch.search_all(
            center, radius, catalog_db=['BROKEN', self.url],
            pedantic=self.pedantic, verbose=self.verbose)

        assert len(all_results) == 1

        tab_1 = all_results[self.url]

        assert tab_1.array.size > 0

    def test_async(self):
        async_search = conesearch.AsyncConeSearch(
            SCS_CENTER, SCS_RADIUS, pedantic=self.pedantic)

        tab = async_search.get(timeout=data.conf.remote_timeout)

        assert async_search.done()
        assert tab.array.size > 0

    def test_async_all(self):
        async_search_all = conesearch.AsyncSearchAll(
            SCS_CENTER, SCS_RADIUS, pedantic=self.pedantic)

        all_results = async_search_all.get(timeout=(data.conf.remote_timeout * 3))

        assert async_search_all.done()
        for tab in all_results.values():
            assert tab.array.size > 0

    @pytest.mark.parametrize(('center', 'radius'),
                             [((SCS_RA, SCS_DEC), 0.8),
                              (SCS_CENTER, 0.8 * u.degree)])
    def test_prediction(self,  center, radius):
        """Prediction tests are not very accurate but will have to do."""
        t_1, tab_1 = conesearch.conesearch_timer(
            center, radius, catalog_db=self.url,
            pedantic=self.pedantic, verbose=self.verbose)
        n_1 = tab_1.array.size

        t_2, n_2 = conesearch.predict_search(
            self.url, center, radius,
            pedantic=self.pedantic, verbose=self.verbose)

        assert n_2 > 0 and n_2 <= n_1 * 1.5
        assert t_2 > 0 and t_2 <= t_1 * 1.5

    def test_prediction_neg_radius(self):
        with pytest.raises(ConeSearchError):
            t, n = conesearch.predict_search(
                self.url, SCS_CENTER, -1, pedantic=self.pedantic,
                verbose=self.verbose)

    def teardown_class(self):
        conf.reset('conesearch_dbname')
        data.conf.reset('remote_timeout')


class TestErrorResponse(object):
    """Test Cone Search error response handling.

    This is defined in Section 2.3 of Simple Cone Search Version 1.03,
    IVOA Recommendation, 22 February 2008.

    Also see https://github.com/astropy/astropy/issues/1001

    """
    def setup_class(self):
        self.datadir = 'data'
        self.pedantic = False
        self.conesearch_errmsg = {
            'conesearch_error1.xml': 'Error in input RA value: as3f',
            'conesearch_error2.xml': 'Error in input RA value: as3f',
            'conesearch_error3.xml': 'Invalid data type: text/html',
            'conesearch_error4.xml': 'Invalid data type: text/html'}

    def conesearch_compare(self, xmlfile, msg):
        """Bypassing Cone Search query and just imitating the reply,
        then check if appropriate error message is caught.

        """
        # conesearch_error4.xml is a wont-fix for now
        if xmlfile == 'conesearch_error4.xml':
            pytest.xfail('Currently not supported, '
                         'see astropy.io.votable.exceptions.W22')

        url = get_pkg_data_filename(os.path.join(self.datadir, xmlfile))
        try:
            r = vos_catalog._vo_service_request(url, self.pedantic, {})
        except VOSError as e:
            assert msg in str(e)

    @pytest.mark.parametrize(('id'), [1, 2, 3, 4])
    def test_conesearch_response(self, id):
        xml = 'conesearch_error{0}.xml'.format(id)
        msg = self.conesearch_errmsg[xml]
        self.conesearch_compare(xml, msg)


def test_validate_coord():
    """Valid coordinates should not raise an error."""
    result = conesearch._validate_coord(
        ICRS(6.02233 * u.degree, -72.08144 * u.degree))
    np.testing.assert_allclose(result, [6.022330000000011, -72.08144])

    result = conesearch._validate_coord(
        SkyCoord(6.02233 * u.degree, -72.08144 * u.degree, frame='icrs'))
    np.testing.assert_allclose(result, [6.022330000000011, -72.08144])

    result = conesearch._validate_coord((-0, 0))
    np.testing.assert_allclose(result, [0, 0])

    result = conesearch._validate_coord((-1, -1))
    np.testing.assert_allclose(result, [359, -1])
