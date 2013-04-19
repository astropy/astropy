# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for `astropy.vo.client`

Examples
--------
Running inside Python:

>>> import astropy
>>> astropy.test('vo.client', remote_data=True)

Running from top level via command line::

    python setup.py test -P vo.client --remote-data

Running from ``astropy/vo/client/tests`` directory::

    setenv ASTROPY_USE_SYSTEM_PYTEST 1
    py.test test_vo.py --remote-data

"""
# STDLIB
import os

# THIRD-PARTY
import numpy as np

# LOCAL
from .. import conesearch, vos_catalog
from ....tests.helper import pytest, remote_data
from ....utils.data import get_pkg_data_filename
from ....utils.data import REMOTE_TIMEOUT


@remote_data
def test_basic_db():
    """Read dummy ``basic.json`` database to test underlying database
    functionality.

    """
    basic_db = vos_catalog.get_remote_catalog_db('basic')
    assert sorted(basic_db.keys()) == ['__version__', 'catalogs', 'content']
    assert basic_db['content'] == ['A', 'B', 'C']

    assert basic_db.list_catalogs() == ['foo']
    assert basic_db.list_catalogs(pattern='whatever') == []

    foo_cat1 = basic_db.get_catalog('foo')
    for k, v in basic_db.get_catalogs():
        assert k == 'foo'
        assert v._tree == foo_cat1._tree == {'title': 'bar', 'url': 'bar.foo'}

    foo_cat2 = basic_db.get_catalog_by_url('bar.foo')
    for k, v in basic_db.get_catalogs_by_url('bar.foo'):
        assert k == 'foo'
        assert v._tree == foo_cat2._tree == {'title': 'bar', 'url': 'bar.foo'}

    try:
        x = basic_db.get_catalog('not_there')
    except vos_catalog.VOSError:
        pass

    assert vos_catalog.list_catalogs('basic') == ['foo']


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

        # Search to perform
        self.ra = 0
        self.dec = 0
        self.sr = 0.1

        # Avoid downloading the full database
        conesearch.CONESEARCH_DBNAME.set('conesearch_simple')

        self.verbose = False
        self.pedantic = False

    def test_cat_listing(self):
        assert (conesearch.list_catalogs() ==
                ['BROKEN', 'USNO ACT', 'USNO NOMAD', 'USNO-A2', 'USNO-B1'])
        assert (conesearch.list_catalogs(pattern='usno*a') ==
                ['USNO ACT', 'USNO NOMAD', 'USNO-A2'])

    def test_one_search(self):
        """This does not necessarily uses ``self.url`` because of
        unordered dict in JSON tree.

        """
        tab_1 = conesearch.conesearch(
            self.ra, self.dec, self.sr,
            pedantic=self.pedantic, verbose=self.verbose)

        assert tab_1.array.size > 0

    def test_searches(self):
        tab_2 = conesearch.conesearch(
            self.ra, self.dec, self.sr, catalog_db=self.url,
            pedantic=self.pedantic, verbose=self.verbose)

        tab_3 = conesearch.conesearch(
            self.ra, self.dec, self.sr, catalog_db=[self.catname, self.url],
            pedantic=self.pedantic, verbose=self.verbose)

        tab_4 = conesearch.conesearch(
            self.ra, self.dec, self.sr,
            catalog_db=vos_catalog.get_remote_catalog_db(
                conesearch.CONESEARCH_DBNAME()),
            pedantic=self.pedantic, verbose=self.verbose)

        assert tab_2.url == tab_3.url
        np.testing.assert_array_equal(tab_2.array, tab_3.array)

        # If this fails, it is because of dict hashing, no big deal.
        if tab_2.url == tab_4.url:
            np.testing.assert_array_equal(tab_2.array, tab_4.array)
        else:
            pytest.xfail('conesearch_simple.json used a different URL')

    def test_async(self):
        async_search = conesearch.AsyncConeSearch(
            self.ra, self.dec, self.sr, pedantic=self.pedantic)

        tab = async_search.get(timeout=REMOTE_TIMEOUT())

        assert async_search.done()
        assert tab.array.size > 0

    def test_prediction(self):
        """Prediction tests are not very accurate but will have to do."""
        t_1, tab_1 = conesearch.conesearch_timer(
            self.ra, self.dec, self.sr, catalog_db=self.url,
            pedantic=self.pedantic, verbose=self.verbose)
        n_1 = tab_1.array.size

        t_2, n_2 = conesearch.predict_search(
            self.url, self.ra, self.dec, self.sr,
            pedantic=self.pedantic, verbose=self.verbose)

        assert n_2 > 0 and n_2 <= n_1 * 1.5
        assert t_2 > 0 and t_2 <= t_1 * 1.5

    def teardown_class(self):
        conesearch.CONESEARCH_DBNAME.set(
            conesearch.CONESEARCH_DBNAME.defaultvalue)


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
        except vos_catalog.VOSError as e:
            assert msg in str(e)

    @pytest.mark.parametrize(('id'), [1, 2, 3, 4])
    def test_conesearch_response(self, id):
        xml = 'conesearch_error{0}.xml'.format(id)
        msg = self.conesearch_errmsg[xml]
        self.conesearch_compare(xml, msg)
