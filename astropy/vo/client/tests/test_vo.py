# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for `astropy.vo.client`

Examples
--------
Running from top level via command line::

    python setup.py test -P vo.client --remote-data

Running from `astropy/vo/client/tests` directory::

    setenv ASTROPY_USE_SYSTEM_PYTEST 1
    py.test test_vo.py --remote-data

"""
# THIRD-PARTY
import numpy

# LOCAL
from .. import conesearch, vos_catalog
from ....tests.helper import remote_data


@remote_data
def test_basic_db():
    """
    Read dummy `basic.json` database to test underlying database
    functionality.

    """
    basic_db = vos_catalog.get_remote_catalog_db('basic')
    assert sorted(basic_db.keys()) == ['__version__', 'catalogs', 'content']
    assert basic_db['content'] == ['A', 'B', 'C']

    assert basic_db.list_catalogs() == ['foo']
    assert basic_db.list_catalogs(match_string='whatever', sort=True) == []

    for k, v in basic_db.get_catalogs():
        assert k == 'foo'
        assert v._tree == 'bar'

    foo_cat = basic_db.get_catalog('foo')
    assert foo_cat._tree == 'bar'

    try:
        x = basic_db.get_catalog('not_there')
    except vos_catalog.VOSError as e:
        pass

    assert vos_catalog.list_catalogs('basic') == ['foo']


@remote_data
class TestConeSearch():
    """
    Test Cone Search on a pre-defined access URL.

    .. note::

        This test will fail if the URL becomes inaccessible,
        which is beyond AstroPy's control. When this happens,
        change the test to use a different URL.

        At the time this was written, `pedantic=True` will
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
        assert conesearch.list_catalogs(sort=True) == \
            ['BROKEN', 'USNO ACT', 'USNO NOMAD', 'USNO-A2', 'USNO-B1']

    def test_one_search(self):
        """
        This does not necessarily uses `self.url` because of
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
        numpy.testing.assert_array_equal(tab_2.array, tab_3.array)

        assert tab_2.url == tab_4.url
        numpy.testing.assert_array_equal(tab_2.array, tab_4.array)

    def test_async(self):
        async_search = conesearch.AsyncConeSearch(
            self.ra, self.dec, self.sr, pedantic=self.pedantic)

        tab = async_search.get(timeout=vos_catalog.TIMEOUT())

        assert async_search.is_done()
        assert tab.array.size > 0

    def test_prediction(self):
        """Prediction tests are not very accurate but will have to do."""
        t_1, n_1, url_1 = conesearch.conesearch_timer(
            self.ra, self.dec, self.sr, catalog_db=self.url,
            pedantic=self.pedantic, verbose=self.verbose)

        t_2, n_2 = conesearch.predict_search(
            self.ra, self.dec, self.sr, catalog_db=self.url,
            pedantic=self.pedantic, verbose=self.verbose)

        assert n_2 > 0 and n_2 <= n_1 * 1.5
        assert t_2 > 0 and t_2 <= t_1 * 1.5
