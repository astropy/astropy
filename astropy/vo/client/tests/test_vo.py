"""Tests for `astropy.vo.client.vos_catalog`"""

from .. import conesearch, vos_catalog


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
    except VOSError as e:
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
        self.url = 'http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=NOMAD&'

        # Search to perform
        self.ra = 0
        self.dec = 0
        self.sr = 0.3

        # Avoid downloading the full database
        self.old_stype = conesearch._SERVICE_TYPE
        conesearch._SERVICE_TYPE = 'conesearch_simple'

    def test(self):
        #self.url, self.ra, self.dec, self.sr
        #catdb_none = conesearch.conesearch(self.ra, self.dec, self.sr, pedantic=False)
        pass

    def teardown_class(self):
        """Undo local change before exiting."""
        conesearch._SERVICE_TYPE = self.old_stype
