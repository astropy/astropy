# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for `astropy.vo.client.vos_catalog`.

.. note::

    ``vos_catalog`` convenience functions are indirectly tested in
    `astropy.vo.client.tests.test_conesearch`.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# STDLIB
import os
import shutil
import tempfile

# LOCAL
from .. import vos_catalog
from ..exceptions import (VOSError, MissingCatalog, DuplicateCatalogName,
                          DuplicateCatalogURL)
from ....tests.helper import pytest, remote_data, catch_warnings
from ....utils.exceptions import AstropyDeprecationWarning
from ....utils.data import get_pkg_data_filename


__doctest_skip__ = ['*']

DB_FILE = get_pkg_data_filename(os.path.join('data', 'basic.json'))


class TestCatalog(object):
    """Test VOSCatalog class."""
    def setup_class(self):
        self.cat = vos_catalog.VOSCatalog.create(
            'foo', 'bar.foo', description='test', my_num=100)

    def test_set_attr(self):
        """See if fields are set properly."""
        assert self.cat['title'] == 'foo'
        assert self.cat['url'] == 'bar.foo'
        assert self.cat['description'] == 'test'
        assert self.cat['my_num'] == 100

    def test_legal_deletion(self):
        """Test deletion."""
        self.cat.delete_attribute('description')
        assert sorted(self.cat) == ['my_num', 'title', 'url']

    def test_illegal_deletion(self):
        """Deletion of compulsory key is now allowed."""
        with pytest.raises(VOSError):
            self.cat.delete_attribute('title')


def test_db_missing_catalog_key():
    """Database must have 'catalogs' key."""
    with pytest.raises(VOSError):
        db = vos_catalog.VOSDatabase({})


def test_db_illegal_catalog():
    """Database must have legal catalog, if given.

    This also tests VOSCatalog init check and
    database ``create_empty()`` method.

    """
    db = vos_catalog.VOSDatabase.create_empty()
    db._catalogs['foo'] = {'foo': 'bar'}
    with pytest.raises(VOSError):
        db = vos_catalog.VOSDatabase(db._tree)


class TestDatabase(object):
    """Test VOSDatabase class."""
    def setup_class(self):
        """Use ``from_json()`` method to init."""
        self.db = vos_catalog.VOSDatabase.from_json(DB_FILE)

    def test_set_attr(self):
        """See if database is set up properly."""
        assert sorted(self.db) == ['__version__', 'catalogs', 'content']
        assert self.db['content'] == ['A', 'B', 'C']
        assert self.db._url_keys['bar.foo'] == ['foo']
        assert self.db.version == 1
        assert len(self.db) == 1

    @pytest.mark.parametrize(
        ('methodname1', 'methodname2', 'arg'),
        [('get_catalog', 'get_catalogs', 'foo'),
         ('get_catalog_by_url', 'get_catalogs_by_url', 'bar.foo')])
    def test_get_catalog(self, methodname1, methodname2, arg):
        """Test catalog retrieval."""
        method1 = getattr(self.db, methodname1)
        method2 = getattr(self.db, methodname2)

        foo_cat = method1(arg)

        if methodname2.endswith('url'):
            foo_iter = method2(arg)
        else:
            foo_iter = method2()

        for k, v in foo_iter:
            assert k == 'foo'
            assert (v._tree == foo_cat._tree ==
                    {'title': 'bar', 'url': 'bar.foo'})

        with pytest.raises(MissingCatalog):
            method1('foofoo')

    def test_changes(self):
        """Test catalog addition, listing, deletion, and merge."""
        # Addition and listing
        self.db.add_catalog(
            'new_cat', vos_catalog.VOSCatalog.create('new_cat', 'new_url'))
        self.db.add_catalog_by_url('bar', 'bar.foo', allow_duplicate_url=True)
        assert self.db.list_catalogs() == ['bar', 'foo', 'new_cat']
        assert self.db.list_catalogs(pattern='f') == ['foo']
        assert self.db.list_catalogs_by_url() == ['bar.foo', 'new_url']
        assert self.db.list_catalogs_by_url(pattern='r.f') == ['bar.foo']

        with pytest.raises(DuplicateCatalogName):
            self.db.add_catalog(
                'foo', vos_catalog.VOSCatalog.create('my_title', 'my_url'))

        with pytest.raises(DuplicateCatalogURL):
            self.db.add_catalog(
                'foo3', vos_catalog.VOSCatalog.create('foo3', 'bar.foo'))

        # Deletion
        # Note: Deletion by URL calls delete_catalog() method.
        self.db.delete_catalog_by_url('bar.foo')
        assert self.db.list_catalogs() == ['new_cat']
        assert self.db.list_catalogs_by_url() == ['new_url']
        assert self.db._url_keys['bar.foo'] == []

        # Merge
        other_db = vos_catalog.VOSDatabase.create_empty()
        other_db.add_catalog(
            'o_cat', vos_catalog.VOSCatalog.create('o_title', 'o_url'))
        other_db.add_catalog(
            'o_cat2', vos_catalog.VOSCatalog.create('o_title2', 'o_url2'))
        new_db = self.db.merge(other_db)
        assert new_db.list_catalogs() == ['new_cat', 'o_cat', 'o_cat2']
        assert new_db.list_catalogs_by_url() == ['new_url', 'o_url', 'o_url2']

        # Make sure inputs are unchanged
        assert self.db.list_catalogs() == ['new_cat']
        assert self.db.list_catalogs_by_url() == ['new_url']
        assert other_db.list_catalogs() == ['o_cat', 'o_cat2']
        assert other_db.list_catalogs_by_url() == ['o_url', 'o_url2']

    def test_illegal_addition(self):
        """Test illegal catalog addition."""
        with pytest.raises(VOSError):
            self.db.add_catalog('foo2', 'not_a_cat')

    @pytest.mark.parametrize(
        'methodname', ['delete_catalog', 'delete_catalog_by_url'])
    def test_illegal_deletion(self, methodname):
        """Test illegal catalog deletion."""
        method = getattr(self.db, methodname)

        with pytest.raises(MissingCatalog):
            method('not_there')

    def test_illegal_merge(self):
        """Test illegal database merger."""
        with pytest.raises(VOSError):
            new_db = self.db.merge('not_a_db')

        other_db = vos_catalog.VOSDatabase.create_empty()
        other_db['__version__'] = 99
        with pytest.raises(VOSError):
            new_db = self.db.merge(other_db)


class TestWriteJson(object):
    """Test writing database to JSON file."""
    def setup_class(self):
        self.outdir = tempfile.mkdtemp()

    def test_write(self):
        outfile = os.path.join(self.outdir, 'test_1.json')
        db = vos_catalog.VOSDatabase.from_json(DB_FILE)
        db.to_json(outfile)
        with pytest.raises(OSError) as exc:
            db.to_json(outfile)
        assert str(exc.value).endswith("test_1.json exists.")
        with catch_warnings(AstropyDeprecationWarning) as warning_lines:
            db.to_json(outfile, clobber=True)
            assert warning_lines[0].category == AstropyDeprecationWarning
            assert (str(warning_lines[0].message) == '"clobber" was '
                    'deprecated in version 1.3 and will be removed in a '
                    'future version. Use argument "overwrite" instead.')
        db.to_json(outfile, overwrite=True)

        # Read it back in
        db2 = vos_catalog.VOSDatabase.from_json(outfile)
        assert db.list_catalogs() == db2.list_catalogs()
        assert db.list_catalogs_by_url() == db2.list_catalogs_by_url()

    def teardown_class(self):
        shutil.rmtree(self.outdir)


@remote_data
def test_db_from_registry():
    """Test database created from VO registry.

    .. note::

        We have no control of the remote registry.
        This test just makes sure it does not crash,
        but does not check for quality of data.

    """
    from ...validator import conf

    db = vos_catalog.VOSDatabase.from_registry(
        conf.conesearch_master_list, encoding='binary', show_progress=False)

    # Should have over 9k catalogs; Update test if this changes.
    assert len(db) > 9000
