# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for `astropy.vo.validator`

.. note::

    This test will fail if external URL query status
    changes. This is beyond the control of AstroPy.
    When this happens, rerun or update the test.

Examples
--------
Running inside Python::

    >>> import astropy
    >>> astropy.test('vo.validator', remote_data=True)

Running from top level via command line::

    $ python setup.py test -P vo.validator --remote-data

Running from ``astropy/vo/validator/tests`` directory::

    $ setenv ASTROPY_USE_SYSTEM_PYTEST 1
    $ py.test test_validate.py --remote-data

"""
# STDLIB
import json
import os
import shutil
import sys
import tempfile

# LOCAL
from .. import inspect, validate
from ...client.vos_catalog import BASEURL
from ....tests.helper import pytest, remote_data
from ....utils.data import _find_pkg_data_path, get_pkg_data_filename
from ....utils.data import REMOTE_TIMEOUT


__doctest_skip__ = ['*']


@remote_data
class TestConeSearchValidation(object):
    """Validation on a small subset of Cone Search sites."""
    def setup_class(self):
        self.datadir = 'data'
        self.out_dir = tempfile.mkdtemp()
        self.filenames = {'good': 'conesearch_good.json',
                          'warn': 'conesearch_warn.json',
                          'excp': 'conesearch_exception.json',
                          'nerr': 'conesearch_error.json'}

        validate.CS_MSTR_LIST.set(get_pkg_data_filename(os.path.join(
            self.datadir, 'vao_conesearch_sites_121107_subset.xml')))

        REMOTE_TIMEOUT.set(30)

    @pytest.mark.parametrize(('parallel'), [True, False])
    def test_validation(self, parallel):
        if os.path.exists(self.out_dir):
            shutil.rmtree(self.out_dir)

        # For some reason, Python 3.3 does not work
        # using multiprocessing.Pool and callback function.
        # See http://bugs.python.org/issue16307
        try:
            validate.check_conesearch_sites(
                destdir=self.out_dir, parallel=parallel, url_list=None)
        except AssertionError as e:
            if parallel and sys.version_info >= (3, 3):
                pytest.xfail('See http://bugs.python.org/issue16307')
            else:
                raise

        for val in self.filenames.values():
            _compare_catnames(get_pkg_data_filename(
                os.path.join(self.datadir, val)),
                os.path.join(self.out_dir, val))

    def test_url_list(self):
        local_outdir = os.path.join(self.out_dir, 'subtmp1')
        local_list = [
            'http://www.google.com/foo&',
            'http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/252/out&']
        # Same multiprocessing problem in Python 3.3 as above
        validate.check_conesearch_sites(destdir=local_outdir, parallel=False,
                                        url_list=local_list)
        _compare_catnames(get_pkg_data_filename(
            os.path.join(self.datadir, 'conesearch_good_subset.json')),
            os.path.join(local_outdir, 'conesearch_good.json'))

    def teardown_class(self):
        validate.CS_MSTR_LIST.set(validate.CS_MSTR_LIST.defaultvalue)
        REMOTE_TIMEOUT.set(REMOTE_TIMEOUT.defaultvalue)
        shutil.rmtree(self.out_dir)


class TestConeSearchResults(object):
    """Inspection of `TestConeSearchValidation` results."""
    def setup_class(self):
        self.datadir = 'data'
        self.out_dir = tempfile.mkdtemp()
        BASEURL.set(_find_pkg_data_path(self.datadir) + os.sep)
        self.r = inspect.ConeSearchResults()

    def test_catkeys(self):
        assert (self.r.catkeys['good'] ==
                ['HST Guide Star Catalog 2.3 1',
                 'The USNO-A2.0 Catalogue (Monet+ 1998) 1'])
        assert self.r.catkeys['warn'] == []
        assert self.r.catkeys['exception'] == []
        assert self.r.catkeys['error'] == []

    def gen_cmp(self, func, oname, *args, **kwargs):
        dat_file = get_pkg_data_filename(os.path.join(self.datadir, oname))
        out_file = os.path.join(self.out_dir, oname)
        with open(out_file, 'w') as fout:
            func(fout=fout, *args, **kwargs)

        with open(dat_file) as f1:
            with open(out_file) as f2:
                assert f1.read() == f2.read()

    def test_tally(self):
        self.gen_cmp(self.r.tally, 'tally.out')

    def test_listcats(self):
        self.gen_cmp(self.r.list_cats, 'listcats1.out', 'good')
        self.gen_cmp(self.r.list_cats, 'listcats2.out', 'good',
                     ignore_noncrit=True)

    def test_printcat(self):
        self.gen_cmp(self.r.print_cat, 'printcat.out',
                     'The USNO-A2.0 Catalogue (Monet+ 1998) 1')

    def teardown_class(self):
        BASEURL.set(BASEURL.defaultvalue)
        shutil.rmtree(self.out_dir)


def _load_catnames(fname):
    with open(fname, 'r') as fd:
        js = json.load(fd)
        cats = sorted(js['catalogs'])
    return cats


def _compare_catnames(fname1, fname2):
    cat1 = _load_catnames(fname1)
    cat2 = _load_catnames(fname2)
    assert cat1 == cat2
