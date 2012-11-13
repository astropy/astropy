# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for `astropy.vo.server`

Examples
--------
Running from top level via command line::

    python setup.py test -P vo.server --remote-data

Running from `astropy/vo/server/tests` directory::

    setenv ASTROPY_USE_SYSTEM_PYTEST 1
    py.test test_validate.py --remote-data

"""
# STDLIB
import json
import os
import shutil
import tempfile

# LOCAL
from .. import validate
from ....config import get_data_filename
from ....tests.helper import pytest, remote_data


@remote_data
class TestConeSearchValidation():
    """
    Validation on a small subset of Cone Search sites.

    .. note::

        This test will fail if external URL query status
        changes. This is beyond the control of AstroPy.
        When this happens, update the test.

    """
    def setup_class(self):
        self.datadir = 'data'
        self.out_dir = tempfile.mkdtemp()
        self.filenames = {'good': 'conesearch_good.json',
                          'warn': 'conesearch_warn.json',
                          'excp': 'conesearch_exception.json',
                          'nerr': 'conesearch_error.json'}

        validate.CS_MSTR_LIST.set(get_data_filename(os.path.join(
            self.datadir, 'vao_conesearch_sites_121107_subset.xml')))

    @pytest.mark.parametrize(('multiproc'), [True, False])
    def test_validation(self, multiproc):
        if os.path.exists(self.out_dir):
            shutil.rmtree(self.out_dir)

        validate.check_conesearch_sites(
            destdir=self.out_dir, multiproc=multiproc)

        for val in self.filenames.values():
            _compare_catnames(get_data_filename(
                os.path.join(self.datadir, val)),
                os.path.join(self.out_dir, val))

        # Symbolic link
        _compare_catnames(get_data_filename(
            os.path.join(self.datadir, self.filenames['good'])),
            os.path.join(self.out_dir, 'conesearch.json'))

    def test_url_list(self):
        local_outdir = os.path.join(self.out_dir, 'subtmp1')
        local_list = [
            'http://heasarc.gsfc.nasa.gov/cgi-bin/vo/cone/coneGet.pl?'
            'table=batse4b&',
            'http://vizier.u-strasbg.fr/viz-bin/votable/-A?'
            '-source=J/AJ/130/2212/table3&']
        validate.check_conesearch_sites(destdir=local_outdir,
                                        url_list=local_list)
        _compare_catnames(get_data_filename(
            os.path.join(self.datadir, 'conesearch_warn_subset.json')),
            os.path.join(local_outdir, 'conesearch_warn.json'))

    def teardown_class(self):
        shutil.rmtree(self.out_dir)


def _load_catnames(fname):
    with open(fname, 'r') as fd:
        js = json.load(fd)
        cats = sorted(js['catalogs'].keys())
    return cats


def _compare_catnames(fname1, fname2):
    cat1 = _load_catnames(fname1)
    cat2 = _load_catnames(fname2)
    assert cat1 == cat2
