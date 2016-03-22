# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for `astropy.vo.validator.validate`.

.. note::

    This test will fail if external URL query status
    changes. This is beyond the control of AstroPy.
    When this happens, rerun or update the test.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

# STDLIB
import os
import shutil
import tempfile

# LOCAL
from .. import conf
from .. import validate
from ...client.vos_catalog import VOSDatabase
from ....tests.helper import pytest, remote_data
from ....utils.data import get_pkg_data_filename
from ....utils import data


__doctest_skip__ = ['*']


@remote_data
class TestConeSearchValidation(object):
    """Validation on a small subset of Cone Search sites."""
    def setup_class(self):
        self.datadir = 'data'
        self.out_dir = tempfile.mkdtemp()
        self.filenames = {
            'good': 'conesearch_good.json',
            'warn': 'conesearch_warn.json',
            'excp': 'conesearch_exception.json',
            'nerr': 'conesearch_error.json'}

        conf.conesearch_master_list = get_pkg_data_filename(os.path.join(
            self.datadir, 'vao_conesearch_sites_121107_subset.xml'))

        data.conf.remote_timeout = 30

    @staticmethod
    def _compare_catnames(fname1, fname2):
        db1 = VOSDatabase.from_json(fname1)
        db2 = VOSDatabase.from_json(fname2)
        assert db1.list_catalogs() == db2.list_catalogs()

    @pytest.mark.parametrize(('parallel'), [True, False])
    def test_validation(self, parallel):
        if os.path.exists(self.out_dir):
            shutil.rmtree(self.out_dir)

        validate.check_conesearch_sites(
            destdir=self.out_dir, parallel=parallel, url_list=None)

        for val in self.filenames.values():
            self._compare_catnames(get_pkg_data_filename(
                os.path.join(self.datadir, val)),
                os.path.join(self.out_dir, val))

    @pytest.mark.parametrize(('parallel'), [True, False])
    def test_url_list(self, parallel):
        local_outdir = os.path.join(self.out_dir, 'subtmp1')
        local_list = [
            'http://www.google.com/foo&',
            'http://vizier.u-strasbg.fr/viz-bin/votable/-A?-out.all&-source=I/252/out&']
        validate.check_conesearch_sites(destdir=local_outdir,
                                        parallel=parallel,
                                        url_list=local_list)
        self._compare_catnames(get_pkg_data_filename(
            os.path.join(self.datadir, self.filenames['good'])),
            os.path.join(local_outdir, self.filenames['good']))

    def teardown_class(self):
        conf.reset('conesearch_master_list')
        data.conf.reset('remote_timeout')
        shutil.rmtree(self.out_dir)
