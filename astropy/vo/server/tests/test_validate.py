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
import filecmp
import os
import shutil
import tempfile

# LOCAL
from .. import validate
from ....config import get_data_filename
from ....tests.helper import remote_data


@remote_data
class TestConeSearchValidation():
    """Validation on a small subset of Cone Search sites."""
    def setup_class(self):
        self.datadir = 'data' + os.sep
        self.out_dir = tempfile.mkdtemp()
        self.filenames = {'good': 'conesearch_good.json',
                          'warn': 'conesearch_warn.json',
                          'excp': 'conesearch_exception.json',
                          'nerr': 'conesearch_error.json'}

        validate.CS_MSTR_LIST.set(get_data_filename(
            self.datadir + 'vao_conesearch_sites_121107_subset.xml'))

    def test_validation(self):
        validate.check_conesearch_sites(destdir=self.out_dir)

        for val in self.filenames.values():
            assert filecmp.cmp(get_data_filename(self.datadir + val),
                               self.out_dir + os.sep + val)

        # Symbolic link
        assert filecmp.cmp(
            get_data_filename(self.datadir + self.filenames['warn']),
            self.out_dir + os.sep + 'conesearch.json')

    def teardown_class(self):
        shutil.rmtree(self.out_dir)
