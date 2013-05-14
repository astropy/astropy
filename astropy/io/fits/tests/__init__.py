# Licensed under a 3-clause BSD style license - see PYFITS.rst
from __future__ import division  # confidence high

import os
import shutil
import stat
import tempfile
import time
import warnings

from ... import fits


class FitsTestCase(object):
    def setup(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), 'data')
        self.temp_dir = tempfile.mkdtemp(prefix='fits-test-')

        # Restore global settings to defaults
        # TODO: Replace this when there's a better way to in the config API to
        # force config values to their defaults
        fits.ENABLE_RECORD_VALUED_KEYWORD_CARDS.set(True)
        fits.EXTENSION_NAME_CASE_SENSITIVE.set(False)
        fits.STRIP_HEADER_WHITESPACE.set(True)
        fits.USE_MEMMAP.set(True)

        # Ignore deprecation warnings--this only affects Python 2.5 and 2.6,
        # since deprecation warnings are ignored by defualt on 2.7
        warnings.simplefilter('ignore')
        warnings.simplefilter('always', UserWarning)

    def teardown(self):
        warnings.resetwarnings()
        if hasattr(self, 'temp_dir') and os.path.exists(self.temp_dir):
            tries = 3
            while tries:
                try:
                    shutil.rmtree(self.temp_dir)
                    break
                except OSError:
                    # Probably couldn't delete the file because for whatever
                    # reason a handle to it is still open/hasn't been
                    # garbage-collected
                    time.sleep(0.5)
                    tries -= 1

    def copy_file(self, filename):
        """Copies a backup of a test data file to the temp dir and sets its
        mode to writeable.
        """

        shutil.copy(self.data(filename), self.temp(filename))
        os.chmod(self.temp(filename), stat.S_IREAD | stat.S_IWRITE)

    def data(self, filename):
        """Returns the path to a test data file."""

        return os.path.join(self.data_dir, filename)

    def temp(self, filename):
        """ Returns the full path to a file in the test temp dir."""

        return os.path.join(self.temp_dir, filename)
