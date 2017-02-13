# Licensed under a 3-clause BSD style license - see PYFITS.rst
from __future__ import division  # confidence high

import os
import shutil
import stat
import tempfile
import time

from ... import fits


class FitsTestCase(object):
    def setup(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), 'data')
        self.temp_dir = tempfile.mkdtemp(prefix='fits-test-')

        # Restore global settings to defaults
        # TODO: Replace this when there's a better way to in the config API to
        # force config values to their defaults
        fits.conf.enable_record_valued_keyword_cards = True
        fits.conf.extension_name_case_sensitive = False
        fits.conf.strip_header_whitespace = True
        fits.conf.use_memmap = True

    def teardown(self):
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

        fits.conf.reset('enable_record_valued_keyword_cards')
        fits.conf.reset('extension_name_case_sensitive')
        fits.conf.reset('strip_header_whitespace')
        fits.conf.reset('use_memmap')

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
