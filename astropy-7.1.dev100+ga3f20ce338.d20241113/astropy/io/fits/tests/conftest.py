# Licensed under a 3-clause BSD style license - see PYFITS.rst

import os
import pathlib
import shutil
import stat
import tempfile
import time

import pytest

from astropy.io import fits


@pytest.fixture(
    params=[False, "str", "pathlib"], ids=["", "home_is_data", "home_is_data, pathlib"]
)
def home_is_data(request, monkeypatch):
    """
    Pytest fixture to run a test case both with and without tilde paths.

    In the tilde-path case, calls like self.data('filename.fits') will
    produce '~/filename.fits', and environment variables will be temporarily
    modified so that '~' resolves to the data directory.
    """
    # This checks the value specified in the fixture annotation
    if request.param:
        # `request.instance` refers to the test case that's using this fixture.
        request.instance.monkeypatch = monkeypatch
        request.instance.set_home_as_data()

        request.instance.set_paths_via_pathlib(request.param == "pathlib")


@pytest.fixture(
    params=[False, "str", "pathlib"], ids=["", "home_is_data", "home_is_data, pathlib"]
)
def home_is_temp(request, monkeypatch):
    """
    Pytest fixture to run a test case both with and without tilde paths.

    In the tilde-path case, calls like self.temp('filename.fits') will
    produce '~/filename.fits', and environment variables will be temporarily
    modified so that '~' resolves to the temp directory. These files will also
    be tracked so that, after the test case, we can verify no files were written
    to a literal tilde path.
    """
    # This checks the value specified in the fixture annotation
    if request.param:
        # `request.instance` refers to the test case that's using this fixture.
        request.instance.monkeypatch = monkeypatch
        request.instance.set_home_as_temp()

        request.instance.set_paths_via_pathlib(request.param == "pathlib")


class FitsTestCase:
    def setup_method(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data")
        self.temp_dir = tempfile.mkdtemp(prefix="fits-test-")

        self.home_is_data = False
        self.home_is_temp = False
        self.temp_files_used = set()
        self.use_pathlib = False

        # Restore global settings to defaults
        # TODO: Replace this when there's a better way to in the config API to
        # force config values to their defaults
        fits.conf.enable_record_valued_keyword_cards = True
        fits.conf.extension_name_case_sensitive = False
        fits.conf.strip_header_whitespace = True
        fits.conf.use_memmap = True

    def teardown_method(self):
        if self.home_is_temp:
            # Verify that no files were written to a literal tilde path
            for temp_file, temp_file_no_tilde in self.temp_files_used:
                assert not os.path.exists(temp_file)
                assert os.path.exists(temp_file_no_tilde)

        if hasattr(self, "temp_dir") and os.path.exists(self.temp_dir):
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

        fits.conf.reset("enable_record_valued_keyword_cards")
        fits.conf.reset("extension_name_case_sensitive")
        fits.conf.reset("strip_header_whitespace")
        fits.conf.reset("use_memmap")

    def copy_file(self, filename):
        """Copies a backup of a test data file to the temp dir and sets its
        mode to writeable.
        """
        shutil.copy(
            os.path.expanduser(self.data(filename)),
            os.path.expanduser(self.temp(filename)),
        )
        os.chmod(os.path.expanduser(self.temp(filename)), stat.S_IREAD | stat.S_IWRITE)

    def data(self, filename):
        """Returns the path to a test data file."""
        if self.home_is_data:
            prefix = "~"
        else:
            prefix = self.data_dir

        if self.use_pathlib:
            return pathlib.Path(prefix, filename)
        return os.path.join(prefix, filename)

    def temp(self, filename):
        """Returns the full path to a file in the test temp dir."""
        real_target = os.path.join(self.temp_dir, filename)
        if self.home_is_temp:
            prefix = "~"
            # Record the '~' path and the intended path, for use
            # in `home_is_temp`
            self.temp_files_used.add((os.path.join(prefix, filename), real_target))
        else:
            prefix = self.temp_dir

        if self.use_pathlib:
            return pathlib.Path(prefix, filename)
        return os.path.join(prefix, filename)

    def set_home_as_data(self):
        """
        This overrides the HOME environment variable, so that paths beginning
        with '~/' expand to the data directory. Used by the `home_is_data`
        fixture.
        """
        self.home_is_data = True
        # For Unix
        self.monkeypatch.setenv("HOME", self.data_dir)
        # For Windows
        self.monkeypatch.setenv("USERPROFILE", self.data_dir)

    def set_home_as_temp(self):
        """
        This overrides the HOME environment variable, so that paths beginning
        with '~/' expand to the temp directory. In conjunction with
        self.temp(), temporary files are tracked as they are created, so we can
        verify they end up in the temporary directory and not unexpected places
        in the filesystem. Used by the `home_is_temp` fixture.
        """
        self.home_is_temp = True
        # For Unix
        self.monkeypatch.setenv("HOME", self.temp_dir)
        # For Windows
        self.monkeypatch.setenv("USERPROFILE", self.temp_dir)

    def set_paths_via_pathlib(self, use_pathlib):
        self.use_pathlib = use_pathlib
