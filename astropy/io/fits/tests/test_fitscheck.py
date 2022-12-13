# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re

import pytest

from astropy import __version__ as version
from astropy.io import fits
from astropy.io.fits.scripts import fitscheck
from astropy.utils.exceptions import AstropyUserWarning

from .conftest import FitsTestCase


class TestFitscheck(FitsTestCase):
    def test_help(self):
        with pytest.raises(SystemExit) as e:
            fitscheck.main(["-h"])
        assert e.value.code == 0

    def test_version(self, capsys):
        with pytest.raises(SystemExit) as e:
            fitscheck.main(["--version"])
            out = capsys.readouterr()[0]
            assert out == f"fitscheck {version}"
        assert e.value.code == 0

    def test_missing_file(self, capsys):
        assert fitscheck.main(["missing.fits"]) == 1
        stdout, stderr = capsys.readouterr()
        assert "No such file or directory" in stderr

    def test_valid_file(self, capsys):
        testfile = self.data("checksum.fits")

        assert fitscheck.main([testfile]) == 0
        assert fitscheck.main([testfile, "--compliance"]) == 0

        assert fitscheck.main([testfile, "-v"]) == 0
        stdout, stderr = capsys.readouterr()
        assert "OK" in stderr

    def test_remove_checksums(self, capsys):
        self.copy_file("checksum.fits")
        testfile = self.temp("checksum.fits")
        assert fitscheck.main([testfile, "--checksum", "remove"]) == 1
        assert fitscheck.main([testfile]) == 1
        stdout, stderr = capsys.readouterr()
        assert "MISSING" in stderr

    def test_no_checksums(self, capsys):
        testfile = self.data("arange.fits")

        assert fitscheck.main([testfile]) == 1
        stdout, stderr = capsys.readouterr()
        assert "Checksum not found" in stderr

        assert fitscheck.main([testfile, "--ignore-missing"]) == 0
        stdout, stderr = capsys.readouterr()
        assert stderr == ""

    def test_overwrite_invalid(self, caplog):
        """
        Tests that invalid checksum or datasum are overwritten when the file is
        saved.
        """
        reffile = self.temp("ref.fits")
        with fits.open(self.data("tb.fits")) as hdul:
            hdul.writeto(reffile, checksum=True)

        # replace checksums with wrong ones
        testfile = self.temp("test.fits")
        with fits.open(self.data("tb.fits")) as hdul:
            hdul[0].header["DATASUM"] = "1       "
            hdul[0].header["CHECKSUM"] = "8UgqATfo7TfoATfo"
            hdul[1].header["DATASUM"] = "2349680925"
            hdul[1].header["CHECKSUM"] = "11daD8bX98baA8bU"
            hdul.writeto(testfile)

        assert fitscheck.main([testfile]) == 1
        assert re.match(
            r"BAD.*Checksum verification failed for HDU", caplog.records[0].message
        )
        caplog.clear()

        with pytest.warns(AstropyUserWarning):
            assert fitscheck.main([testfile, "--write", "--force"]) == 1
        assert re.match(
            r"BAD.*Checksum verification failed for HDU", caplog.records[0].message
        )
        caplog.clear()

        # check that the file was fixed
        assert fitscheck.main([testfile]) == 0
