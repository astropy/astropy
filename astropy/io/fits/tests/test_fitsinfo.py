# Licensed under a 3-clause BSD style license - see LICENSE.rst

import subprocess

import pytest

from astropy import __version__ as version
from astropy.io.fits.scripts import fitsinfo

from .conftest import FitsTestCase


class TestFitsinfo(FitsTestCase):
    def test_help(self):
        with pytest.raises(SystemExit) as e:
            fitsinfo.main(["-h"])
        assert e.value.code == 0

    def test_version(self, capsys):
        script = "fitsinfo"
        script_result = subprocess.run(
            [script, "--version"], capture_output=True, text=True, check=False
        )

        assert script_result.returncode == 0
        assert script_result.stdout.strip() == f"{script} {version}"

    def test_onefile(self, capsys):
        fitsinfo.main([self.data("arange.fits")])
        out, err = capsys.readouterr()
        out = out.splitlines()
        assert len(out) == 3
        assert out[1].startswith(
            "No.    Name      Ver    Type      Cards   Dimensions   Format"
        )
        assert out[2].startswith(
            "  0  PRIMARY       1 PrimaryHDU       7   (11, 10, 7)   int32"
        )

    def test_multiplefiles(self, capsys):
        fitsinfo.main([self.data("arange.fits"), self.data("ascii.fits")])
        out, err = capsys.readouterr()
        out = out.splitlines()
        assert len(out) == 8
        assert out[1].startswith(
            "No.    Name      Ver    Type      Cards   Dimensions   Format"
        )
        assert out[2].startswith(
            "  0  PRIMARY       1 PrimaryHDU       7   (11, 10, 7)   int32"
        )
        assert out[3] == ""
        assert out[7].startswith(
            "  1                1 TableHDU        20   5R x 2C   [E10.4, I5]"
        )
