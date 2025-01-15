# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

import numpy as np
import pytest

from astropy import __version__ as version
from astropy.io import fits
from astropy.io.fits import FITSDiff, HDUList, Header, ImageHDU
from astropy.io.fits.convenience import writeto
from astropy.io.fits.hdu import PrimaryHDU, hdulist
from astropy.io.fits.scripts import fitsdiff
from astropy.utils.misc import _NOT_OVERWRITING_MSG_MATCH

from .conftest import FitsTestCase


class TestFITSDiff_script(FitsTestCase):
    def test_help(self):
        with pytest.raises(SystemExit) as e:
            fitsdiff.main(["-h"])
        assert e.value.code == 0

    def test_version(self, capsys):
        with pytest.raises(SystemExit) as e:
            fitsdiff.main(["--version"])
            out = capsys.readouterr()[0]
            assert out == f"fitsdiff {version}"
        assert e.value.code == 0

    def test_noargs(self):
        with pytest.raises(SystemExit) as e:
            fitsdiff.main([""])
        assert e.value.code == 2

    def test_oneargargs(self):
        with pytest.raises(SystemExit) as e:
            fitsdiff.main(["file1"])
        assert e.value.code == 2

    def test_nodiff(self):
        a = np.arange(100).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp("testa.fits")
        tmp_b = self.temp("testb.fits")
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        numdiff = fitsdiff.main([tmp_a, tmp_b])
        assert numdiff == 0

    def test_onediff(self):
        a = np.arange(100).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        b[1, 0] = 12
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp("testa.fits")
        tmp_b = self.temp("testb.fits")
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        numdiff = fitsdiff.main([tmp_a, tmp_b])
        assert numdiff == 1

    def test_manydiff(self, capsys):
        a = np.arange(100).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a + 1
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp("testa.fits")
        tmp_b = self.temp("testb.fits")
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)

        numdiff = fitsdiff.main([tmp_a, tmp_b])
        out, err = capsys.readouterr()
        assert numdiff == 1
        assert out.splitlines()[-4:] == [
            "        a> 9",
            "        b> 10",
            "     ...",
            "     100 different pixels found (100.00% different).",
        ]

        numdiff = fitsdiff.main(["-n", "1", tmp_a, tmp_b])
        out, err = capsys.readouterr()
        assert numdiff == 1
        assert out.splitlines()[-4:] == [
            "        a> 0",
            "        b> 1",
            "     ...",
            "     100 different pixels found (100.00% different).",
        ]

    def test_outputfile(self):
        a = np.arange(100).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        b[1, 0] = 12
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp("testa.fits")
        tmp_b = self.temp("testb.fits")
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)

        numdiff = fitsdiff.main(["-o", self.temp("diff.txt"), tmp_a, tmp_b])
        assert numdiff == 1
        with open(self.temp("diff.txt")) as f:
            out = f.read()
        assert out.splitlines()[-4:] == [
            "     Data differs at [1, 2]:",
            "        a> 10",
            "        b> 12",
            "     1 different pixels found (1.00% different).",
        ]

    def test_atol(self):
        a = np.arange(100, dtype=float).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        b[1, 0] = 11
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp("testa.fits")
        tmp_b = self.temp("testb.fits")
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)

        numdiff = fitsdiff.main(["-a", "1", tmp_a, tmp_b])
        assert numdiff == 0

        numdiff = fitsdiff.main(["--exact", "-a", "1", tmp_a, tmp_b])
        assert numdiff == 1

    def test_rtol(self):
        a = np.arange(100, dtype=float).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        b[1, 0] = 11
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp("testa.fits")
        tmp_b = self.temp("testb.fits")
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        numdiff = fitsdiff.main(["-r", "1e-1", tmp_a, tmp_b])
        assert numdiff == 0

    def test_rtol_diff(self, capsys):
        a = np.arange(100, dtype=float).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        b[1, 0] = 11
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp("testa.fits")
        tmp_b = self.temp("testb.fits")
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        numdiff = fitsdiff.main(["-r", "1e-2", tmp_a, tmp_b])
        assert numdiff == 1
        out, err = capsys.readouterr()
        assert (
            out
            == f"""
 fitsdiff: {version}
 a: {tmp_a}
 b: {tmp_b}
 Maximum number of different data values to be reported: 10
 Relative tolerance: 0.01, Absolute tolerance: 0.0

Primary HDU:

   Data contains differences:
     Data differs at [1, 2]:
        a> 10.0
         ?  ^
        b> 11.0
         ?  ^
     1 different pixels found (1.00% different).
"""
        )
        assert err == ""

    def test_wildcard(self):
        tmp1 = self.temp("tmp_file1")
        with pytest.raises(SystemExit) as e:
            fitsdiff.main([tmp1 + "*", "ACME"])
        assert e.value.code == 2

    def test_not_quiet(self, capsys):
        a = np.arange(100).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp("testa.fits")
        tmp_b = self.temp("testb.fits")
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        numdiff = fitsdiff.main([tmp_a, tmp_b])
        assert numdiff == 0
        out, err = capsys.readouterr()
        assert (
            out
            == f"""
 fitsdiff: {version}
 a: {tmp_a}
 b: {tmp_b}
 Maximum number of different data values to be reported: 10
 Relative tolerance: 0.0, Absolute tolerance: 0.0

No differences found.
"""
        )
        assert err == ""

    def test_quiet(self, capsys):
        a = np.arange(100).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp("testa.fits")
        tmp_b = self.temp("testb.fits")
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        numdiff = fitsdiff.main(["-q", tmp_a, tmp_b])
        assert numdiff == 0
        out, err = capsys.readouterr()
        assert out == ""
        assert err == ""

    @pytest.mark.slow
    def test_path(self, capsys):
        os.mkdir(self.temp("sub/"))
        tmp_b = self.temp("sub/ascii.fits")

        tmp_g = self.temp("sub/group.fits")
        tmp_h = self.data("group.fits")
        with hdulist.fitsopen(tmp_h) as hdu_b:
            hdu_b.writeto(tmp_g)

        writeto(tmp_b, np.arange(100).reshape(10, 10))

        # one modified file and a directory
        assert fitsdiff.main(["-q", self.data_dir, tmp_b]) == 1
        assert fitsdiff.main(["-q", tmp_b, self.data_dir]) == 1

        # two directories
        tmp_d = self.temp("sub/")
        assert fitsdiff.main(["-q", self.data_dir, tmp_d]) == 1
        assert fitsdiff.main(["-q", tmp_d, self.data_dir]) == 1
        with pytest.warns(
            UserWarning,
            match=r"Field 'ORBPARM' has a repeat " r"count of 0 in its format code",
        ):
            assert fitsdiff.main(["-q", self.data_dir, self.data_dir]) == 0

        # no match
        tmp_c = self.data("arange.fits")
        fitsdiff.main([tmp_c, tmp_d])
        out, err = capsys.readouterr()
        assert "'arange.fits' has no match in" in err

        # globbing
        with pytest.warns(
            UserWarning,
            match=r"Field 'ORBPARM' has a repeat " r"count of 0 in its format code",
        ):
            assert fitsdiff.main(["-q", self.data_dir + "/*.fits", self.data_dir]) == 0
        assert fitsdiff.main(["-q", self.data_dir + "/g*.fits", tmp_d]) == 0

        # one file and a directory
        tmp_f = self.data("tb.fits")
        assert fitsdiff.main(["-q", tmp_f, self.data_dir]) == 0
        assert fitsdiff.main(["-q", self.data_dir, tmp_f]) == 0

    def test_ignore_hdus(self):
        a = np.arange(100).reshape(10, 10)
        b = a.copy() + 1
        ha = Header([("A", 1), ("B", 2), ("C", 3)])
        phdu_a = PrimaryHDU(header=ha)
        phdu_b = PrimaryHDU(header=ha)
        ihdu_a = ImageHDU(data=a, name="SCI")
        ihdu_b = ImageHDU(data=b, name="SCI")
        hdulist_a = HDUList([phdu_a, ihdu_a])
        hdulist_b = HDUList([phdu_b, ihdu_b])
        tmp_a = self.temp("testa.fits")
        tmp_b = self.temp("testb.fits")
        hdulist_a.writeto(tmp_a)
        hdulist_b.writeto(tmp_b)

        numdiff = fitsdiff.main([tmp_a, tmp_b])
        assert numdiff == 1

        numdiff = fitsdiff.main([tmp_a, tmp_b, "-u", "SCI"])
        assert numdiff == 0

    def test_ignore_hdus_report(self, capsys):
        a = np.arange(100).reshape(10, 10)
        b = a.copy() + 1
        ha = Header([("A", 1), ("B", 2), ("C", 3)])
        phdu_a = PrimaryHDU(header=ha)
        phdu_b = PrimaryHDU(header=ha)
        ihdu_a = ImageHDU(data=a, name="SCI")
        ihdu_b = ImageHDU(data=b, name="SCI")
        hdulist_a = HDUList([phdu_a, ihdu_a])
        hdulist_b = HDUList([phdu_b, ihdu_b])
        tmp_a = self.temp("testa.fits")
        tmp_b = self.temp("testb.fits")
        hdulist_a.writeto(tmp_a)
        hdulist_b.writeto(tmp_b)

        numdiff = fitsdiff.main([tmp_a, tmp_b, "-u", "SCI"])
        assert numdiff == 0
        out, err = capsys.readouterr()
        assert "testa.fits" in out
        assert "testb.fits" in out


@pytest.mark.skip(reason="fails intentionally to show open files (see PR #10159)")
def test_fitsdiff_openfile(tmp_path):
    """Make sure that failing FITSDiff doesn't leave open files."""
    path1 = tmp_path / "file1.fits"
    path2 = tmp_path / "file2.fits"

    hdulist = HDUList([PrimaryHDU(), ImageHDU(data=np.zeros(5))])
    hdulist.writeto(path1)
    hdulist[1].data[0] = 1
    hdulist.writeto(path2)

    diff = FITSDiff(path1, path2)
    assert diff.identical, diff.report()


class Test_FITSDiff(FitsTestCase):
    def test_FITSDiff_report(self, home_is_temp):
        self.copy_file("test0.fits")
        fits.setval(self.temp("test0.fits"), "TESTKEY", value="testval")
        d = FITSDiff(self.data("test0.fits"), self.temp("test0.fits"))
        assert not d.identical
        d.report(self.temp("diff_report.txt"))

        with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
            d.report(self.temp("diff_report.txt"), overwrite=False)
        d.report(self.temp("diff_report.txt"), overwrite=True)

        with open(os.path.expanduser(self.temp("diff_report.txt"))) as f:
            assert "Extra keyword 'TESTKEY' in b: 'testval'" in f.read()
