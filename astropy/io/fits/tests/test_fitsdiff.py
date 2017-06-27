# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys

import numpy as np

from . import FitsTestCase
from ..hdu import PrimaryHDU
from ..scripts import fitsdiff
from ....tests.helper import catch_warnings
from ....tests.helper import pytest
from ....utils.exceptions import AstropyDeprecationWarning
from ....version import version


class TestFITSDiff_script(FitsTestCase):
    def setup_method(self, method):
        self.sys_argv_orig = sys.argv
        sys.argv = ["fitsdiff"]

    def teardown_method(self, method):
        sys.argv = self.sys_argv_orig

    def test_noargs(self):
        with pytest.raises(SystemExit) as e:
            fitsdiff.main()
        assert e.value.code == 2

    def test_oneargargs(self):
        testargs = ["file1"]
        sys.argv += testargs
        with pytest.raises(SystemExit) as e:
            fitsdiff.main()
        assert e.value.code == 2

    def test_nodiff(self):
        a = np.arange(100).reshape((10, 10))
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        testargs = [tmp_a, tmp_b]
        sys.argv += testargs
        numdiff = fitsdiff.main()
        assert numdiff == 0

    def test_onediff(self):
        a = np.arange(100).reshape((10, 10))
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        b[1, 0] = 12
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        testargs = [tmp_a, tmp_b]
        sys.argv += testargs
        numdiff = fitsdiff.main()
        assert numdiff == 1

    def test_rtol(self):
        a = np.arange(100, dtype=np.float).reshape((10, 10))
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        b[1, 0] = 11
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        testargs = ["-r", "1e-1", tmp_a, tmp_b]
        sys.argv += testargs
        numdiff = fitsdiff.main()
        assert numdiff == 0

    def test_rtol_diff(self, capsys):
        a = np.arange(100, dtype=np.float).reshape((10, 10))
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        b[1, 0] = 11
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        testargs = ["-r", "1e-2", tmp_a, tmp_b]
        sys.argv += testargs
        numdiff = fitsdiff.main()
        assert numdiff == 1
        out, err = capsys.readouterr()
        assert out == """
 fitsdiff: {}
 a: {}
 b: {}
 Maximum number of different data values to be reported: 10
 Relative tolerance: 0.01, Absolute tolerance: 0.0

Primary HDU:\n\n   Data contains differences:
     Data differs at [1, 2]:
        a> 10.0
         ?  ^
        b> 11.0
         ?  ^
     1 different pixels found (1.00% different).\n""".format(version, tmp_a, tmp_b)
        assert err == ""

    def test_fitsdiff_script_both_d_and_r(self, capsys):
        a = np.arange(100).reshape((10, 10))
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        testargs = ["-r", "1e-4", "-d", "1e-2", tmp_a, tmp_b]
        sys.argv += testargs
        with catch_warnings(AstropyDeprecationWarning) as warning_lines:
            fitsdiff.main()
            # `rtol` is always ignored when `tolerance` is provided
            assert warning_lines[0].category == AstropyDeprecationWarning
            assert (str(warning_lines[0].message) ==
                    '"-d" ("--difference-tolerance") was deprecated in version 2.0 '
                    'and will be removed in a future version. '
                    'Use "-r" ("--relative-tolerance") instead.')
        out, err = capsys.readouterr()
        assert out == """
 fitsdiff: {}
 a: {}
 b: {}
 Maximum number of different data values to be reported: 10
 Relative tolerance: 0.01, Absolute tolerance: 0.0

No differences found.\n""".format(version, tmp_a, tmp_b)

    def test_wildcard(self):
        tmp1 = self.temp("tmp_file1")
        testargs = [tmp1+"*", "ACME"]
        sys.argv += testargs
        with pytest.raises(SystemExit) as e:
            fitsdiff.main()
        assert e.value.code == 2

    def test_not_quiet(self, capsys):
        a = np.arange(100).reshape((10, 10))
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        testargs = [tmp_a, tmp_b]
        sys.argv += testargs
        numdiff = fitsdiff.main()
        assert numdiff == 0
        out, err = capsys.readouterr()
        assert out == """
 fitsdiff: {}
 a: {}
 b: {}
 Maximum number of different data values to be reported: 10
 Relative tolerance: 0.0, Absolute tolerance: 0.0

No differences found.\n""".format(version, tmp_a, tmp_b)
        assert err == ""

    def test_quiet(self, capsys):
        a = np.arange(100).reshape((10, 10))
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        testargs = ["-q", tmp_a, tmp_b]
        sys.argv += testargs
        numdiff = fitsdiff.main()
        assert numdiff == 0
        out, err = capsys.readouterr()
        assert out == ""
        assert err == ""
