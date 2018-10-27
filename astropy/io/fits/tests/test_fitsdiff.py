# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
import os

from . import FitsTestCase
from ..convenience import writeto
from ..hdu import PrimaryHDU, hdulist
from .. import Header, ImageHDU, HDUList
from ..scripts import fitsdiff
from ....tests.helper import catch_warnings
from ....utils.exceptions import AstropyDeprecationWarning
from ....version import version


class TestFITSDiff_script(FitsTestCase):

    def test_noargs(self):
        with pytest.raises(SystemExit) as e:
            fitsdiff.main()
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
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
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
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        numdiff = fitsdiff.main([tmp_a, tmp_b])
        assert numdiff == 1

    def test_manydiff(self, capsys):
        a = np.arange(100).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a + 1
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)

        numdiff = fitsdiff.main([tmp_a, tmp_b])
        out, err = capsys.readouterr()
        assert numdiff == 1
        assert out.splitlines()[-4:] == [
            '        a> 9',
            '        b> 10',
            '     ...',
            '     100 different pixels found (100.00% different).']

        numdiff = fitsdiff.main(['-n', '1', tmp_a, tmp_b])
        out, err = capsys.readouterr()
        assert numdiff == 1
        assert out.splitlines()[-4:] == [
            '        a> 0',
            '        b> 1',
            '     ...',
            '     100 different pixels found (100.00% different).']

    def test_outputfile(self):
        a = np.arange(100).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        b[1, 0] = 12
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)

        numdiff = fitsdiff.main(['-o', self.temp('diff.txt'), tmp_a, tmp_b])
        assert numdiff == 1
        with open(self.temp('diff.txt')) as f:
            out = f.read()
        assert out.splitlines()[-4:] == [
            '     Data differs at [1, 2]:',
            '        a> 10',
            '        b> 12',
            '     1 different pixels found (1.00% different).']

    def test_atol(self):
        a = np.arange(100, dtype=float).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        b[1, 0] = 11
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
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
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
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
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        numdiff = fitsdiff.main(["-r", "1e-2", tmp_a, tmp_b])
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
        a = np.arange(100).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        with catch_warnings(AstropyDeprecationWarning) as warning_lines:
            fitsdiff.main(["-r", "1e-4", "-d", "1e-2", tmp_a, tmp_b])
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
        with pytest.raises(SystemExit) as e:
            fitsdiff.main([tmp1+"*", "ACME"])
        assert e.value.code == 2

    def test_not_quiet(self, capsys):
        a = np.arange(100).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        numdiff = fitsdiff.main([tmp_a, tmp_b])
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
        a = np.arange(100).reshape(10, 10)
        hdu_a = PrimaryHDU(data=a)
        b = a.copy()
        hdu_b = PrimaryHDU(data=b)
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdu_a.writeto(tmp_a)
        hdu_b.writeto(tmp_b)
        numdiff = fitsdiff.main(["-q", tmp_a, tmp_b])
        assert numdiff == 0
        out, err = capsys.readouterr()
        assert out == ""
        assert err == ""

    def test_path(self, capsys):
        os.mkdir(self.temp('sub/'))
        tmp_b = self.temp('sub/ascii.fits')

        tmp_g = self.temp('sub/group.fits')
        tmp_h = self.data('group.fits')
        with hdulist.fitsopen(tmp_h) as hdu_b:
            hdu_b.writeto(tmp_g)

        writeto(tmp_b, np.arange(100).reshape(10, 10))

        # one modified file and a directory
        assert fitsdiff.main(["-q", self.data_dir, tmp_b]) == 1
        assert fitsdiff.main(["-q", tmp_b, self.data_dir]) == 1

        # two directories
        tmp_d = self.temp('sub/')
        assert fitsdiff.main(["-q", self.data_dir, tmp_d]) == 1
        assert fitsdiff.main(["-q", tmp_d, self.data_dir]) == 1
        with pytest.warns(UserWarning, match="Field 'ORBPARM' has a repeat "
                          "count of 0 in its format code"):
            assert fitsdiff.main(["-q", self.data_dir, self.data_dir]) == 0

        # no match
        tmp_c = self.data('arange.fits')
        fitsdiff.main([tmp_c, tmp_d])
        out, err = capsys.readouterr()
        assert "'arange.fits' has no match in" in err

        # globbing
        with pytest.warns(UserWarning, match="Field 'ORBPARM' has a repeat "
                          "count of 0 in its format code"):
            assert fitsdiff.main(["-q", self.data_dir+'/*.fits',
                                  self.data_dir]) == 0
        assert fitsdiff.main(["-q", self.data_dir+'/g*.fits', tmp_d]) == 0

        # one file and a directory
        tmp_f = self.data('tb.fits')
        assert fitsdiff.main(["-q", tmp_f, self.data_dir]) == 0
        assert fitsdiff.main(["-q", self.data_dir, tmp_f]) == 0

    def test_ignore_hdus(self):
        a = np.arange(100).reshape(10, 10)
        b = a.copy() + 1
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        phdu_a = PrimaryHDU(header=ha)
        phdu_b = PrimaryHDU(header=ha)
        ihdu_a = ImageHDU(data=a, name='SCI')
        ihdu_b = ImageHDU(data=b, name='SCI')
        hdulist_a = HDUList([phdu_a, ihdu_a])
        hdulist_b = HDUList([phdu_b, ihdu_b])
        tmp_a = self.temp('testa.fits')
        tmp_b = self.temp('testb.fits')
        hdulist_a.writeto(tmp_a)
        hdulist_b.writeto(tmp_b)

        numdiff = fitsdiff.main([tmp_a, tmp_b])
        assert numdiff == 1

        numdiff = fitsdiff.main([tmp_a, tmp_b, "-u", "SCI"])
        assert numdiff == 0
