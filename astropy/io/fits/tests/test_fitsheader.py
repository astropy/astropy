# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from . import FitsTestCase
from astropy.io.fits.scripts import fitsheader


class TestFITSheader_script(FitsTestCase):

    def test_noargs(self):
        with pytest.raises(SystemExit) as e:
            fitsheader.main(['-h'])
        assert e.value.code == 0

    def test_file_exists(self, capsys):
        fitsheader.main([self.data('arange.fits')])
        out, err = capsys.readouterr()
        assert out.splitlines()[1].startswith(
            'SIMPLE  =                    T / conforms to FITS standard')
        assert err == ''

    def test_by_keyword(self, capsys):
        fitsheader.main(['-k', 'NAXIS', self.data('arange.fits')])
        out, err = capsys.readouterr()
        assert out.splitlines()[1].startswith(
            'NAXIS   =                    3 / number of array dimensions')

        fitsheader.main(['-k', 'NAXIS*', self.data('arange.fits')])
        out, err = capsys.readouterr()
        out = out.splitlines()
        assert len(out) == 5
        assert out[1].startswith('NAXIS')
        assert out[2].startswith('NAXIS1')
        assert out[3].startswith('NAXIS2')
        assert out[4].startswith('NAXIS3')

        fitsheader.main(['-k', 'RANDOMKEY', self.data('arange.fits')])
        out, err = capsys.readouterr()
        assert err.startswith('WARNING') and 'RANDOMKEY' in err
        assert not err.startswith('ERROR')

    def test_by_extension(self, capsys):
        fitsheader.main(['-e', '1', self.data('test0.fits')])
        out, err = capsys.readouterr()
        assert len(out.splitlines()) == 62

        fitsheader.main(['-e', '3', '-k', 'BACKGRND', self.data('test0.fits')])
        out, err = capsys.readouterr()
        assert out.splitlines()[1].startswith('BACKGRND=                 312.')

        fitsheader.main(['-e', '0', '-k', 'BACKGRND', self.data('test0.fits')])
        out, err = capsys.readouterr()
        assert err.startswith('WARNING')

        fitsheader.main(['-e', '3', '-k', 'FOO', self.data('test0.fits')])
        out, err = capsys.readouterr()
        assert err.startswith('WARNING')

    def test_table(self, capsys):
        fitsheader.main(['-t', '-k', 'BACKGRND', self.data('test0.fits')])
        out, err = capsys.readouterr()
        out = out.splitlines()
        assert len(out) == 5
        assert out[1].endswith('|   1 | BACKGRND | 316.0 |')
        assert out[2].endswith('|   2 | BACKGRND | 351.0 |')
        assert out[3].endswith('|   3 | BACKGRND | 312.0 |')
        assert out[4].endswith('|   4 | BACKGRND | 323.0 |')

        fitsheader.main(['-t', '-e', '0', '-k', 'NAXIS',
                         self.data('arange.fits'),
                         self.data('ascii.fits'),
                         self.data('blank.fits')])
        out, err = capsys.readouterr()
        out = out.splitlines()
        assert len(out) == 4
        assert out[1].endswith('|   0 |   NAXIS |     3 |')
        assert out[2].endswith('|   0 |   NAXIS |     0 |')
        assert out[3].endswith('|   0 |   NAXIS |     2 |')

    def test_fitsort(self, capsys):
        fitsheader.main(['-e', '0', '-f', '-k', 'EXPSTART', '-k', 'EXPTIME',
                         self.data('test0.fits'), self.data('test1.fits')])
        out, err = capsys.readouterr()
        out = out.splitlines()
        assert len(out) == 4
        assert out[2].endswith('test0.fits 49491.65366175    0.23')
        assert out[3].endswith('test1.fits 49492.65366175    0.22')

        fitsheader.main(['-e', '0', '-f', '-k', 'EXPSTART', '-k', 'EXPTIME',
                         self.data('test0.fits')])
        out, err = capsys.readouterr()
        out = out.splitlines()
        assert len(out) == 3
        assert out[2].endswith('test0.fits 49491.65366175    0.23')

        fitsheader.main(['-f', '-k', 'NAXIS',
                         self.data('tdim.fits'), self.data('test1.fits')])
        out, err = capsys.readouterr()
        out = out.splitlines()
        assert len(out) == 4
        assert out[0].endswith('0:NAXIS 1:NAXIS 2:NAXIS 3:NAXIS 4:NAXIS')
        assert out[2].endswith('tdim.fits       0       2      --      --      --')
        assert out[3].endswith('test1.fits       0       2       2       2       2')

        # check that files without required keyword are present
        fitsheader.main(['-f', '-k', 'DATE-OBS',
                         self.data('table.fits'), self.data('test0.fits')])
        out, err = capsys.readouterr()
        out = out.splitlines()
        assert len(out) == 4
        assert out[2].endswith('table.fits       --')
        assert out[3].endswith('test0.fits 19/05/94')

        # check that COMMENT and HISTORY are excluded
        fitsheader.main(['-e', '0', '-f', self.data('tb.fits')])
        out, err = capsys.readouterr()
        out = out.splitlines()
        assert len(out) == 3
        assert out[2].endswith('tb.fits   True     16     0   True '
                               'STScI-STSDAS/TABLES  tb.fits       1')

    def test_dotkeyword(self, capsys):
        fitsheader.main(['-e', '0', '-k', 'ESO DET ID',
                         self.data('fixed-1890.fits')])
        out, err = capsys.readouterr()
        out = out.splitlines()
        assert len(out) == 2
        assert out[1].strip().endswith("HIERARCH ESO DET ID = 'DV13' / Detector system Id")

        fitsheader.main(['-e', '0', '-k', 'ESO.DET.ID',
                         self.data('fixed-1890.fits')])
        out, err = capsys.readouterr()
        out = out.splitlines()
        assert len(out) == 2
        assert out[1].strip().endswith("HIERARCH ESO DET ID = 'DV13' / Detector system Id")
