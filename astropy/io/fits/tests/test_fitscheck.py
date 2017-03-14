import sys

import numpy as np
import re

from . import FitsTestCase
from ..hdu import PrimaryHDU
import astropy.io.fits as fits
from ..scripts import fitscheck
from ....tests.helper import catch_warnings
from ....tests.helper import pytest
from ....utils.exceptions import AstropyDeprecationWarning
from ....version import version


class TestFITSCheck_script(FitsTestCase):

    script_name = 'fitscheck'

    def test_forced_write(self, capsys, monkeypatch):

        self.copy_file('test0.fits')
        args = [self.script_name, '--force', '--write', self.temp('test0.fits')]
        
        with fits.open(self.data('test0.fits')) as hdulist:
            with pytest.raises(KeyError):
                hdulist[0].header['CHECKSUM']
            with pytest.raises(KeyError): 
                hdulist[0].header['DATASUM']

        monkeypatch.setattr(sys, 'argv', args)

        fitscheck.main()
        # out, err = capsys.readouterr()
        # print(out)

        with fits.open(self.temp('test0.fits')) as hdulist:
            for hdu in hdulist:
                assert hdu.header['CHECKSUM']
                assert hdu.header['DATASUM']

    def test_write_nonexistent_without_force(self, capsys, monkeypatch):
        '''
        User should not be able to add/update CHECKSUM and DATASUM if :

        CHECKSUM and DATASUM are non-existent/erroneous AND only --write flag is used
          --force flag is necessary in such a case
        '''
        
        self.copy_file('test0.fits')
        args = [self.script_name, '--write', self.temp('test0.fits')]

        monkeypatch.setattr(sys, 'argv', args)

        fitscheck.main()
        out, err = capsys.readouterr()

        assert '1 errors' in err
        assert not out

        with fits.open(self.temp('test0.fits')) as hdulist:
            with pytest.raises(KeyError):
                hdulist[0].header['CHECKSUM']
            with pytest.raises(KeyError): 
                hdulist[0].header['DATASUM']

    def test_delete_existing_keywords(self, capsys, monkeypatch):
        '''
        Regression test for https://github.com/astropy/astropy/issues/5874

        User should be able to delete CHECKSUM and DATASUM keywords from all HDUs 

        fitscheck --checksum none --write filename
        '''

        self.copy_file('checksum.fits')
        args = [self.script_name, '-k', 'none', '--write', self.temp('checksum.fits')]

        monkeypatch.setattr(sys, 'argv', args)

        with fits.open(self.temp('checksum.fits')) as hdulist:
            assert hdulist[0].header['CHECKSUM']
            assert hdulist[0].header['DATASUM']

        fitscheck.main()
        out, err = capsys.readouterr()
        assert not err
        with fits.open(self.temp('checksum.fits')) as hdulist:
            with pytest.raises(KeyError):
                hdulist[0].header['CHECKSUM']
            with pytest.raises(KeyError): 
                hdulist[0].header['DATASUM']

    def test_check_compliance_only(self, capsys, monkeypatch):
        pass