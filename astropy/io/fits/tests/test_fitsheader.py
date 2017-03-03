# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys

import numpy as np
import re

from . import FitsTestCase
from ..hdu import PrimaryHDU
from ..scripts import fitsheader
from ....tests.helper import catch_warnings
from ....tests.helper import pytest
from ....utils.exceptions import AstropyDeprecationWarning
from ....version import version
 

class TestFITSheader_script(FitsTestCase):

    script_name = 'fitsheader'
    def test_noargs(self, monkeypatch):
        monkeypatch.setattr(sys, 'argv', [self.script_name])
        with pytest.raises(SystemExit) as e:
            fitsheader.main()
        assert e.value.code == 2

    def test_file_exists(self, capsys, monkeypatch):
        args = [self.script_name, self.data('test0.fits')]
        monkeypatch.setattr(sys, 'argv', args)
        fitsheader.main()
        out, err = capsys.readouterr()

        assert not err.startswith('WARNING')
        assert not err.startswith('ERROR')

    @pytest.mark.parametrize('test_kw,expected', [
        ('BSCALE', 'BSCALE  =           1.000000E0 / REAL = TAPE*BSCALE + BZERO')
    ])
    def test_by_existing_keyword(self, test_kw, expected, capsys, monkeypatch):
        args = [self.script_name, '-k', test_kw, self.data('test0.fits')]
        monkeypatch.setattr(sys, 'argv', args)
        fitsheader.main()
        out, err = capsys.readouterr()

        assert expected in out
    
    # Different test because stderr will not be none in case of testing by existing keyword
    # fitsheader searches for a keyword in all the HDU's
  
    def test_by_non_existing_keyword(self, capsys, monkeypatch):
        args = [self.script_name, '-k', 'RANDOMKEY', self.data('test0.fits')]
        monkeypatch.setattr(sys, 'argv', args)
        fitsheader.main()
        out, err = capsys.readouterr()

        assert err.startswith('WARNING') and 'RANDOMKEY' in err
        assert not err.startswith('ERROR')

    @pytest.mark.parametrize('test_ext', [
        ('0'),
        ('PRIMARY')
    ])
    def test_by_extension(self, test_ext, capsys, monkeypatch):

        args = [self.script_name, '-e', test_ext, self.data('test0.fits')]
        monkeypatch.setattr(sys, 'argv', args)
        fitsheader.main()
        out, err = capsys.readouterr()

        # Something like -
        # "# HDU 0 in ./io/fits/tests/data/test0.fits:"
        # preceds output if extension is correct
        # Clipping that and checking if there is valid ASCII data below that
        out = out.split('\n', 1)[1]

        assert re.match('[a-zA-Z0-9_/=.]+', out) is not None
        assert not err.startswith('WARNING')
        assert not err.startswith('ERROR')

    @pytest.mark.parametrize('test_ext,test_kw', [
        ('0' ,'BSCALE'),
        ('SCI', 'EXPNAME')
    ])
    def test_by_correct_extension_and_keyword(self, capsys, test_ext, test_kw, monkeypatch):
        args = [self.script_name, '-e', test_ext, '-k', test_kw, self.data('test0.fits')]
        monkeypatch.setattr(sys, 'argv', args)
        fitsheader.main()
        out, err = capsys.readouterr()
        out = out.split('\n', 1)[1]

        assert re.match('[a-zA-Z0-9_/=.]+', out) is not None
        assert not err.startswith('WARNING')
        assert not err.startswith('ERROR')

    @pytest.mark.parametrize('test_ext,test_kw', [
        ('0' ,'RANDOM'),
        ('RANDOM_EXT', 'EXPNAME'),
        ('9', 'EXPNAME')
    ])
    def test_by_incorrect_extension_and_keyword(self, capsys, test_ext, test_kw, monkeypatch):
        args = [self.script_name, '-e', test_ext, '-k', test_kw, self.data('test0.fits')]
        monkeypatch.setattr(sys, 'argv', args)
        print(args)
        fitsheader.main()
        out, err = capsys.readouterr()

        # Something like -
        # "# HDU 0 in ./io/fits/tests/data/test0.fits:"
        # preceds output if keyword is wrong and extension is correct
        # Clipping that
        try:
            out = out.split('\n', 1)[1]
        except IndexError:
            out = ''

        assert re.match('[a-zA-Z0-9_/=.]+', out) is None
        assert err.startswith('WARNING')
        assert not err.startswith('ERROR')

    def test_table_formatting(self):
        pass


