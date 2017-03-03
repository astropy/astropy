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

# sys.argv = ['fitsheader']
# sys.argv += ['./data/test0.fits']
# sys.stdout = open('randomfile', 'w') #fixture for each method, teardown deletes it. open in self.temp directory
# fitsheader.main()
# sys.stdout.close()
# f = open('randomfile', 'r')
# assert 

# @pytest.fixture(scope=module , params=['./tests/test0.fits'])
# def resource_testfile_setup(request):
#     with open(request.params[0]) as test_fits:
#         yield test_fits

@pytest.fixture(scope='module')
def warning_regex():
    regex = re.compile('^WARNING.*')
    return regex

@pytest.fixture(scope='module')
def error_regex():
    regex = re.compile('^ERROR.*')
    return regex

class TestFITSheader_script(FitsTestCase):

    def setup_method(self,  method):
        self.sys_argv_orig = sys.argv
        sys.argv = ["fitsheader"]

    def teardown_method(self, method):
        sys.argv = self.sys_argv_orig

    def test_noargs(self):
        with pytest.raises(SystemExit) as e:
            fitsheader.main()
        assert e.value.code == 2

    @pytest.mark.parametrize('test_file', [
        ('test0.fits'),
        pytest.mark.xfail(('random.fits'), reason = 'bad file')
        
    ])
    def test_file_exists(self, test_file, capsys, warning_regex, error_regex):
        sys.argv += [self.data(test_file)]

        fitsheader.main()
        out, err = capsys.readouterr()

        assert error_regex.match(err) is None
        assert warning_regex.match(err) is None

    @pytest.mark.parametrize('test_kw,expected', [
        ('BSCALE', 'BSCALE  =           1.000000E0 / REAL = TAPE*BSCALE + BZERO'),
        pytest.mark.xfail(('LRFWAVE', 'LRFWAVE = 0.0 / linear ramp filter wavelength'), reason = 'bad format')
        
    ])
    def test_by_existing_keyword(self, test_kw, expected, capsys, warning_regex):
        sys.argv += ['-k']
        sys.argv += [test_kw]
        sys.argv += [self.data('test0.fits')]

        fitsheader.main()
        out, err = capsys.readouterr()

        assert expected in out
    
    # Different test because stderr will not be none in case of testing by existing keyword
    # fitsheader searches for a keyword in all the HDU's
    @pytest.mark.parametrize('test_kw', [
        ('RANDOMKEY')
    ])    
    def test_by_non_existing_keyword(self, test_kw, capsys, warning_regex):
        sys.argv += ['-k']
        sys.argv += [test_kw]
        sys.argv += [self.data('test0.fits')]

        fitsheader.main()
        out, err = capsys.readouterr()

        assert warning_regex.match(err) is not None

    @pytest.mark.parametrize('test_ext', [
        ('0'),
        pytest.mark.xfail(('9')),
        ('PRIMARY'),
        pytest.mark.xfail(('RANDOM'))
    ])
    def test_by_extension(self, test_ext, capsys, warning_regex, error_regex):
        sys.argv += ['-e']
        sys.argv += [test_ext]
        sys.argv += [self.data('test0.fits')]

        fitsheader.main()
        out, err = capsys.readouterr()

        # Something like -
        # "# HDU 0 in ./io/fits/tests/data/test0.fits:"
        # preceds output if extension is correct
        # Clipping that and checking if there is valid ASCII data below that
        out = out.split('\n', 1)[1]

        assert re.match('[a-zA-Z0-9_/=.]+', out) is not None
        assert warning_regex.match(err) is None
        assert error_regex.match(err) is None

    @pytest.mark.parametrize('test_ext,test_kw', [
        ('0' ,'BSCALE'),
        ('SCI', 'EXPNAME')
    ])
    def test_by_correct_extension_and_keyword(self, capsys, test_ext, test_kw, warning_regex, error_regex):
        sys.argv += ['-e']
        sys.argv += [test_ext]
        sys.argv += ['-k']
        sys.argv += [test_kw]
        sys.argv += [self.data('test0.fits')]

        fitsheader.main()
        out, err = capsys.readouterr()
        out = out.split('\n', 1)[1]

        assert re.match('[a-zA-Z0-9_/=.]+', out) is not None
        assert warning_regex.match(err) is None
        assert error_regex.match(err) is None

    @pytest.mark.parametrize('test_ext,test_kw', [
        ('0' ,'RANDOM'),
        ('RANDOM_EXT', 'EXPNAME'),
        ('9', 'EXPNAME')
    ])
    def test_by_incorrect_extension_and_keyword(self, capsys, test_ext, test_kw, warning_regex, error_regex):
        sys.argv += ['-e']
        sys.argv += [test_ext]
        sys.argv += ['-k']
        sys.argv += [test_kw]
        sys.argv += [self.data('test0.fits')]
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
        assert warning_regex.match(err) is not None
        assert error_regex.match(err) is None


