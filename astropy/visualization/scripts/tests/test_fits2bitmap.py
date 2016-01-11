# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from ....tests.helper import pytest
from ....io import fits

try:
    import matplotlib
    HAS_MATPLOTLIB = True
    from ..fits2bitmap import fits2bitmap, main
except:
    HAS_MATPLOTLIB = False


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_fits2bitmap_function(tmpdir):
    filename = tmpdir.join('test.fits').strpath
    fits.writeto(filename, np.ones((128, 128)))
    fits2bitmap(filename)


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_fits2bitmap_script(tmpdir):
    filename = tmpdir.join('test.fits').strpath
    fits.writeto(filename, np.ones((128, 128)))
    main([filename, '-e', '0'])


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_exten_num(tmpdir):
    filename = tmpdir.join('test.fits').strpath
    data = np.ones((100, 100))
    hdu1 = fits.PrimaryHDU()
    hdu2 = fits.ImageHDU(data)
    hdulist = fits.HDUList([hdu1, hdu2])
    hdulist.writeto(filename)
    main([filename, '-e', '1'])


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_exten_name(tmpdir):
    filename = tmpdir.join('test.fits').strpath
    data = np.ones((100, 100))
    hdu1 = fits.PrimaryHDU()
    extname = 'SCI'
    hdu2 = fits.ImageHDU(data)
    hdu2.header['EXTNAME'] = extname
    hdulist = fits.HDUList([hdu1, hdu2])
    hdulist.writeto(filename)
    main([filename, '-e', extname])


@pytest.mark.skipif('not HAS_MATPLOTLIB')
@pytest.mark.parametrize('file_exten', ['.gz', '.bz2'])
def test_compressed_fits(tmpdir, file_exten):
    filename = tmpdir.join('test.fits' + file_exten).strpath
    fits.writeto(filename, np.ones((128, 128)))
    main([filename, '-e', '0'])
