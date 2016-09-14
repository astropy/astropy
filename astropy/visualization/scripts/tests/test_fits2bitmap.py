# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from ....tests.helper import pytest
from ....io import fits

try:
    import matplotlib  # pylint: disable=W0611
    HAS_MATPLOTLIB = True
    from ..fits2bitmap import fits2bitmap, main
except ImportError:
    HAS_MATPLOTLIB = False


@pytest.mark.skipif('not HAS_MATPLOTLIB')
class TestFits2Bitmap(object):
    def setup_class(self):
        self.filename = 'test.fits'

    def test_function(self, tmpdir):
        filename = tmpdir.join(self.filename).strpath
        fits.writeto(filename, np.ones((128, 128)))
        fits2bitmap(filename)

    def test_script(self, tmpdir):
        filename = tmpdir.join(self.filename).strpath
        fits.writeto(filename, np.ones((128, 128)))
        main([filename, '-e', '0'])

    def test_exten_num(self, tmpdir):
        filename = tmpdir.join(self.filename).strpath
        data = np.ones((100, 100))
        hdu1 = fits.PrimaryHDU()
        hdu2 = fits.ImageHDU(data)
        hdulist = fits.HDUList([hdu1, hdu2])
        hdulist.writeto(filename)
        main([filename, '-e', '1'])

    def test_exten_name(self, tmpdir):
        filename = tmpdir.join(self.filename).strpath
        data = np.ones((100, 100))
        hdu1 = fits.PrimaryHDU()
        extname = 'SCI'
        hdu2 = fits.ImageHDU(data)
        hdu2.header['EXTNAME'] = extname
        hdulist = fits.HDUList([hdu1, hdu2])
        hdulist.writeto(filename)
        main([filename, '-e', extname])

    @pytest.mark.parametrize('file_exten', ['.gz', '.bz2'])
    def test_compressed_fits(self, tmpdir, file_exten):
        filename = tmpdir.join('test.fits' + file_exten).strpath
        fits.writeto(filename, np.ones((128, 128)))
        main([filename, '-e', '0'])
