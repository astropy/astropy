# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np

from astropy.io import fits

try:
    import matplotlib  # pylint: disable=W0611
    import matplotlib.image as mpimg
    HAS_MATPLOTLIB = True
    from astropy.visualization.scripts.fits2bitmap import fits2bitmap, main
except ImportError:
    HAS_MATPLOTLIB = False


@pytest.mark.skipif('not HAS_MATPLOTLIB')
class TestFits2Bitmap:
    def setup_class(self):
        self.filename = 'test.fits'
        self.array = np.arange(16384).reshape((128, 128))

    def test_function(self, tmpdir):
        filename = tmpdir.join(self.filename).strpath
        fits.writeto(filename, self.array)
        fits2bitmap(filename)

    def test_script(self, tmpdir):
        filename = tmpdir.join(self.filename).strpath
        fits.writeto(filename, self.array)
        main([filename, '-e', '0'])

    def test_exten_num(self, tmpdir):
        filename = tmpdir.join(self.filename).strpath
        hdu1 = fits.PrimaryHDU()
        hdu2 = fits.ImageHDU(self.array)
        hdulist = fits.HDUList([hdu1, hdu2])
        hdulist.writeto(filename)
        main([filename, '-e', '1'])

    def test_exten_name(self, tmpdir):
        filename = tmpdir.join(self.filename).strpath
        hdu1 = fits.PrimaryHDU()
        extname = 'SCI'
        hdu2 = fits.ImageHDU(self.array)
        hdu2.header['EXTNAME'] = extname
        hdulist = fits.HDUList([hdu1, hdu2])
        hdulist.writeto(filename)
        main([filename, '-e', extname])

    @pytest.mark.parametrize('file_exten', ['.gz', '.bz2'])
    def test_compressed_fits(self, tmpdir, file_exten):
        filename = tmpdir.join('test.fits' + file_exten).strpath
        fits.writeto(filename, self.array)
        main([filename, '-e', '0'])

    def test_orientation(self, tmpdir):
        """
        Regression test to check the image vertical orientation/origin.
        """

        filename = tmpdir.join(self.filename).strpath
        out_filename = 'fits2bitmap_test.png'
        out_filename = tmpdir.join(out_filename).strpath
        data = np.zeros((32, 32))
        data[0:16, :] = 1.
        fits.writeto(filename, data)
        main([filename, '-e', '0', '-o', out_filename])

        img = mpimg.imread(out_filename)
        assert img[0, 0, 0] == 0
        assert img[31, 31, 0] == 1
