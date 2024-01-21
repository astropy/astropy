# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

from astropy.io import fits
from astropy.utils.compat.optional_deps import HAS_MATPLOTLIB

if HAS_MATPLOTLIB:
    import matplotlib.image as mpimg

    from astropy.visualization.scripts.fits2bitmap import fits2bitmap, main


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="requires matplotlib")
class TestFits2Bitmap:
    def setup_class(self):
        self.filename = "test.fits"
        self.array = np.arange(16384).reshape((128, 128))

    def test_function(self, tmp_path):
        filename = tmp_path / self.filename
        fits.writeto(filename, self.array)
        fits2bitmap(filename)

    def test_script(self, tmp_path):
        filename = str(tmp_path / self.filename)
        fits.writeto(filename, self.array)
        main([filename, "-e", "0"])

    def test_exten_num(self, tmp_path):
        filename = str(tmp_path / self.filename)
        hdu1 = fits.PrimaryHDU()
        hdu2 = fits.ImageHDU(self.array)
        hdulist = fits.HDUList([hdu1, hdu2])
        hdulist.writeto(filename)
        main([filename, "-e", "1"])

    def test_exten_name(self, tmp_path):
        filename = str(tmp_path / self.filename)
        hdu1 = fits.PrimaryHDU()
        extname = "SCI"
        hdu2 = fits.ImageHDU(self.array)
        hdu2.header["EXTNAME"] = extname
        hdulist = fits.HDUList([hdu1, hdu2])
        hdulist.writeto(filename)
        main([filename, "-e", extname])

    @pytest.mark.parametrize("file_exten", [".gz", ".bz2"])
    def test_compressed_fits(self, tmp_path, file_exten):
        filename = str(tmp_path / f"test.fits{file_exten}")
        fits.writeto(filename, self.array)
        main([filename, "-e", "0"])

    def test_orientation(self, tmp_path):
        """
        Regression test to check the image vertical orientation/origin.
        """

        filename = str(tmp_path / self.filename)
        out_filename = "fits2bitmap_test.png"
        out_filename = str(tmp_path / out_filename)
        data = np.zeros((32, 32))
        data[0:16, :] = 1.0
        fits.writeto(filename, data)
        main([filename, "-e", "0", "-o", out_filename])

        img = mpimg.imread(out_filename)
        assert img[0, 0, 0] == 0
        assert img[31, 31, 0] == 1

    def test_min_max_cut_deprecations(self, tmp_path):
        filename = str(tmp_path / self.filename)
        fits.writeto(filename, self.array)
        with pytest.raises(SystemExit):
            main([filename, "--min_cut=0.1", "--vmin=0.1"])
        with pytest.raises(SystemExit):
            main([filename, "--max_cut=0.1", "--vmax=0.1"])
