# Licensed under a 3-clause BSD style license - see PYFITS.rst

import numpy as np

from astropy.io import fits
from astropy.io.fits.tests.test_checksum import BaseChecksumTests


class TestChecksumFunctions(BaseChecksumTests):
    # All checksums have been verified against CFITSIO

    def test_compressed_image_data(self):
        with fits.open(self.data("comp.fits")) as h1:
            h1.writeto(self.temp("tmp.fits"), overwrite=True, checksum=True)

            with fits.open(
                self.temp("tmp.fits"), checksum=True, disable_image_compression=True
            ) as h2:
                assert "CHECKSUM" in h2[0].header
                assert h2[0].header["CHECKSUM"] == "D8iBD6ZAD6fAD6ZA"
                assert "DATASUM" in h2[0].header
                assert h2[0].header["DATASUM"] == "0"
                assert "CHECKSUM" in h2[1].header
                assert h2[1].header["CHECKSUM"] == "ZeAbdb8aZbAabb7a"
                assert "DATASUM" in h2[1].header
                assert h2[1].header["DATASUM"] == "113055149"
                assert h2[1].header._countblanks() == 78
                assert "ZHECKSUM" not in h2[1].header
                assert "ZDATASUM" not in h2[1].header

            with fits.open(self.temp("tmp.fits"), checksum=True) as h2:
                assert np.all(h1[1].data == h2[1].data)
                assert "CHECKSUM" in h2[0].header
                assert h2[0].header["CHECKSUM"] == "D8iBD6ZAD6fAD6ZA"
                assert "DATASUM" in h2[0].header
                assert h2[0].header["DATASUM"] == "0"
                assert "CHECKSUM" not in h2[1].header
                assert "DATASUM" not in h2[1].header
                assert h2[1].header._countblanks() == 78

    def test_failing_compressed_datasum(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/4587
        """
        n = np.ones((10, 10), dtype="float32")
        comp_hdu = fits.CompImageHDU(n)
        comp_hdu.writeto(self.temp("tmp.fits"), checksum=True)

        with fits.open(self.temp("tmp.fits"), checksum=True) as hdul:
            assert np.all(hdul[1].data == comp_hdu.data)

    def test_compressed_image_data_int16(self):
        n = np.arange(100, dtype="int16")
        hdu = fits.ImageHDU(n)
        comp_hdu = fits.CompImageHDU(hdu.data, hdu.header)
        comp_hdu.writeto(self.temp("tmp.fits"), checksum=True)
        hdu.writeto(self.temp("uncomp.fits"), checksum=True)

        with (
            fits.open(self.temp("tmp.fits"), checksum=True) as hdul,
            fits.open(
                self.temp("tmp.fits"), checksum=True, disable_image_compression=True
            ) as hdul_comp,
            fits.open(self.temp("uncomp.fits"), checksum=True) as hdul_uncomp,
        ):
            # Check that the data was decompressed properly
            assert np.all(hdul[1].data == comp_hdu.data)
            assert np.all(hdul[1].data == hdu.data)

            # Check primary HDU checksums
            assert "CHECKSUM" in hdul[0].header
            assert hdul[0].header["CHECKSUM"] == "D8iBD6ZAD6fAD6ZA"
            assert "DATASUM" in hdul[0].header
            assert hdul[0].header["DATASUM"] == "0"

            # Check first HDU checksum in compressed data
            assert "CHECKSUM" in hdul_comp[1].header
            assert hdul_comp[1].header["CHECKSUM"] == "g4NEg2LEg2LEg2LE"
            assert "DATASUM" in hdul_comp[1].header
            assert hdul_comp[1].header["DATASUM"] == "2453673070"
            assert "CHECKSUM" in hdul_comp[1].header
            assert "ZHECKSUM" not in hdul_comp[1].header
            assert "ZDATASUM" not in hdul_comp[1].header

            # Check no checksums in decompressed data
            assert "CHECKSUM" not in hdul[1].header
            assert "DATASUM" not in hdul[1].header

            # Check checksums in original image data
            assert "CHECKSUM" in hdul_uncomp[1].header
            assert hdul_uncomp[1].header["CHECKSUM"] == "ZE94eE91ZE91bE91"
            assert "DATASUM" in hdul_uncomp[1].header
            assert hdul_uncomp[1].header["DATASUM"] == "160565700"

            # If in future we make it easy to compute the decompressed image
            # checksum in CompImageHDU.header, the following test should pass
            # assert header_comp["ZHECKSUM"] == header_uncomp["CHECKSUM"]
            # assert header_comp["ZDATASUM"] == header_uncomp["DATASUM"]

    def test_compressed_image_data_float32(self):
        n = np.arange(100, dtype="float32")
        hdu = fits.ImageHDU(n)
        comp_hdu = fits.CompImageHDU(hdu.data, hdu.header)
        comp_hdu.writeto(self.temp("tmp.fits"), checksum=True)
        hdu.writeto(self.temp("uncomp.fits"), checksum=True)

        with (
            fits.open(self.temp("tmp.fits"), checksum=True) as hdul,
            fits.open(
                self.temp("tmp.fits"), checksum=True, disable_image_compression=True
            ) as hdul_comp,
            fits.open(self.temp("uncomp.fits"), checksum=True) as hdul_uncomp,
        ):
            # Check that the data was decompressed properly
            assert np.all(hdul[1].data == comp_hdu.data)
            assert np.all(hdul[1].data == hdu.data)

            # Check primary HDU checksums
            assert "CHECKSUM" in hdul[0].header
            assert hdul[0].header["CHECKSUM"] == "D8iBD6ZAD6fAD6ZA"
            assert "DATASUM" in hdul[0].header
            assert hdul[0].header["DATASUM"] == "0"

            # Check first HDU checksum in compressed data
            assert "CHECKSUM" in hdul_comp[1].header
            assert "DATASUM" in hdul_comp[1].header

            # The checksum ends up being different on Windows and s390/bigendian,
            # possibly due to slight floating point differences? See gh-10921.
            # TODO fix these so they work on all platforms; otherwise pointless.
            # assert hdul_comp[1].header['CHECKSUM'] == 'eATIf3SHe9SHe9SH'
            # assert hdul_comp[1].header['DATASUM'] == '1277667818'

            # Check no checksums in decompressed data
            assert "CHECKSUM" not in hdul[1].header
            assert "DATASUM" not in hdul[1].header

            # Check checksums in original image data
            assert "CHECKSUM" in hdul_uncomp[1].header
            assert hdul_uncomp[1].header["CHECKSUM"] == "Cgr5FZo2Cdo2CZo2"
            assert "DATASUM" in hdul_uncomp[1].header
            assert hdul_uncomp[1].header["DATASUM"] == "2393636889"

            # If in future we make it easy to compute the decompressed image
            # checksum in CompImageHDU.header, the following test should pass
            # assert header_comp["ZHECKSUM"] == header_uncomp["CHECKSUM"]
            # assert header_comp["ZDATASUM"] == header_uncomp["DATASUM"]
