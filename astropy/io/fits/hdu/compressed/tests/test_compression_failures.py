# Licensed under a 3-clause BSD style license - see PYFITS.rst

import numpy as np
import pytest

from astropy.io import fits
from astropy.io.fits.hdu.compressed._tiled_compression import compress_image_data
from astropy.io.fits.tests.conftest import FitsTestCase

MAX_INT = np.iinfo(np.intc).max
MAX_LONG = np.iinfo(int).max
MAX_LONGLONG = np.iinfo(np.longlong).max


class TestCompressionFunction(FitsTestCase):
    def test_wrong_argument_number(self):
        with pytest.raises(TypeError):
            compress_image_data(1, 2)

    def test_unknown_compression_type(self):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        bintable = hdu._get_bintable_without_data()
        bintable.header["ZCMPTYPE"] = "fun"
        with pytest.raises(ValueError) as exc:
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )
        assert "Unrecognized compression type: fun" in str(exc.value)

    def test_zbitpix_unknown(self):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        bintable = hdu._get_bintable_without_data()
        bintable.header["ZBITPIX"] = 13
        with pytest.raises(ValueError) as exc:
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )
        assert "Invalid value for BITPIX: 13" in str(exc.value)

    def test_data_none(self):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        hdu.data = None
        hdu.tile_shape = None
        bintable = hdu._get_bintable_without_data()
        with pytest.raises(TypeError) as exc:
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )
        assert "Image data must be a numpy.ndarray" in str(exc.value)

    def test_invalid_tform(self):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        bintable = hdu._get_bintable_without_data()
        bintable.header["TFORM1"] = "TX"
        with pytest.raises(RuntimeError) as exc:
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )
        assert "TX" in str(exc.value) and "TFORM" in str(exc.value)

    def test_invalid_zdither(self):
        hdu = fits.CompImageHDU(np.ones((10, 10)), quantize_method=1)
        bintable = hdu._get_bintable_without_data()
        bintable.header["ZDITHER0"] = "a"
        with pytest.raises(TypeError):
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )

    @pytest.mark.parametrize("kw", ["ZNAXIS", "ZBITPIX"])
    def test_header_missing_keyword(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        bintable = hdu._get_bintable_without_data()
        del bintable.header[kw]
        with pytest.raises(KeyError) as exc:
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )
        assert kw in str(exc.value)

    @pytest.mark.parametrize("kw", ["ZNAXIS", "ZVAL1", "ZVAL2", "ZBLANK", "BLANK"])
    def test_header_value_int_overflow(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        bintable = hdu._get_bintable_without_data()
        bintable.header[kw] = MAX_INT + 1
        with pytest.raises(OverflowError):
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )

    @pytest.mark.parametrize("kw", ["ZTILE1", "ZNAXIS1"])
    def test_header_value_long_overflow(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        bintable = hdu._get_bintable_without_data()
        bintable.header[kw] = MAX_LONG + 1
        with pytest.raises(OverflowError):
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )

    @pytest.mark.parametrize("kw", ["NAXIS1", "NAXIS2", "TNULL1", "PCOUNT", "THEAP"])
    def test_header_value_longlong_overflow(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        bintable = hdu._get_bintable_without_data()
        bintable.header[kw] = MAX_LONGLONG + 1
        with pytest.raises(OverflowError):
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )

    @pytest.mark.parametrize("kw", ["ZVAL3"])
    def test_header_value_float_overflow(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        bintable = hdu._get_bintable_without_data()
        bintable.header[kw] = 1e300
        with pytest.raises(OverflowError):
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )

    @pytest.mark.parametrize("kw", ["NAXIS1", "NAXIS2", "TFIELDS", "PCOUNT"])
    def test_header_value_negative(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        bintable = hdu._get_bintable_without_data()
        bintable.header[kw] = -1
        with pytest.raises(ValueError) as exc:
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )
        assert f"{kw} should not be negative." in str(exc.value)

    @pytest.mark.parametrize(("kw", "limit"), [("ZNAXIS", 999), ("TFIELDS", 999)])
    def test_header_value_exceeds_custom_limit(self, kw, limit):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        bintable = hdu._get_bintable_without_data()
        bintable.header[kw] = limit + 1
        with pytest.raises(ValueError) as exc:
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )
        assert kw in str(exc.value)

    @pytest.mark.parametrize(
        "kw", ["TTYPE1", "TFORM1", "ZCMPTYPE", "ZNAME1", "ZQUANTIZ"]
    )
    def test_header_value_no_string(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        bintable = hdu._get_bintable_without_data()
        bintable.header[kw] = 1
        with pytest.raises(TypeError):
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )

    @pytest.mark.parametrize("kw", ["TZERO1", "TSCAL1"])
    def test_header_value_no_double(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        bintable = hdu._get_bintable_without_data()
        bintable.header[kw] = "1"
        with pytest.raises(TypeError):
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )

    @pytest.mark.parametrize("kw", ["ZSCALE", "ZZERO"])
    def test_header_value_no_double_int_image(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10), dtype=np.int32))
        bintable = hdu._get_bintable_without_data()
        bintable.header[kw] = "1"
        with pytest.raises(TypeError):
            compress_image_data(
                hdu.data, hdu.compression_type, bintable.header, bintable.columns
            )
