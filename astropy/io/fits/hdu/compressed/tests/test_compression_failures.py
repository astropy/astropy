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

    @pytest.mark.parametrize(
        ("numpy_dtype", "fits_format", "big_endian_dtype"),
        [
            (np.int16, "1PI", ">i2"),
            (np.float32, "1PE", ">f4"),
        ],
    )
    def test_uncompressed_data_column(
        self, tmp_path, numpy_dtype, fits_format, big_endian_dtype
    ):
        """
        Regression test for https://github.com/astropy/astropy/issues/15477

        Test that data stored in UNCOMPRESSED_DATA column is correctly read.
        UNCOMPRESSED_DATA columns can have integer (PI, PJ) or floating-point
        (PE, PD) formats. This tests the case where some tiles have their data
        in UNCOMPRESSED_DATA instead of COMPRESSED_DATA (which can happen when
        compression fails for certain tiles).
        """
        tile_size = 5
        original_data = np.arange(100, dtype=numpy_dtype).reshape((10, 10))
        hdu = fits.CompImageHDU(
            data=original_data,
            compression_type="RICE_1",
            tile_shape=(tile_size, tile_size),
        )

        filename = tmp_path / "test_uncompressed_tiles.fits"
        hdu.writeto(filename)

        # Read and modify: put first tile's data in UNCOMPRESSED_DATA instead.
        # This simulates a file created by another tool (like cfitsio directly).
        with fits.open(filename, disable_image_compression=True) as hdul:
            orig_table = hdul[1].data
            orig_header = hdul[1].header.copy()
            orig_cols = hdul[1].columns
            n_tiles = len(orig_table)

            compressed_col_data = np.empty(n_tiles, dtype=object)
            uncompressed_col_data = np.empty(n_tiles, dtype=object)

            for i in range(n_tiles):
                if i == 0:
                    # First tile: put in UNCOMPRESSED_DATA, leave COMPRESSED_DATA empty
                    compressed_col_data[i] = np.array([], dtype=np.uint8)
                    tile_data = original_data[:tile_size, :tile_size].flatten()
                    uncompressed_col_data[i] = tile_data.astype(big_endian_dtype)
                else:
                    # Other tiles: keep in COMPRESSED_DATA
                    compressed_col_data[i] = orig_table["COMPRESSED_DATA"][i]
                    uncompressed_col_data[i] = np.array([], dtype=numpy_dtype)

            # Build column list: COMPRESSED_DATA, UNCOMPRESSED_DATA, plus any
            # other columns from the original (ZSCALE, ZZERO for float data)
            new_col_list = [
                fits.Column(
                    name="COMPRESSED_DATA",
                    format=orig_cols["COMPRESSED_DATA"].format,
                    array=compressed_col_data,
                ),
                fits.Column(
                    name="UNCOMPRESSED_DATA",
                    format=f"{fits_format}({tile_size * tile_size})",
                    array=uncompressed_col_data,
                ),
            ]

            # Preserve other columns (ZSCALE, ZZERO, etc.) needed for quantized data
            for col in orig_cols:
                if col.name not in ("COMPRESSED_DATA", "GZIP_COMPRESSED_DATA"):
                    new_col_list.append(
                        fits.Column(
                            name=col.name,
                            format=col.format,
                            array=orig_table[col.name],
                        )
                    )

            new_cols = fits.ColDefs(new_col_list)
            new_table = fits.BinTableHDU.from_columns(new_cols, header=orig_header)
            modified_filename = tmp_path / "test_uncompressed_tiles_modified.fits"
            fits.HDUList([hdul[0].copy(), new_table]).writeto(modified_filename)

        # Verify the data is correctly read. Use atol=0.1 to account for
        # quantization noise in float data (ZSCALE is typically ~0.13 for this data)
        with fits.open(modified_filename) as hdul:
            np.testing.assert_allclose(hdul[1].data, original_data, atol=0.1)
