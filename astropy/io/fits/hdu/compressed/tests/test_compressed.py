# Licensed under a 3-clause BSD style license - see PYFITS.rst

import io
import math
import os
import pickle
import re
import time
from io import BytesIO
from itertools import product

import numpy as np
import pytest
from hypothesis import given
from hypothesis.extra.numpy import basic_indices
from numpy.testing import assert_allclose, assert_equal

from astropy.io import fits
from astropy.io.fits.hdu.compressed import (
    COMPRESSION_TYPES,
    DITHER_SEED_CHECKSUM,
    SUBTRACTIVE_DITHER_1,
)
from astropy.io.fits.tests.conftest import FitsTestCase
from astropy.io.fits.tests.test_table import comparerecords
from astropy.io.fits.verify import VerifyWarning
from astropy.utils.data import download_file
from astropy.utils.misc import NumpyRNGContext


class TestCompressedImage(FitsTestCase):
    def test_empty(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/2595
        """

        hdu = fits.CompImageHDU()
        assert hdu.data is None
        hdu.writeto(self.temp("test.fits"))

        with fits.open(self.temp("test.fits"), mode="update") as hdul:
            assert len(hdul) == 2
            assert isinstance(hdul[1], fits.CompImageHDU)
            assert hdul[1].data is None

            # Now test replacing the empty data with an array and see what
            # happens
            hdul[1].data = np.arange(100, dtype=np.int32)

        with fits.open(self.temp("test.fits")) as hdul:
            assert len(hdul) == 2
            assert isinstance(hdul[1], fits.CompImageHDU)
            assert np.all(hdul[1].data == np.arange(100, dtype=np.int32))

    @pytest.mark.parametrize(
        ("data", "compression_type", "quantize_level"),
        [
            (np.zeros((2, 10, 10), dtype=np.float32), "RICE_1", 16),
            (np.zeros((2, 10, 10), dtype=np.float32), "GZIP_1", -0.01),
            (np.zeros((2, 10, 10), dtype=np.float32), "GZIP_2", -0.01),
            (np.zeros((100, 100)) + 1, "HCOMPRESS_1", 16),
            (np.zeros((10, 10)), "PLIO_1", 16),
        ],
    )
    @pytest.mark.parametrize("byte_order", ["<", ">"])
    def test_comp_image(self, data, compression_type, quantize_level, byte_order):
        data = data.view(data.dtype.newbyteorder(byte_order))
        primary_hdu = fits.PrimaryHDU()
        ofd = fits.HDUList(primary_hdu)
        chdu = fits.CompImageHDU(
            data,
            name="SCI",
            compression_type=compression_type,
            quantize_level=quantize_level,
        )
        ofd.append(chdu)
        ofd.writeto(self.temp("test_new.fits"), overwrite=True)
        ofd.close()
        with fits.open(self.temp("test_new.fits")) as fd:
            assert (fd[1].data == data).all()
            assert fd[1].header["NAXIS"] == chdu.header["NAXIS"]
            assert fd[1].header["NAXIS1"] == chdu.header["NAXIS1"]
            assert fd[1].header["NAXIS2"] == chdu.header["NAXIS2"]
            assert fd[1].header["BITPIX"] == chdu.header["BITPIX"]

    @pytest.mark.remote_data
    def test_comp_image_quantize_level(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/5969

        Test that quantize_level is used.

        """

        np.random.seed(42)

        # Basically what scipy.datasets.ascent() does.
        fname = download_file(
            "https://github.com/scipy/dataset-ascent/blob/main/ascent.dat?raw=true"
        )
        with open(fname, "rb") as f:
            scipy_data = np.array(pickle.load(f))

        data = scipy_data + np.random.randn(512, 512) * 10

        fits.ImageHDU(data).writeto(self.temp("im1.fits"))
        fits.CompImageHDU(
            data,
            compression_type="RICE_1",
            quantize_method=1,
            quantize_level=-1,
            dither_seed=5,
        ).writeto(self.temp("im2.fits"))
        fits.CompImageHDU(
            data,
            compression_type="RICE_1",
            quantize_method=1,
            quantize_level=-100,
            dither_seed=5,
        ).writeto(self.temp("im3.fits"))

        im1 = fits.getdata(self.temp("im1.fits"))
        im2 = fits.getdata(self.temp("im2.fits"))
        im3 = fits.getdata(self.temp("im3.fits"))

        assert not np.array_equal(im2, im3)
        assert np.isclose(np.min(im1 - im2), -0.5, atol=1e-3)
        assert np.isclose(np.max(im1 - im2), 0.5, atol=1e-3)
        assert np.isclose(np.min(im1 - im3), -50, atol=1e-1)
        assert np.isclose(np.max(im1 - im3), 50, atol=1e-1)

    def test_comp_image_hcompression_1_invalid_data(self):
        """
        Tests compression with the HCOMPRESS_1 algorithm with data that is
        not 2D and has a non-2D tile size.
        """

        pytest.raises(
            ValueError,
            fits.CompImageHDU,
            np.zeros((2, 10, 10), dtype=np.float32),
            name="SCI",
            compression_type="HCOMPRESS_1",
            quantize_level=16,
            tile_shape=(2, 10, 10),
        )

    def test_comp_image_hcompress_image_stack(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/171

        Tests that data containing more than two dimensions can be
        compressed with HCOMPRESS_1 so long as the user-supplied tile size can
        be flattened to two dimensions.
        """

        cube = np.arange(300, dtype=np.float32).reshape(3, 10, 10)
        hdu = fits.CompImageHDU(
            data=cube,
            name="SCI",
            compression_type="HCOMPRESS_1",
            quantize_level=16,
            tile_shape=(1, 5, 5),
        )
        hdu.writeto(self.temp("test.fits"))

        with fits.open(self.temp("test.fits")) as hdul:
            # HCOMPRESSed images are allowed to deviate from the original by
            # about 1/quantize_level of the RMS in each tile.
            assert np.abs(hdul["SCI"].data - cube).max() < 1.0 / 15.0

    def test_subtractive_dither_seed(self):
        """
        Regression test for https://github.com/spacetelescope/PyFITS/issues/32

        Ensure that when floating point data is compressed with the
        SUBTRACTIVE_DITHER_1 quantization method that the correct ZDITHER0 seed
        is added to the header, and that the data can be correctly
        decompressed.
        """

        array = np.arange(100.0).reshape(10, 10)
        csum = (array[0].view("uint8").sum() % 10000) + 1
        hdu = fits.CompImageHDU(
            data=array,
            quantize_method=SUBTRACTIVE_DITHER_1,
            dither_seed=DITHER_SEED_CHECKSUM,
        )
        hdu.writeto(self.temp("test.fits"))

        with fits.open(self.temp("test.fits")) as hdul:
            comp_header = hdul[1]._bintable.header
            assert "ZQUANTIZ" in comp_header
            assert comp_header["ZQUANTIZ"] == "SUBTRACTIVE_DITHER_1"
            assert "ZDITHER0" in comp_header
            assert comp_header["ZDITHER0"] == csum
            assert np.all(hdul[1].data == array)

    def test_disable_image_compression(self):
        with fits.open(self.data("comp.fits"), disable_image_compression=True) as hdul:
            # The compressed image HDU should show up as a BinTableHDU, but
            # *not* a CompImageHDU
            assert isinstance(hdul[1], fits.BinTableHDU)
            assert not isinstance(hdul[1], fits.CompImageHDU)

        with fits.open(self.data("comp.fits")) as hdul:
            assert isinstance(hdul[1], fits.CompImageHDU)

    def test_open_comp_image_in_update_mode(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/167

        Similar to test_open_scaled_in_update_mode(), but specifically for
        compressed images.
        """

        # Copy the original file before making any possible changes to it
        self.copy_file("comp.fits")
        mtime = os.stat(self.temp("comp.fits")).st_mtime

        time.sleep(1)

        fits.open(self.temp("comp.fits"), mode="update").close()

        # Ensure that no changes were made to the file merely by immediately
        # opening and closing it.
        assert mtime == os.stat(self.temp("comp.fits")).st_mtime

    @pytest.mark.slow
    def test_open_scaled_in_update_mode_compressed(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/88 2

        Identical to test_open_scaled_in_update_mode() but with a compressed
        version of the scaled image.
        """

        # Copy+compress the original file before making any possible changes to
        # it
        with fits.open(self.data("scale.fits"), do_not_scale_image_data=True) as hdul:
            chdu = fits.CompImageHDU(data=hdul[0].data, header=hdul[0].header)
            chdu.header["BZERO"] = hdul[0].header["BZERO"]
            chdu.header["BSCALE"] = hdul[0].header["BSCALE"]
            chdu.writeto(self.temp("scale.fits"))
        mtime = os.stat(self.temp("scale.fits")).st_mtime

        time.sleep(1)

        # Now open the file in update mode and close immediately. Note that we
        # need to set do_not_scale_image_data otherwise the data is scaled upon
        # being opened.
        fits.open(
            self.temp("scale.fits"),
            mode="update",
            do_not_scale_image_data=True,
        ).close()

        # Ensure that no changes were made to the file merely by immediately
        # opening and closing it.
        assert mtime == os.stat(self.temp("scale.fits")).st_mtime

        # Insert a slight delay to ensure the mtime does change when the file
        # is changed
        time.sleep(1)

        hdul = fits.open(self.temp("scale.fits"), "update")
        hdul[1].data
        hdul.close()

        # Now the file should be updated with the rescaled data
        assert mtime != os.stat(self.temp("scale.fits")).st_mtime
        hdul = fits.open(self.temp("scale.fits"), mode="update")
        assert hdul[1].data.dtype == np.dtype("float32")
        assert hdul[1].header["BITPIX"] == -32
        assert "BZERO" not in hdul[1].header
        assert "BSCALE" not in hdul[1].header

        # Try reshaping the data, then closing and reopening the file; let's
        # see if all the changes are preserved properly
        hdul[1].data.shape = (42, 10)
        hdul.close()

        hdul = fits.open(self.temp("scale.fits"))
        assert hdul[1].shape == (42, 10)
        assert hdul[1].data.dtype == np.dtype("float32")
        assert hdul[1].header["BITPIX"] == -32
        assert "BZERO" not in hdul[1].header
        assert "BSCALE" not in hdul[1].header
        hdul.close()

    def test_write_comp_hdu_direct_from_existing(self):
        with fits.open(self.data("comp.fits")) as hdul:
            hdul[1].writeto(self.temp("test.fits"))

        with fits.open(self.data("comp.fits")) as hdul1:
            with fits.open(self.temp("test.fits")) as hdul2:
                assert np.all(hdul1[1].data == hdul2[1].data)
                assert comparerecords(
                    hdul1[1].compressed_data, hdul2[1].compressed_data
                )

    def test_rewriting_large_scaled_image_compressed(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/88 1

        Identical to test_rewriting_large_scaled_image() but with a compressed
        image.
        """

        with fits.open(
            self.data("fixed-1890.fits"), do_not_scale_image_data=True
        ) as hdul:
            chdu = fits.CompImageHDU(data=hdul[0].data, header=hdul[0].header)
            chdu.writeto(self.temp("fixed-1890-z.fits"))

        hdul = fits.open(self.temp("fixed-1890-z.fits"))
        orig_data = hdul[1].data
        hdul.writeto(self.temp("test_new.fits"), overwrite=True)
        hdul.close()
        hdul = fits.open(self.temp("test_new.fits"))
        assert (hdul[1].data == orig_data).all()
        hdul.close()

        # Just as before, but this time don't touch hdul[0].data before writing
        # back out--this is the case that failed in
        # https://aeon.stsci.edu/ssb/trac/pyfits/ticket/84
        hdul = fits.open(self.temp("fixed-1890-z.fits"))
        hdul.writeto(self.temp("test_new.fits"), overwrite=True)
        hdul.close()
        hdul = fits.open(self.temp("test_new.fits"))
        assert (hdul[1].data == orig_data).all()
        hdul.close()

        # Test opening/closing/reopening a scaled file in update mode
        hdul = fits.open(self.temp("fixed-1890-z.fits"), do_not_scale_image_data=True)
        hdul.writeto(
            self.temp("test_new.fits"), overwrite=True, output_verify="silentfix"
        )
        hdul.close()
        hdul = fits.open(self.temp("test_new.fits"))
        orig_data = hdul[1].data
        hdul.close()
        hdul = fits.open(self.temp("test_new.fits"), mode="update")
        hdul.close()
        hdul = fits.open(self.temp("test_new.fits"))
        assert (hdul[1].data == orig_data).all()
        hdul.close()

    def test_scale_back_compressed(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/88 3

        Identical to test_scale_back() but uses a compressed image.
        """

        # Create a compressed version of the scaled image
        with fits.open(self.data("scale.fits"), do_not_scale_image_data=True) as hdul:
            chdu = fits.CompImageHDU(data=hdul[0].data, header=hdul[0].header)
            chdu.header["BZERO"] = hdul[0].header["BZERO"]
            chdu.header["BSCALE"] = hdul[0].header["BSCALE"]
            chdu.writeto(self.temp("scale.fits"))

        with fits.open(self.temp("scale.fits"), mode="update", scale_back=True) as hdul:
            orig_bitpix = hdul[1].header["BITPIX"]
            orig_bzero = hdul[1].header["BZERO"]
            orig_bscale = hdul[1].header["BSCALE"]
            orig_data = hdul[1].data.copy()
            hdul[1].data[0] = 0

        with fits.open(self.temp("scale.fits"), do_not_scale_image_data=True) as hdul:
            assert hdul[1].header["BITPIX"] == orig_bitpix
            assert hdul[1].header["BZERO"] == orig_bzero
            assert hdul[1].header["BSCALE"] == orig_bscale

            zero_point = math.floor(-orig_bzero / orig_bscale)
            assert (hdul[1].data[0] == zero_point).all()

        with fits.open(self.temp("scale.fits")) as hdul:
            assert_equal(hdul[1].data[1:], orig_data[1:])
            # Extra test to ensure that after everything the data is still the
            # same as in the original uncompressed version of the image
            with fits.open(self.data("scale.fits")) as hdul2:
                # Recall we made the same modification to the data in hdul
                # above
                hdul2[0].data[0] = 0
                assert_equal(hdul[1].data, hdul2[0].data)

    def test_lossless_gzip_compression(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/198"""

        rng = np.random.default_rng(42)
        noise = rng.normal(size=(20, 20))

        chdu1 = fits.CompImageHDU(data=noise, compression_type="GZIP_1")
        # First make a test image with lossy compression and make sure it
        # wasn't compressed perfectly.  This shouldn't happen ever, but just to
        # make sure the test non-trivial.
        chdu1.writeto(self.temp("test.fits"))

        with fits.open(self.temp("test.fits")) as h:
            assert np.abs(noise - h[1].data).max() > 0.0

        del h

        chdu2 = fits.CompImageHDU(
            data=noise, compression_type="GZIP_1", quantize_level=0.0
        )  # No quantization
        chdu2.writeto(self.temp("test.fits"), overwrite=True)

        with fits.open(self.temp("test.fits")) as h:
            assert (noise == h[1].data).all()

    def test_compression_column_tforms(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/199"""

        # Some interestingly tiled data so that some of it is quantized and
        # some of it ends up just getting gzip-compressed
        data2 = (np.arange(1, 8, dtype=np.float32) * 10)[:, np.newaxis] + np.arange(
            1, 7
        )
        np.random.seed(1337)
        data1 = np.random.uniform(size=(6 * 4, 7 * 4))
        data1[: data2.shape[0], : data2.shape[1]] = data2
        chdu = fits.CompImageHDU(data1, compression_type="RICE_1", tile_shape=(6, 7))
        chdu.writeto(self.temp("test.fits"))

        with fits.open(self.temp("test.fits"), disable_image_compression=True) as h:
            assert re.match(r"^1PB\(\d+\)$", h[1].header["TFORM1"])
            assert re.match(r"^1PB\(\d+\)$", h[1].header["TFORM2"])

    def test_compression_update_header(self):
        """Regression test for
        https://github.com/spacetelescope/PyFITS/issues/23
        """

        self.copy_file("comp.fits")
        with fits.open(self.temp("comp.fits"), mode="update") as hdul:
            assert isinstance(hdul[1], fits.CompImageHDU)
            hdul[1].header["test1"] = "test"
            hdul[1]._header["test2"] = "test2"

        with fits.open(self.temp("comp.fits")) as hdul:
            assert "test1" in hdul[1].header
            assert hdul[1].header["test1"] == "test"
            assert "test2" in hdul[1].header
            assert hdul[1].header["test2"] == "test2"

        # Test update via index now:
        with fits.open(self.temp("comp.fits"), mode="update") as hdul:
            hdr = hdul[1].header
            hdr[hdr.index("TEST1")] = "foo"

        with fits.open(self.temp("comp.fits")) as hdul:
            assert hdul[1].header["TEST1"] == "foo"

        # Test slice updates
        with fits.open(self.temp("comp.fits"), mode="update") as hdul:
            hdul[1].header["TEST*"] = "qux"

        with fits.open(self.temp("comp.fits")) as hdul:
            assert list(hdul[1].header["TEST*"].values()) == ["qux", "qux"]

        with fits.open(self.temp("comp.fits"), mode="update") as hdul:
            hdr = hdul[1].header
            idx = hdr.index("TEST1")
            hdr[idx : idx + 2] = "bar"

        with fits.open(self.temp("comp.fits")) as hdul:
            assert list(hdul[1].header["TEST*"].values()) == ["bar", "bar"]

        # Test updating a specific COMMENT card duplicate
        with fits.open(self.temp("comp.fits"), mode="update") as hdul:
            hdul[1].header[("COMMENT", 1)] = "I am fire. I am death!"

        with fits.open(self.temp("comp.fits")) as hdul:
            assert hdul[1].header["COMMENT"][1] == "I am fire. I am death!"
            assert hdul[1]._header["COMMENT"][1] == "I am fire. I am death!"

        # Test deleting by keyword and by slice
        with fits.open(self.temp("comp.fits"), mode="update") as hdul:
            hdr = hdul[1].header
            del hdr["COMMENT"]
            idx = hdr.index("TEST1")
            del hdr[idx : idx + 2]

        with fits.open(self.temp("comp.fits")) as hdul:
            assert "COMMENT" not in hdul[1].header
            assert "COMMENT" not in hdul[1]._header
            assert "TEST1" not in hdul[1].header
            assert "TEST1" not in hdul[1]._header
            assert "TEST2" not in hdul[1].header
            assert "TEST2" not in hdul[1]._header

    def test_compression_update_header_with_reserved(self, tmp_path):
        """
        Ensure that setting reserved keywords related to the table data
        structure on CompImageHDU image headers fails.
        """

        def test_set_keyword(hdr, keyword, value):
            with pytest.warns(UserWarning) as w:
                hdr[keyword] = value
            assert len(w) == 1
            assert str(w[0].message).startswith(f"Keyword {keyword!r} is reserved")
            assert keyword not in hdr

        with fits.open(self.data("comp.fits")) as hdul:
            hdr = hdul[1].header
            hdr["TFIELDS"] = 8
            hdr["THEAP"] = 1000
            hdr["TTYPE1"] = "Foo"
            hdr["ZCMPTYPE"] = "ASDF"
            hdr["ZVAL1"] = "Foo"
            with pytest.warns() as record:
                hdul.writeto(tmp_path / "test.fits")
            assert len(record) == 5
            for i, keyword in enumerate(
                ("TFIELDS", "THEAP", "TTYPE1", "ZCMPTYPE", "ZVAL1")
            ):
                assert f"Keyword {keyword!r} is reserved" in record[i].message.args[0]

    def test_compression_header_append(self, tmp_path):
        with fits.open(self.data("comp.fits")) as hdul:
            imghdr = hdul[1].header

            imghdr.append("TFIELDS")

            imghdr.append(("FOO", "bar", "qux"), end=True)
            assert "FOO" in imghdr
            assert imghdr[-1] == "bar"

            imghdr.append(("CHECKSUM", "abcd1234"))
            assert "CHECKSUM" in imghdr
            assert imghdr["CHECKSUM"] == "abcd1234"

            with pytest.warns(
                VerifyWarning, match="Keyword 'TFIELDS' is reserved"
            ) as w:
                hdul.writeto(tmp_path / "updated.fits")

        with fits.open(
            tmp_path / "updated.fits", disable_image_compression=True
        ) as hdulc:
            tblhdr = hdulc[1].header

            assert "FOO" in tblhdr
            assert tblhdr["FOO"] == "bar"

            assert "CHECKSUM" not in tblhdr
            assert "ZHECKSUM" in tblhdr
            assert tblhdr["ZHECKSUM"] == "abcd1234"

    def test_compression_header_append2(self):
        """
        Regression test for issue https://github.com/astropy/astropy/issues/5827
        """
        with fits.open(self.data("comp.fits")) as hdul:
            header = hdul[1].header
            while len(header) < 1000:
                header.append()  # pad with grow room

            # Append stats to header:
            header.append(("Q1_OSAVG", 1, "[adu] quadrant 1 overscan mean"))
            header.append(("Q1_OSSTD", 1, "[adu] quadrant 1 overscan stddev"))
            header.append(("Q1_OSMED", 1, "[adu] quadrant 1 overscan median"))

    def test_compression_header_insert(self, tmp_path):
        with fits.open(self.data("comp.fits")) as hdul:
            imghdr = hdul[1].header

            # First try inserting a restricted keyword
            imghdr.insert(1000, "TFIELDS")

            # First try keyword-relative insert
            imghdr.insert("TELESCOP", ("OBSERVER", "Phil Plait"))
            assert "OBSERVER" in imghdr
            assert imghdr.index("OBSERVER") == imghdr.index("TELESCOP") - 1

            # Next let's see if an index-relative insert winds up being
            # sensible
            idx = imghdr.index("OBSERVER")
            imghdr.insert("OBSERVER", ("FOO",))
            assert "FOO" in imghdr
            assert imghdr.index("FOO") == idx

            with pytest.warns(
                VerifyWarning, match="Keyword 'TFIELDS' is reserved"
            ) as w:
                hdul.writeto(tmp_path / "updated.fits")

        with fits.open(
            tmp_path / "updated.fits", disable_image_compression=True
        ) as hdulc:
            tblhdr = hdulc[1].header

            assert tblhdr.count("TFIELDS") == 1

            assert "OBSERVER" in tblhdr
            assert tblhdr.index("OBSERVER") == tblhdr.index("TELESCOP") - 1

            assert "FOO" in tblhdr
            assert tblhdr.index("FOO") == tblhdr.index("OBSERVER") - 1

    def test_compression_header_set_before_after(self, tmp_path):
        with fits.open(self.data("comp.fits")) as hdul:
            imghdr = hdul[1].header

            imghdr.set("ZBITPIX", 77, "asdf", after="GCOUNT")

            with pytest.warns(UserWarning, match="Keyword 'ZBITPIX' is reserved"):
                hdul.writeto(tmp_path / "updated1.fits")

        with fits.open(
            tmp_path / "updated1.fits", disable_image_compression=True
        ) as hdulc:
            tblhdr = hdulc[1].header

            assert tblhdr.count("ZBITPIX") == 1
            assert tblhdr["ZBITPIX"] != 77

        with fits.open(self.data("comp.fits")) as hdul:
            imghdr = hdul[1].header

            imghdr.set("FOO", 2, before="OBJECT")
            imghdr.set("BAR", 3, after="OBJECT")
            assert imghdr.index("FOO") == imghdr.index("OBJECT") - 1
            assert imghdr.index("BAR") == imghdr.index("OBJECT") + 1
            assert imghdr["FOO"] == 2
            assert imghdr["BAR"] == 3

            hdul.writeto(tmp_path / "updated2.fits")

        with fits.open(
            tmp_path / "updated2.fits", disable_image_compression=True
        ) as hdulc:
            tblhdr = hdulc[1].header

            assert tblhdr.index("FOO") == tblhdr.index("OBJECT") - 1
            assert tblhdr.index("BAR") == tblhdr.index("OBJECT") + 1
            assert tblhdr["FOO"] == 2
            assert tblhdr["BAR"] == 3

    def test_compression_header_append_commentary(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/2363
        """

        hdu = fits.CompImageHDU(np.array([0], dtype=np.int32))
        hdu.header["COMMENT"] = "hello world"
        assert hdu.header["COMMENT"] == ["hello world"]
        hdu.writeto(self.temp("test.fits"))

        with fits.open(self.temp("test.fits")) as hdul:
            assert hdul[1].header["COMMENT"] == ["hello world"]

    def test_compression_with_gzip_column(self):
        """
        Regression test for https://github.com/spacetelescope/PyFITS/issues/71
        """

        arr = np.zeros((2, 7000), dtype="float32")

        # The first row (which will be the first compressed tile) has a very
        # wide range of values that will be difficult to quantize, and should
        # result in use of a GZIP_COMPRESSED_DATA column
        arr[0] = np.linspace(0, 1, 7000)
        arr[1] = np.random.normal(size=7000)

        hdu = fits.CompImageHDU(data=arr)
        hdu.writeto(self.temp("test.fits"))

        with fits.open(self.temp("test.fits")) as hdul:
            comp_hdu = hdul[1]

            # GZIP-compressed tile should compare exactly
            assert np.all(comp_hdu.data[0] == arr[0])
            # The second tile uses lossy compression and may be somewhat off,
            # so we don't bother comparing it exactly

    def test_duplicate_compression_header_keywords(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/2750

        Tests that the fake header (for the compressed image) can still be read
        even if the real header contained a duplicate ZTENSION keyword (the
        issue applies to any keyword specific to the compression convention,
        however).
        """

        arr = np.arange(100, dtype=np.int32)
        hdu = fits.CompImageHDU(data=arr)
        hdu.writeto(self.temp("test1.fits"))

        # append the duplicate keyword
        with fits.open(
            self.temp("test1.fits"), disable_image_compression=True
        ) as hdulc:
            hdulc[1].header.append(("ZTENSION", "IMAGE"))
            header = hdulc[1].header
            hdulc.writeto(self.temp("test2.fits"))

        with fits.open(self.temp("test2.fits")) as hdul:
            assert header == hdul[1]._bintable.header
            # There's no good reason to have a duplicate keyword, but
            # technically it isn't invalid either :/
            assert hdul[1]._bintable.header.count("ZTENSION") == 2

    def test_scale_bzero_with_compressed_int_data(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/4600
        and https://github.com/astropy/astropy/issues/4588

        Identical to test_scale_bzero_with_int_data() but uses a compressed
        image.
        """

        a = np.arange(100, 200, dtype=np.int16)

        hdu1 = fits.CompImageHDU(data=a.copy())
        hdu2 = fits.CompImageHDU(data=a.copy())
        # Previously the following line would throw a TypeError,
        # now it should be identical to the integer bzero case
        hdu1.scale("int16", bzero=99.0)
        hdu2.scale("int16", bzero=99)
        assert np.allclose(hdu1.data, hdu2.data)

    def test_scale_back_compressed_uint_assignment(self):
        """
        Extend fix for #4600 to assignment to data

        Identical to test_scale_back_uint_assignment() but uses a compressed
        image.

        Suggested by:
        https://github.com/astropy/astropy/pull/4602#issuecomment-208713748
        """

        a = np.arange(100, 200, dtype=np.uint16)
        fits.CompImageHDU(a).writeto(self.temp("test.fits"))
        with fits.open(self.temp("test.fits"), mode="update", scale_back=True) as hdul:
            hdul[1].data[:] = 0
            assert np.allclose(hdul[1].data, 0)

    def test_compressed_header_double_extname(self):
        """Test that a double EXTNAME with one default value does not
        mask the non-default value."""
        with fits.open(self.data("double_ext.fits")) as hdul:
            hdu = hdul[1]

            # The non-default name should be returned.
            assert hdu.name == "ccd00"
            assert "EXTNAME" in hdu.header
            assert hdu.name == hdu.header["EXTNAME"]

            # There should be 1 non-default EXTNAME entries.
            indices = hdu.header._keyword_indices["EXTNAME"]
            assert len(indices) == 1

            # Test header sync from property set.
            new_name = "NEW_NAME"
            hdu.name = new_name
            assert hdu.name == new_name
            assert hdu.header["EXTNAME"] == new_name

            # Check that setting the header will change the name property.
            hdu.header["EXTNAME"] = "NEW2"
            assert hdu.name == "NEW2"

            hdul.writeto(self.temp("tmp.fits"), overwrite=True)
            with fits.open(self.temp("tmp.fits")) as hdul1:
                hdu1 = hdul1[1]
                assert len(hdu.header._keyword_indices["EXTNAME"]) == 1
                assert hdu1.name == "NEW2"

            # Check that deleting EXTNAME will and setting the name will
            # work properly.
            del hdu.header["EXTNAME"]
            hdu.name = "RE-ADDED"
            assert hdu.name == "RE-ADDED"

            with pytest.raises(TypeError):
                hdu.name = 42

    def test_compressed_header_extname(self):
        """Test consistent EXTNAME / hdu name interaction."""
        name = "FOO"
        hdu = fits.CompImageHDU(data=np.arange(10), name=name)
        assert hdu._header["EXTNAME"] == name
        assert hdu.header["EXTNAME"] == name
        assert hdu.name == name

        name = "BAR"
        hdu.name = name
        assert hdu._header["EXTNAME"] == name
        assert hdu.header["EXTNAME"] == name
        assert hdu.name == name

        assert len(hdu._header._keyword_indices["EXTNAME"]) == 1

    def test_compressed_header_minimal(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/11694

        Tests that CompImageHDU can be initialized with a Header that
        contains few or no cards, and doesn't require specific cards
        such as 'BITPIX' or 'NAXIS'.
        """
        fits.CompImageHDU(data=np.arange(10), header=fits.Header())
        header = fits.Header({"HELLO": "world"})
        hdu = fits.CompImageHDU(data=np.arange(10), header=header)
        assert hdu.header["HELLO"] == "world"

    @pytest.mark.parametrize(
        ("keyword", "dtype", "expected"),
        [
            ("BSCALE", np.uint8, np.float32),
            ("BSCALE", np.int16, np.float32),
            ("BSCALE", np.int32, np.float64),
            ("BZERO", np.uint8, np.float32),
            ("BZERO", np.int16, np.float32),
            ("BZERO", np.int32, np.float64),
        ],
    )
    def test_compressed_scaled_float(self, keyword, dtype, expected):
        """
        If BSCALE,BZERO is set to floating point values, the image
        should be floating-point.

        https://github.com/astropy/astropy/pull/6492

        Parameters
        ----------
        keyword : `str`
            Keyword to set to a floating-point value to trigger
            floating-point pixels.
        dtype : `numpy.dtype`
            Type of original array.
        expected : `numpy.dtype`
            Expected type of uncompressed array.
        """
        value = 1.23345  # A floating-point value
        hdu = fits.CompImageHDU(np.arange(0, 10, dtype=dtype))
        hdu.header[keyword] = value
        hdu.writeto(self.temp("test.fits"))
        del hdu
        with fits.open(self.temp("test.fits")) as hdu:
            assert hdu[1].header[keyword] == value
            assert hdu[1].data.dtype == expected

    @pytest.mark.parametrize(
        "dtype", (np.uint8, np.int16, np.uint16, np.int32, np.uint32)
    )
    def test_compressed_integers(self, dtype):
        """Test that the various integer dtypes are correctly written and read.

        Regression test for https://github.com/astropy/astropy/issues/9072

        """
        mid = np.iinfo(dtype).max // 2
        data = np.arange(mid - 50, mid + 50, dtype=dtype)
        testfile = self.temp("test.fits")
        hdu = fits.CompImageHDU(data=data)
        hdu.writeto(testfile, overwrite=True)
        new = fits.getdata(testfile)
        np.testing.assert_array_equal(data, new)

    @pytest.mark.parametrize(
        ("dtype", "compression_type"), product(("f", "i4"), COMPRESSION_TYPES)
    )
    def test_write_non_contiguous_data(self, dtype, compression_type):
        """
        Regression test for https://github.com/astropy/astropy/issues/2150

        This used to require changing the whole array to be C-contiguous before
        passing to CFITSIO, but we no longer need this - our explicit conversion
        to bytes in the compression codecs returns contiguous bytes for each
        tile on-the-fly.
        """

        orig = np.arange(400, dtype=dtype).reshape((20, 20), order="f")[::2, ::2]
        assert not orig.flags.contiguous
        primary = fits.PrimaryHDU()
        hdu = fits.CompImageHDU(orig, compression_type=compression_type)
        hdulist = fits.HDUList([primary, hdu])
        hdulist.writeto(self.temp("test.fits"))

        actual = fits.getdata(self.temp("test.fits"))
        assert_equal(orig, actual)

    def test_slice_and_write_comp_hdu(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/9955
        """
        with fits.open(self.data("comp.fits")) as hdul:
            hdul[1].data = hdul[1].data[:200, :100]
            assert not hdul[1].data.flags.contiguous
            hdul[1].writeto(self.temp("test.fits"))

        with fits.open(self.data("comp.fits")) as hdul1:
            with fits.open(self.temp("test.fits")) as hdul2:
                assert_equal(hdul1[1].data[:200, :100], hdul2[1].data)

    def test_comp_image_properties_default(self):
        chdu = fits.CompImageHDU(np.zeros((3, 4, 5)))
        assert chdu.tile_shape == (1, 1, 5)
        assert chdu.compression_type == "RICE_1"

    def test_comp_image_properties_set(self):
        chdu = fits.CompImageHDU(
            np.zeros((3, 4, 5)), compression_type="PLIO_1", tile_shape=(2, 3, 4)
        )
        assert chdu.tile_shape == (2, 3, 4)
        assert chdu.compression_type == "PLIO_1"

    def test_compressed_optional_prefix_tform(self, tmp_path):
        # Regression test for a bug that caused an error if a
        # compressed file had TFORM missing the optional 1 prefix

        data = np.zeros((3, 4, 5))

        hdu1 = fits.CompImageHDU(data=data)
        hdu1.writeto(tmp_path / "compressed.fits")

        with fits.open(
            tmp_path / "compressed.fits", disable_image_compression=True, mode="update"
        ) as hdul:
            assert hdul[1].header["TFORM1"] == "1PB(0)"
            assert hdul[1].header["TFORM2"] == "1PB(24)"
            hdul[1].header["TFORM1"] = "PB(0)"
            hdul[1].header["TFORM2"] = "PB(24)"

        with fits.open(
            tmp_path / "compressed.fits", disable_image_compression=True
        ) as hdul:
            assert hdul[1].header["TFORM1"] == "PB(0)"
            assert hdul[1].header["TFORM2"] == "PB(24)"

        with fits.open(tmp_path / "compressed.fits") as hdul:
            assert_equal(hdul[1].data, data)

    def test_info(self):
        """
        Make sure .info() works correctly when CompImageHDUs are present
        """
        output = io.StringIO()
        with fits.open(self.data("comp.fits")) as hdul:
            hdul.info(output=output)
        output.seek(0)

        # Note: ignore the first line which just gives the filename
        actual = output.read().splitlines()[1:]

        expected = [
            "No.    Name      Ver    Type      Cards   Dimensions   Format",
            "0  PRIMARY       1 PrimaryHDU       4   ()",
            "1  COMPRESSED_IMAGE    1 CompImageHDU    105   (440, 300)   int16",
        ]

        for line_actual, line_expected in zip(actual, expected, strict=True):
            assert line_actual.strip() == line_expected.strip()

    def test_shape(self):
        with fits.open(self.data("comp.fits")) as hdul:
            assert hdul[1].header["NAXIS1"] == 440
            assert hdul[1].header["NAXIS2"] == 300
            assert hdul[1].shape == (300, 440)
            hdul[1].data = np.ones((120, 150))
            assert hdul[1].header["NAXIS1"] == 150
            assert hdul[1].header["NAXIS2"] == 120
            assert hdul[1].shape == (120, 150)

    def test_inplace_data_modify(self, tmp_path):
        self.copy_file("comp.fits")

        with fits.open(self.temp("comp.fits"), mode="update") as hdul:
            data = hdul[1].data
            data[0] = 0

        with fits.open(self.temp("comp.fits")) as hdul:
            assert hdul[1].data[0, 0] == 0

    def test_summary_noload(self):
        # Make sure that calling info() (and therefore CompImageHDU.summary)
        # does not cause the data to be loaded.
        with fits.open(self.data("comp.fits")) as hdul:
            summary = hdul.info(output=False)
            assert summary == [
                (0, "PRIMARY", 1, "PrimaryHDU", 4, (), "", ""),
                (
                    1,
                    "COMPRESSED_IMAGE",
                    1,
                    "CompImageHDU",
                    105,
                    (440, 300),
                    "int16",
                    "",
                ),
            ]
            assert not hdul[1]._data_loaded

    def test_fileinfo(self):
        with fits.open(self.data("comp.fits")) as hdul:
            res = hdul.fileinfo(1)

        assert res["datLoc"] == 14400
        assert res["datSpan"] == 72000
        assert res["filemode"] == "readonly"
        assert res["filename"] == self.data("comp.fits")
        assert res["hdrLoc"] == 2880
        assert not res["resized"]


class TestCompHDUSections:
    @pytest.fixture(autouse=True)
    def setup_method(self, tmp_path):
        shape = (13, 17, 25)
        self.data = np.arange(np.prod(shape)).reshape(shape).astype(np.int32)

        header1 = fits.Header()
        hdu1 = fits.CompImageHDU(
            self.data, header1, compression_type="RICE_1", tile_shape=(5, 4, 5)
        )

        header2 = fits.Header()
        hdu2 = fits.CompImageHDU(
            self.data, header2, compression_type="RICE_1", tile_shape=(5, 4, 5)
        )
        hdu2.header["BSCALE"] = 2
        hdu2.header["BZERO"] = 100
        hdulist = fits.HDUList([fits.PrimaryHDU(), hdu1, hdu2])
        hdulist.writeto(tmp_path / "sections.fits")

        self.hdul = fits.open(tmp_path / "sections.fits")
        self.hdul2 = fits.open(tmp_path / "sections.fits")

    def teardown_method(self):
        self.hdul.close()
        self.hdul = None

    @given(basic_indices((13, 17, 25)))
    def test_section_slicing(self, index):
        assert_equal(self.hdul[1].section[index], self.hdul[1].data[index])
        assert_equal(self.hdul[1].section[index], self.data[index])

    @given(basic_indices((13, 17, 25)))
    def test_section_slicing_scaling(self, index):
        assert_equal(self.hdul[2].section[index], self.hdul[2].data[index])
        assert_equal(self.hdul[2].section[index], self.data[index] * 2 + 100)

    def test_section_properties(self):
        assert self.hdul[1].section.dtype is np.dtype("int32")
        assert self.hdul[1].section.ndim == 3


def test_comphdu_fileobj():
    # Regression test for a bug that caused an error to happen
    # internally when reading the data if requested data shapes
    # were not plain integers - this was triggered when accessing
    # sections on data backed by certain kinds of objects such as
    # BytesIO (but not regular file handles)

    data = np.arange(6).reshape((2, 3)).astype(np.int32)

    byte_buffer = BytesIO()

    header = fits.Header()
    hdu = fits.CompImageHDU(data, header, compression_type="RICE_1")
    hdu.writeto(byte_buffer)

    byte_buffer.seek(0)

    hdu2 = fits.open(byte_buffer, mode="readonly")[1]
    assert hdu2.section[1, 2] == 5


def test_comphdu_bscale(tmp_path):
    """
    Regression test for a bug that caused extensions that used BZERO and BSCALE
    that got turned into CompImageHDU to end up with BZERO/BSCALE before the
    TFIELDS.
    """

    filename1 = tmp_path / "3hdus.fits"
    filename2 = tmp_path / "3hdus_comp.fits"

    x = np.random.random((100, 100)) * 100

    x0 = fits.PrimaryHDU()
    x1 = fits.ImageHDU(np.array(x - 50, dtype=int), uint=True)
    x1.header["BZERO"] = 20331
    x1.header["BSCALE"] = 2.3
    hdus = fits.HDUList([x0, x1])
    hdus.writeto(filename1)

    # fitsverify (based on cfitsio) should fail on this file, only seeing the
    # first HDU.
    with fits.open(filename1) as hdus:
        hdus[1] = fits.CompImageHDU(
            data=hdus[1].data.astype(np.uint32), header=hdus[1].header
        )
        hdus.writeto(filename2)

    # open again and verify
    with fits.open(filename2) as hdus:
        hdus[1].verify("exception")


def test_image_write_readonly(tmp_path):
    # Regression test to make sure that we can write out read-only arrays (#5512)

    x = np.array([1.0, 2.0, 3.0])
    x.setflags(write=False)
    ghdu = fits.CompImageHDU(data=x)
    # add_datasum does not work for CompImageHDU
    # ghdu.add_datasum()

    filename = tmp_path / "test2.fits"

    ghdu.writeto(filename)

    with fits.open(filename) as hdulist:
        assert_equal(hdulist[1].data, [1.0, 2.0, 3.0])


def test_uint_option(tmp_path):
    """
    Check that the uint= option works correctly
    """

    filename = tmp_path / "uint_test.fits"

    data = (2 ** (1 + np.arange(16).reshape((4, 4))) - 1).astype(np.uint16)

    hdu = fits.CompImageHDU(data)
    hdu.writeto(filename)

    with fits.open(filename) as hdulist:
        assert hdulist[1].data.dtype == np.dtype("uint16")
        assert_equal(hdulist[1].data, data)

    with fits.open(filename, uint=False) as hdulist:
        assert hdulist[1].data.dtype == np.dtype("float32")
        assert_allclose(hdulist[1].data, data)


def test_incorrect_bzero(tmp_path):
    """
    Regression test for https://github.com/astropy/astropy/issues/5999 which is
    a bug that caused BZERO to be incorrectly set to a value in a compressed
    FITS HDU if a header was passed in with a value even though the data was
    specified as a floating-point value.
    """

    data = np.arange(0, 100, dtype=np.uint16)
    hdu = fits.ImageHDU(data=data)
    data = hdu.data.astype(np.float64)
    header = hdu.header

    hdulist = fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(data=data, header=header),
            fits.CompImageHDU(data=data, header=header),
        ]
    )

    hdulist.writeto(tmp_path / "test_bzero.fits")

    for hdu in hdulist:
        assert hdu.header.get("BZERO") is None
        assert hdu.header.get("BSCALE") is None

    with fits.open(tmp_path / "test_bzero.fits") as hdulist_read:
        for hdu in hdulist_read:
            assert hdu.header.get("BZERO") is None
            assert hdu.header.get("BSCALE") is None
        np.testing.assert_array_equal(hdulist_read[1].data, hdulist_read[2].data)


def test_custom_extname():
    """
    Regression test for a bug that caused specifying a custom name for a
    CompImageHDU to not work if an existing image header was used.
    """
    data = np.arange(0, 100, dtype=np.uint16)
    hdu = fits.ImageHDU(data=data)
    header = hdu.header
    fits.CompImageHDU(header=header, name="compressed")


def test_pickle():
    """
    Regression test for https://github.com/astropy/astropy/issues/10512 which
    was a bug that caused CompImageHDU to not be picklable.
    """
    data = np.array([1, 2, 3])
    hdu = fits.CompImageHDU(data)
    a = pickle.loads(pickle.dumps(hdu))
    assert_equal(hdu.data, data)


def test_compression_settings_delayed_data(tmp_path):
    """
    Regression test for https://github.com/astropy/astropy/issues/12216 which
    was a bug that caused compression settings to be ignored if the data was
    only set later.
    """

    hdu1 = fits.CompImageHDU(quantize_level=-32)
    hdu2 = fits.CompImageHDU(quantize_level=-32, data=np.array([0.0]))

    with NumpyRNGContext(42):
        data = np.random.random((16, 16))

    hdu1.data = data
    hdu2.data = data

    hdu1.writeto(tmp_path / "data1.fits")
    hdu2.writeto(tmp_path / "data2.fits")

    with fits.open(tmp_path / "data1.fits") as hdulist1_read:
        with fits.open(tmp_path / "data2.fits") as hdulist2_read:
            assert_equal(hdulist1_read[1].data, hdulist2_read[1].data)


def test_header_assignment_issue(tmp_path):
    """
    Regression test for https://github.com/astropy/astropy/issues/14081 which
    was a bug that caused entries in header to not be preserved under certain
    conditions if copied from another header.
    """

    ih = fits.ImageHDU()
    ih.header["test"] = "right"

    ch = fits.CompImageHDU()
    ch.header["test"] = "wrong"
    ch.header = ih.header
    assert ch.header["test"] == "right"

    ch.writeto(tmp_path / "test_header.fits")
    with fits.open(tmp_path / "test_header.fits") as chdl:
        assert chdl[1].header["test"] == "right"


def test_section_unwritten():
    """
    Regression test for https://github.com/astropy/astropy/issues/14611 which
    was a bug that caused CompImageHDU.section to not work correctly if the
    file was not written out to disk first.
    """

    data = np.arange(21 * 33).reshape((21, 33)).astype(np.int32)
    header = fits.Header()
    hdu = fits.CompImageHDU(data, header, compression_type="RICE_1", tile_shape=(5, 6))
    assert_equal(hdu.section[...], data)
    assert hdu.section[3, 4] == data[3, 4]


EXPECTED_HEADER = """
XTENSION= 'IMAGE   '           / Image extension
BITPIX  =                   16 / data type of original image
NAXIS   =                    2 / dimension of original image
NAXIS1  =                   10 / length of original image axis
NAXIS2  =                   10 / length of original image axis
PCOUNT  =                    0 / number of parameters
GCOUNT  =                    1 / number of groups
END
""".lstrip()


def test_header():
    """
    Check that the header is correct when reading in compressed images and
    correctly shows the image dimensions.
    """

    filename = os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "..",
        "tests",
        "data",
        "compressed_image.fits",
    )

    with fits.open(filename) as hdulist:
        assert hdulist[1].header == fits.Header.fromstring(EXPECTED_HEADER, sep="\n")


def test_rice_one_alias():
    # Regression test for a bug that caused RICE_ONE (which we document as an
    # acceptable alias) to no longer be recognized.
    chdu = fits.CompImageHDU(np.zeros((3, 4, 5)))
    chdu.compression_type = "RICE_ONE"


def test_header_order(tmp_path):
    # Make sure that comments and history entries appear in the same order and
    # same location after compression and decompression.

    data = np.random.random((128, 128))
    header = fits.Header()

    header["a"] = 1
    header["b"] = 2
    header["c"] = 3
    header["d"] = 4

    header.add_history("history one", before="a")
    header.add_comment("comment one", before="a")
    header.add_comment("comment two", after="a")
    header.add_comment("comment three", before="b")
    header.add_history("history two", before="b")
    header.add_comment("comment four", after="b")
    header.add_comment("comment five", after="d")
    header.add_history("history three")
    header.add_blank()
    header.add_blank()
    header.add_blank()

    hdulist = fits.HDUList([fits.PrimaryHDU(), fits.CompImageHDU(data, header)])

    hdulist.writeto(tmp_path / "test.fits", overwrite=True)

    with fits.open(tmp_path / "test.fits") as hdulist2:
        assert hdulist[1].header.tostring("\n") == hdulist2[1].header.tostring("\n")


def test_hdu_lazy_loading(tmp_path):
    # Lazy loading of HDUs relies on parameters such as _data_offset and so on,
    # so we need to make sure these correctly point to the internal BinTableHDU.

    hdulist = fits.HDUList([fits.PrimaryHDU()])

    for idx in range(6):
        data = np.ones((2**idx, 2**idx)) * idx
        hdulist.append(fits.CompImageHDU(data))

    hdulist.writeto(tmp_path / "multi_hdu.fits")

    with open(tmp_path / "multi_hdu.fits", "rb") as fileobj:
        with fits.open(fileobj) as hdulist:
            for idx in range(6):
                assert_equal(hdulist[idx + 1].data, idx)
                assert hdulist[idx + 1].data.shape == (2**idx, 2**idx)
                fileobj.seek(0)


@pytest.mark.parametrize(
    "init_kwargs, expected_zflags",
    [
        pytest.param(
            dict(quantize_level=-32, data=np.array([0.0])),
            {
                "ZNAME1": "BLOCKSIZE",
                "ZVAL1": 32,
                "ZNAME2": "BYTEPIX",
                "ZVAL2": 4,
                "ZNAME3": "NOISEBIT",
                "ZVAL3": -32,
                "ZCMPTYPE": "RICE_1",
                "ZQUANTIZ": "NO_DITHER",
            },
            id="quantize_level_w_data",
        ),
        pytest.param(
            dict(quantize_level=-32),
            {
                "ZNAME1": "BLOCKSIZE",
                "ZVAL1": 32,
                "ZNAME2": "BYTEPIX",
                "ZVAL2": 4,
                "ZNAME3": "NOISEBIT",
                "ZVAL3": -32,
                "ZCMPTYPE": "RICE_1",
                "ZQUANTIZ": "NO_DITHER",
            },
            id="quantize_level_wo_data",
        ),
        pytest.param(
            dict(quantize_method=2, dither_seed=1, data=np.array([0.0])),
            {
                "ZNAME1": "BLOCKSIZE",
                "ZVAL1": 32,
                "ZNAME2": "BYTEPIX",
                "ZVAL2": 4,
                "ZNAME3": "NOISEBIT",
                "ZVAL3": 16,
                "ZCMPTYPE": "RICE_1",
                "ZQUANTIZ": "SUBTRACTIVE_DITHER_2",
            },
            id="quantize_level_dither_seed_w_data",
        ),
        pytest.param(
            dict(quantize_method=2, dither_seed=1),
            {
                "ZNAME1": "BLOCKSIZE",
                "ZVAL1": 32,
                "ZNAME2": "BYTEPIX",
                "ZVAL2": 4,
                "ZNAME3": "NOISEBIT",
                "ZVAL3": 16,
                "ZCMPTYPE": "RICE_1",
                "ZQUANTIZ": "SUBTRACTIVE_DITHER_2",
            },
            id="quantize_level_dither_seed_wo_data",
        ),
    ],
)
def test_compression_options_with_mutated_data(
    init_kwargs,
    expected_zflags,
    tmp_path,
):
    # see https://github.com/astropy/astropy/issues/12216
    prng = np.random.default_rng(seed=0)
    data = prng.random((50, 30)).astype("float32")

    hdu = fits.CompImageHDU(**init_kwargs)
    hdu.data = data

    hdu.writeto(tmp_path / "t.fits")

    hdr = fits.getheader(tmp_path / "t.fits", ext=1, disable_image_compression=True)

    zflags = {k: hdr[k] for k in expected_zflags}
    assert zflags == expected_zflags


def test_reserved_keywords_stripped(tmp_path):
    # Regression test for a bug that caused THEAP, ZBLANK, ZSCALE and ZZERO to
    # not be correctly stripped from the compressed header when decompressing
    #
    # See also https://github.com/astropy/astropy/issues/18067

    data = np.arange(6).reshape((2, 3))

    hdu = fits.CompImageHDU(data)
    hdu.writeto(tmp_path / "compressed.fits")

    with fits.open(
        tmp_path / "compressed.fits", disable_image_compression=True
    ) as hduc:
        hduc[1].header["THEAP"] = hduc[1].header["NAXIS1"] * hduc[1].header["NAXIS2"]
        hduc[1].header["ZBLANK"] = 1231212
        hduc[1].header["ZSCALE"] = 2
        hduc[1].header["ZZERO"] = 10
        hduc[1].writeto(tmp_path / "compressed_with_extra.fits")

    with fits.open(tmp_path / "compressed_with_extra.fits") as hdud:
        assert "THEAP" not in hdud[1].header
        assert "ZBLANK" not in hdud[1].header
        assert "ZSCALE" not in hdud[1].header
        assert "ZZERO" not in hdud[1].header
