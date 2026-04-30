# Licensed under a 3-clause BSD style license - see PYFITS.rst

import io
import os
import pathlib
import shutil
import subprocess
import warnings

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from astropy import units as u
from astropy.io import fits
from astropy.io.fits import printdiff
from astropy.io.fits.connect import REMOVE_KEYWORDS
from astropy.io.fits.tests.test_table import _assert_attr_col
from astropy.table import Table
from astropy.utils.exceptions import AstropyUserWarning

from .conftest import FitsTestCase

HAS_FUNPACK = shutil.which("funpack") is not None
HAS_FPACK = shutil.which("fpack") is not None


class TestConvenience(FitsTestCase):
    def test_resource_warning(self):
        warnings.simplefilter("always", ResourceWarning)
        _ = fits.getdata(self.data("test0.fits"))
        _ = fits.getheader(self.data("test0.fits"))

    def test_fileobj_not_closed(self):
        """
        Tests that file-like objects are not closed after being passed
        to convenience functions.

        Regression test for https://github.com/astropy/astropy/issues/5063
        """

        f = open(self.data("test0.fits"), "rb")
        _ = fits.getdata(f)
        assert not f.closed

        f.seek(0)
        _ = fits.getheader(f)
        assert not f.closed

        f.close()  # Close it now

    def test_table_to_hdu(self):
        table = Table(
            [[1, 2, 3], ["a", "b", "c"], [2.3, 4.5, 6.7]],
            names=["a", "b", "c"],
            dtype=["i", "U1", "f"],
        )
        table["a"].unit = "m/s"
        table["b"].unit = "not-a-unit"

        with pytest.warns(
            u.UnitsWarning, match="'not-a-unit' did not parse as fits unit"
        ) as w:
            hdu = fits.table_to_hdu(table, name="MYTABLE")

        assert len(w) == 1
        assert hdu.header["EXTNAME"] == "MYTABLE"

        # Check that TUNITn cards appear in the correct order
        # (https://github.com/astropy/astropy/pull/5720)
        assert hdu.header.index("TUNIT1") < hdu.header.index("TTYPE2")

        assert isinstance(hdu, fits.BinTableHDU)
        filename = self.temp("test_table_to_hdu.fits")
        hdu.writeto(filename)

        assert fits.getval(filename, "EXTNAME", ext=1) == "MYTABLE"

    def test_masked_table_to_hdu(self):
        i = np.ma.MaskedArray([1, 2, 3], mask=[True, False, False])
        s = np.ma.MaskedArray(["a", "b", "c"], mask=[False, True, True])
        c = np.ma.MaskedArray([2.3 + 1j, 4.5 + 0j, 6.7 - 1j], mask=[True, False, True])
        f = np.ma.MaskedArray([2.3, 4.5, 6.7], mask=[True, False, True])
        table = Table([i, s, c, f], names=["i", "s", "c", "f"])
        # Check that FITS standard is used in replacing masked values.
        hdu = fits.table_to_hdu(table)
        assert isinstance(hdu, fits.BinTableHDU)
        assert hdu.header["TNULL1"] == i.fill_value
        assert_array_equal(hdu.data["i"], i.filled())
        assert_array_equal(hdu.data["s"], s.filled(""))
        assert_array_equal(hdu.data["c"], c.filled(np.nan))
        assert_array_equal(hdu.data["c"].real, c.real.filled(np.nan))
        assert_array_equal(hdu.data["c"].imag, c.imag.filled(np.nan))
        assert_array_equal(hdu.data["c"], c.filled(complex(np.nan, np.nan)))
        assert_array_equal(hdu.data["f"], f.filled(np.nan))
        filename = self.temp("test_table_to_hdu.fits")
        hdu.writeto(filename, overwrite=True)

    def test_masked_integer_arrays(self):
        # Regression test for #18817
        testfile = self.temp("test_masked_integer_arrays.fits")
        t_w = Table(
            rows=[
                [[np.ma.masked, np.ma.masked]],
                [[1, 2]],
                [[1, np.ma.masked]],
            ],
            names=["a"],
        )
        t_w.write(testfile, overwrite=True)
        t_r = Table.read(testfile)
        assert repr(t_w) == repr(t_r)

    def test_table_non_stringifyable_unit_to_hdu(self):
        table = Table(
            [[1, 2, 3], ["a", "b", "c"], [2.3, 4.5, 6.7]],
            names=["a", "b", "c"],
            dtype=["i", "U1", "f"],
        )
        table["a"].unit = u.core.IrreducibleUnit("test")

        with pytest.warns(
            AstropyUserWarning, match="The unit 'test' could not be saved"
        ) as w:
            fits.table_to_hdu(table)
        assert len(w) == 1

    def test_table_to_hdu_convert_comment_convention(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/6079
        """
        table = Table(
            [[1, 2, 3], ["a", "b", "c"], [2.3, 4.5, 6.7]],
            names=["a", "b", "c"],
            dtype=["i", "U1", "f"],
        )
        table.meta["comments"] = ["This", "is", "a", "comment"]
        hdu = fits.table_to_hdu(table)

        assert hdu.header.get("comment") == ["This", "is", "a", "comment"]
        with pytest.raises(ValueError):
            hdu.header.index("comments")

    def test_table_to_hdu_filter_reserved(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/9387
        """
        diag = "be ignored since it conflicts with a FITS reserved keyword"
        ins_cards = {
            "EXPTIME": 32.1,
            "XTENSION": "NEWTABLE",
            "NAXIS": 1,
            "NAXIS1": 3,
            "NAXIS2": 9,
            "PCOUNT": 42,
            "OBSERVER": "Adams",
        }
        table = Table(
            [[1, 2, 3], ["a", "b", "c"], [2.3, 4.5, 6.7]],
            names=["a", "b", "c"],
            dtype=["i4", "U1", "f8"],
        )
        table.meta.update(ins_cards)

        with pytest.warns(
            AstropyUserWarning, match=rf"Meta-data keyword \w+ will {diag}"
        ) as w:
            hdu = fits.table_to_hdu(table)

        # This relies on the warnings being raised in the order of the
        # meta dict (note that the first and last card are legitimate keys)
        assert len(w) == len(ins_cards) - 2
        for i, key in enumerate(list(ins_cards)[1:-1]):
            assert f"Meta-data keyword {key}" in str(w[i].message)

        assert hdu.header.get("XTENSION") == "BINTABLE"
        assert hdu.header.get("NAXIS") == 2
        assert hdu.header.get("NAXIS1") == 13
        assert hdu.header.get("NAXIS2") == 3
        assert hdu.header.get("PCOUNT") == 0
        np.testing.assert_almost_equal(hdu.header.get("EXPTIME"), 3.21e1)

    @pytest.mark.parametrize("card", REMOVE_KEYWORDS)
    def test_table_to_hdu_warn_reserved(self, card):
        """
        Test warning for each keyword in ..connect.REMOVE_KEYWORDS, 1 by 1
        """
        diag = "be ignored since it conflicts with a FITS reserved keyword"
        res_cards = {
            "XTENSION": "BINTABLE",
            "BITPIX": 8,
            "NAXIS": 2,
            "NAXIS1": 12,
            "NAXIS2": 3,
            "PCOUNT": 0,
            "GCOUNT": 1,
            "TFIELDS": 2,
            "THEAP": None,
        }
        ins_cards = {
            "XTENSION": "TABLE",
            "BITPIX": 16,
            "NAXIS": 1,
            "NAXIS1": 2,
            "NAXIS2": 6,
            "PCOUNT": 2,
            "GCOUNT": 2,
            "TFIELDS": 4,
            "THEAP": 36,
        }

        table = Table(
            [[1.0, 2.0, 3.0], [2.3, 4.5, 6.7]],
            names=["wavelength", "flux"],
            dtype=["f8", "f4"],
        )
        table.meta["ORIGIN"] = "Min.Silly Walks"
        table.meta[card] = ins_cards[card]
        assert table.meta.get(card) != res_cards[card]

        with pytest.warns(
            AstropyUserWarning, match=f"Meta-data keyword {card} will {diag}"
        ):
            hdu = fits.table_to_hdu(table)

        assert hdu.header.get(card) == res_cards[card]
        assert hdu.header.get("ORIGIN") == "Min.Silly Walks"

    def test_table_to_hdu_filter_incompatible(self):
        """
        Test removal of unsupported data types from header
        """
        table = Table(
            [[1, 2, 3], ["a", "b", "c"], [2.3, 4.5, 6.7]],
            names=["a", "b", "c"],
            dtype=["i4", "U1", "f8"],
        )
        table.meta.update(
            {
                "OBSDATE": "2001-05-26",
                "RAMP": np.arange(5),
                "TARGETS": {"PRIMARY": 1, "SECONDAR": 3},
            }
        )
        with pytest.warns(
            AstropyUserWarning,
            match=r"Attribute \S+ of type "
            r".+ cannot be added to FITS Header - skipping",
        ):
            hdu = fits.table_to_hdu(table)

        assert hdu.header.get("OBSDATE") == "2001-05-26"
        assert "RAMP" not in hdu.header
        assert "TARGETS" not in hdu.header

    def test_table_writeto_header(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/5988
        """
        data = np.zeros((5,), dtype=[("x", float), ("y", int)])
        h_in = fits.Header()
        h_in["ANSWER"] = (42.0, "LTU&E")
        filename = self.temp("tabhdr42.fits")
        fits.writeto(filename, data=data, header=h_in, overwrite=True)
        h_out = fits.getheader(filename, ext=1)
        assert h_out["ANSWER"] == 42

    def test_image_extension_update_header(self, home_is_temp):
        """
        Test that _makehdu correctly includes the header. For example in the
        fits.update convenience function.
        """
        filename = self.temp("twoextension.fits")

        hdus = [fits.PrimaryHDU(np.zeros((10, 10))), fits.ImageHDU(np.zeros((10, 10)))]

        # Try to update a non-existent file
        with pytest.raises(FileNotFoundError, match="No such file"):
            fits.update(
                filename, np.zeros((10, 10)), header=fits.Header([("WHAT", 100)]), ext=1
            )

        fits.HDUList(hdus).writeto(filename)

        fits.update(
            filename, np.zeros((10, 10)), header=fits.Header([("WHAT", 100)]), ext=1
        )
        h_out = fits.getheader(filename, ext=1)
        assert h_out["WHAT"] == 100

    def test_printdiff(self):
        """
        Test that FITSDiff can run the different inputs without crashing.
        """

        # Testing different string input options
        assert printdiff(self.data("arange.fits"), self.data("blank.fits")) is None
        assert (
            printdiff(self.data("arange.fits"), self.data("blank.fits"), ext=0) is None
        )
        assert (
            printdiff(
                self.data("o4sp040b0_raw.fits"),
                self.data("o4sp040b0_raw.fits"),
                extname="sci",
            )
            is None
        )

        # This may seem weird, but check printdiff to see, need to test
        # incorrect second file
        with pytest.raises(OSError):
            printdiff("o4sp040b0_raw.fits", "fakefile.fits", extname="sci")

        # Test HDU object inputs
        with fits.open(self.data("stddata.fits"), mode="readonly") as in1:
            with fits.open(self.data("checksum.fits"), mode="readonly") as in2:
                assert printdiff(in1[0], in2[0]) is None

                with pytest.raises(ValueError):
                    printdiff(in1[0], in2[0], ext=0)

                assert printdiff(in1, in2) is None

                with pytest.raises(NotImplementedError):
                    printdiff(in1, in2, 0)

    def test_tabledump(self):
        """
        A simple test of the dump method.
        Also regression test for https://github.com/astropy/astropy/issues/6937
        """
        datastr = (
            '"                    1" "abc" "     3.70000007152557" "                    0"\n'
            '"                    2" "xy " "     6.69999971389771" "                    1"\n'
        )
        cdstr = (
            'c1 1J I11              ""               ""'
            '               -2147483647      ""               ""              \n'
            'c2 3A A3               ""               ""'
            '               ""               ""               ""              \n'
            'c3 1E G15.7            ""               ""'
            '               ""               3                0.4             \n'
            'c4 1L L6               ""               ""'
            '               ""               ""               ""              \n'
        )
        # copy fits file to the temp directory
        testfile = self.copy_file("tb.fits")

        # test without datafile
        fits.tabledump(testfile)
        assert os.path.isfile(self.temp("tb_1.txt"))

        # test with datafile
        fits.tabledump(testfile, datafile=self.temp("test_tb.txt"))
        assert os.path.isfile(self.temp("test_tb.txt"))

        # test with datafile and cdfile
        datafile = self.temp("data.txt")
        cdfile = self.temp("coldefs.txt")
        fits.tabledump(testfile, datafile, cdfile)
        assert os.path.isfile(datafile)
        with open(datafile) as data:
            assert data.read() == datastr
        with open(cdfile) as coldefs:
            assert coldefs.read() == cdstr

    @pytest.mark.parametrize("tablename", ["table.fits", "tb.fits"])
    def test_dump_load_round_trip(self, tablename):
        """
        A simple test of the dump/load methods; dump the data, column, and
        header files and try to reload the table from them.
        """

        # copy fits file to the temp directory
        testfile = self.copy_file(tablename)

        datafile = self.temp("data.txt")
        cdfile = self.temp("coldefs.txt")
        hfile = self.temp("header.txt")
        fits.tabledump(testfile, datafile, cdfile, hfile)

        new_tbhdu = fits.tableload(datafile, cdfile, hfile)

        with fits.open(testfile) as hdul:
            _assert_attr_col(new_tbhdu, hdul[1])

    def test_append_filename(self, home_is_temp):
        """
        Test fits.append with a filename argument.
        """
        data = np.arange(6)
        testfile = self.temp("test_append_1.fits")

        # Test case 1: creation of file
        fits.append(testfile, data=data, checksum=True)

        # Test case 2: append to existing file, with verify=True
        # Also test that additional keyword can be passed to fitsopen
        fits.append(testfile, data=data * 2, checksum=True, ignore_blank=True)

        # Test case 3: append to existing file, with verify=False
        fits.append(testfile, data=data * 3, checksum=True, verify=False)

        with fits.open(testfile, checksum=True) as hdu1:
            np.testing.assert_array_equal(hdu1[0].data, data)
            np.testing.assert_array_equal(hdu1[1].data, data * 2)
            np.testing.assert_array_equal(hdu1[2].data, data * 3)

    @pytest.mark.parametrize("mode", ["wb", "wb+", "ab", "ab+"])
    def test_append_filehandle(self, tmp_path, mode):
        """
        Test fits.append with a file handle argument.
        """
        append_file = tmp_path / "append.fits"
        with append_file.open(mode) as handle:
            fits.append(filename=handle, data=np.ones((4, 4)))

    def test_append_with_header(self):
        """
        Test fits.append with a fits Header, which triggers detection of the
        HDU class. Regression test for
        https://github.com/astropy/astropy/issues/8660
        """
        testfile = self.temp("test_append_1.fits")
        with fits.open(self.data("test0.fits")) as hdus:
            for hdu in hdus:
                fits.append(testfile, hdu.data, hdu.header, checksum=True)

        with fits.open(testfile, checksum=True) as hdus:
            assert len(hdus) == 5

    def test_pathlib(self):
        testfile = pathlib.Path(self.temp("test.fits"))
        data = np.arange(10)
        hdulist = fits.HDUList([fits.PrimaryHDU(data)])
        hdulist.writeto(testfile)

        with fits.open(testfile) as hdul:
            np.testing.assert_array_equal(hdul[0].data, data)

    def test_getdata_ext_given(self):
        prihdu = fits.PrimaryHDU(data=np.zeros((5, 5), dtype=int))
        exthdu1 = fits.ImageHDU(data=np.ones((5, 5), dtype=int))
        exthdu2 = fits.ImageHDU(data=2 * np.ones((5, 5), dtype=int))
        hdulist = fits.HDUList([prihdu, exthdu1, exthdu2])
        buf = io.BytesIO()
        hdulist.writeto(buf)

        for ext in [0, 1, 2]:
            buf.seek(0)
            data = fits.getdata(buf, ext=ext)
            assert data[0, 0] == ext

    def test_getdata_ext_given_nodata(self):
        prihdu = fits.PrimaryHDU(data=np.zeros((5, 5), dtype=int))
        exthdu1 = fits.ImageHDU(data=np.ones((5, 5), dtype=int))
        exthdu2 = fits.ImageHDU(data=None)
        hdulist = fits.HDUList([prihdu, exthdu1, exthdu2])
        buf = io.BytesIO()
        hdulist.writeto(buf)
        buf.seek(0)

        with pytest.raises(IndexError, match="No data in HDU #2."):
            fits.getdata(buf, ext=2)

    def test_getdata_ext_not_given_with_data_in_primary(self):
        prihdu = fits.PrimaryHDU(data=np.zeros((5, 5), dtype=int))
        exthdu1 = fits.ImageHDU(data=None)
        exthdu2 = fits.ImageHDU(data=None)
        hdulist = fits.HDUList([prihdu, exthdu1, exthdu2])
        buf = io.BytesIO()
        hdulist.writeto(buf)
        buf.seek(0)

        data = fits.getdata(buf)
        assert data[0, 0] == 0

    def test_getdata_ext_not_given_with_data_in_ext(self):
        # tests fallback mechanism
        prihdu = fits.PrimaryHDU(data=None)
        exthdu1 = fits.ImageHDU(data=np.ones((5, 5), dtype=int))
        exthdu2 = fits.ImageHDU(data=None)
        hdulist = fits.HDUList([prihdu, exthdu1, exthdu2])
        buf = io.BytesIO()
        hdulist.writeto(buf)
        buf.seek(0)

        data = fits.getdata(buf)
        assert data[0, 0] == 1

    def test_getdata_ext_not_given_nodata_any(self):
        # tests exception raised when there is no data in either
        # Primary HDU or first extension HDU
        prihdu = fits.PrimaryHDU(data=None)
        exthdu1 = fits.ImageHDU(data=None)
        exthdu2 = fits.ImageHDU(data=np.ones((5, 5), dtype=int))
        hdulist = fits.HDUList([prihdu, exthdu1, exthdu2])
        buf = io.BytesIO()
        hdulist.writeto(buf)
        buf.seek(0)

        with pytest.raises(
            IndexError, match="No data in either Primary or first extension HDUs."
        ):
            fits.getdata(buf)

    def test_getdata_ext_not_given_nodata_noext(self):
        # tests exception raised when there is no data in the
        # Primary HDU and there are no extension HDUs
        prihdu = fits.PrimaryHDU(data=None)
        hdulist = fits.HDUList([prihdu])
        buf = io.BytesIO()
        hdulist.writeto(buf)
        buf.seek(0)

        with pytest.raises(
            IndexError, match="No data in Primary HDU and no extension HDU found."
        ):
            fits.getdata(buf)

    def test_getdata_lower_upper(self, tmp_path):
        t = Table([np.zeros(5), np.ones(5), np.arange(5)], names=["a", "b", "c"])
        t.write(tmp_path / "lower.fits")
        t = Table([np.zeros(5), np.ones(5), np.arange(5)], names=["A", "B", "C"])
        t.write(tmp_path / "upper.fits")

        ref = np.zeros(5)
        lowup = fits.getdata(tmp_path / "lower.fits", 1, upper=True)
        assert lowup.dtype.names == ("A", "B", "C")
        assert_array_equal(lowup["A"], ref)
        assert_array_equal(lowup["a"], ref)
        assert_array_equal(lowup.A, ref)

        uplow = fits.getdata(tmp_path / "upper.fits", 1, lower=True)
        assert uplow.dtype.names == ("a", "b", "c")
        assert_array_equal(uplow["A"], ref)
        assert_array_equal(uplow["a"], ref)
        assert_array_equal(uplow.a, ref)

    @pytest.fixture
    def make_files_for_packing(self):
        # Make random but reproducible data
        np.random.seed(2305235)
        data_header = fits.Header()
        data_header["a"] = "b"
        data_header["c"] = "d"
        nx, ny = 512, 517
        data_hdu = fits.PrimaryHDU(
            data=np.random.poisson(1000, size=(ny, nx)).astype(np.int32),
            header=data_header,
        )
        data_hdu.add_checksum()

        table = Table({"col1": ["foo", "bar"], "col2": ["test1", "test2"]})
        table_hdu = fits.BinTableHDU(table)
        table_hdu.add_checksum()

        # Define 2 cases. One with the primary HDU having data,
        # and the second will only have a header
        self.hdulist_pack1 = fits.HDUList([data_hdu, table_hdu])

        # Case 2
        primary_header = fits.Header()
        primary_header["e"] = "f"
        primary_header["g"] = "h"
        primary_header["i"] = "j"

        primary_hdu = fits.PrimaryHDU(header=primary_header)
        primary_hdu.add_checksum()

        data_header = fits.Header()
        data_header["k"] = "l"
        data_header["m"] = "n"
        nx, ny = 479, 492
        data_hdu = fits.ImageHDU(
            data=np.random.poisson(1500, size=(ny, nx)).astype(np.int32),
            header=data_header,
        )
        data_hdu.add_checksum()

        table = Table({"col1": ["foo", "bar"], "col2": ["test1", "test2"]})
        table_hdu = fits.BinTableHDU(table)
        table_hdu.add_checksum()

        # This file should stay as the primary header, even in the
        # compressed file
        self.hdulist_pack2 = fits.HDUList([primary_hdu, data_hdu, table_hdu])

    def test_pack_unpack_round_trip(self, make_files_for_packing):
        # Test if fits files can round trip pack and unpack
        packed_hdulist1 = fits.pack(self.hdulist_pack1)
        unpacked_hdulist1 = fits.unpack(packed_hdulist1)
        assert fits.FITSDiff(unpacked_hdulist1, self.hdulist_pack1).identical

        packed_hdulist2 = fits.pack(self.hdulist_pack2)
        unpacked_hdulist2 = fits.unpack(packed_hdulist2)
        assert fits.FITSDiff(unpacked_hdulist2, self.hdulist_pack2).identical

    @pytest.mark.skipif(not HAS_FUNPACK, reason="requires fpack binary from cfitsio")
    def test_packed_funpack_compatability(self, make_files_for_packing):
        # Test that pack creates files cfitsio funpack can unpack
        for i in ["1", "2"]:
            packed_hdulist = fits.pack(getattr(self, f"hdulist_pack{i}"))
            packed_filename = self.temp(f"packed{i}.fits.fz")
            packed_hdulist.writeto(packed_filename, overwrite=True)
            funpack_result = subprocess.run(
                ["funpack", packed_filename],
                capture_output=True,
                text=True,
                check=False,
            )
            assert funpack_result.returncode == 0, (
                f"funpack failed:\n{funpack_result.stderr}"
            )
            with fits.open(packed_filename.replace(".fz", "")) as unpacked_hdulist:
                assert fits.FITSDiff(
                    unpacked_hdulist,
                    getattr(self, f"hdulist_pack{i}"),
                    ignore_keywords=["CHECKSUM"],
                    ignore_comments=["SIMPLE"],
                ).identical

    @pytest.mark.skipif(not HAS_FPACK, reason="requires funpack binary from cfitsio")
    def test_fpacked_unpack_compatability(self, make_files_for_packing):
        # Test that unpack works on a cfitsio fpacked file
        for i in ["1", "2"]:
            packed_filename = self.temp(f"packed{i}.fits")
            getattr(self, f"hdulist_pack{i}").writeto(packed_filename, overwrite=True)
            fpack_result = subprocess.run(
                ["fpack", packed_filename],
                capture_output=True,
                text=True,
                check=False,
            )
            assert fpack_result.returncode == 0, f"fpack failed:\n{fpack_result.stderr}"
            with fits.open(packed_filename + ".fz") as packed_hdulist:
                unpacked_hdulist = fits.unpack(packed_hdulist)
                assert fits.FITSDiff(
                    unpacked_hdulist,
                    getattr(self, f"hdulist_pack{i}"),
                    ignore_keywords=["CHECKSUM"],
                    ignore_comments=["SIMPLE"],
                ).identical

    def test_packed_headers(self, make_files_for_packing):
        # Test that the final packed headers are what we expect
        expected_primary_headers = {
            "1": [
                (b"SIMPLE", b"T"),
                (b"BITPIX", b"8"),
                (b"NAXIS", b"0"),
                (b"EXTEND", b"T"),
            ],
            "2": [
                (b"SIMPLE", b"T"),
                (b"BITPIX", b"8"),
                (b"NAXIS", b"0"),
                (b"EXTEND", b"T"),
                (b"E", b"'f       '"),
                (b"G", b"'h       '"),
                (b"I", b"'j       '"),
            ],
        }
        expected_image_headers = {
            "1": [
                (b"XTENSION", b"'BINTABLE'"),
                (b"BITPIX", b"8"),
                (b"NAXIS", b"2"),
                (b"NAXIS1", b"8"),
                (b"NAXIS2", b"517"),
                (b"PCOUNT", b"259822"),
                (b"GCOUNT", b"1"),
                (b"TFIELDS", b"1"),
                (b"TTYPE1", b"'COMPRESSED_DATA'"),
                (b"TFORM1", b"'1PB(511)'"),
                (b"ZIMAGE", b"T"),
                (b"ZSIMPLE", b"T"),
                (b"ZBITPIX", b"32"),
                (b"ZNAXIS", b"2"),
                (b"ZNAXIS1", b"512"),
                (b"ZNAXIS2", b"517"),
                (b"ZTILE1", b"512"),
                (b"ZTILE2", b"1"),
                (b"ZCMPTYPE", b"'RICE_1  '"),
                (b"ZNAME1", b"'BLOCKSIZE'"),
                (b"ZVAL1", b"32"),
                (b"ZNAME2", b"'BYTEPIX '"),
                (b"ZVAL2", b"4"),
                (b"EXTNAME", b"'COMPRESSED_IMAGE'"),
                (b"A", b"'b       '"),
                (b"C", b"'d       '"),
            ],
            "2": [
                (b"XTENSION", b"'BINTABLE'"),
                (b"BITPIX", b"8"),
                (b"NAXIS", b"2"),
                (b"NAXIS1", b"8"),
                (b"NAXIS2", b"492"),
                (b"PCOUNT", b"239433"),
                (b"GCOUNT", b"1"),
                (b"TFIELDS", b"1"),
                (b"TTYPE1", b"'COMPRESSED_DATA'"),
                (b"TFORM1", b"'1PB(495)'"),
                (b"ZIMAGE", b"T"),
                (b"ZTENSION", b"'IMAGE   '"),
                (b"ZBITPIX", b"32"),
                (b"ZNAXIS", b"2"),
                (b"ZNAXIS1", b"479"),
                (b"ZNAXIS2", b"492"),
                (b"ZPCOUNT", b"0"),
                (b"ZGCOUNT", b"1"),
                (b"ZTILE1", b"479"),
                (b"ZTILE2", b"1"),
                (b"ZCMPTYPE", b"'RICE_1  '"),
                (b"ZNAME1", b"'BLOCKSIZE'"),
                (b"ZVAL1", b"32"),
                (b"ZNAME2", b"'BYTEPIX '"),
                (b"ZVAL2", b"4"),
                (b"EXTNAME", b"'COMPRESSED_IMAGE'"),
                (b"K", b"'l       '"),
                (b"M", b"'n       '"),
            ],
        }

        for i in ["1", "2"]:
            packed_hdulist = fits.pack(getattr(self, f"hdulist_pack{i}"))
            packed_buffer = io.BytesIO()
            packed_hdulist.writeto(packed_buffer)
            packed_buffer.seek(0)
            primary_header = _simple_str_header_parser(packed_buffer.read(2880))
            assert primary_header == expected_primary_headers[i]
            image_header = _simple_str_header_parser(packed_buffer.read(2880))
            assert image_header == expected_image_headers[i]


def _simple_str_header_parser(header_str):
    """
    A simple parser to check that the header on disk makes sense.
    We throw away comments, checksum, and datasum keywords.
    """
    chunk_size = 80
    header = []
    for i in range(0, len(header_str), chunk_size):
        chunk = header_str[i : i + chunk_size]
        bytes_chunk = bytes(chunk)
        if len(bytes_chunk.strip()) == 0:
            continue
        if bytes_chunk.startswith(b"COMMENT"):
            continue
        if bytes_chunk.strip() == b"END":
            break
        key_value = bytes_chunk.split(b"/")[0]
        key, value = key_value.split(b"=", 1)
        key = key.strip()
        value = value.strip()
        if key == b"CHECKSUM" or key == b"DATASUM":
            continue
        if key == b"ZHECKSUM" or key == b"ZDATASUM":
            continue
        header.append((key, value))

    return header
