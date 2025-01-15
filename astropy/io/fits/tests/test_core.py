# Licensed under a 3-clause BSD style license - see PYFITS.rst

import errno
import gzip
import io
import mmap
import os
import pathlib
import shutil
import sys
import urllib.request
import zipfile
from unittest.mock import patch

import numpy as np
import pytest

from astropy.io import fits
from astropy.io.fits.convenience import _getext
from astropy.io.fits.diff import FITSDiff
from astropy.io.fits.file import GZIP_MAGIC, _File
from astropy.io.tests import safeio
from astropy.utils import data

# NOTE: Python can be built without bz2.
from astropy.utils.compat.optional_deps import HAS_BZ2
from astropy.utils.data import conf
from astropy.utils.exceptions import AstropyUserWarning
from astropy.utils.misc import _NOT_OVERWRITING_MSG_MATCH

from .conftest import FitsTestCase

if HAS_BZ2:
    import bz2


class TestCore(FitsTestCase):
    def test_missing_file(self):
        with pytest.raises(OSError):
            fits.open(self.temp("does-not-exist.fits"))

    def test_naxisj_check(self):
        with fits.open(self.data("o4sp040b0_raw.fits")) as hdulist:
            hdulist[1].header["NAXIS3"] = 500

            assert "NAXIS3" in hdulist[1].header
            hdulist.verify("silentfix")
            assert "NAXIS3" not in hdulist[1].header

    def test_byteswap(self):
        p = fits.PrimaryHDU()
        lst = fits.HDUList()

        n = np.array([1, 60000, 0], dtype="u2").astype("i2")

        c = fits.Column(name="foo", format="i2", bscale=1, bzero=32768, array=n)
        t = fits.BinTableHDU.from_columns([c])

        lst.append(p)
        lst.append(t)

        lst.writeto(self.temp("test.fits"), overwrite=True)

        with fits.open(self.temp("test.fits")) as p:
            assert p[1].data[1]["foo"] == 60000.0

    def test_fits_file_path_object(self):
        """
        Testing when fits file is passed as pathlib.Path object #4412.
        """
        fpath = pathlib.Path(self.data("tdim.fits"))
        with fits.open(fpath) as hdulist:
            assert hdulist[0].filebytes() == 2880
            assert hdulist[1].filebytes() == 5760

            with fits.open(self.data("tdim.fits")) as hdulist2:
                assert FITSDiff(hdulist2, hdulist).identical is True

    def test_fits_pathlike_object(self):
        """
        Testing when fits file is passed as os.PathLike object #11579.
        """

        class TPath(os.PathLike):
            def __init__(self, path):
                self.path = path

            def __fspath__(self):
                return str(self.path)

        fpath = TPath(self.data("tdim.fits"))
        with fits.open(fpath) as hdulist:
            assert hdulist[0].filebytes() == 2880
            assert hdulist[1].filebytes() == 5760

            with fits.open(self.data("tdim.fits")) as hdulist2:
                assert FITSDiff(hdulist2, hdulist).identical is True

    def test_fits_file_bytes_object(self):
        """
        Testing when fits file is passed as bytes.
        """
        with fits.open(self.data("tdim.fits").encode()) as hdulist:
            assert hdulist[0].filebytes() == 2880
            assert hdulist[1].filebytes() == 5760

            with fits.open(self.data("tdim.fits")) as hdulist2:
                assert FITSDiff(hdulist2, hdulist).identical is True

    def test_add_del_columns(self):
        p = fits.ColDefs([])
        p.add_col(fits.Column(name="FOO", format="3J"))
        p.add_col(fits.Column(name="BAR", format="1I"))
        assert p.names == ["FOO", "BAR"]
        p.del_col("FOO")
        assert p.names == ["BAR"]

    def test_add_del_columns2(self):
        hdulist = fits.open(self.data("tb.fits"))
        table = hdulist[1]
        assert table.data.dtype.names == ("c1", "c2", "c3", "c4")
        assert table.columns.names == ["c1", "c2", "c3", "c4"]
        table.columns.del_col("c1")
        assert table.data.dtype.names == ("c2", "c3", "c4")
        assert table.columns.names == ["c2", "c3", "c4"]

        table.columns.del_col("c3")
        assert table.data.dtype.names == ("c2", "c4")
        assert table.columns.names == ["c2", "c4"]

        table.columns.add_col(fits.Column("foo", "3J"))
        assert table.data.dtype.names == ("c2", "c4", "foo")
        assert table.columns.names == ["c2", "c4", "foo"]

        hdulist.writeto(self.temp("test.fits"), overwrite=True)
        hdulist.close()
        # NOTE: If you see a warning, might be related to
        # https://github.com/spacetelescope/PyFITS/issues/44
        with fits.open(self.temp("test.fits")) as hdulist:
            table = hdulist[1]
            assert table.data.dtype.names == ("c2", "c4", "foo")
            assert table.columns.names == ["c2", "c4", "foo"]

    def test_update_header_card(self):
        """A very basic test for the Header.update method--I'd like to add a
        few more cases to this at some point.
        """

        header = fits.Header()
        comment = "number of bits per data pixel"
        header["BITPIX"] = (16, comment)
        assert "BITPIX" in header
        assert header["BITPIX"] == 16
        assert header.comments["BITPIX"] == comment

        header.update(BITPIX=32)
        assert header["BITPIX"] == 32
        assert header.comments["BITPIX"] == ""

    def test_set_card_value(self):
        """Similar to test_update_header_card(), but tests the the
        `header['FOO'] = 'bar'` method of updating card values.
        """

        header = fits.Header()
        comment = "number of bits per data pixel"
        card = fits.Card.fromstring(f"BITPIX  = 32 / {comment}")
        header.append(card)

        header["BITPIX"] = 32

        assert "BITPIX" in header
        assert header["BITPIX"] == 32
        assert header.cards[0].keyword == "BITPIX"
        assert header.cards[0].value == 32
        assert header.cards[0].comment == comment

    def test_uint(self):
        filename = self.data("o4sp040b0_raw.fits")
        with fits.open(filename, uint=False) as hdulist_f:
            with fits.open(filename, uint=True) as hdulist_i:
                assert hdulist_f[1].data.dtype == np.float32
                assert hdulist_i[1].data.dtype == np.uint16
                assert np.all(hdulist_f[1].data == hdulist_i[1].data)

    def test_fix_missing_card_append(self):
        hdu = fits.ImageHDU()
        errs = hdu.req_cards("TESTKW", None, None, "foo", "silentfix", [])
        assert len(errs) == 1
        assert "TESTKW" in hdu.header
        assert hdu.header["TESTKW"] == "foo"
        assert hdu.header.cards[-1].keyword == "TESTKW"

    def test_fix_invalid_keyword_value(self):
        hdu = fits.ImageHDU()
        hdu.header["TESTKW"] = "foo"
        errs = hdu.req_cards("TESTKW", None, lambda v: v == "foo", "foo", "ignore", [])
        assert len(errs) == 0

        # Now try a test that will fail, and ensure that an error will be
        # raised in 'exception' mode
        errs = hdu.req_cards(
            "TESTKW", None, lambda v: v == "bar", "bar", "exception", []
        )
        assert len(errs) == 1
        assert errs[0][1] == "'TESTKW' card has invalid value 'foo'."

        # See if fixing will work
        hdu.req_cards("TESTKW", None, lambda v: v == "bar", "bar", "silentfix", [])
        assert hdu.header["TESTKW"] == "bar"

    def test_unfixable_missing_card(self):
        class TestHDU(fits.hdu.base.NonstandardExtHDU):
            def _verify(self, option="warn"):
                errs = super()._verify(option)
                hdu.req_cards("TESTKW", None, None, None, "fix", errs)
                return errs

            @classmethod
            def match_header(cls, header):
                # Since creating this HDU class adds it to the registry we
                # don't want the file reader to possibly think any actual
                # HDU from a file should be handled by this class
                return False

        hdu = TestHDU(header=fits.Header())
        with pytest.raises(fits.VerifyError):
            hdu.verify("fix")

    def test_exception_on_verification_error(self):
        hdu = fits.ImageHDU()
        del hdu.header["XTENSION"]
        with pytest.raises(fits.VerifyError):
            hdu.verify("exception")

    def test_ignore_verification_error(self):
        hdu = fits.ImageHDU()
        del hdu.header["NAXIS"]
        # The default here would be to issue a warning; ensure that no warnings
        # or exceptions are raised
        hdu.verify("ignore")

        # Make sure the error wasn't fixed either, silently or otherwise
        assert "NAXIS" not in hdu.header

    def test_unrecognized_verify_option(self):
        hdu = fits.ImageHDU()
        with pytest.raises(ValueError):
            hdu.verify("foobarbaz")

    def test_errlist_basic(self):
        # Just some tests to make sure that _ErrList is setup correctly.
        # No arguments
        error_list = fits.verify._ErrList()
        assert error_list == []
        # Some contents - this is not actually working, it just makes sure they
        # are kept.
        error_list = fits.verify._ErrList([1, 2, 3])
        assert error_list == [1, 2, 3]

    def test_combined_verify_options(self):
        """
        Test verify options like fix+ignore.
        """

        def make_invalid_hdu():
            hdu = fits.ImageHDU()
            # Add one keyword to the header that contains a fixable defect, and one
            # with an unfixable defect
            c1 = fits.Card.fromstring("test    = '    test'")
            c2 = fits.Card.fromstring("P.I.    = '  Hubble'")
            hdu.header.append(c1)
            hdu.header.append(c2)
            return hdu

        # silentfix+ignore should be completely silent
        hdu = make_invalid_hdu()
        hdu.verify("silentfix+ignore")

        # silentfix+warn should be quiet about the fixed HDU and only warn
        # about the unfixable one
        hdu = make_invalid_hdu()
        with pytest.warns(AstropyUserWarning, match="Illegal keyword name") as w:
            hdu.verify("silentfix+warn")
        assert len(w) == 4

        # silentfix+exception should only mention the unfixable error in the
        # exception
        hdu = make_invalid_hdu()
        with pytest.raises(fits.VerifyError, match=r"Illegal keyword name") as excinfo:
            hdu.verify("silentfix+exception")
        assert "not upper case" not in str(excinfo.value)

        # fix+ignore is not too useful, but it should warn about the fixed
        # problems while saying nothing about the unfixable problems
        hdu = make_invalid_hdu()
        with pytest.warns(AstropyUserWarning, match="not upper case") as w:
            hdu.verify("fix+ignore")
        assert len(w) == 4

        # fix+warn
        hdu = make_invalid_hdu()
        with pytest.warns(AstropyUserWarning) as w:
            hdu.verify("fix+warn")
        assert len(w) == 6
        assert "not upper case" in str(w[2].message)
        assert "Illegal keyword name" in str(w[4].message)

        # fix+exception
        hdu = make_invalid_hdu()
        with pytest.raises(fits.VerifyError, match=r"Illegal keyword name") as excinfo:
            hdu.verify("fix+exception")
        assert "not upper case" in str(excinfo.value)

    def test_getext(self):
        """
        Test the various different ways of specifying an extension header in
        the convenience functions.
        """
        filename = self.data("test0.fits")

        hl, ext = _getext(filename, "readonly", 1)
        assert ext == 1
        hl.close()

        pytest.raises(ValueError, _getext, filename, "readonly", 1, 2)
        pytest.raises(ValueError, _getext, filename, "readonly", (1, 2))
        pytest.raises(ValueError, _getext, filename, "readonly", "sci", "sci")
        pytest.raises(TypeError, _getext, filename, "readonly", 1, 2, 3)

        hl, ext = _getext(filename, "readonly", ext=1)
        assert ext == 1
        hl.close()

        hl, ext = _getext(filename, "readonly", ext=("sci", 2))
        assert ext == ("sci", 2)
        hl.close()

        pytest.raises(
            TypeError, _getext, filename, "readonly", 1, ext=("sci", 2), extver=3
        )
        pytest.raises(
            TypeError, _getext, filename, "readonly", ext=("sci", 2), extver=3
        )

        hl, ext = _getext(filename, "readonly", "sci")
        assert ext == ("sci", 1)
        hl.close()

        hl, ext = _getext(filename, "readonly", "sci", 1)
        assert ext == ("sci", 1)
        hl.close()

        hl, ext = _getext(filename, "readonly", ("sci", 1))
        assert ext == ("sci", 1)
        hl.close()

        hl, ext = _getext(
            filename, "readonly", "sci", extver=1, do_not_scale_image_data=True
        )
        assert ext == ("sci", 1)
        hl.close()

        pytest.raises(TypeError, _getext, filename, "readonly", "sci", ext=1)
        pytest.raises(TypeError, _getext, filename, "readonly", "sci", 1, extver=2)

        hl, ext = _getext(filename, "readonly", extname="sci")
        assert ext == ("sci", 1)
        hl.close()

        hl, ext = _getext(filename, "readonly", extname="sci", extver=1)
        assert ext == ("sci", 1)
        hl.close()

        pytest.raises(TypeError, _getext, filename, "readonly", extver=1)

    def test_extension_name_case_sensitive(self):
        """
        Tests that setting fits.conf.extension_name_case_sensitive at
        runtime works.
        """

        hdu = fits.ImageHDU()
        hdu.name = "sCi"
        assert hdu.name == "SCI"
        assert hdu.header["EXTNAME"] == "SCI"

        with fits.conf.set_temp("extension_name_case_sensitive", True):
            hdu = fits.ImageHDU()
            hdu.name = "sCi"
            assert hdu.name == "sCi"
            assert hdu.header["EXTNAME"] == "sCi"

        hdu.name = "sCi"
        assert hdu.name == "SCI"
        assert hdu.header["EXTNAME"] == "SCI"

    def test_hdu_fromstring(self):
        """
        Tests creating a fully-formed HDU object from a string containing the
        bytes of the HDU.
        """
        infile = self.data("test0.fits")
        outfile = self.temp("test.fits")

        with open(infile, "rb") as fin:
            dat = fin.read()

        offset = 0
        with fits.open(infile) as hdul:
            hdulen = hdul[0]._data_offset + hdul[0]._data_size
            hdu = fits.PrimaryHDU.fromstring(dat[:hdulen])
            assert isinstance(hdu, fits.PrimaryHDU)
            assert hdul[0].header == hdu.header
            assert hdu.data is None

        hdu.header["TEST"] = "TEST"
        hdu.writeto(outfile)
        with fits.open(outfile) as hdul:
            assert isinstance(hdu, fits.PrimaryHDU)
            assert hdul[0].header[:-1] == hdu.header[:-1]
            assert hdul[0].header["TEST"] == "TEST"
            assert hdu.data is None

        with fits.open(infile) as hdul:
            for ext_hdu in hdul[1:]:
                offset += hdulen
                hdulen = len(str(ext_hdu.header)) + ext_hdu._data_size
                hdu = fits.ImageHDU.fromstring(dat[offset : offset + hdulen])
                assert isinstance(hdu, fits.ImageHDU)
                assert ext_hdu.header == hdu.header
                assert (ext_hdu.data == hdu.data).all()

    def test_nonstandard_hdu(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/157

        Tests that "Nonstandard" HDUs with SIMPLE = F are read and written
        without prepending a superfluous and unwanted standard primary HDU.
        """

        data = np.arange(100, dtype=np.uint8)
        hdu = fits.PrimaryHDU(data=data)
        hdu.header["SIMPLE"] = False
        hdu.writeto(self.temp("test.fits"))

        info = [(0, "", 1, "NonstandardHDU", 5, (), "", "")]
        with fits.open(self.temp("test.fits")) as hdul:
            assert hdul.info(output=False) == info
            # NonstandardHDUs just treat the data as an unspecified array of
            # bytes.  The first 100 bytes should match the original data we
            # passed in...the rest should be zeros padding out the rest of the
            # FITS block
            assert (hdul[0].data[:100] == data).all()
            assert (hdul[0].data[100:] == 0).all()

    def test_extname(self):
        """Test getting/setting the EXTNAME of an HDU."""

        h1 = fits.PrimaryHDU()
        assert h1.name == "PRIMARY"
        # Normally a PRIMARY HDU should not have an EXTNAME, though it should
        # have a default .name attribute
        assert "EXTNAME" not in h1.header

        # The current version of the FITS standard does allow PRIMARY HDUs to
        # have an EXTNAME, however.
        h1.name = "NOTREAL"
        assert h1.name == "NOTREAL"
        assert h1.header.get("EXTNAME") == "NOTREAL"

        # Updating the EXTNAME in the header should update the .name
        h1.header["EXTNAME"] = "TOOREAL"
        assert h1.name == "TOOREAL"

        # If we delete an EXTNAME keyword from a PRIMARY HDU it should go back
        # to the default
        del h1.header["EXTNAME"]
        assert h1.name == "PRIMARY"

        # For extension HDUs the situation is a bit simpler:
        h2 = fits.ImageHDU()
        assert h2.name == ""
        assert "EXTNAME" not in h2.header
        h2.name = "HELLO"
        assert h2.name == "HELLO"
        assert h2.header.get("EXTNAME") == "HELLO"
        h2.header["EXTNAME"] = "GOODBYE"
        assert h2.name == "GOODBYE"

    def test_extver_extlevel(self):
        """Test getting/setting the EXTVER and EXTLEVEL of and HDU."""

        # EXTVER and EXTNAME work exactly the same; their semantics are, for
        # now, to be inferred by the user.  Although they should never be less
        # than 1, the standard does not explicitly forbid any value so long as
        # it's an integer
        h1 = fits.PrimaryHDU()
        assert h1.ver == 1
        assert h1.level == 1
        assert "EXTVER" not in h1.header
        assert "EXTLEVEL" not in h1.header

        h1.ver = 2
        assert h1.header.get("EXTVER") == 2
        h1.header["EXTVER"] = 3
        assert h1.ver == 3
        del h1.header["EXTVER"]
        h1.ver == 1

        h1.level = 2
        assert h1.header.get("EXTLEVEL") == 2
        h1.header["EXTLEVEL"] = 3
        assert h1.level == 3
        del h1.header["EXTLEVEL"]
        assert h1.level == 1

        pytest.raises(TypeError, setattr, h1, "ver", "FOO")
        pytest.raises(TypeError, setattr, h1, "level", "BAR")

    def test_consecutive_writeto(self):
        """
        Regression test for an issue where calling writeto twice on the same
        HDUList could write a corrupted file.

        https://github.com/spacetelescope/PyFITS/issues/40 is actually a
        particular instance of this problem, though isn't unique to sys.stdout.
        """

        with fits.open(self.data("test0.fits")) as hdul1:
            # Add a bunch of header keywords so that the data will be forced to
            # new offsets within the file:
            for idx in range(40):
                hdul1[1].header[f"TEST{idx}"] = "test"

            hdul1.writeto(self.temp("test1.fits"))
            hdul1.writeto(self.temp("test2.fits"))

            # Open a second handle to the original file and compare it to hdul1
            # (We only compare part of the one header that was modified)
            # Compare also with the second writeto output
            with fits.open(self.data("test0.fits")) as hdul2:
                with fits.open(self.temp("test2.fits")) as hdul3:
                    for hdul in (hdul1, hdul3):
                        for idx, hdus in enumerate(zip(hdul2, hdul)):
                            hdu2, hdu = hdus
                            if idx != 1:
                                assert hdu.header == hdu2.header
                            else:
                                assert hdu2.header == hdu.header[: len(hdu2.header)]
                            assert np.all(hdu.data == hdu2.data)


class TestConvenienceFunctions(FitsTestCase):
    def test_writeto(self, home_is_temp):
        """
        Simple test for writing a trivial header and some data to a file
        with the `writeto()` convenience function.
        """
        filename = self.temp("array.fits")
        data = np.zeros((100, 100))
        header = fits.Header()
        fits.writeto(filename, data, header=header, overwrite=True)
        with fits.open(filename) as hdul:
            assert len(hdul) == 1
            assert (data == hdul[0].data).all()

    def test_writeto_2(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/107

        Test of `writeto()` with a trivial header containing a single keyword.
        """
        filename = self.temp("array.fits")
        data = np.zeros((100, 100))
        header = fits.Header()
        header.set("CRPIX1", 1.0)
        fits.writeto(
            filename, data, header=header, overwrite=True, output_verify="silentfix"
        )
        with fits.open(filename) as hdul:
            assert len(hdul) == 1
            assert (data == hdul[0].data).all()
            assert "CRPIX1" in hdul[0].header
            assert hdul[0].header["CRPIX1"] == 1.0

    def test_writeto_overwrite(self, home_is_temp):
        """
        Ensure the `overwrite` keyword works as it should
        """
        filename = self.temp("array.fits")
        data = np.zeros((100, 100))
        header = fits.Header()
        fits.writeto(filename, data, header=header)

        with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
            fits.writeto(filename, data, header=header, overwrite=False)

        fits.writeto(filename, data, header=header, overwrite=True)

        with fits.open(filename) as hdul:
            assert len(hdul) == 1
            assert (data == hdul[0].data).all()


class TestFileFunctions(FitsTestCase):
    """
    Tests various basic I/O operations, specifically in the
    astropy.io.fits.file._File class.
    """

    def test_open_nonexistent(self):
        """Test that trying to open a non-existent file results in an
        OSError (and not some other arbitrary exception).
        """

        with pytest.raises(OSError, match=r"No such file or directory"):
            fits.open(self.temp("foobar.fits"))

        # But opening in ostream or append mode should be okay, since they
        # allow writing new files
        for mode in ("ostream", "append"):
            with fits.open(self.temp("foobar.fits"), mode=mode) as _:
                pass

            assert os.path.exists(self.temp("foobar.fits"))
            os.remove(self.temp("foobar.fits"))

    def test_open_file_handle(self):
        # Make sure we can open a FITS file from an open file handle
        with open(self.data("test0.fits"), "rb") as handle:
            with fits.open(handle) as _:
                pass

        with open(self.temp("temp.fits"), "wb") as handle:
            with fits.open(handle, mode="ostream") as _:
                pass

        # Opening without explicitly specifying binary mode should fail
        with pytest.raises(ValueError):
            with open(self.data("test0.fits")) as handle:
                with fits.open(handle) as _:
                    pass

        # All of these read modes should fail
        for mode in ["r", "rt"]:
            with pytest.raises(ValueError):
                with open(self.data("test0.fits"), mode=mode) as handle:
                    with fits.open(handle) as _:
                        pass

        # These update or write modes should fail as well
        for mode in ["w", "wt", "w+", "wt+", "r+", "rt+", "a", "at", "a+", "at+"]:
            with pytest.raises(ValueError):
                with open(self.temp("temp.fits"), mode=mode) as handle:
                    with fits.open(handle) as _:
                        pass

    def test_fits_file_handle_mode_combo(self):
        # This should work fine since no mode is given
        with open(self.data("test0.fits"), "rb") as handle:
            with fits.open(handle) as _:
                pass

        # This should work fine since the modes are compatible
        with open(self.data("test0.fits"), "rb") as handle:
            with fits.open(handle, mode="readonly") as _:
                pass

        # This should not work since the modes conflict
        with pytest.raises(ValueError):
            with open(self.data("test0.fits"), "rb") as handle:
                with fits.open(handle, mode="ostream") as _:
                    pass

    def test_open_from_url(self):
        file_url = "file:///" + self.data("test0.fits").lstrip("/")
        with urllib.request.urlopen(file_url) as urlobj:
            with fits.open(urlobj) as _:
                pass

        # It will not be possible to write to a file that is from a URL object
        for mode in ("ostream", "append", "update"):
            with pytest.raises(ValueError):
                with urllib.request.urlopen(file_url) as urlobj:
                    with fits.open(urlobj, mode=mode) as _:
                        pass

    @pytest.mark.remote_data(source="astropy")
    def test_open_from_remote_url(self):
        for dataurl in (conf.dataurl, conf.dataurl_mirror):
            remote_url = f"{dataurl}/allsky/allsky_rosat.fits"
            try:
                with urllib.request.urlopen(remote_url) as urlobj:
                    with fits.open(urlobj) as fits_handle:
                        assert len(fits_handle) == 1

                for mode in ("ostream", "append", "update"):
                    with pytest.raises(ValueError):
                        with urllib.request.urlopen(remote_url) as urlobj:
                            with fits.open(urlobj, mode=mode) as fits_handle:
                                assert len(fits_handle) == 1
            except (urllib.error.HTTPError, urllib.error.URLError):
                continue
            else:
                break
        else:
            raise Exception("Could not download file")

    def test_open_gzipped(self):
        gzip_file = self._make_gzip_file()
        with fits.open(gzip_file) as fits_handle:
            assert fits_handle._file.compression == "gzip"
            assert len(fits_handle) == 5
        with fits.open(gzip.GzipFile(gzip_file)) as fits_handle:
            assert fits_handle._file.compression == "gzip"
            assert len(fits_handle) == 5

    def test_open_gzipped_from_handle(self):
        with open(self._make_gzip_file(), "rb") as handle:
            with fits.open(handle) as fits_handle:
                assert fits_handle._file.compression == "gzip"

    def test_detect_gzipped(self):
        """Test detection of a gzip file when the extension is not .gz."""
        with fits.open(self._make_gzip_file("test0.fz")) as fits_handle:
            assert fits_handle._file.compression == "gzip"
            assert len(fits_handle) == 5

    def test_writeto_append_mode_gzip(self):
        """Regression test for
        https://github.com/spacetelescope/PyFITS/issues/33

        Check that a new GzipFile opened in append mode can be used to write
        out a new FITS file.
        """

        # Note: when opening a GzipFile the 'b+' is superfluous, but this was
        # still how the original test case looked
        # Note: with statement not supported on GzipFile in older Python
        # versions
        fileobj = gzip.GzipFile(self.temp("test.fits.gz"), "ab+")
        h = fits.PrimaryHDU()
        try:
            h.writeto(fileobj)
        finally:
            fileobj.close()

        with fits.open(self.temp("test.fits.gz")) as hdul:
            assert hdul[0].header == h.header

    def test_fits_update_mode_gzip(self):
        """Test updating a GZipped FITS file"""

        with fits.open(self._make_gzip_file("update.gz"), mode="update") as fits_handle:
            hdu = fits.ImageHDU(data=[x for x in range(100)])
            fits_handle.append(hdu)

        with fits.open(self.temp("update.gz")) as new_handle:
            assert len(new_handle) == 6
            assert (new_handle[-1].data == [x for x in range(100)]).all()

    def test_fits_append_mode_gzip(self):
        """Make sure that attempting to open an existing GZipped FITS file in
        'append' mode raises an error"""

        with pytest.raises(OSError):
            with fits.open(self._make_gzip_file("append.gz"), mode="append") as _:
                pass

    @pytest.mark.skipif(not HAS_BZ2, reason="Python built without bz2 module")
    def test_open_bzipped(self):
        bzip_file = self._make_bzip2_file()
        with fits.open(bzip_file) as fits_handle:
            assert fits_handle._file.compression == "bzip2"
            assert len(fits_handle) == 5

        with fits.open(bz2.BZ2File(bzip_file)) as fits_handle:
            assert fits_handle._file.compression == "bzip2"
            assert len(fits_handle) == 5

    @pytest.mark.skipif(not HAS_BZ2, reason="Python built without bz2 module")
    def test_open_bzipped_from_handle(self):
        with open(self._make_bzip2_file(), "rb") as handle:
            with fits.open(handle) as fits_handle:
                assert fits_handle._file.compression == "bzip2"
                assert len(fits_handle) == 5

    @pytest.mark.skipif(not HAS_BZ2, reason="Python built without bz2 module")
    def test_detect_bzipped(self):
        """Test detection of a bzip2 file when the extension is not .bz2."""
        with fits.open(self._make_bzip2_file("test0.xx")) as fits_handle:
            assert fits_handle._file.compression == "bzip2"
            assert len(fits_handle) == 5

    @pytest.mark.skipif(not HAS_BZ2, reason="Python built without bz2 module")
    def test_writeto_bzip2_fileobj(self):
        """Test writing to a bz2.BZ2File file like object"""
        fileobj = bz2.BZ2File(self.temp("test.fits.bz2"), "w")
        h = fits.PrimaryHDU()
        try:
            h.writeto(fileobj)
        finally:
            fileobj.close()

        with fits.open(self.temp("test.fits.bz2")) as hdul:
            assert hdul[0].header == h.header

    @pytest.mark.skipif(not HAS_BZ2, reason="Python built without bz2 module")
    def test_writeto_bzip2_filename(self):
        """Test writing to a bzip2 file by name"""
        filename = self.temp("testname.fits.bz2")
        h = fits.PrimaryHDU()
        h.writeto(filename)

        with fits.open(self.temp("testname.fits.bz2")) as hdul:
            assert hdul[0].header == h.header

    def test_open_zipped(self):
        zip_file = self._make_zip_file()
        with fits.open(zip_file) as fits_handle:
            assert fits_handle._file.compression == "zip"
            assert len(fits_handle) == 5
        with fits.open(zipfile.ZipFile(zip_file)) as fits_handle:
            assert fits_handle._file.compression == "zip"
            assert len(fits_handle) == 5

    def test_open_zipped_from_handle(self):
        with open(self._make_zip_file(), "rb") as handle:
            with fits.open(handle) as fits_handle:
                assert fits_handle._file.compression == "zip"
                assert len(fits_handle) == 5

    def test_detect_zipped(self):
        """Test detection of a zip file when the extension is not .zip."""

        zf = self._make_zip_file(filename="test0.fz")
        with fits.open(zf) as fits_handle:
            assert len(fits_handle) == 5

    def test_open_zipped_writeable(self):
        """Opening zipped files in a writeable mode should fail."""

        zf = self._make_zip_file()
        pytest.raises(OSError, fits.open, zf, "update")
        pytest.raises(OSError, fits.open, zf, "append")

        zf = zipfile.ZipFile(zf, "a")
        pytest.raises(OSError, fits.open, zf, "update")
        pytest.raises(OSError, fits.open, zf, "append")

    def test_read_open_astropy_gzip_file(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/2774

        This tests reading from a ``GzipFile`` object from Astropy's
        compatibility copy of the ``gzip`` module.
        """
        gf = gzip.GzipFile(self._make_gzip_file())
        try:
            assert len(fits.open(gf)) == 5
        finally:
            gf.close()

    def test_open_multiple_member_zipfile(self):
        """
        Opening zip files containing more than one member files should fail
        as there's no obvious way to specify which file is the FITS file to
        read.
        """

        zfile = zipfile.ZipFile(self.temp("test0.zip"), "w")
        zfile.write(self.data("test0.fits"))
        zfile.writestr("foo", "bar")
        zfile.close()

        with pytest.raises(OSError):
            fits.open(zfile.filename)

    def test_read_open_file(self):
        """Read from an existing file object."""

        with open(self.data("test0.fits"), "rb") as f:
            assert len(fits.open(f)) == 5

    def test_read_closed_file(self):
        """Read from an existing file object that's been closed."""

        f = open(self.data("test0.fits"), "rb")
        f.close()
        with fits.open(f) as f2:
            assert len(f2) == 5

    def test_read_open_gzip_file(self):
        """Read from an open gzip file object."""

        gf = gzip.GzipFile(self._make_gzip_file())
        try:
            assert len(fits.open(gf)) == 5
        finally:
            gf.close()

    def test_open_gzip_file_for_writing(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/195."""

        gf = self._make_gzip_file()
        with fits.open(gf, mode="update") as h:
            h[0].header["EXPFLAG"] = "ABNORMAL"
            h[1].data[0, 0] = 1
        with fits.open(gf) as h:
            # Just to make sure the update worked; if updates work
            # normal writes should work too...
            assert h[0].header["EXPFLAG"] == "ABNORMAL"
            assert h[1].data[0, 0] == 1

    def test_write_read_gzip_file(self, home_is_temp):
        """
        Regression test for https://github.com/astropy/astropy/issues/2794

        Ensure files written through gzip are readable.
        """

        data = np.arange(100)
        hdu = fits.PrimaryHDU(data=data)
        hdu.writeto(self.temp("test.fits.gz"))

        with open(os.path.expanduser(self.temp("test.fits.gz")), "rb") as f:
            assert f.read(3) == GZIP_MAGIC

        with fits.open(self.temp("test.fits.gz")) as hdul:
            assert np.all(hdul[0].data == data)

    @pytest.mark.parametrize("ext", ["gz", "bz2", "zip"])
    def test_compressed_ext_but_not_compressed(self, ext):
        testfile = self.temp(f"test0.fits.{ext}")
        shutil.copy(self.data("test0.fits"), testfile)

        with fits.open(testfile) as hdul:
            assert len(hdul) == 5

        fits.append(testfile, np.arange(5))

        with fits.open(testfile) as hdul:
            assert len(hdul) == 6

    def test_read_file_like_object(self):
        """Test reading a FITS file from a file-like object."""

        filelike = io.BytesIO()
        with open(self.data("test0.fits"), "rb") as f:
            filelike.write(f.read())
        filelike.seek(0)
        assert len(fits.open(filelike)) == 5

    def test_updated_file_permissions(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/79

        Tests that when a FITS file is modified in update mode, the file
        permissions are preserved.
        """

        filename = self.temp("test.fits")
        hdul = [fits.PrimaryHDU(), fits.ImageHDU()]
        hdul = fits.HDUList(hdul)
        hdul.writeto(filename)

        old_mode = os.stat(filename).st_mode

        hdul = fits.open(filename, mode="update")
        hdul.insert(1, fits.ImageHDU())
        hdul.flush()
        hdul.close()

        assert old_mode == os.stat(filename).st_mode

    def test_fileobj_mode_guessing(self):
        """Tests whether a file opened without a specified io.fits mode
        ('readonly', etc.) is opened in a mode appropriate for the given file
        object.
        """

        self.copy_file("test0.fits")

        # Opening in text mode should outright fail
        for mode in ("r", "w", "a"):
            with open(self.temp("test0.fits"), mode) as f:
                pytest.raises(ValueError, fits.HDUList.fromfile, f)

        # Need to re-copy the file since opening it in 'w' mode blew it away
        self.copy_file("test0.fits")

        with open(self.temp("test0.fits"), "rb") as f:
            with fits.HDUList.fromfile(f) as h:
                assert h.fileinfo(0)["filemode"] == "readonly"

        for mode in ("wb", "ab"):
            with open(self.temp("test0.fits"), mode) as f:
                with fits.HDUList.fromfile(f) as h:
                    # Basically opening empty files for output streaming
                    assert len(h) == 0

        # Need to re-copy the file since opening it in 'w' mode blew it away
        self.copy_file("test0.fits")

        with open(self.temp("test0.fits"), "wb+") as f:
            with fits.HDUList.fromfile(f) as h:
                # wb+ still causes an existing file to be overwritten so there
                # are no HDUs
                assert len(h) == 0

        # Need to re-copy the file since opening it in 'w' mode blew it away
        self.copy_file("test0.fits")

        with open(self.temp("test0.fits"), "rb+") as f:
            with fits.HDUList.fromfile(f) as h:
                assert h.fileinfo(0)["filemode"] == "update"

        with open(self.temp("test0.fits"), "ab+") as f:
            with fits.HDUList.fromfile(f) as h:
                assert h.fileinfo(0)["filemode"] == "append"

    def test_mmap_unwriteable(self):
        """Regression test for https://github.com/astropy/astropy/issues/968

        Temporarily patches mmap.mmap to exhibit platform-specific bad
        behavior.
        """

        class MockMmap(mmap.mmap):
            def flush(self):
                raise OSError("flush is broken on this platform")

        old_mmap = mmap.mmap
        mmap.mmap = MockMmap

        # Force the mmap test to be rerun
        _File.__dict__["_mmap_available"]._cache.clear()

        try:
            self.copy_file("test0.fits")
            with pytest.warns(
                AstropyUserWarning, match=r"mmap\.flush is unavailable"
            ) as w:
                with fits.open(
                    self.temp("test0.fits"), mode="update", memmap=True
                ) as h:
                    h[1].data[0, 0] = 999

            assert len(w) == 1

            # Double check that writing without mmap still worked
            with fits.open(self.temp("test0.fits")) as h:
                assert h[1].data[0, 0] == 999
        finally:
            mmap.mmap = old_mmap
            _File.__dict__["_mmap_available"]._cache.clear()

    @pytest.mark.openfiles_ignore
    def test_mmap_allocate_error(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/1380

        Temporarily patches mmap.mmap to raise an OSError if mode is ACCESS_COPY.
        """

        mmap_original = mmap.mmap

        # We patch mmap here to raise an error if access=mmap.ACCESS_COPY, which
        # emulates an issue that an OSError is raised if the available address
        # space is less than the size of the file even if memory mapping is used.

        def mmap_patched(*args, **kwargs):
            if kwargs.get("access") == mmap.ACCESS_COPY:
                exc = OSError()
                exc.errno = errno.ENOMEM
                raise exc
            else:
                return mmap_original(*args, **kwargs)

        with fits.open(self.data("test0.fits"), memmap=True) as hdulist:
            with patch.object(mmap, "mmap", side_effect=mmap_patched) as p:
                with pytest.warns(
                    AstropyUserWarning,
                    match=r"Could not memory " r"map array with mode='readonly'",
                ):
                    data = hdulist[1].data
                p.reset_mock()
            assert not data.flags.writeable

    def test_mmap_closing(self):
        """
        Tests that the mmap reference is closed/removed when there aren't any
        HDU data references left.
        """

        if not _File._mmap_available:
            pytest.xfail("not expected to work on platforms without mmap support")

        with fits.open(self.data("test0.fits"), memmap=True) as hdul:
            assert hdul._file._mmap is None

            hdul[1].data
            assert hdul._file._mmap is not None

            del hdul[1].data
            # Should be no more references to data in the file so close the
            # mmap
            assert hdul._file._mmap is None

            hdul[1].data
            hdul[2].data
            del hdul[1].data
            # hdul[2].data is still references so keep the mmap open
            assert hdul._file._mmap is not None
            del hdul[2].data
            assert hdul._file._mmap is None

        assert hdul._file._mmap is None

        with fits.open(self.data("test0.fits"), memmap=True) as hdul:
            hdul[1].data

        # When the only reference to the data is on the hdu object, and the
        # hdulist it belongs to has been closed, the mmap should be closed as
        # well
        assert hdul._file._mmap is None

        with fits.open(self.data("test0.fits"), memmap=True) as hdul:
            data = hdul[1].data
            # also make a copy
            data_copy = data.copy()

        # The HDUList is closed; in fact, get rid of it completely
        del hdul

        # The data array should still work though...
        assert np.all(data == data_copy)

    def test_uncloseable_file(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/2356

        Demonstrates that FITS files can still be read from file-like objects
        that don't have an obvious "open" or "closed" state.
        """

        class MyFileLike:
            def __init__(self, foobar):
                self._foobar = foobar

            def read(self, n):
                return self._foobar.read(n)

            def seek(self, offset, whence=os.SEEK_SET):
                self._foobar.seek(offset, whence)

            def tell(self):
                return self._foobar.tell()

        with open(self.data("test0.fits"), "rb") as f:
            fileobj = MyFileLike(f)

            with fits.open(fileobj) as hdul1:
                with fits.open(self.data("test0.fits")) as hdul2:
                    assert hdul1.info(output=False) == hdul2.info(output=False)
                    for hdu1, hdu2 in zip(hdul1, hdul2):
                        assert hdu1.header == hdu2.header
                        if hdu1.data is not None and hdu2.data is not None:
                            assert np.all(hdu1.data == hdu2.data)

    def test_write_bytesio_discontiguous(self):
        """
        Regression test related to
        https://github.com/astropy/astropy/issues/2794#issuecomment-55441539

        Demonstrates that writing an HDU containing a discontiguous Numpy array
        should work properly.
        """

        data = np.arange(100)[::3]
        hdu = fits.PrimaryHDU(data=data)
        fileobj = io.BytesIO()
        hdu.writeto(fileobj)

        fileobj.seek(0)

        with fits.open(fileobj) as h:
            assert np.all(h[0].data == data)

    def test_write_bytesio(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/2463

        Test against `io.BytesIO`.  `io.StringIO` is not supported.
        """

        self._test_write_string_bytes_io(io.BytesIO())

    @pytest.mark.skipif(
        sys.platform.startswith("win32"), reason="Cannot test on Windows"
    )
    def test_filename_with_colon(self):
        """
        Test reading and writing a file with a colon in the filename.

        Regression test for https://github.com/astropy/astropy/issues/3122
        """

        # Skip on Windows since colons in filenames makes NTFS sad.

        filename = "APEXHET.2014-04-01T15:18:01.000.fits"
        hdu = fits.PrimaryHDU(data=np.arange(10))
        hdu.writeto(self.temp(filename))

        with fits.open(self.temp(filename)) as hdul:
            assert np.all(hdul[0].data == hdu.data)

    def test_writeto_full_disk(self, monkeypatch):
        """
        Test that it gives a readable error when trying to write an hdulist
        to a full disk.
        """

        def _writeto(self, array):
            raise OSError("Fake error raised when writing file.")

        def get_free_space_in_dir(path):
            return 0

        with pytest.raises(OSError) as exc:
            monkeypatch.setattr(fits.hdu.base._BaseHDU, "_writeto", _writeto)
            monkeypatch.setattr(data, "get_free_space_in_dir", get_free_space_in_dir)

            n = np.arange(0, 1000, dtype="int64")
            hdu = fits.PrimaryHDU(n)
            hdulist = fits.HDUList(hdu)
            filename = self.temp("test.fits")

            with open(filename, mode="wb") as fileobj:
                hdulist.writeto(fileobj)

        assert (
            "Not enough space on disk: requested 8000, available 0. "
            "Fake error raised when writing file." == exc.value.args[0]
        )

    def test_flush_full_disk(self, monkeypatch):
        """
        Test that it gives a readable error when trying to update an hdulist
        to a full disk.
        """
        filename = self.temp("test.fits")
        hdul = [fits.PrimaryHDU(), fits.ImageHDU()]
        hdul = fits.HDUList(hdul)
        hdul[0].data = np.arange(0, 1000, dtype="int64")
        hdul.writeto(filename)

        def _writedata(self, fileobj):
            raise OSError("Fake error raised when writing file.")

        def get_free_space_in_dir(path):
            return 0

        monkeypatch.setattr(fits.hdu.base._BaseHDU, "_writedata", _writedata)
        monkeypatch.setattr(data, "get_free_space_in_dir", get_free_space_in_dir)

        with pytest.raises(OSError) as exc:
            with fits.open(filename, mode="update") as hdul:
                hdul[0].data = np.arange(0, 1000, dtype="int64")
                hdul.insert(1, fits.ImageHDU())
                hdul.flush()

        assert (
            "Not enough space on disk: requested 8000, available 0. "
            "Fake error raised when writing file." == exc.value.args[0]
        )

    def _test_write_string_bytes_io(self, fileobj):
        """
        Implemented for both test_write_stringio and test_write_bytesio.
        """

        with fits.open(self.data("test0.fits")) as hdul:
            hdul.writeto(fileobj)
            hdul2 = fits.HDUList.fromstring(fileobj.getvalue())
            assert FITSDiff(hdul, hdul2).identical

    def _make_gzip_file(self, filename="test0.fits.gz"):
        gzfile = self.temp(filename)
        with open(self.data("test0.fits"), "rb") as f:
            gz = gzip.open(gzfile, "wb")
            gz.write(f.read())
            gz.close()

        return gzfile

    def test_write_overwrite(self, home_is_temp):
        filename = self.temp("test_overwrite.fits")
        hdu = fits.PrimaryHDU(data=np.arange(10))
        hdu.writeto(filename)
        with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
            hdu.writeto(filename)
        hdu.writeto(filename, overwrite=True)

    def _make_zip_file(self, mode="copyonwrite", filename="test0.fits.zip"):
        zfile = zipfile.ZipFile(self.temp(filename), "w")
        zfile.write(self.data("test0.fits"))
        zfile.close()

        return zfile.filename

    def _make_bzip2_file(self, filename="test0.fits.bz2"):
        bzfile = self.temp(filename)
        with open(self.data("test0.fits"), "rb") as f:
            bz = bz2.BZ2File(bzfile, "w")
            bz.write(f.read())
            bz.close()

        return bzfile

    def test_simulateonly(self):
        """Write to None simulates writing."""

        with fits.open(self.data("test0.fits")) as hdul:
            hdul.writeto(None)
            hdul[0].writeto(None)
            hdul[0].header.tofile(None)

    def test_bintablehdu_zero_bytes(self):
        """Make sure we don't have any zero-byte writes in BinTableHDU"""

        bright = np.rec.array(
            [
                (1, "Sirius", -1.45, "A1V"),
                (2, "Canopus", -0.73, "F0Ib"),
                (3, "Rigil Kent", -0.1, "G2V"),
            ],
            formats="int16,a20,float32,a10",
            names="order,name,mag,Sp",
        )

        hdu_non_zero = fits.BinTableHDU(bright)
        # use safeio, a special file handler meant to fail on zero-byte writes
        fh = safeio.CatchZeroByteWriter(open(self.temp("bright.fits"), mode="wb"))
        hdu_non_zero.writeto(fh)
        fh.close()

    def test_primaryhdu_zero_bytes(self):
        """
        Make sure we don't have any zero-byte writes from an ImageHDU
        (or other) of `size % BLOCK_SIZE == 0`
        """

        hdu_img_2880 = fits.PrimaryHDU(data=np.arange(720, dtype="i4"))
        # use safeio, a special file handler meant to fail on zero-byte writes
        fh = safeio.CatchZeroByteWriter(open(self.temp("image.fits"), mode="wb"))
        hdu_img_2880.writeto(fh)
        fh.close()


class TestStreamingFunctions(FitsTestCase):
    """Test functionality of the StreamingHDU class."""

    def test_streaming_hdu(self, home_is_temp):
        shdu = self._make_streaming_hdu(self.temp("new.fits"))
        assert isinstance(shdu.size, int)
        assert shdu.size == 100

        arr = np.arange(25, dtype=np.int32).reshape((5, 5))
        shdu.write(arr)
        assert shdu.writecomplete
        shdu.close()

        with fits.open(self.temp("new.fits")) as hdul:
            assert len(hdul) == 1
            assert (hdul[0].data == arr).all()

    def test_streaming_hdu_file_wrong_mode(self):
        """
        Test that streaming an HDU to a file opened in the wrong mode fails as
        expected.
        """
        with pytest.raises(ValueError):
            with open(self.temp("new.fits"), "wb") as f:
                header = fits.Header()
                fits.StreamingHDU(f, header)

    def test_streaming_hdu_write_file(self):
        """Test streaming an HDU to an open file object."""

        arr = np.zeros((5, 5), dtype=np.int32)
        with open(self.temp("new.fits"), "ab+") as f:
            shdu = self._make_streaming_hdu(f)
            shdu.write(arr)
            assert shdu.writecomplete
            assert shdu.size == 100
        with fits.open(self.temp("new.fits")) as hdul:
            assert len(hdul) == 1
            assert (hdul[0].data == arr).all()

    def test_streaming_hdu_write_file_like(self):
        """Test streaming an HDU to an open file-like object."""

        arr = np.zeros((5, 5), dtype=np.int32)
        # The file-like object underlying a StreamingHDU must be in binary mode
        sf = io.BytesIO()
        shdu = self._make_streaming_hdu(sf)
        shdu.write(arr)
        assert shdu.writecomplete
        assert shdu.size == 100

        sf.seek(0)
        hdul = fits.open(sf)
        assert len(hdul) == 1
        assert (hdul[0].data == arr).all()

    def test_streaming_hdu_append_extension(self):
        arr = np.zeros((5, 5), dtype=np.int32)
        with open(self.temp("new.fits"), "ab+") as f:
            shdu = self._make_streaming_hdu(f)
            shdu.write(arr)
        # Doing this again should update the file with an extension
        with open(self.temp("new.fits"), "ab+") as f:
            shdu = self._make_streaming_hdu(f)
            shdu.write(arr)

    def test_fix_invalid_extname(self, capsys):
        phdu = fits.PrimaryHDU()
        ihdu = fits.ImageHDU()
        ihdu.header["EXTNAME"] = 12345678
        hdul = fits.HDUList([phdu, ihdu])
        filename = self.temp("temp.fits")

        pytest.raises(
            fits.VerifyError, hdul.writeto, filename, output_verify="exception"
        )
        with pytest.warns(
            fits.verify.VerifyWarning, match=r"Verification reported errors"
        ):
            hdul.writeto(filename, output_verify="fix")
        with fits.open(filename):
            assert hdul[1].name == "12345678"
            assert hdul[1].header["EXTNAME"] == "12345678"

        hdul.close()

    def _make_streaming_hdu(self, fileobj):
        hd = fits.Header()
        hd["SIMPLE"] = (True, "conforms to FITS standard")
        hd["BITPIX"] = (32, "array data type")
        hd["NAXIS"] = (2, "number of array dimensions")
        hd["NAXIS1"] = 5
        hd["NAXIS2"] = 5
        hd["EXTEND"] = True
        return fits.StreamingHDU(fileobj, hd)

    def test_blank_ignore(self):
        with fits.open(self.data("blank.fits"), ignore_blank=True) as f:
            assert f[0].data.flat[0] == 2

    def test_error_if_memmap_impossible(self):
        pth = self.data("blank.fits")
        with fits.open(pth, memmap=True) as hdul:
            with pytest.raises(ValueError):
                hdul[0].data

        # However, it should not fail if do_not_scale_image_data was used:
        # See https://github.com/astropy/astropy/issues/3766
        with fits.open(pth, memmap=True, do_not_scale_image_data=True) as hdul:
            hdul[0].data  # Just make sure it doesn't crash
