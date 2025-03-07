# Licensed under a 3-clause BSD style license - see PYFITS.rst

import warnings

import numpy as np
import pytest

from astropy.io import fits
from astropy.io.fits.hdu.base import _ValidHDU

from .conftest import FitsTestCase
from .test_table import comparerecords


class BaseChecksumTests(FitsTestCase):
    def setup_method(self):
        super().setup_method()
        self._oldfilters = warnings.filters[:]
        warnings.filterwarnings("error", message="Checksum verification failed")
        warnings.filterwarnings("error", message="Datasum verification failed")

        # Monkey-patch the _get_timestamp method so that the checksum
        # timestamps (and hence the checksum themselves) are always the same
        self._old_get_timestamp = _ValidHDU._get_timestamp
        _ValidHDU._get_timestamp = lambda self: "2013-12-20T13:36:10"

    def teardown_method(self):
        super().teardown_method()
        warnings.filters = self._oldfilters
        _ValidHDU._get_timestamp = self._old_get_timestamp


class TestChecksumFunctions(BaseChecksumTests):
    # All checksums have been verified against CFITSIO

    def test_sample_file(self):
        hdul = fits.open(self.data("checksum.fits"), checksum=True)
        assert hdul._read_all
        hdul.close()

    def test_image_create(self):
        n = np.arange(100, dtype=np.int64)
        hdu = fits.PrimaryHDU(n)
        hdu.writeto(self.temp("tmp.fits"), overwrite=True, checksum=True)
        with fits.open(self.temp("tmp.fits"), checksum=True) as hdul:
            assert (hdu.data == hdul[0].data).all()
            assert "CHECKSUM" in hdul[0].header
            assert "DATASUM" in hdul[0].header

            assert hdul[0].header["CHECKSUM"] == "ZHMkeGKjZGKjbGKj"
            assert hdul[0].header["DATASUM"] == "4950"

    def test_scaled_data(self):
        with fits.open(self.data("scale.fits")) as hdul:
            orig_data = hdul[0].data.copy()
            hdul[0].scale("int16", "old")
            hdul.writeto(self.temp("tmp.fits"), overwrite=True, checksum=True)
            with fits.open(self.temp("tmp.fits"), checksum=True) as hdul1:
                assert (hdul1[0].data == orig_data).all()
                assert "CHECKSUM" in hdul1[0].header
                assert hdul1[0].header["CHECKSUM"] == "cUmaeUjZcUjacUjW"
                assert "DATASUM" in hdul1[0].header
                assert hdul1[0].header["DATASUM"] == "1891563534"

    def test_scaled_data_auto_rescale(self):
        """
        Regression test for
        https://github.com/astropy/astropy/issues/3883#issuecomment-115122647

        Ensure that when scaled data is automatically rescaled on
        opening/writing a file that the checksum and datasum are computed for
        the rescaled array.
        """

        with fits.open(self.data("scale.fits")) as hdul:
            # Write out a copy of the data with the rescaling applied
            hdul.writeto(self.temp("rescaled.fits"))

        # Reopen the new file and save it back again with a checksum
        with fits.open(self.temp("rescaled.fits")) as hdul:
            hdul.writeto(self.temp("rescaled2.fits"), overwrite=True, checksum=True)

        # Now do like in the first writeto but use checksum immediately
        with fits.open(self.data("scale.fits")) as hdul:
            hdul.writeto(self.temp("rescaled3.fits"), checksum=True)

        # Also don't rescale the data but add a checksum
        with fits.open(self.data("scale.fits"), do_not_scale_image_data=True) as hdul:
            hdul.writeto(self.temp("scaled.fits"), checksum=True)

        # Must used nested with statements to support older Python versions
        # (but contextlib.nested is not available in newer Pythons :(
        with fits.open(self.temp("rescaled2.fits")) as hdul1:
            with fits.open(self.temp("rescaled3.fits")) as hdul2:
                with fits.open(self.temp("scaled.fits")) as hdul3:
                    hdr1 = hdul1[0].header
                    hdr2 = hdul2[0].header
                    hdr3 = hdul3[0].header
                    assert hdr1["DATASUM"] == hdr2["DATASUM"]
                    assert hdr1["CHECKSUM"] == hdr2["CHECKSUM"]
                    assert hdr1["DATASUM"] != hdr3["DATASUM"]
                    assert hdr1["CHECKSUM"] != hdr3["CHECKSUM"]

    def test_uint16_data(self):
        checksums = [
            ("aDcXaCcXaCcXaCcX", "0"),
            ("oYiGqXi9oXiEoXi9", "1746888714"),
            ("VhqQWZoQVfoQVZoQ", "0"),
            ("4cPp5aOn4aOn4aOn", "0"),
            ("8aCN8X9N8aAN8W9N", "1756785133"),
            ("UhqdUZnbUfnbUZnb", "0"),
            ("4cQJ5aN94aNG4aN9", "0"),
        ]
        with fits.open(self.data("o4sp040b0_raw.fits"), uint=True) as hdul:
            hdul.writeto(self.temp("tmp.fits"), overwrite=True, checksum=True)
            with fits.open(self.temp("tmp.fits"), uint=True, checksum=True) as hdul1:
                for idx, (hdu_a, hdu_b) in enumerate(zip(hdul, hdul1)):
                    if hdu_a.data is None or hdu_b.data is None:
                        assert hdu_a.data is hdu_b.data
                    else:
                        assert (hdu_a.data == hdu_b.data).all()

                    assert "CHECKSUM" in hdul[idx].header
                    assert hdul[idx].header["CHECKSUM"] == checksums[idx][0]
                    assert "DATASUM" in hdul[idx].header
                    assert hdul[idx].header["DATASUM"] == checksums[idx][1]

    def test_groups_hdu_data(self):
        imdata = np.arange(100.0)
        imdata.shape = (10, 1, 1, 2, 5)
        pdata1 = np.arange(10) + 0.1
        pdata2 = 42
        x = fits.hdu.groups.GroupData(
            imdata, parnames=["abc", "xyz"], pardata=[pdata1, pdata2], bitpix=-32
        )
        hdu = fits.GroupsHDU(x)
        hdu.writeto(self.temp("tmp.fits"), overwrite=True, checksum=True)
        with fits.open(self.temp("tmp.fits"), checksum=True) as hdul:
            assert comparerecords(hdul[0].data, hdu.data)
            assert "CHECKSUM" in hdul[0].header
            assert hdul[0].header["CHECKSUM"] == "3eDQAZDO4dDOAZDO"
            assert "DATASUM" in hdul[0].header
            assert hdul[0].header["DATASUM"] == "2797758084"

    def test_binary_table_data(self):
        a1 = np.array(["NGC1001", "NGC1002", "NGC1003"])
        a2 = np.array([11.1, 12.3, 15.2])
        col1 = fits.Column(name="target", format="20A", array=a1)
        col2 = fits.Column(name="V_mag", format="E", array=a2)
        cols = fits.ColDefs([col1, col2])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.writeto(self.temp("tmp.fits"), overwrite=True, checksum=True)
        with fits.open(self.temp("tmp.fits"), checksum=True) as hdul:
            assert comparerecords(tbhdu.data, hdul[1].data)
            assert "CHECKSUM" in hdul[0].header
            assert hdul[0].header["CHECKSUM"] == "D8iBD6ZAD6fAD6ZA"
            assert "DATASUM" in hdul[0].header
            assert hdul[0].header["DATASUM"] == "0"
            assert "CHECKSUM" in hdul[1].header
            assert hdul[1].header["CHECKSUM"] == "aD1Oa90MaC0Ma90M"
            assert "DATASUM" in hdul[1].header
            assert hdul[1].header["DATASUM"] == "1062205743"

    def test_variable_length_table_data(self):
        c1 = fits.Column(
            name="var",
            format="PJ()",
            array=np.array([[45.0, 56], np.array([11, 12, 13])], "O"),
        )
        c2 = fits.Column(name="xyz", format="2I", array=[[11, 3], [12, 4]])
        tbhdu = fits.BinTableHDU.from_columns([c1, c2])
        tbhdu.writeto(self.temp("tmp.fits"), overwrite=True, checksum=True)
        with fits.open(self.temp("tmp.fits"), checksum=True) as hdul:
            assert comparerecords(tbhdu.data, hdul[1].data)
            assert "CHECKSUM" in hdul[0].header
            assert hdul[0].header["CHECKSUM"] == "D8iBD6ZAD6fAD6ZA"
            assert "DATASUM" in hdul[0].header
            assert hdul[0].header["DATASUM"] == "0"
            assert "CHECKSUM" in hdul[1].header
            assert hdul[1].header["CHECKSUM"] == "YIGoaIEmZIEmaIEm"
            assert "DATASUM" in hdul[1].header
            assert hdul[1].header["DATASUM"] == "1507485"

    def test_variable_length_table_data2(self):
        """regression test for #12119"""

        time_data = [
            np.array([2021, 1, 5, 10, 5, 30], dtype=np.uint16),
            np.array([2021, 2, 19, 11, 19, 56], dtype=np.uint16),
            np.array([2021, 4, 21, 16, 10, 24], dtype=np.uint16),
            np.array([2021, 7, 22, 14, 42, 20], dtype=np.uint16),
        ]
        time_col = fits.Column(name="time", format="6I", array=time_data)

        version_data = ["5.45.70", "5.45.71", "5.45.102", "5.50.109"]
        version_col = fits.Column(name="Version", format="PA(8)", array=version_data)
        columns = [time_col, version_col]

        testfile = self.temp("tmp.fits")
        tbl = fits.BinTableHDU.from_columns(columns, name="DemoBinTable")
        hdul = fits.HDUList([fits.PrimaryHDU(), tbl])

        # here checksum is computed from in-memory data, which was producing
        # a wrong checksum and warnings when reading back the file
        hdul.writeto(testfile, checksum=True)

        testfile2 = self.temp("tmp2.fits")
        with fits.open(testfile, checksum=True) as hdul:
            checksum = hdul[1]._checksum
            datasum = hdul[1]._datasum
            # so write again the file but here data was not loaded so checksum
            # is computed directly from the file bytes, which was producing
            # a correct checksum. Below we compare both to make sure they are
            # consistent.
            hdul.writeto(testfile2, checksum=True)

        with fits.open(testfile2, checksum=True) as hdul:
            assert checksum == hdul[1]._checksum
            assert datasum == hdul[1]._datasum

    def test_ascii_table_data(self):
        a1 = np.array(["abc", "def"])
        r1 = np.array([11.0, 12.0])
        c1 = fits.Column(name="abc", format="A3", array=a1)
        # This column used to be E format, but the single-precision float lost
        # too much precision when scaling so it was changed to a D
        c2 = fits.Column(name="def", format="D", array=r1, bscale=2.3, bzero=0.6)
        c3 = fits.Column(name="t1", format="I", array=[91, 92, 93])
        x = fits.ColDefs([c1, c2, c3])
        hdu = fits.TableHDU.from_columns(x)
        hdu.writeto(self.temp("tmp.fits"), overwrite=True, checksum=True)
        with fits.open(self.temp("tmp.fits"), checksum=True) as hdul:
            assert comparerecords(hdu.data, hdul[1].data)
            assert "CHECKSUM" in hdul[0].header
            assert hdul[0].header["CHECKSUM"] == "D8iBD6ZAD6fAD6ZA"
            assert "DATASUM" in hdul[0].header
            assert hdul[0].header["DATASUM"] == "0"

            assert "CHECKSUM" in hdul[1].header
            assert hdul[1].header["CHECKSUM"] == "3rKFAoI94oICAoI9"
            assert "DATASUM" in hdul[1].header
            assert hdul[1].header["DATASUM"] == "1914653725"

    def test_open_with_no_keywords(self):
        hdul = fits.open(self.data("arange.fits"), checksum=True)
        hdul.close()

    def test_append(self):
        hdul = fits.open(self.data("tb.fits"))
        hdul.writeto(self.temp("tmp.fits"), overwrite=True)
        n = np.arange(100)
        fits.append(self.temp("tmp.fits"), n, checksum=True)
        hdul.close()
        hdul = fits.open(self.temp("tmp.fits"), checksum=True)
        assert hdul[0]._checksum is None
        hdul.close()

    def test_writeto_convenience(self):
        n = np.arange(100)
        fits.writeto(self.temp("tmp.fits"), n, overwrite=True, checksum=True)
        hdul = fits.open(self.temp("tmp.fits"), checksum=True)
        self._check_checksums(hdul[0])
        hdul.close()

    def test_hdu_writeto(self):
        n = np.arange(100, dtype="int16")
        hdu = fits.ImageHDU(n)
        hdu.writeto(self.temp("tmp.fits"), checksum=True)
        hdul = fits.open(self.temp("tmp.fits"), checksum=True)
        self._check_checksums(hdul[0])
        hdul.close()

    def test_hdu_writeto_existing(self):
        """
        Tests that when using writeto with checksum=True, a checksum and
        datasum are added to HDUs that did not previously have one.

        Regression test for https://github.com/spacetelescope/PyFITS/issues/8
        """

        with fits.open(self.data("tb.fits")) as hdul:
            hdul.writeto(self.temp("test.fits"), checksum=True)

        with fits.open(self.temp("test.fits")) as hdul:
            assert "CHECKSUM" in hdul[0].header
            # These checksums were verified against CFITSIO
            assert hdul[0].header["CHECKSUM"] == "7UgqATfo7TfoATfo"
            assert "DATASUM" in hdul[0].header
            assert hdul[0].header["DATASUM"] == "0"
            assert "CHECKSUM" in hdul[1].header
            assert hdul[1].header["CHECKSUM"] == "99daD8bX98baA8bU"
            assert "DATASUM" in hdul[1].header
            assert hdul[1].header["DATASUM"] == "1829680925"

    def test_datasum_only(self):
        n = np.arange(100, dtype="int16")
        hdu = fits.ImageHDU(n)
        hdu.writeto(self.temp("tmp.fits"), overwrite=True, checksum="datasum")
        with fits.open(self.temp("tmp.fits"), checksum=True) as hdul:
            if not (hasattr(hdul[0], "_datasum") and hdul[0]._datasum):
                pytest.fail("Missing DATASUM keyword")

            if not (hasattr(hdul[0], "_checksum") and not hdul[0]._checksum):
                pytest.fail("Non-empty CHECKSUM keyword")

    def test_open_update_mode_preserve_checksum(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/148 where
        checksums are being removed from headers when a file is opened in
        update mode, even though no changes were made to the file.
        """

        self.copy_file("checksum.fits")

        with fits.open(self.temp("checksum.fits")) as hdul:
            data = hdul[1].data.copy()

        hdul = fits.open(self.temp("checksum.fits"), mode="update")
        hdul.close()

        with fits.open(self.temp("checksum.fits")) as hdul:
            assert "CHECKSUM" in hdul[1].header
            assert "DATASUM" in hdul[1].header
            assert comparerecords(data, hdul[1].data)

    def test_open_update_mode_update_checksum(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/148, part
        2.  This ensures that if a file contains a checksum, the checksum is
        updated when changes are saved to the file, even if the file was opened
        with the default of checksum=False.

        An existing checksum and/or datasum are only stripped if the file is
        opened with checksum='remove'.
        """

        self.copy_file("checksum.fits")
        with fits.open(self.temp("checksum.fits")) as hdul:
            header = hdul[1].header.copy()
            data = hdul[1].data.copy()

        with fits.open(self.temp("checksum.fits"), mode="update") as hdul:
            hdul[1].header["FOO"] = "BAR"
            hdul[1].data[0]["TIME"] = 42

        with fits.open(self.temp("checksum.fits")) as hdul:
            header2 = hdul[1].header
            data2 = hdul[1].data
            assert header2[:-3] == header[:-2]
            assert "CHECKSUM" in header2
            assert "DATASUM" in header2
            assert header2["FOO"] == "BAR"
            assert (data2["TIME"][1:] == data["TIME"][1:]).all()
            assert data2["TIME"][0] == 42

        with fits.open(
            self.temp("checksum.fits"), mode="update", checksum="remove"
        ) as hdul:
            pass

        with fits.open(self.temp("checksum.fits")) as hdul:
            header2 = hdul[1].header
            data2 = hdul[1].data
            assert header2[:-1] == header[:-2]
            assert "CHECKSUM" not in header2
            assert "DATASUM" not in header2
            assert header2["FOO"] == "BAR"
            assert (data2["TIME"][1:] == data["TIME"][1:]).all()
            assert data2["TIME"][0] == 42

    def test_overwrite_invalid(self):
        """
        Tests that invalid checksum or datasum are overwritten when the file is
        saved.
        """

        reffile = self.temp("ref.fits")
        with fits.open(self.data("tb.fits")) as hdul:
            hdul.writeto(reffile, checksum=True)

        testfile = self.temp("test.fits")
        with fits.open(self.data("tb.fits")) as hdul:
            hdul[0].header["DATASUM"] = "1       "
            hdul[0].header["CHECKSUM"] = "8UgqATfo7TfoATfo"
            hdul[1].header["DATASUM"] = "2349680925"
            hdul[1].header["CHECKSUM"] = "11daD8bX98baA8bU"
            hdul.writeto(testfile)

        with fits.open(testfile) as hdul:
            hdul.writeto(self.temp("test2.fits"), checksum=True)

        with fits.open(self.temp("test2.fits")) as hdul:
            with fits.open(reffile) as ref:
                assert "CHECKSUM" in hdul[0].header
                # These checksums were verified against CFITSIO
                assert hdul[0].header["CHECKSUM"] == ref[0].header["CHECKSUM"]
                assert "DATASUM" in hdul[0].header
                assert hdul[0].header["DATASUM"] == "0"
                assert "CHECKSUM" in hdul[1].header
                assert hdul[1].header["CHECKSUM"] == ref[1].header["CHECKSUM"]
                assert "DATASUM" in hdul[1].header
                assert hdul[1].header["DATASUM"] == ref[1].header["DATASUM"]

    def _check_checksums(self, hdu):
        if not (hasattr(hdu, "_datasum") and hdu._datasum):
            pytest.fail("Missing DATASUM keyword")

        if not (hasattr(hdu, "_checksum") and hdu._checksum):
            pytest.fail("Missing CHECKSUM keyword")
