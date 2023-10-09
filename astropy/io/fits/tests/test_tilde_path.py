# Licensed under a 3-clause BSD style license - see PYFITS.rst

import os

import numpy as np
import pytest

from astropy.io import fits
from astropy.utils.misc import _NOT_OVERWRITING_MSG_MATCH

from .conftest import FitsTestCase


class TestTildePaths(FitsTestCase):
    """
    Exercises a few functions, just to ensure they run with tilde paths (i.e.
    paths like '~/filename.fits'). Also exercises a few subclasses. Most of the
    rest of the testing of tilde path handling is done by adding `home_is_data`
    and `home_is_temp` fixtures (defined and explained in __init__.py) to
    appropriate test cases, so that they are run both with and without tilde
    paths.
    """

    def test_fits_info(self, home_is_data):
        fits.info(self.data("tb.fits"), output=False)

    def test_fits_printdiff(self, home_is_data):
        fits.printdiff(self.data("test0.fits"), self.data("tb.fits"))

    def test_fits_get_data(self, home_is_data):
        fits.getdata(self.data("test0.fits"))
        fits.getdata(self.data("test0.fits"), header=True)

    def test_fits_get_header(self, home_is_data):
        fits.getheader(self.data("test0.fits"))
        fits.getheader(self.data("tb.fits"), ext=1)

    def test_fits_get_set_del_val(self, home_is_temp):
        self.copy_file("test0.fits")
        filename = self.temp("test0.fits")
        assert fits.getval(filename, "shutter") == "B"

        fits.setval(filename, "shutter", value="C")
        assert fits.getval(filename, "shutter") == "C"

        with pytest.raises(KeyError):
            fits.getval(filename, "missing")

        fits.setval(filename, "missing", value="C")
        assert fits.getval(filename, "missing") == "C"

        fits.delval(filename, "missing")
        with pytest.raises(KeyError):
            fits.getval(filename, "missing")

    def test_header_formatter(self, home_is_data):
        from astropy.io.fits.scripts import fitsheader

        hf = fitsheader.HeaderFormatter(self.data("zerowidth.fits"))
        hf.close()

        thf = fitsheader.TableHeaderFormatter(self.data("zerowidth.fits"))
        thf.close()

    def test_BinTableHDU_dump_load(self, home_is_temp):
        bright = np.rec.array(
            [
                (1, "Serius", -1.45, "A1V"),
                (2, "Canopys", -0.73, "F0Ib"),
                (3, "Rigil Kent", -0.1, "G2V"),
            ],
            formats="int16,S20,float32,S10",
            names="order,name,mag,Sp",
        )
        hdu = fits.BinTableHDU(bright)

        hdu.dump(self.temp("data.txt"), self.temp("cdfile.txt"), self.temp("hfile.txt"))
        with pytest.raises(OSError, match="already exists"):
            hdu.dump(
                self.temp("data.txt"),
                self.temp("cdfile.txt"),
                self.temp("hfile.txt"),
                overwrite=False,
            )
        hdu.dump(
            self.temp("data.txt"),
            self.temp("cdfile.txt"),
            self.temp("hfile.txt"),
            overwrite=True,
        )

        fits.BinTableHDU.load(
            self.temp("data.txt"), self.temp("cdfile.txt"), self.temp("hfile.txt")
        )

        self.copy_file("tb.fits")
        with fits.open(self.temp("tb.fits")) as hdul:
            hdu = hdul[1]
            hdu.dump()
        assert os.path.exists(os.path.expanduser(self.temp("tb.txt")))

    def test_BinTableHDU_writeto(self, home_is_temp):
        bright = np.rec.array(
            [
                (1, "Serius", -1.45, "A1V"),
                (2, "Canopys", -0.73, "F0Ib"),
                (3, "Rigil Kent", -0.1, "G2V"),
            ],
            formats="int16,S20,float32,S10",
            names="order,name,mag,Sp",
        )
        hdu = fits.BinTableHDU(bright)

        hdu.writeto(self.temp("table.fits"))
        with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
            hdu.writeto(self.temp("table.fits"), overwrite=False)
        hdu.writeto(self.temp("table.fits"), overwrite=True)

    def test_TableHDU_writeto(self, home_is_temp):
        bright = np.rec.array(
            [
                (1, "Serius", -1.45, "A1V"),
                (2, "Canopys", -0.73, "F0Ib"),
                (3, "Rigil Kent", -0.1, "G2V"),
            ],
            formats="int16,S20,float32,S10",
            names="order,name,mag,Sp",
        )
        hdu = fits.TableHDU.from_columns(bright, nrows=2)

        hdu.writeto(self.temp("table.fits"))
        with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
            hdu.writeto(self.temp("table.fits"), overwrite=False)
        hdu.writeto(self.temp("table.fits"), overwrite=True)

    def fits_tabledump(self, home_is_temp):
        self.copy_file("tb.fits")

        fits.tabledump(
            self.temp("tb.fits"),
            self.temp("data.txt"),
            self.temp("cdfile.txt"),
            self.temp("hfile.txt"),
        )
        with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
            fits.tabledump(
                self.temp("tb.fits"),
                self.temp("data.txt"),
                self.temp("cdfile.txt"),
                self.temp("hfile.txt"),
                overwrite=False,
            )
        fits.tabledump(
            self.temp("tb.fits"),
            self.temp("data.txt"),
            self.temp("cdfile.txt"),
            self.temp("hfile.txt"),
            overwrite=True,
        )

        fits.tableload(
            self.temp("data.txt"), self.temp("cdfile.txt"), self.temp("hfile.txt")
        )

    def test_ImageHDU_writeto(self, home_is_temp):
        hdu = fits.ImageHDU(np.arange(100).reshape((10, 10)))
        hdu.writeto(self.temp("image.fits"))
        with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
            hdu.writeto(self.temp("image.fits"), overwrite=False)
        hdu.writeto(self.temp("image.fits"), overwrite=True)

    def test_CompImageHDU_writeto(self, home_is_temp):
        hdu = fits.CompImageHDU(np.arange(100).reshape((10, 10)).astype(np.int32))
        hdu.writeto(self.temp("image.fits"))
        with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
            hdu.writeto(self.temp("image.fits"), overwrite=False)
        hdu.writeto(self.temp("image.fits"), overwrite=True)
