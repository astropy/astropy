# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Tests of NDDataBase

import numpy as np

from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.nddata.nddata_base import NDDataBase


class MinimalSubclass(NDDataBase):
    def __init__(self):
        super().__init__()

    @property
    def data(self):
        return None

    @property
    def mask(self):
        return super().mask

    @property
    def unit(self):
        return super().unit

    @property
    def wcs(self):
        return super().wcs

    @property
    def meta(self):
        return super().meta

    @property
    def uncertainty(self):
        return super().uncertainty

    @property
    def psf(self):
        return super().psf


class MinimalSubclassNoPSF(NDDataBase):
    def __init__(self):
        super().__init__()

    @property
    def data(self):
        return None

    @property
    def mask(self):
        return super().mask

    @property
    def unit(self):
        return super().unit

    @property
    def wcs(self):
        return super().wcs

    @property
    def meta(self):
        return super().meta

    @property
    def uncertainty(self):
        return super().uncertainty


def test_nddata_base_subclass():
    a = MinimalSubclass()
    assert a.meta is None
    assert a.data is None
    assert a.mask is None
    assert a.unit is None
    assert a.wcs is None
    assert a.uncertainty is None
    assert a.psf is None


def test_omitting_psf_is_ok():
    # Make sure that psf does not need to be overridden when creating a subclass
    b = MinimalSubclassNoPSF()
    assert b.psf is None


def test_read_with_matching_units(tmp_path, caplog):
    """Test reading a FITS file with a matching unit in BUNIT."""
    # Create a temporary FITS file
    fits_file = tmp_path / "test.fits"
    hdu = fits.PrimaryHDU(data=np.zeros((10, 10)))
    hdu.header["BUNIT"] = "adu"
    hdu.writeto(fits_file, overwrite=True)

    # Read the FITS file and capture log messages
    with caplog.at_level("WARNING"):  # Set the logging level to WARNING
        ccd = CCDData.read(fits_file, unit="adu")

    # Assert no warnings were logged
    assert len(caplog.records) == 0

    # Assert that the CCDData object has the correct unit
    assert isinstance(ccd, CCDData)
    assert ccd.unit == u.adu
    assert ccd.data.shape == (10, 10)
