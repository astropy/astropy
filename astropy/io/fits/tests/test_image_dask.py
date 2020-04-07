# Tests related to writing dask arrays to FITS files in an efficient way

import pytest

import numpy as np
from astropy.io import fits
from astropy.io.fits import ImageHDU, PrimaryHDU

da = pytest.importorskip("dask.array")


@pytest.fixture
def dask_array_in_mem():
    return da.from_array(np.random.random((1322, 755))).rechunk((59, 55))


def test_construct_image_hdu(dask_array_in_mem):
    hdu = ImageHDU(data=dask_array_in_mem)
    assert isinstance(hdu.data, da.Array)


def test_construct_hdulist(dask_array_in_mem):
    hdu = ImageHDU(data=dask_array_in_mem)
    hdulist = fits.HDUList([hdu])
    assert isinstance(hdulist[0].data, da.Array)


def test_save_primary_hdu(dask_array_in_mem, tmp_path):

    # Saving a Primary HDU directly

    filename = tmp_path / 'test.fits'

    hdu = PrimaryHDU(data=dask_array_in_mem)
    hdu.writeto(filename)

    with fits.open(filename) as hdulist_new:
        assert isinstance(hdulist_new[0].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[0].data, dask_array_in_mem.compute())


def test_save_image_hdu(dask_array_in_mem, tmp_path):

    # Saving an image HDU directly

    filename = tmp_path / 'test.fits'

    hdu = ImageHDU(data=dask_array_in_mem)
    hdu.writeto(filename)

    with fits.open(filename) as hdulist_new:
        assert isinstance(hdulist_new[1].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[1].data, dask_array_in_mem.compute())


def test_save_hdulist(dask_array_in_mem, tmp_path):

    # Saving an HDUList

    filename = tmp_path / 'test.fits'

    hdu = PrimaryHDU(data=dask_array_in_mem)
    hdulist = fits.HDUList([hdu])
    assert isinstance(hdulist[0].data, da.Array)
    hdulist.writeto(filename)

    with fits.open(filename) as hdulist_new:
        assert isinstance(hdulist_new[0].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[0].data, dask_array_in_mem.compute())


def test_long_header(dask_array_in_mem, tmp_path):

    # Make sure things work correctly if there is a long header in the HDU.

    filename = tmp_path / 'test.fits'

    # NOTE: we deliberately set up a long header here rather than add the
    # keys one by one to hdu.header as adding the header in one go used to
    # cause issues, so this acts as a regression test.
    header = fits.Header()
    for index in range(2048):
        header[f'KEY{index:x}'] = 0.

    hdu = PrimaryHDU(data=dask_array_in_mem, header=header)
    hdu.writeto(filename)

    with fits.open(filename) as hdulist_new:
        assert len(hdulist_new[0].header) == 2053
        assert isinstance(hdulist_new[0].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[0].data, dask_array_in_mem.compute())
