# Tests related to writing dask arrays to FITS files in an efficient way

import numpy as np
import pytest

from astropy.io import fits
from astropy.io.fits import ImageHDU, PrimaryHDU

da = pytest.importorskip("dask.array")


@pytest.fixture
def dask_array_in_mem():
    return da.random.uniform(-1000, 1000, (1322, 755)).rechunk((59, 55))


def test_construct_image_hdu(dask_array_in_mem):
    hdu = ImageHDU(data=dask_array_in_mem)
    assert isinstance(hdu.data, da.Array)


def test_construct_hdulist(dask_array_in_mem):
    hdu = ImageHDU(data=dask_array_in_mem)
    hdulist = fits.HDUList([hdu])
    assert isinstance(hdulist[0].data, da.Array)


def test_save_primary_hdu(dask_array_in_mem, tmp_path):
    # Saving a Primary HDU directly

    filename = tmp_path / "test.fits"

    hdu = PrimaryHDU(data=dask_array_in_mem)
    hdu.writeto(filename)

    with fits.open(filename) as hdulist_new:
        assert isinstance(hdulist_new[0].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[0].data, dask_array_in_mem.compute())


def test_save_image_hdu(dask_array_in_mem, tmp_path):
    # Saving an image HDU directly

    filename = tmp_path / "test.fits"

    hdu = ImageHDU(data=dask_array_in_mem)
    hdu.writeto(filename)

    with fits.open(filename) as hdulist_new:
        assert isinstance(hdulist_new[1].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[1].data, dask_array_in_mem.compute())


def test_save_hdulist(dask_array_in_mem, tmp_path):
    # Saving an HDUList

    filename = tmp_path / "test.fits"

    hdu1 = PrimaryHDU(data=dask_array_in_mem)
    hdu2 = ImageHDU(data=np.random.random((128, 128)))
    hdu3 = ImageHDU(data=dask_array_in_mem * 2)
    hdulist = fits.HDUList([hdu1, hdu2, hdu3])
    assert isinstance(hdulist[0].data, da.Array)
    hdulist.writeto(filename)

    with fits.open(filename) as hdulist_new:
        assert isinstance(hdulist_new[0].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[0].data, dask_array_in_mem.compute())
        assert isinstance(hdulist_new[1].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[1].data, hdu2.data)
        assert isinstance(hdulist_new[2].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[2].data, dask_array_in_mem.compute() * 2)


def test_long_header(dask_array_in_mem, tmp_path):
    # Make sure things work correctly if there is a long header in the HDU.

    filename = tmp_path / "test.fits"

    # NOTE: we deliberately set up a long header here rather than add the
    # keys one by one to hdu.header as adding the header in one go used to
    # cause issues, so this acts as a regression test.
    header = fits.Header()
    for index in range(2048):
        header[f"KEY{index:x}"] = 0.0

    hdu = PrimaryHDU(data=dask_array_in_mem, header=header)
    hdu.writeto(filename)

    with fits.open(filename) as hdulist_new:
        assert len(hdulist_new[0].header) == 2053
        assert isinstance(hdulist_new[0].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[0].data, dask_array_in_mem.compute())


VALID_DTYPES = (">i2", "<i2", ">i4", "<i4", ">i8", "<i8", ">f4", "<f4", ">f8", "<f8")


@pytest.mark.parametrize("dtype", VALID_DTYPES)
def test_dtypes(dask_array_in_mem, tmp_path, dtype):
    filename = tmp_path / "test.fits"

    array = dask_array_in_mem.astype(dtype)

    hdu = PrimaryHDU(data=array)
    hdu.writeto(filename)

    with fits.open(filename) as hdulist_new:
        assert isinstance(hdulist_new[0].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[0].data, array.compute())


def test_scaled(dask_array_in_mem, tmp_path):
    filename = tmp_path / "test.fits"

    hdu = PrimaryHDU(data=dask_array_in_mem)
    hdu.scale("int32", bzero=-1000, bscale=1e-6)
    hdu.writeto(filename)

    with fits.open(filename) as hdulist_new:
        assert isinstance(hdulist_new[0].data, np.ndarray)
        np.testing.assert_allclose(
            hdulist_new[0].data, dask_array_in_mem.compute(), atol=1e-5
        )


def test_scaled_minmax(dask_array_in_mem, tmp_path):
    filename = tmp_path / "test.fits"

    hdu = PrimaryHDU(data=dask_array_in_mem)
    hdu.scale("int32", option="minmax")
    hdu.writeto(filename)

    with fits.open(filename) as hdulist_new:
        assert isinstance(hdulist_new[0].data, np.ndarray)
        np.testing.assert_allclose(
            hdulist_new[0].data, dask_array_in_mem.compute(), atol=1e-5
        )


def test_append(dask_array_in_mem, tmp_path):
    # Test append mode

    filename = tmp_path / "test.fits"

    fits.append(filename, dask_array_in_mem)
    fits.append(filename, np.arange(10))

    with fits.open(filename) as hdulist_new:
        assert isinstance(hdulist_new[0].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[0].data, dask_array_in_mem.compute())
        assert isinstance(hdulist_new[1].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[1].data, np.arange(10))


# @pytest.mark.parametrize('mode', ['rb+', 'ab', 'ab+', 'wb', 'wb+'])
@pytest.mark.parametrize("mode", ["wb", "wb+"])
def test_file_handle(mode, dask_array_in_mem, tmp_path):
    filename = tmp_path / "test.fits"
    hdu1 = PrimaryHDU(data=dask_array_in_mem)
    hdu2 = ImageHDU(data=np.arange(10))
    hdulist = fits.HDUList([hdu1, hdu2])

    with filename.open(mode=mode) as fp:
        hdulist.writeto(fp)

    with fits.open(filename) as hdulist_new:
        assert isinstance(hdulist_new[0].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[0].data, dask_array_in_mem.compute())
        assert isinstance(hdulist_new[1].data, np.ndarray)
        np.testing.assert_allclose(hdulist_new[1].data, np.arange(10))
