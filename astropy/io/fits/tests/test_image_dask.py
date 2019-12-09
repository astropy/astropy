import pytest

pytest.importorskip("dask.array")

import numpy as np
import dask.array
dask.config.set(scheduler='synchronous')

from astropy.io import fits
from astropy.io.fits import ImageHDU


@pytest.fixture
def dask_array_in_mem():
    return dask.array.from_array(np.random.random((10, 10)))


def test_construct_image_hdu(dask_array_in_mem):
    hdu = ImageHDU(data=dask_array_in_mem)
    assert isinstance(hdu.data, dask.array.Array)


def test_construct_hdul(dask_array_in_mem):
    hdu = ImageHDU(data=dask_array_in_mem)
    hdul = fits.HDUList([hdu])
    assert isinstance(hdul[0].data, dask.array.Array)


def test_save_hdul(dask_array_in_mem, tmp_path):
    fname = "/tmp/test_save_hdul.fits"
    hdu = ImageHDU(data=dask_array_in_mem)
    hdul = fits.HDUList([hdu])
    assert isinstance(hdul[0].data, dask.array.Array)

    with open(fname, "wb+") as fobj:
        hdu.writeto(fobj)


    with fits.open(fname) as newhdul:
        assert isinstance(newhdul[1].data, np.ndarray)
        assert np.allclose(newhdul[1].data, dask_array_in_mem.compute())
