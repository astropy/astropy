from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from astropy.io import fits


@pytest.fixture
def canonical_data_base_path():
    return Path(__file__).parent / "data"


@pytest.fixture(
    params=(Path(__file__).parent / "data").glob("m13_*.fits"), ids=lambda x: x.name
)
def canonical_int_hdus(request):
    """
    This fixture provides 4 files downloaded from https://fits.gsfc.nasa.gov/registry/tilecompression.html

    Which are used as canonical tests of data not compressed by Astropy.
    """
    with fits.open(request.param) as hdul:
        yield hdul[1]


@pytest.fixture
def original_int_hdu(canonical_data_base_path):
    with fits.open(canonical_data_base_path / "m13.fits") as hdul:
        yield hdul[0]


# pytest-openfiles does not correctly check for open files when the files are
# opened in a fixture, so we skip the check here.
# https://github.com/astropy/pytest-openfiles/issues/32
@pytest.mark.openfiles_ignore
def test_canonical_data(original_int_hdu, canonical_int_hdus):
    assert_allclose(original_int_hdu.data, canonical_int_hdus.data)


def test_zblank_support(canonical_data_base_path, tmp_path):

    # This uses a test 12x12 image which contains a NaN value in the [1, 1]
    # pixel - it was compressed using fpack which automatically added a ZBLANK
    # header keyword

    reference = np.arange(144).reshape((12, 12)).astype(float)
    reference[1, 1] = np.nan

    with fits.open(canonical_data_base_path / "compressed_with_nan.fits") as hdul:
        assert_equal(np.round(hdul[1].data), reference)

    # Now generate a file ourselves and check that the output has the ZBLANK
    # keyword set automatically

    hdu = fits.CompImageHDU(data=reference, compression_type="RICE_1", tile_size=(6, 6))

    hdu.writeto(tmp_path / "test_zblank.fits")

    with fits.open(tmp_path / "test_zblank.fits") as hdul:
        assert "ZBLANK" in hdul[1].header
        assert_equal(np.round(hdul[1].data), reference)
