from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from astropy.io import fits

from .conftest import fitsio_param_to_astropy_param


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


@pytest.mark.parametrize(
    ("shape", "tile_dim"),
    (
        ([5, 5, 5], [5, 5, 5]),
        ([5, 5, 5], [5, 5, 1]),  # something for HCOMPRESS
        ([10, 15, 20], [5, 5, 5]),
        ([10, 5, 12], [5, 5, 5]),
        # TODO: There's a stupid bit of code in CompImageHDU which stops this working.
        # ([2, 3, 4, 5], [1, 1, 2, 3]),
        ([2, 3, 4, 5], [5, 5, 1, 1]),
    ),
)
def test_roundtrip_high_D(
    numpy_rng, compression_type, compression_param, tmp_path, dtype, shape, tile_dim
):
    if compression_type == "HCOMPRESS_1" and (
        len(shape) < 2 or np.count_nonzero(np.array(tile_dim) != 1) != 2
    ):
        pytest.xfail("HCOMPRESS requires 2D tiles.")
    random = numpy_rng.uniform(high=255, size=shape)
    # Set first value to be exactly zero as zero values require special treatment
    # for SUBTRACTIVE_DITHER_2
    random.ravel()[0] = 0.0
    original_data = random.astype(dtype)

    filename = tmp_path / f"{compression_type}_{dtype}.fits"

    param = fitsio_param_to_astropy_param(compression_param)
    hdu = fits.CompImageHDU(
        data=original_data,
        compression_type=compression_type,
        tile_size=tile_dim,
        **param,
    )
    hdu.writeto(filename)

    with fits.open(filename) as hdul:
        a = hdul[1].data
    #     np.testing.assert_allclose(original_data, hdul[1].data, atol=1)
