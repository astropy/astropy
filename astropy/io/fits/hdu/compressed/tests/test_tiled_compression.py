import hashlib
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from astropy.io import fits
from astropy.io.fits.hdu.compressed._codecs import PLIO1
from astropy.io.fits.hdu.compressed._compression import CfitsioException

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

    hdu = fits.CompImageHDU(
        data=reference, compression_type="RICE_1", tile_shape=(6, 6)
    )

    hdu.writeto(tmp_path / "test_zblank.fits")

    with fits.open(tmp_path / "test_zblank.fits") as hdul:
        assert "ZBLANK" in hdul[1].header
        assert_equal(np.round(hdul[1].data), reference)


@pytest.mark.parametrize(
    ("shape", "tile_shape"),
    (
        ([10, 10], [5, 5]),  # something for HCOMPRESS
        ([5, 5, 5], [5, 5, 5]),
        # ([5, 5, 5], [5, 5, 1]),  # something for HCOMPRESS
        ([10, 15, 20], [5, 5, 5]),
        ([10, 5, 12], [5, 5, 5]),
        # TODO: There's a stupid bit of code in CompImageHDU which stops this working.
        # ([2, 3, 4, 5], [1, 1, 2, 3]),
        ([2, 3, 4, 5], [5, 5, 1, 1]),
    ),
)
def test_roundtrip_high_D(
    numpy_rng, compression_type, compression_param, tmp_path, dtype, shape, tile_shape
):
    if compression_type == "HCOMPRESS_1" and (
        # We don't have at least a 2D image
        len(shape) < 2
        or
        # We don't have 2D tiles
        np.count_nonzero(np.array(tile_shape) != 1) != 2
        or
        # TODO: The following restrictions can be lifted with some extra work.
        # The tile is not the first two dimensions of the data
        tile_shape[0] == 1
        or tile_shape[1] == 1
        or
        # The tile dimensions not an integer multiple of the array dims
        np.count_nonzero(np.array(shape[:2]) % tile_shape[:2]) != 0
    ):
        pytest.xfail("HCOMPRESS requires 2D tiles.")
    random = numpy_rng.uniform(high=255, size=shape)
    # Set first value to be exactly zero as zero values require special treatment
    # for SUBTRACTIVE_DITHER_2
    random.ravel()[0] = 0.0
    original_data = random.astype(dtype)

    dtype_sanitizer = {
        ">": "big",
        "<": "little",
        "=": "native",
    }
    filename = (
        tmp_path / f"{compression_type}_{dtype[1:]}_{dtype_sanitizer[dtype[0]]}.fits"
    )

    param = fitsio_param_to_astropy_param(compression_param)
    hdu = fits.CompImageHDU(
        data=original_data,
        compression_type=compression_type,
        tile_shape=tile_shape,
        **param,
    )
    hdu.writeto(filename)

    atol = 0
    if compression_param.get("qmethod", None) is not None:
        # This is a horrific hack We are comparing quantized data to unquantized
        # data here, so there can be pretty large differences.  What this test
        # is really checking for is arrays which are *completely* different,
        # which would indicate the compression has not worked.
        atol = 17

    with fits.open(filename) as hdul:
        a = hdul[1].data
        np.testing.assert_allclose(original_data, hdul[1].data, atol=atol)


def test_plio_1_out_of_range():
    pc = PLIO1(tilesize=10)
    data = np.arange(-10, 0).astype(np.int32)

    with pytest.raises(ValueError):
        pc.encode(data)


def test_invalid_tile(tmp_path):
    # Regression test for a bug that caused a segmentation fault if the data
    # for a tile would have resulted in a warning

    m13_rice_path = Path(__file__).parent / "data" / "m13_rice.fits"

    # For this test, we will change the length of the first chunk so that it is
    # invalid. To do this, we just change the bytes directly so we first check
    # a checksum of the file to make sure it is as we expect to start with.

    with open(m13_rice_path, "rb") as f:
        content = f.read()

    assert hashlib.sha256(content).hexdigest()[:8] == "de6d2f69"

    # We change bytes 8640 to 8643 which are the length of the first tile:

    assert content[8640:8644] == b"\x00\x00\x00\x96"

    with open(tmp_path / "m13_corrupted.fits", "wb") as f:
        f.write(content[:8640])
        f.write(b"\x00\x00\x00\x95")
        f.write(content[8644:])

    with fits.open(tmp_path / "m13_corrupted.fits") as hdulist:
        # Access the data to make sure we decompress it
        with pytest.raises(
            CfitsioException,
            match="decompression error: hit end of compressed byte stream",
        ):
            hdulist[1].data.sum()
