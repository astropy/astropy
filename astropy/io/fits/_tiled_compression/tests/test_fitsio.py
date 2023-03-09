"""
This test file uses the https://github.com/esheldon/fitsio package to verify
our compression and decompression routines against the implementation in
cfitsio.

*Note*: The fitsio library is GPL licensed, therefore it could be interpreted
 that so is this test file. Given that this test file isn't imported anywhere
 else in the code this shouldn't cause us any issues. Please bear this in mind
 when editing this file.
"""
import os

import numpy as np
import pytest

from astropy.io import fits

from .conftest import _expand, fitsio_param_to_astropy_param

# This is so that tox can force this file to be run, and not be silently
# skipped on CI, but in all other test runs it's skipped if fitsio isn't present.
if "ASTROPY_ALWAYS_TEST_FITSIO" in os.environ:
    import fitsio
else:
    fitsio = pytest.importorskip("fitsio")


@pytest.fixture(
    scope="module",
    params=_expand(
        [((10,),), ((5,), (1,), (3,))],
        [((12, 12),), ((1, 12), (4, 5), (6, 6), None)],
        [((15, 15),), ((1, 15), (5, 1), (5, 5))],
        [
            ((15, 15, 15),),
            ((5, 5, 1), (5, 7, 1), (1, 5, 4), (1, 1, 15), (15, 1, 5)),
        ],
        # Test the situation where the tile shape is passed larger than the
        # array shape
        [
            (
                (4, 4, 5),
                (5, 5, 5),
            ),
            (
                (5, 5, 1),
                None,
            ),
        ],
        # Test shapes which caused errors
        # This one we can't test here as it causes segfaults in cfitsio
        # It is tested in test_roundtrip_high_D though.
        # [
        #     ((3, 4, 5),),
        #     ((1, 2, 3),),
        # ],
        # >3D Data are not currently supported by cfitsio
    ),
    ids=lambda x: f"shape: {x[0]} tile_dims: {x[1]}",
)
def array_shapes_tile_dims(request, compression_type):
    shape, tile_dims = request.param
    # H_COMPRESS needs >=2D data and always 2D tiles
    if compression_type == "HCOMPRESS_1":
        if (
            # We don't have at least a 2D image
            len(shape) < 2
            or
            # We don't have 2D tiles
            np.count_nonzero(np.array(tile_dims) != 1) != 2
            or
            # TODO: The following restrictions can be lifted with some extra work.
            # The tile is not the first two dimensions of the data
            tile_dims[0] == 1
            or tile_dims[1] == 1
            or
            # The tile dimensions not an integer multiple of the array dims
            np.count_nonzero(np.array(shape[:2]) % tile_dims[:2]) != 0
        ):
            pytest.xfail(
                "HCOMPRESS requires 2D tiles, from the first two"
                "dimensions, and an integer number of tiles along the first two"
                "axes."
            )
    return shape, tile_dims


@pytest.fixture(scope="module")
def tile_dims(array_shapes_tile_dims):
    return array_shapes_tile_dims[1]


@pytest.fixture(scope="module")
def data_shape(array_shapes_tile_dims):
    return array_shapes_tile_dims[0]


@pytest.fixture(scope="module")
def base_original_data(data_shape, dtype, numpy_rng, compression_type):
    random = numpy_rng.uniform(high=255, size=data_shape)
    # Set first value to be exactly zero as zero values require special treatment
    # for SUBTRACTIVE_DITHER_2
    random.ravel()[0] = 0.0
    # There seems to be a bug with the fitsio library where HCOMPRESS doesn't
    # work with int16 random data, so use a bit for structured test data.
    if compression_type.startswith("HCOMPRESS") and "i2" in dtype or "u1" in dtype:
        random = np.arange(np.prod(data_shape)).reshape(data_shape)
    return random.astype(dtype)


@pytest.fixture(scope="module")
def fitsio_compressed_file_path(
    tmp_path_factory,
    comp_param_dtype,
    base_original_data,
    data_shape,  # For debugging
    tile_dims,
):
    compression_type, param, dtype = comp_param_dtype

    if (
        base_original_data.ndim > 2
        and "u1" in dtype
        and compression_type == "HCOMPRESS_1"
    ):
        pytest.xfail("fitsio won't write these")

    if compression_type == "PLIO_1" and "f" in dtype:
        # fitsio fails with a compression error
        pytest.xfail("fitsio fails to write these")

    if compression_type == "NOCOMPRESS":
        pytest.xfail("fitsio does not support NOCOMPRESS")

    if (
        compression_type == "HCOMPRESS_1"
        and "f" in dtype
        and param.get("qmethod", None) == 2
    ):
        # fitsio writes these files with very large/incorrect zzero values, whereas
        # qmethod == 1 works (and the two methods should be identical except for the
        # treatment of zeros)
        pytest.xfail("fitsio writes these files with very large/incorrect zzero values")

    tmp_path = tmp_path_factory.mktemp("fitsio")
    original_data = base_original_data.astype(dtype)

    filename = tmp_path / f"{compression_type}_{dtype}.fits"
    fits = fitsio.FITS(filename, "rw")
    fits.write(original_data, compress=compression_type, tile_dims=tile_dims, **param)

    return filename


@pytest.fixture(scope="module")
def astropy_compressed_file_path(
    comp_param_dtype,
    tmp_path_factory,
    base_original_data,
    data_shape,  # For debugging
    tile_dims,
):
    compression_type, param, dtype = comp_param_dtype
    original_data = base_original_data.astype(dtype)

    tmp_path = tmp_path_factory.mktemp("astropy")
    filename = tmp_path / f"{compression_type}_{dtype}.fits"

    param = fitsio_param_to_astropy_param(param)
    hdu = fits.CompImageHDU(
        data=original_data,
        compression_type=compression_type,
        tile_shape=None if tile_dims is None else tile_dims[::-1],
        **param,
    )
    hdu.writeto(filename)

    return filename


def test_decompress(
    fitsio_compressed_file_path,
    comp_param_dtype,
):
    compression_type, param, dtype = comp_param_dtype

    with fits.open(fitsio_compressed_file_path) as hdul:
        data = hdul[1].data

        assert hdul[1]._header["ZCMPTYPE"].replace("ONE", "1") == compression_type
        assert hdul[1].data.dtype.kind == np.dtype(dtype).kind
        assert hdul[1].data.dtype.itemsize == np.dtype(dtype).itemsize

    # The data might not always match the original data exactly in the case of
    # lossy compression so instead of comparing the array read by astropy to the
    # original data, we compare it to the data read in by fitsio (as those
    # should match)

    fts = fitsio.FITS(fitsio_compressed_file_path)
    data2 = fts[1].read()
    np.testing.assert_allclose(data, data2)

    # The first value should be exactly equal to zero when using SUBTRACTIVE_DITHER_2
    if param.get("qmethod", None) == 2:
        assert data.ravel()[0] == 0.0


def test_compress(
    astropy_compressed_file_path,
    compression_type,
    dtype,
):
    if compression_type == "NOCOMPRESS":
        pytest.xfail("fitsio does not support NOCOMPRESS")

    fts = fitsio.FITS(astropy_compressed_file_path, "r")
    header = fts[1].read_header()
    data = fts[1].read()

    assert header["ZCMPTYPE"] == compression_type
    assert data.dtype.kind == np.dtype(dtype).kind
    assert data.dtype.itemsize == np.dtype(dtype).itemsize

    # The data might not always match the original data exactly in the case of
    # lossy compression so instead of comparing the array read by fitsio to the
    # original data, we compare it to the data read in by astropy (as those
    # should match)

    with fits.open(astropy_compressed_file_path) as hdul:
        np.testing.assert_allclose(data, hdul[1].data)
