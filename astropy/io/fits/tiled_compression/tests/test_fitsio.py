"""
This test file uses the https://github.com/esheldon/fitsio package to verify
our compression and decompression routines against the implementation in
cfitsio.

*Note*: The fitsio library is GPL licensed, therefore it could be interpreted
 that so is this test file. Given that this test file isn't imported anywhere
 else in the code this shouldn't cause us any issues. Please bear this in mind
 when editing this file.
"""
import numpy as np
import pytest

from astropy.io import fits

fitsio = pytest.importorskip("fitsio")


@pytest.fixture(
    scope="module",
    params=[
        (10,),
        (12, 12),
        (15, 15),
        (15, 15, 15),
        # >3D Data are not currently supported with astropy
        # (15, 15, 15, 15),
    ],
    ids=lambda x: f"shape: {x}",
)
def base_original_data(request, compression_type_dtype):
    if compression_type_dtype[0] == "HCOMPRESS_1" and len(request.param) != 2:
        pytest.xfail("HCOMPRESS is 2D only apparently")
    size = np.product(request.param)
    return np.arange(size).reshape(request.param)


@pytest.fixture(
    scope="module",
    ids=lambda x: f"tiles per axis: {x}",
    params=[3, 2, 1, (1, 3, -99), (1, 2, 3), (1, -99, -99)],
)
def tile_dims(request, base_original_data, compression_type_dtype):
    compression_type = compression_type_dtype[0]
    tile_scale_factor = np.array(request.param)
    if tile_scale_factor.ndim > 0:
        tile_scale_factor = tile_scale_factor[: base_original_data.ndim]

    if compression_type == "HCOMPRESS_1":
        # If HCOMPRESS then we can only have 2D tiles so all other dims must
        # have tile length of 1
        if tile_scale_factor.ndim > 2:
            tile_scale_factor[2:] = -99

        if np.count_nonzero(tile_scale_factor != -99) != 2:
            pytest.xfail("HCOMPRESS needs two non-unity tile dimensions.")

    tile_scale_factor = np.where(
        tile_scale_factor == -99, base_original_data.shape, tile_scale_factor
    )

    return tuple(
        np.array(np.array(base_original_data.shape) / tile_scale_factor, dtype=int)
    )


@pytest.fixture(scope="module")
def fitsio_compressed_file_path(
    tmp_path_factory,
    compression_type_dtype,
    base_original_data,
    tile_dims,
):
    compression_type, dtype = compression_type_dtype
    if base_original_data.ndim > 2 and "u1" in dtype:
        pytest.xfail("These don't work")
    tmp_path = tmp_path_factory.mktemp("fitsio")
    original_data = base_original_data.astype(dtype)

    filename = tmp_path / f"{compression_type}_{dtype}.fits"
    fits = fitsio.FITS(filename, "rw")
    fits.write(original_data, compress=compression_type, tile_dims=tile_dims)

    return filename


@pytest.fixture(scope="module")
def astropy_compressed_file_path(
    tmp_path_factory, compression_type_dtype, base_original_data
):
    compression_type, dtype = compression_type_dtype
    if base_original_data.ndim > 2 and "u1" in dtype:
        pytest.xfail("These don't work")
    original_data = base_original_data.astype(dtype)

    tmp_path = tmp_path_factory.mktemp("astropy")
    filename = tmp_path / f"{compression_type}_{dtype}.fits"
    hdu = fits.CompImageHDU(data=original_data, compression_type=compression_type)
    hdu.writeto(filename)

    return filename


def test_decompress(
    fitsio_compressed_file_path, base_original_data, compression_type_dtype
):
    compression_type, dtype = compression_type_dtype

    with fits.open(fitsio_compressed_file_path) as hdul:
        data = hdul[1].data

        assert hdul[1]._header["ZCMPTYPE"] == compression_type
        assert hdul[1].data.dtype.kind == np.dtype(dtype).kind
        assert hdul[1].data.dtype.itemsize == np.dtype(dtype).itemsize
        # assert hdul[1].data.dtype.byteorder == np.dtype(dtype).byteorder
    np.testing.assert_allclose(data, base_original_data)


def test_compress(
    astropy_compressed_file_path, base_original_data, compression_type_dtype
):
    compression_type, dtype = compression_type_dtype

    fits = fitsio.FITS(astropy_compressed_file_path, "r")
    header = fits[1].read_header()
    data = fits[1].read()

    assert header["ZCMPTYPE"] == compression_type
    assert data.dtype.kind == np.dtype(dtype).kind
    assert data.dtype.itemsize == np.dtype(dtype).itemsize
    # assert data.dtype.byteorder == np.dtype(dtype).byteorder
    np.testing.assert_allclose(data, base_original_data)
