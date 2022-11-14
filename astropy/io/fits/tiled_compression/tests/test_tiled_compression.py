from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_equal

from astropy.io import fits
from astropy.io.fits.tiled_compression import (
    compress_hdu,
    compress_tile,
    decompress_hdu,
    decompress_tile,
)
from astropy.io.fits.tiled_compression.tiled_compression import (
    _buffer_to_array,
    _header_to_settings,
)

COMPRESSION_TYPES = [
    "GZIP_1",
    "GZIP_2",
    "RICE_1",
    "HCOMPRESS_1",
    "PLIO_1",
]

parameters = []
for compression_type in COMPRESSION_TYPES:
    # io.fits doesn't seem able to compress 64-bit data, even though e.g. GZIP_?
    # and HCOMPRESS_1 should be able to handle it.
    for itemsize in [1, 2, 4]:
        for endian in ["<", ">"]:
            format = "u" if itemsize == 1 else "i"
            parameters.append((compression_type, f"{endian}{format}{itemsize}"))


@pytest.mark.parametrize(("compression_type", "dtype"), parameters)
def test_basic(tmp_path, compression_type, dtype):

    # Generate compressed file dynamically

    original_data = np.arange(144).reshape((12, 12)).astype(dtype)

    header = fits.Header()

    hdu = fits.CompImageHDU(
        original_data, header, compression_type=compression_type, tile_size=(4, 4)
    )

    hdu.writeto(tmp_path / "test.fits")

    # Load in raw compressed data
    hdulist = fits.open(tmp_path / "test.fits", disable_image_compression=True)

    settings = settings_from_hdu(hdulist[1], original_data)

    # Test decompression of the first tile

    compressed_tile_bytes = hdulist[1].data["COMPRESSED_DATA"][0].tobytes()

    tile_data_buffer = decompress_tile(
        compressed_tile_bytes, algorithm=compression_type, **settings
    )

    # TODO: determine whether we are happy with having to interpret the returned bytes from
    # the GZip codec or whether we want to set the dtype as a setting to the codec.
    if compression_type.startswith("GZIP"):
        # NOTE: It looks like the data is stored as big endian data even if it was
        # originally little-endian.
        tile_data = (
            np.asarray(tile_data_buffer)
            .view(original_data.dtype.newbyteorder(">"))
            .reshape(tile_shape)
        )
    else:
        tile_data = np.asarray(tile_data_buffer).reshape(tile_shape)

    assert_equal(tile_data, original_data[:4, :4])

    # Now compress the original data and compare to compressed bytes. Since
    # the exact compressed bytes might not match (e.g. for GZIP it will depend
    # on the compression level) we instead put the compressed bytes into the
    # original BinTableHDU, then read it in as a normal compressed HDU and make
    # sure the final data match.

    if compression_type.startswith("GZIP"):
        # NOTE: It looks like the data is stored as big endian data even if it was
        # originally little-endian.
        tile_data_buffer = (
            original_data[:4, :4].astype(original_data.dtype.newbyteorder(">")).data
        )
    else:
        tile_data_buffer = original_data[:4, :4].data

    compressed_tile_bytes = compress_tile(
        tile_data_buffer, algorithm=compression_type, **settings
    )

    # Then check that it also round-trips if we go through fits.open
    if compression_type == "PLIO_1":
        hdulist[1].data["COMPRESSED_DATA"][0] = np.frombuffer(
            compressed_tile_bytes, dtype=np.int16
        )
    else:
        hdulist[1].data["COMPRESSED_DATA"][0] = np.frombuffer(
            compressed_tile_bytes, dtype=np.uint8
        )
    hdulist[1].writeto(tmp_path / "updated.fits")
    hdulist.close()
    hdulist_new = fits.open(tmp_path / "updated.fits")
    assert_equal(hdulist_new[1].data, original_data)
    hdulist_new.close()


@pytest.mark.parametrize(("compression_type", "dtype"), parameters)
def test_decompress_hdu(tmp_path, compression_type, dtype):

    # NOTE: for now this test is designed to compare the Python implementation
    # of decompress_hdu with the C implementation - once we get rid of the C
    # implementation we should update this test.

    original_data = np.arange(144).reshape((12, 12)).astype(dtype)

    header = fits.Header()

    hdu = fits.CompImageHDU(
        original_data, header, compression_type=compression_type, tile_size=(4, 4)
    )

    hdu.writeto(tmp_path / "test.fits")

    # Load in CompImageHDU
    hdulist = fits.open(tmp_path / "test.fits")
    hdu = hdulist[1]

    data = decompress_hdu(hdu)

    assert_equal(data, original_data)

    # NOTE: the following block can be removed once we remove the C implementation
    from astropy.io.fits.compression import decompress_hdu as decompress_hdu_c

    data_c = decompress_hdu_c(hdu)
    assert_equal(data, data_c)

    hdulist.close()


def _assert_heap_same(heap_a, heap_b, n_tiles):
    # The exact length and data might not always match because they might be
    # equivalent but store the tiles in a different order or with padding, so
    # we need to interpret the heap data and actually check chunk by chunk
    heap_header_a = heap_a[: n_tiles * 2 * 4].view(">i4")
    heap_header_b = heap_b[: n_tiles * 2 * 4].view(">i4")
    for itile in range(n_tiles):
        len_a = heap_header_a[itile * 2]
        off_a = heap_header_a[itile * 2 + 1] + n_tiles * 2 * 4
        len_b = heap_header_b[itile * 2]
        off_b = heap_header_b[itile * 2 + 1] + n_tiles * 2 * 4
        data_a = heap_a[off_a : off_a + len_a]
        data_b = heap_b[off_b : off_b + len_b]
        if np.all(data_a[:4] == np.array([31, 139, 8, 0], dtype=np.uint8)):
            # gzip compressed - need to compare decompressed bytes as compression
            # is non-deterministic due to e.g. compression level and so on
            from gzip import decompress

            data_a = decompress(data_a.tobytes())
            data_b = decompress(data_b.tobytes())
            assert data_a == data_b
        else:
            assert_equal(data_a, data_b)


@pytest.mark.parametrize(("compression_type", "dtype"), parameters)
def test_compress_hdu(tmp_path, compression_type, dtype):

    # NOTE: for now this test is designed to compare the Python implementation
    # of compress_hdu with the C implementation - once we get rid of the C
    # implementation we should update this test.

    if compression_type == "PLIO_1":
        pytest.xfail()

    original_data = np.arange(144).reshape((12, 12)).astype(dtype)

    header = fits.Header()

    hdu = fits.CompImageHDU(
        original_data, header, compression_type=compression_type, tile_size=(4, 4)
    )

    heap_length, heap_data = compress_hdu(hdu)

    # NOTE: the following block can be removed once we remove the C implementation
    from astropy.io.fits.compression import compress_hdu as compress_hdu_c

    # Note that when CompImageHDU calls compress_hdu it converts to native byte
    # order first, so we have to do this here.
    hdu.data = hdu.data.astype(hdu.data.dtype.newbyteorder("="))
    heap_length_c, heap_data_c = compress_hdu_c(hdu)

    _assert_heap_same(heap_data, heap_data_c, 9)


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
    with fits.open(request.param, disable_image_compression=True) as hdul:
        yield hdul[1]


@pytest.fixture
def original_int_hdu(canonical_data_base_path):
    with fits.open(canonical_data_base_path / "m13.fits") as hdul:
        yield hdul[0]


def test_canonical_data(original_int_hdu, canonical_int_hdus):
    hdr = canonical_int_hdus.header
    tile_size = (hdr["ZTILE2"], hdr["ZTILE1"])
    compression_type = hdr["ZCMPTYPE"]
    original_tile_1 = original_int_hdu.data[: tile_size[0], : tile_size[1]]
    compressed_tile_bytes = canonical_int_hdus.data["COMPRESSED_DATA"][0].tobytes()

    settings = _header_to_settings(canonical_int_hdus.header)
    tile_data_buffer = decompress_tile(
        compressed_tile_bytes, algorithm=compression_type, **settings
    )
    tile_data = _buffer_to_array(tile_data_buffer, hdr)
    np.testing.assert_allclose(original_tile_1, tile_data)
