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
from contextlib import nullcontext

import numpy as np
import pytest

from astropy.io import fits
from astropy.utils import NumpyRNGContext
from astropy.utils.exceptions import AstropyUserWarning

from .conftest import COMPRESSION_TYPES, _expand, fitsio_param_to_astropy_param

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
    if (compression_type.startswith("HCOMPRESS") and "i2" in dtype) or "u1" in dtype:
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

    if compression_type == "HCOMPRESS_1" and "u2" in dtype:
        # cfitsio's HCOMPRESS encoder underestimates the output buffer size
        # for uint16 data on certain tile configurations, raising "encode:
        # output buffer too small". astropy's encoder handles the same combos.
        pytest.xfail("cfitsio HCOMPRESS encoder buffer-size bug for uint16 input")

    if compression_type == "PLIO_1" and "f" in dtype:
        # fitsio fails with a compression error
        pytest.xfail("fitsio fails to write these")

    if (
        compression_type == "PLIO_1"
        and np.dtype(dtype).kind == "u"
        and np.dtype(dtype).itemsize >= 2
    ):
        # PLIO can't represent the BZERO-shifted unsigned values; cfitsio
        # rejects 4/8-byte cases outright and segfaults on 2-byte ones.
        pytest.xfail("PLIO_1 cannot encode unsigned multi-byte integers")

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

    if (
        compression_type == "PLIO_1"
        and np.dtype(dtype).kind == "u"
        and np.dtype(dtype).itemsize >= 2
    ):
        pytest.xfail("PLIO_1 cannot encode unsigned multi-byte integers")

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

    with fits.open(fitsio_compressed_file_path, disable_image_compression=True) as hdul:
        assert hdul[1].header["ZCMPTYPE"].replace("ONE", "1") == compression_type

    with fits.open(fitsio_compressed_file_path) as hdul:
        data = hdul[1].data
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

    # Use relaxed tolerances for HCOMPRESS_1 with quantization
    if compression_type == "HCOMPRESS_1" and "ZQUANTIZ" in header:
        # SUBTRACTIVE_DITHER_2 and other quantization methods can have small precision differences
        rtol = 1e-5
        atol = 1e-5
    else:
        # Use default tolerances for other cases
        rtol = 1e-7
        atol = 0

    with fits.open(astropy_compressed_file_path) as hdul:
        np.testing.assert_allclose(data, hdul[1].data, rtol=rtol, atol=atol)


@pytest.mark.parametrize("kind", ["i", "u"])
@pytest.mark.parametrize(
    "nbytes,overflow", [(2, False), (4, False), (8, False), (8, True)]
)
@pytest.mark.parametrize("compression_type", COMPRESSION_TYPES)
def test_decompress_integers(nbytes, overflow, compression_type, kind, tmp_path):
    if kind == "u" and compression_type == "PLIO_1" and nbytes >= 2:
        pytest.skip(
            "PLIO_1 cannot encode unsigned multi-byte integers (covered elsewhere)"
        )

    testfile = tmp_path / "test.fits.fz"
    data = np.random.poisson(1000, size=(52, 57)).astype(f"{kind}{nbytes}")
    if overflow:
        # push past the 32-bit limit so the conversion fallback fails
        data += np.iinfo(np.int32).max if kind == "i" else np.uint64(2**32)
    data_hdu = fits.PrimaryHDU(data=data)
    compressed_hdu = fits.CompImageHDU(data=data, compression_type=compression_type)

    if compression_type in ("RICE_1", "PLIO_1", "HCOMPRESS_1") and nbytes == 8:
        test_func = pytest.raises if overflow else pytest.warns
        ctx = test_func(
            ValueError if overflow else AstropyUserWarning,
            match=f"{compression_type} compression doesn't support 64 integers.*",
        )
        nbytes = 4
    else:
        ctx = nullcontext()

    with ctx:
        compressed_hdu.writeto(testfile)

    if overflow:
        # nothing more to test
        return

    with fits.open(testfile, disable_image_compression=True) as hdul:
        assert hdul[1].header["ZCMPTYPE"] == compression_type
        if compression_type == "RICE_1":
            assert hdul[1].header["ZNAME2"] == "BYTEPIX"
            assert hdul[1].header["ZVAL2"] == nbytes

    with fits.open(testfile) as hdul:
        np.testing.assert_array_equal(data, hdul[1].data)
        assert hdul[1].data.dtype.kind == np.dtype(data.dtype).kind
        assert hdul[1].data.dtype.itemsize == nbytes

    if compression_type != "NOCOMPRESS" and nbytes != 8:
        # fitsio does not support NOCOMPRESS or 64-bit data
        fts = fitsio.FITS(testfile)
        data2 = fts[1].read()
        np.testing.assert_array_equal(data, data2)


INTEGER_DTYPES_FULL_RANGE = [
    "<i2",
    ">i2",
    "<i4",
    ">i4",
    "<i8",
    ">i8",
    "<u1",
    ">u1",
    "<u2",
    ">u2",
    "<u4",
    ">u4",
    "<u8",
    ">u8",
]


def _full_range_integer_data(dtype, shape=(32, 32)):
    info = np.iinfo(dtype)
    sentinels = np.array(
        [
            info.min,
            info.min + 1,
            -1 if info.min < 0 else 0,
            0,
            1,
            info.max - 1,
            info.max,
        ],
        dtype=dtype,
    )
    with NumpyRNGContext(0):
        data = np.random.randint(
            info.min,
            info.max + 1,
            size=shape,
            dtype=dtype.lstrip("<>"),
        ).astype(dtype)
    data.ravel()[: sentinels.size] = sentinels
    return data


@pytest.mark.parametrize("dtype", INTEGER_DTYPES_FULL_RANGE)
@pytest.mark.parametrize("compression_type", COMPRESSION_TYPES)
def test_integer_full_range_roundtrip(compression_type, dtype, tmp_path):
    """Cross-check exact round-tripping of boundary-spanning integer data
    between astropy and fitsio for every standard FITS integer dtype and
    every compression algorithm that supports it."""
    data = _full_range_integer_data(dtype)
    np_dtype = np.dtype(dtype)
    info = np.iinfo(dtype)
    astropy_path = tmp_path / "astropy.fits"
    hdu = fits.CompImageHDU(data=data, compression_type=compression_type)

    # 64-bit integer data with full-range sentinels overflows the 32-bit
    # conversion that RICE_1 and HCOMPRESS_1 fall back to. (PLIO_1 + i8 also
    # lands here; PLIO_1 + u8 is caught by the unsigned-multi-byte rejection
    # below.) For these combos cfitsio also rejects the input outright, so
    # cross-check that fitsio raises too.
    if (
        compression_type in ("RICE_1", "HCOMPRESS_1")
        and np_dtype.kind in ("i", "u")
        and np_dtype.itemsize == 8
    ) or (
        compression_type == "PLIO_1" and np_dtype.kind == "i" and np_dtype.itemsize == 8
    ):
        with pytest.raises(
            ValueError,
            match=(
                f"{compression_type} compression doesn't support 64 integers, "
                "but data cannot be converted to 32 bits without overflow"
            ),
        ):
            hdu.writeto(astropy_path)
        if compression_type == "HCOMPRESS_1":
            with pytest.raises(
                OSError,
                match=r"writing T(U)?LONGLONG to compressed image is not supported",
            ):
                with fitsio.FITS(tmp_path / "fitsio.fits", "rw") as fts:
                    fts.write(data, compress=compression_type)
        return

    # PLIO_1 cannot encode unsigned multi-byte integers because the FITS
    # BZERO=2**(N-1) offset astropy applies to unsigned data produces negative
    # values that PLIO rejects. Astropy raises a clean ValueError for all
    # unsigned itemsize >= 2; cfitsio raises only for itemsize >= 4 and
    # segfaults on itemsize == 2 with full-range data, so the cross-check
    # only runs for the 4/8-byte cases.
    if compression_type == "PLIO_1" and np_dtype.kind == "u" and np_dtype.itemsize >= 2:
        with pytest.raises(
            ValueError,
            match=r"PLIO_1 compression does not support unsigned integers",
        ):
            hdu.writeto(astropy_path)
        if np_dtype.itemsize >= 4:
            with pytest.raises(
                ValueError,
                match=r"Unsigned 4/8-byte integers currently not allowed",
            ):
                with fitsio.FITS(tmp_path / "fitsio.fits", "rw") as fts:
                    fts.write(data, compress=compression_type)
        return

    # PLIO_1 also can't encode signed values outside [0, 2**24 - 1].
    if compression_type == "PLIO_1" and (info.min < 0 or info.max > 2**24 - 1):
        with pytest.raises(
            ValueError,
            match=r"data out of range for PLIO compression",
        ):
            hdu.writeto(astropy_path)
        return

    # 1. astropy writes a compressed file.
    hdu.writeto(astropy_path)

    # 2. astropy reads its own compressed file.
    with fits.open(astropy_path) as hdul:
        rt = hdul[1].data
        assert rt.dtype.kind == np_dtype.kind
        assert rt.dtype.itemsize == np_dtype.itemsize
        np.testing.assert_array_equal(rt, data)

    # fitsio doesn't support NOCOMPRESS, so the cross-checks stop here.
    # cfitsio also refuses to read or write GZIP-compressed 64-bit integer
    # images even though the FITS Tile Compression Convention permits them,
    # so astropy's standard-compliant output cannot be cross-validated for
    # those combinations.
    if compression_type == "NOCOMPRESS" or (
        compression_type in ("GZIP_1", "GZIP_2") and np_dtype.itemsize == 8
    ):
        return

    # 3. fitsio reads astropy's compressed file.
    with fitsio.FITS(astropy_path) as fts:
        rt_fitsio = fts[1].read()
        assert rt_fitsio.dtype.kind == np_dtype.kind
        assert rt_fitsio.dtype.itemsize == np_dtype.itemsize
        np.testing.assert_array_equal(rt_fitsio, data)

    # 4. astropy reads a file fitsio compressed.
    fitsio_path = tmp_path / "fitsio.fits"
    with fitsio.FITS(fitsio_path, "rw") as fts:
        fts.write(data, compress=compression_type)
    with fits.open(fitsio_path) as hdul:
        rt_astropy = hdul[1].data
        assert rt_astropy.dtype.kind == np_dtype.kind
        assert rt_astropy.dtype.itemsize == np_dtype.itemsize
        np.testing.assert_array_equal(rt_astropy, data)
