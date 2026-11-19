# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for astropy.utils.data with compressed files."""

import gzip
import io
import os
import tracemalloc

import pytest

from astropy.utils.compat.optional_deps import HAS_BZ2, HAS_LZMA, HAS_UNCOMPRESSPY
from astropy.utils.data import conf, get_readable_fileobj


@pytest.fixture(
    params=["invalid.dat.bz2", "invalid.dat.xz", "invalid.dat.gz", "invalid.dat.Z"]
)
def bad_compressed(request, tmp_path):
    # These contents have valid headers for their respective file formats, but
    # are otherwise malformed and invalid.
    bz_content = b"BZhinvalid"
    gz_content = b"\x1f\x8b\x08invalid"
    xz_content = b"\xfd7zXZ\x00invalid"
    lzw_content = b"\x1f\x9d\x90invalid"

    datafile = tmp_path / request.param
    filename = str(datafile)

    if filename.endswith(".bz2"):
        contents = bz_content
    elif filename.endswith(".gz"):
        contents = gz_content
    elif filename.endswith(".xz"):
        contents = xz_content
    elif filename.endswith(".Z"):
        contents = lzw_content
    else:
        contents = "invalid"

    datafile.write_bytes(contents)

    return filename


def test_local_data_obj_invalid(bad_compressed):
    is_bz2 = bad_compressed.endswith(".bz2")
    is_xz = bad_compressed.endswith(".xz")
    is_lzw = bad_compressed.endswith(".Z")

    # Note, since these invalid files are created on the fly in order to avoid
    # problems with detection by antivirus software
    # (see https://github.com/astropy/astropy/issues/6520), it is no longer
    # possible to use ``get_pkg_data_fileobj`` to read the files. Technically,
    # they're not local anymore: they just live in a temporary directory
    # created by pytest. However, we can still use get_readable_fileobj for the
    # test.
    if (
        (not HAS_BZ2 and is_bz2)
        or (not HAS_LZMA and is_xz)
        or (not HAS_UNCOMPRESSPY and is_lzw)
    ):
        with pytest.raises(ModuleNotFoundError):
            with get_readable_fileobj(bad_compressed, encoding="binary") as f:
                f.read()
    else:
        with get_readable_fileobj(bad_compressed, encoding="binary") as f:
            assert f.read().rstrip().endswith(b"invalid")


def _compress_lzw(payload):
    # LZW (.Z) can not be written with the standard library or uncompresspy.
    ncompress = pytest.importorskip("ncompress")
    return ncompress.compress(payload)


BAD_COMPRESSORS = [
    pytest.param("gz", lambda b: gzip.compress(b, compresslevel=9), id="gzip"),
]

if HAS_BZ2:
    import bz2

    BAD_COMPRESSORS.append(pytest.param("bz2", bz2.compress, id="bzip2"))
if HAS_LZMA:
    import lzma

    BAD_COMPRESSORS.append(
        pytest.param(
            "xz",
            lambda b: lzma.compress(b, preset=lzma.PRESET_EXTREME),
            id="xz",
        )
    )
if HAS_UNCOMPRESSPY:
    BAD_COMPRESSORS.append(pytest.param("Z", _compress_lzw, id="lzw"))


MiB = 1024**2
# Large enough to be clearly "too much" for the limit, small enough to be
# harmless on any CI machine even if the limit is not enforced.
TOO_BIG_SIZE = 32 * MiB
TOO_BIG_LIMIT = 1 * MiB
READ_CHUNK = 64 * 1024


@pytest.fixture
def bad_factory(tmp_path):
    def make(suffix, compressor):
        compressed = compressor(b"\0" * TOO_BIG_SIZE)
        assert len(compressed) < TOO_BIG_SIZE // 100  # a genuinely huge ratio
        path = tmp_path / f"bad.{suffix}"
        path.write_bytes(compressed)
        return path

    return make


def _read_until_failure(path, encoding):
    """Read in small chunks; return (bytes_delivered, exception, peak_memory)."""
    delivered = 0
    error = None
    tracemalloc.start()
    try:
        with get_readable_fileobj(path, encoding=encoding) as f:
            while chunk := f.read(READ_CHUNK):
                delivered += len(chunk)
    except Exception as e:
        error = e
    finally:
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
    return delivered, error, peak


@pytest.mark.parametrize("encoding", ["binary", "utf-8"])
@pytest.mark.parametrize(("suffix", "compressor"), BAD_COMPRESSORS)
def test_decompressed_size_limit_is_enforced_incrementally(
    bad_factory, suffix, compressor, encoding
):
    """The size limit must trip *while* decompressing, not after the fact.

    So: an error is raised, no more than ~limit bytes were ever delivered, and
    peak Python memory stays far below the size of the full bomb.  The latter
    catches implementations that decompress everything (in memory or to a
    temporary file via a single ``read()``) and only then check the size.
    """
    path = bad_factory(suffix, compressor)
    with (
        conf.set_temp("max_decompressed_size", TOO_BIG_LIMIT),
        conf.set_temp("max_decompression_ratio", 0),
    ):
        delivered, error, peak = _read_until_failure(path, encoding)

    assert isinstance(error, ValueError)
    assert "max_decompressed_size" in str(error)  # tells the user how to override
    assert delivered <= TOO_BIG_LIMIT + READ_CHUNK
    assert peak < TOO_BIG_SIZE // 4


@pytest.mark.parametrize("encoding", ["binary", "utf-8"])
@pytest.mark.parametrize(("suffix", "compressor"), BAD_COMPRESSORS)
def test_decompression_ratio_limit_is_enforced_incrementally(
    bad_factory, suffix, compressor, encoding
):
    """Same as above, but triggered by the compression ratio alone."""
    path = bad_factory(suffix, compressor)
    with (
        conf.set_temp("max_decompressed_size", 0),
        conf.set_temp("max_decompression_ratio", 20),
    ):
        delivered, error, peak = _read_until_failure(path, encoding)

    assert isinstance(error, ValueError)
    assert "max_decompression_ratio" in str(error)
    assert delivered < TOO_BIG_SIZE // 2
    assert peak < TOO_BIG_SIZE // 4


@pytest.mark.parametrize(("suffix", "compressor"), BAD_COMPRESSORS)
def test_decompression_limit_seek_is_counted(bad_factory, suffix, compressor):
    """Seeking forward (or to the end) decompresses data, so it must be limited."""
    path = bad_factory(suffix, compressor)
    with (
        conf.set_temp("max_decompressed_size", TOO_BIG_LIMIT),
        conf.set_temp("max_decompression_ratio", 0),
    ):
        for seek_args in [(TOO_BIG_SIZE // 2,), (0, io.SEEK_END)]:
            with get_readable_fileobj(path, encoding="binary") as f:
                with pytest.raises(ValueError, match="max_decompressed_size"):
                    f.seek(*seek_args)


@pytest.mark.parametrize(("suffix", "compressor"), BAD_COMPRESSORS)
def test_decompression_limit_can_be_disabled_and_rewind_is_free(
    bad_factory, suffix, compressor
):
    path = bad_factory(suffix, compressor)
    with (
        conf.set_temp("max_decompressed_size", 0),
        conf.set_temp("max_decompression_ratio", 0),
    ):
        with get_readable_fileobj(path, encoding="binary") as f:
            assert f.seek(0, io.SEEK_END) == TOO_BIG_SIZE
    # Reading exactly up to the limit, and re-reading after rewinding, is fine.
    with (
        conf.set_temp("max_decompressed_size", TOO_BIG_SIZE),
        conf.set_temp("max_decompression_ratio", 0),
    ):
        with get_readable_fileobj(path, encoding="binary") as f:
            assert len(f.read()) == TOO_BIG_SIZE
            f.seek(0)
            assert len(f.read()) == TOO_BIG_SIZE


def test_rewind_is_free():
    """Rewinding does not count against the limit.

    To test that, we read and rewind a number of times with a file that is otherwise good.
    But if rewind wasn't free, we would exceed the ratio limit after reading a few times.
    """
    payload = os.urandom(1024**2) * 12  # 12 MB, modest ratio
    compressed = gzip.compress(payload)
    with conf.set_temp("max_decompression_ratio", 10):
        with get_readable_fileobj(io.BytesIO(compressed), encoding="binary") as f:
            for i in range(10):
                f.read()
                f.seek(0)


@pytest.mark.parametrize(("suffix", "compressor"), BAD_COMPRESSORS)
def test_legitimate_compressed_file_under_limits(tmp_path, suffix, compressor):
    """Guard against over-eager limits: ordinary data must still read fine."""
    payload = os.urandom(1024) * 64  # 64 KiB, modest ratio
    path = tmp_path / f"ok.{suffix}"
    path.write_bytes(compressor(payload))
    with get_readable_fileobj(path, encoding="binary") as f:
        assert f.read() == payload
