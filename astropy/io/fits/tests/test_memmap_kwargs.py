# Licensed under a 3-clause BSD style license - see PRIMARY_LICENSE.rst

import io

import numpy as np
import pytest

from astropy.io import fits


@pytest.fixture
def image_file(tmp_path):
    filename = tmp_path / "image.fits"
    fits.HDUList(
        [
            fits.PrimaryHDU(np.arange(12, dtype=">f4").reshape(3, 4)),
            fits.ImageHDU(np.arange(6, dtype="u2")),
        ]
    ).writeto(filename)
    return filename


@pytest.mark.parametrize("memmap", [True, False])
@pytest.mark.parametrize("load_data", [True, False])
def test_image(image_file, memmap, load_data):
    with fits.open(image_file, memmap=memmap) as hdul:
        if load_data:
            hdul[0].data
        kwargs = hdul[0].memmap_kwargs()
        assert kwargs == {
            "filename": str(image_file),
            "mode": "r",
            "offset": 2880,
            "shape": (3, 4),
            "dtype": np.dtype(">f4"),
        }
        np.testing.assert_array_equal(np.memmap(**kwargs), hdul[0].data)
        assert hdul[0].memmap_kwargs(raw=True) == kwargs


def test_image_update_mode(image_file):
    with fits.open(image_file, mode="update") as hdul:
        kwargs = hdul[0].memmap_kwargs()
    assert kwargs["mode"] == "r+"
    array = np.memmap(**kwargs)
    array[0, 0] = 42
    array.flush()
    assert fits.getdata(image_file)[0, 0] == 42


def test_image_scaled(image_file):
    with fits.open(image_file) as hdul:
        with pytest.raises(ValueError, match="scaled with BSCALE/BZERO"):
            hdul[1].memmap_kwargs()
        kwargs = hdul[1].memmap_kwargs(raw=True)
        assert kwargs["dtype"] == np.dtype(">i2")
        np.testing.assert_array_equal(
            np.memmap(**kwargs).astype(int) + hdul[1].header["BZERO"],
            hdul[1].data,
        )

    with fits.open(image_file, do_not_scale_image_data=True) as hdul:
        kwargs = hdul[1].memmap_kwargs()
        np.testing.assert_array_equal(np.memmap(**kwargs), hdul[1].data)


def test_image_blank(tmp_path):
    filename = tmp_path / "blank.fits"
    hdu = fits.PrimaryHDU(np.arange(6, dtype="i2"))
    hdu.header["BLANK"] = 3
    hdu.writeto(filename)
    with fits.open(filename) as hdul:
        with pytest.raises(ValueError, match="uses BLANK"):
            hdul[0].memmap_kwargs()
        assert hdul[0].memmap_kwargs(raw=True)["dtype"] == np.dtype(">i2")


def test_image_data_replaced(image_file):
    with fits.open(image_file) as hdul:
        hdul[0].data = hdul[0].data * 2
        for raw in (False, True):
            with pytest.raises(ValueError, match="data has been replaced"):
                hdul[0].memmap_kwargs(raw=raw)


def test_no_data(tmp_path):
    filename = tmp_path / "empty.fits"
    fits.PrimaryHDU().writeto(filename)
    with fits.open(filename) as hdul:
        with pytest.raises(ValueError, match="does not contain any data"):
            hdul[0].memmap_kwargs()


def test_not_in_file():
    with pytest.raises(ValueError, match="not associated with a file"):
        fits.PrimaryHDU(np.ones(3)).memmap_kwargs(raw=True)


def test_file_like_without_name():
    buffer = io.BytesIO()
    fits.PrimaryHDU(np.ones(3)).writeto(buffer)
    buffer.seek(0)
    with fits.open(buffer) as hdul:
        with pytest.raises(ValueError, match="not associated with a file on disk"):
            hdul[0].memmap_kwargs(raw=True)


def test_file_object(image_file):
    with open(image_file, "rb") as fileobj, fits.open(fileobj) as hdul:
        assert hdul[0].memmap_kwargs()["filename"] == str(image_file)


def test_compressed_file(tmp_path):
    filename = tmp_path / "image.fits.gz"
    fits.PrimaryHDU(np.ones(3)).writeto(filename)
    with fits.open(filename) as hdul:
        with pytest.raises(ValueError, match="gzip-compressed"):
            hdul[0].memmap_kwargs(raw=True)


def test_bintable(tmp_path):
    filename = tmp_path / "table.fits"
    fits.BinTableHDU.from_columns(
        [
            fits.Column("x", "D", array=[1.0, 2.0]),
            fits.Column("y", "3J", array=[[1, 2, 3], [4, 5, 6]]),
            fits.Column("s", "4A", array=["ab", "c"]),
        ]
    ).writeto(filename)
    with fits.open(filename) as hdul:
        kwargs = hdul[1].memmap_kwargs()
        assert kwargs["shape"] == (2,)
        assert kwargs["dtype"] == np.dtype(
            [("x", ">f8"), ("y", ">i4", (3,)), ("s", "S4")]
        )
        array = np.memmap(**kwargs)
        for name in ("x", "y"):
            np.testing.assert_array_equal(array[name], hdul[1].data[name])
        np.testing.assert_array_equal(array["s"], [b"ab", b"c"])


@pytest.mark.parametrize(
    "column, reason",
    [
        (fits.Column("v", "PJ()", array=[[1, 2], [3]]), "has variable-length arrays"),
        (fits.Column("v", "L", array=[True, False]), "is logical"),
        (fits.Column("v", "3X", array=[[1, 0, 1], [0, 0, 1]]), "is a bit array"),
        (fits.Column("v", "E", bscale=2.0, array=[1.0, 2.0]), "is scaled"),
        (
            fits.Column("v", "J", bzero=2**31, array=np.array([1, 2], dtype="u4")),
            "is scaled",
        ),
    ],
)
def test_bintable_conversion(tmp_path, column, reason):
    filename = tmp_path / "table.fits"
    fits.BinTableHDU.from_columns(
        [fits.Column("x", "D", array=[1, 2]), column]
    ).writeto(filename)
    with fits.open(filename) as hdul:
        with pytest.raises(ValueError, match=f"Column 'v' {reason}"):
            hdul[1].memmap_kwargs()
        kwargs = hdul[1].memmap_kwargs(raw=True)
        assert kwargs["dtype"].itemsize == hdul[1].header["NAXIS1"]
        np.testing.assert_array_equal(np.memmap(**kwargs)["x"], [1, 2])


def test_bintable_variable_length_descriptors(tmp_path):
    filename = tmp_path / "table.fits"
    fits.BinTableHDU.from_columns(
        [fits.Column("v", "PJ()", array=[[1, 2], [3]])]
    ).writeto(filename)
    with fits.open(filename) as hdul:
        array = np.memmap(**hdul[1].memmap_kwargs(raw=True))
        # Each row holds the (length, heap offset) array descriptor
        np.testing.assert_array_equal(array["v"], [[2, 0], [1, 8]])


def test_ascii_table(tmp_path):
    filename = tmp_path / "table.fits"
    fits.TableHDU.from_columns(
        [
            fits.Column("a", "I5", array=[1, 2]),
            fits.Column("b", "F8.2", array=[1.5, 2.0]),
        ]
    ).writeto(filename)
    with fits.open(filename) as hdul:
        with pytest.raises(ValueError, match="stored as text in an ASCII table"):
            hdul[1].memmap_kwargs()
        array = np.memmap(**hdul[1].memmap_kwargs(raw=True))
        np.testing.assert_array_equal(array["a"], [b"    1", b"    2"])
        np.testing.assert_array_equal(array["b"], [b"    1.50", b"    2.00"])


def test_compressed_image(tmp_path):
    filename = tmp_path / "compressed.fits"
    data = np.arange(12, dtype="f4").reshape(3, 4)
    fits.CompImageHDU(data).writeto(filename)
    with fits.open(filename) as hdul:
        with pytest.raises(ValueError, match="tile-compressed"):
            hdul[1].memmap_kwargs()
        kwargs = hdul[1].memmap_kwargs(raw=True)
        assert kwargs["shape"] == (hdul[1]._bintable.header["NAXIS2"],)
        assert "COMPRESSED_DATA" in kwargs["dtype"].names


def test_random_groups(tmp_path):
    filename = tmp_path / "groups.fits"
    data = fits.GroupData(
        np.arange(20.0).reshape(10, 1, 2),
        parnames=["abc"],
        pardata=[np.arange(10) + 0.5],
        bitpix=-32,
    )
    fits.GroupsHDU(data).writeto(filename)
    with fits.open(filename) as hdul:
        array = np.memmap(**hdul[0].memmap_kwargs())
        assert array.shape == (10,)
        np.testing.assert_array_equal(array["DATA"], hdul[0].data.data)
        np.testing.assert_array_equal(array["abc"], hdul[0].data.par("abc"))
