# Licensed under a 3-clause BSD style license - see PRIMARY_LICENSE.rst

import io

import numpy as np
import pytest

from astropy.io import fits


@pytest.mark.parametrize("mode", ["readonly", "update"])
def test_image(tmp_path, mode):
    filename = tmp_path / "test.fits"
    data = np.arange(12, dtype=">f4").reshape(3, 4)
    fits.PrimaryHDU(data).writeto(filename)
    with fits.open(filename, mode=mode) as hdul:
        info = hdul[0].memmap_info()
    assert info == {
        "filename": str(filename),
        "mode": "r",
        "offset": 2880,
        "shape": (3, 4),
        "dtype": np.dtype(">f4"),
    }
    np.testing.assert_array_equal(np.memmap(**info), data)


def test_bintable(tmp_path):
    filename = tmp_path / "test.fits"
    fits.BinTableHDU.from_columns(
        [
            fits.Column("x", "D", array=[1.0, 2.0]),
            fits.Column("y", "3J", array=[[1, 2, 3], [4, 5, 6]]),
        ]
    ).writeto(filename)
    with fits.open(filename) as hdul:
        array = np.memmap(**hdul[1].memmap_info())
        np.testing.assert_array_equal(array, hdul[1].data)


def test_random_groups(tmp_path):
    filename = tmp_path / "test.fits"
    data = fits.GroupData(
        np.arange(20.0).reshape(10, 1, 2),
        parnames=["abc"],
        pardata=[np.arange(10) + 0.5],
        bitpix=-32,
    )
    fits.GroupsHDU(data).writeto(filename)
    with fits.open(filename) as hdul:
        array = np.memmap(**hdul[0].memmap_info())
        np.testing.assert_array_equal(array["DATA"], hdul[0].data.data)


def _image(**header):
    hdu = fits.PrimaryHDU(np.arange(6, dtype="i2"))
    hdu.header.update(header)
    return hdu


def _table(format, **kwargs):
    return fits.BinTableHDU.from_columns([fits.Column("v", format, **kwargs)])


@pytest.mark.parametrize(
    "hdu, match",
    [
        (_image(BZERO=10), "scaled with BSCALE/BZERO or uses BLANK"),
        (_image(BLANK=3), "scaled with BSCALE/BZERO or uses BLANK"),
        (fits.PrimaryHDU(), "does not contain any data"),
        (fits.CompImageHDU(np.ones((4, 4))), "tile-compressed"),
        (_table("PJ()", array=[[1, 2], [3]]), "Column 'v' is scaled, logical"),
        (_table("L", array=[True, False]), "Column 'v' is scaled, logical"),
        (_table("3X", array=[[1, 0, 1], [0, 0, 1]]), "Column 'v' is scaled, logical"),
        (_table("E", bscale=2.0, array=[1.0, 2.0]), "Column 'v' is scaled, logical"),
        (
            fits.TableHDU.from_columns([fits.Column("v", "I5", array=[1, 2])]),
            "stored as text",
        ),
    ],
)
def test_not_memmappable(tmp_path, hdu, match):
    filename = tmp_path / "test.fits"
    hdu.writeto(filename)
    with fits.open(filename) as hdul:
        with pytest.raises(ValueError, match=match):
            hdul[-1].memmap_info()


def test_not_on_disk(tmp_path):
    with pytest.raises(ValueError, match="not associated with a file"):
        fits.PrimaryHDU(np.ones(3)).memmap_info()

    buffer = io.BytesIO()
    fits.PrimaryHDU(np.ones(3)).writeto(buffer)
    buffer.seek(0)
    with fits.open(buffer) as hdul:
        with pytest.raises(ValueError, match="not associated with a file on disk"):
            hdul[0].memmap_info()

    fits.PrimaryHDU(np.ones(3)).writeto(tmp_path / "test.fits.gz")
    with fits.open(tmp_path / "test.fits.gz") as hdul:
        with pytest.raises(ValueError, match="gzip-compressed"):
            hdul[0].memmap_info()


def test_data_replaced(tmp_path):
    fits.PrimaryHDU(np.ones(3)).writeto(tmp_path / "test.fits")
    with fits.open(tmp_path / "test.fits") as hdul:
        hdul[0].data = hdul[0].data * 2
        with pytest.raises(ValueError, match="data has been replaced"):
            hdul[0].memmap_info()
