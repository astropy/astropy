"""Regression for https://github.com/astropy/astropy/issues/20057.

BinTableHDU.load() and dump() must not leave file handles open when
parsing or writing of the ASCII dump files fails partway through.
"""

import psutil
import pytest

from astropy.io.fits.hdu.table import BinTableHDU


def _open_paths():
    return {f.path for f in psutil.Process().open_files()}


def _assert_no_leak(marker):
    leaked = [p for p in _open_paths() if marker in p]
    assert not leaked, f"leaked handle(s): {leaked}"


def test_load_closes_files_on_parse_error(tmp_path):
    datafile = tmp_path / "data.txt"
    cdfile = tmp_path / "coldefs.txt"

    datafile.write_text("1 2.0\n")
    # Malformed column definition: a single word where five are expected.
    cdfile.write_text("ONLY_ONE_WORD\n")

    with pytest.raises(IndexError):
        BinTableHDU.load(datafile=str(datafile), cdfile=str(cdfile))

    _assert_no_leak("coldefs")


def test_dump_closes_files_on_error(tmp_path):
    import numpy as np

    from astropy.io.fits import Column

    c = Column(name="col", format="J", array=np.array([1, 2]))
    hdu = BinTableHDU.from_columns([c])

    datafile = tmp_path / "data.txt"
    cdfile = tmp_path / "coldefs.txt"

    # Make the file read-only after creation to provoke a write error
    datafile.touch()
    datafile.chmod(0o444)

    with pytest.raises(OSError):
        hdu.dump(datafile=str(datafile), cdfile=str(cdfile))

    _assert_no_leak(str(datafile))
