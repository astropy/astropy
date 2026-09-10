"""Regression for https://github.com/astropy/astropy/issues/20057.

BinTableHDU.load() must not leave file handles open when parsing of the
ASCII dump files fails partway through.
"""

import os

from astropy.io.fits.hdu.table import BinTableHDU


def test_load_closes_files_on_parse_error(tmp_path):
    datafile = tmp_path / "data.txt"
    cdfile = tmp_path / "coldefs.txt"

    datafile.write_text("1 2.0\n")
    # Malformed column definition: a single word where five are expected.
    cdfile.write_text("ONLY_ONE_WORD\n")

    handles_before = (
        set(os.listdir("/proc/self/fd")) if os.path.isdir("/proc/self/fd") else None
    )

    try:
        BinTableHDU.load(datafile=str(datafile), cdfile=str(cdfile))
    except IndexError:
        pass  # expected: words.pop(0) on an exhausted list
    else:
        raise AssertionError("malformed coldefs should fail parsing")

    # The coldefs handle must have been closed even though parsing raised.
    # On POSIX, verify no handle still points at the coldefs file.
    if os.path.isdir("/proc/self/fd"):
        handles_after = set(os.listdir("/proc/self/fd"))
        leaked = handles_after - (handles_before or set())
        for fd in leaked:
            try:
                target = os.readlink(f"/proc/self/fd/{fd}")
            except OSError:
                continue
            assert "coldefs" not in target, f"leaked handle: {fd} -> {target}"


def test_dump_closes_files_on_error(tmp_path):
    """Verify dump() closes file handles even when writing raises."""
    import numpy as np

    from astropy.io.fits import Column
    from astropy.io.fits.hdu.table import BinTableHDU

    # Create a tiny BinTable with one integer column
    c = Column(name="col", format="J", array=np.array([1, 2]))
    hdu = BinTableHDU.from_columns([c])

    datafile = tmp_path / "data.txt"
    cdfile = tmp_path / "coldefs.txt"

    # Make the file read-only after creation to provoke a write error
    datafile.touch()
    datafile.chmod(0o444)

    try:
        hdu.dump(datafile=str(datafile), cdfile=str(cdfile))
    except (PermissionError, OSError):
        pass  # expected: can't write to read-only file
    else:
        raise AssertionError("dump to read-only file should fail")

    # Verify no leaked handles pointing at datafile
    if os.path.isdir("/proc/self/fd"):
        for fd in os.listdir("/proc/self/fd"):
            try:
                target = os.readlink(f"/proc/self/fd/{fd}")
            except OSError:
                continue
            assert str(datafile) not in target, f"leaked handle: {fd} -> {target}"
