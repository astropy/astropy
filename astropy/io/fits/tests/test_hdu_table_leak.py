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
