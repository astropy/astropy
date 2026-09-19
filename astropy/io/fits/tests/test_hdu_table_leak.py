"""Regression for https://github.com/astropy/astropy/issues/20057.

BinTableHDU.load() and dump() must not leave file handles open when
parsing or writing of the ASCII dump files fails partway through.
The leaked handles are detected through the ResourceWarning raised by
the garbage collector for unclosed files.
"""

import gc
import warnings

import pytest

from astropy.io.fits import Column
from astropy.io.fits.column import ColDefs
from astropy.io.fits.hdu.table import BinTableHDU


def _record_resource_warnings(func):
    """Run *func* and return the ``ResourceWarning``s raised while it runs."""
    gc.collect()
    with warnings.catch_warnings(record=True) as catch:
        warnings.simplefilter("always", ResourceWarning)
        func()
        gc.collect()
    return [w for w in catch if issubclass(w.category, ResourceWarning)]


def test_load_closes_files_on_parse_error(tmp_path):
    datafile = tmp_path / "data.txt"
    cdfile = tmp_path / "coldefs.txt"

    datafile.write_text("1 2.0\n")
    # Malformed column definition: a single word where five are expected.
    cdfile.write_text("ONLY_ONE_WORD\n")

    def load():
        with pytest.raises(IndexError):
            BinTableHDU.load(datafile=str(datafile), cdfile=str(cdfile))

    leaked = _record_resource_warnings(load)
    assert not leaked, f"ResourceWarning(s) raised: {[str(w.message) for w in leaked]}"


def test_dump_closes_files_on_error(tmp_path):
    # Make _dump_data fail mid-write: the column definitions reference a
    # column that is not in the data, so writing the table data raises
    # a KeyError after the output files have been opened.
    hdu = BinTableHDU.from_columns([Column(name="col", format="J", array=[1, 2])])
    hdu.columns = ColDefs([hdu.columns[0], Column(name="extra", format="J")])

    datafile = tmp_path / "data.txt"
    cdfile = tmp_path / "coldefs.txt"

    def dump():
        with pytest.raises(KeyError, match="extra"):
            hdu.dump(datafile=str(datafile), cdfile=str(cdfile))

    leaked = _record_resource_warnings(dump)
    assert not leaked, f"ResourceWarning(s) raised: {[str(w.message) for w in leaked]}"
