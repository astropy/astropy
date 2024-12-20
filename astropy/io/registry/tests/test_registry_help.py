# Licensed under a 3-clause BSD style license - see LICENSE.rst

import sys
from io import StringIO

import pytest

from astropy.nddata import CCDData
from astropy.table import Table

SKIPIF_OPTIMIZED_PYTHON = pytest.mark.skipif(
    sys.flags.optimize >= 2, reason="docstrings are not available at runtime"
)
ONLY_OPTIMIZED_PYTHON = pytest.mark.skipif(
    sys.flags.optimize < 2,
    reason="checking behavior specifically if docstrings are not present",
)


@SKIPIF_OPTIMIZED_PYTHON
def test_table_read_help_fits():
    """
    Test dynamically created documentation help via the I/O registry for 'fits'.
    """
    out = StringIO()
    Table.read.help("fits", out)
    doc = out.getvalue()

    # Check a smattering of expected content
    assert "Table.read general documentation" not in doc
    assert "The available built-in formats" not in doc
    assert "Table.read(format='fits') documentation" in doc
    assert "hdu : int or str, optional" in doc


@SKIPIF_OPTIMIZED_PYTHON
def test_table_read_help_ascii():
    """
    Test dynamically created documentation help via the I/O registry for 'ascii'.
    """
    out = StringIO()
    Table.read.help("ascii", out)
    doc = out.getvalue()

    # Check a smattering of expected content
    assert "Table.read general documentation" not in doc
    assert "The available built-in formats" not in doc
    assert "Table.read(format='ascii') documentation" in doc
    assert "delimiter : str" in doc
    assert "ASCII reader 'ascii' details" in doc
    assert "Character-delimited table with a single header line" in doc


@SKIPIF_OPTIMIZED_PYTHON
def test_table_write_help_hdf5():
    """
    Test dynamically created documentation help via the I/O registry for 'hdf5'.
    """
    out = StringIO()
    Table.write.help("hdf5", out)
    doc = out.getvalue()

    # Check a smattering of expected content
    assert "Table.write general documentation" not in doc
    assert "The available built-in formats" not in doc
    assert "Table.write(format='hdf5') documentation" in doc
    assert "Write a Table object to an HDF5 file" in doc
    assert "compression : bool or str or int" in doc


def test_list_formats():
    """
    Test getting list of available formats
    """
    out = StringIO()
    CCDData.write.list_formats(out)
    output = out.getvalue()

    assert (
        output
        == """\
Format Read Write Auto-identify
------ ---- ----- -------------
  fits  Yes   Yes           Yes"""
    )


@SKIPIF_OPTIMIZED_PYTHON
def test_table_write_help_fits():
    """
    Test dynamically created documentation help via the I/O registry for 'fits'.
    """
    out = StringIO()
    Table.write.help("fits", out)
    doc = out.getvalue()

    # Check a smattering of expected content
    assert "Table.write general documentation" not in doc
    assert "The available built-in formats" not in doc
    assert "Table.write(format='fits') documentation" in doc
    assert "Write a Table object to a FITS file" in doc


@SKIPIF_OPTIMIZED_PYTHON
def test_table_write_help_no_format():
    """
    Test dynamically created documentation help via the I/O registry for no
    format provided.
    """
    out = StringIO()
    Table.write.help(out=out)
    doc = out.getvalue()

    # Check a smattering of expected content
    assert "Table.write general documentation" in doc
    assert "The available built-in formats" in doc


@SKIPIF_OPTIMIZED_PYTHON
def test_table_read_help_no_format():
    """
    Test dynamically created documentation help via the I/O registry for not
    format provided.
    """
    out = StringIO()
    Table.read.help(out=out)
    doc = out.getvalue()

    # Check a smattering of expected content
    assert "Table.read general documentation" in doc
    assert "The available built-in formats" in doc


@SKIPIF_OPTIMIZED_PYTHON
def test_ccddata_write_help_fits():
    """
    Test dynamically created documentation help via the I/O registry for 'fits'.
    """
    out = StringIO()
    CCDData.write.help("fits", out)
    doc = out.getvalue()

    # Check a smattering of expected content
    assert "CCDData.write(format='fits') documentation" in doc
    assert "Write CCDData object to FITS file" in doc
    assert "key_uncertainty_type : str, optional" in doc


@SKIPIF_OPTIMIZED_PYTHON
def test_ccddata_read_help_fits():
    """Test dynamically created documentation help via the I/O registry for
    CCDData 'fits'.

    """
    out = StringIO()
    CCDData.read.help("fits", out)
    doc = out.getvalue()

    # Check a smattering of expected content
    assert "CCDData.read(format='fits') documentation" in doc
    assert "Generate a CCDData object from a FITS file" in doc
    assert "hdu_uncertainty : str or None, optional" in doc


@SKIPIF_OPTIMIZED_PYTHON
def test_table_write_help_jsviewer():
    """
    Test dynamically created documentation help via the I/O registry for
    'jsviewer'.
    """
    out = StringIO()
    Table.write.help("jsviewer", out)
    doc = out.getvalue()

    # Check a smattering of expected content
    assert "Table.write general documentation" not in doc
    assert "The available built-in formats" not in doc
    assert "Table.write(format='jsviewer') documentation" in doc


@ONLY_OPTIMIZED_PYTHON
@pytest.mark.parametrize("method", [Table.read, Table.write])
@pytest.mark.parametrize("format", [None, "fits", "jsviewer"])
def test_table_read_help_optimized_mode(method, format):
    out = StringIO()
    with pytest.raises(
        RuntimeError,
        match="The help method is not available under Python's optimized mode.",
    ):
        method.help(format, out=out)
    assert out.getvalue() == ""
