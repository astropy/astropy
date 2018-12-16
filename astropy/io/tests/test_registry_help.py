# Licensed under a 3-clause BSD style license - see LICENSE.rst

from io import StringIO

from astropy.table import Table


def test_table_read_help_fits():
    """
    Test dynamically created documentation help via the I/O registry for 'fits'.
    """
    out = StringIO()
    Table.read.help('fits', out)
    doc = out.getvalue()

    # Check a smattering of expected content
    assert "Table.read general documentation" in doc
    assert "The available built-in formats" in doc
    assert "Table.read(format='fits') documentation" in doc
    assert "hdu : int or str, optional" in doc


def test_table_read_help_ascii():
    """
    Test dynamically created documentation help via the I/O registry for 'ascii'.
    """
    out = StringIO()
    Table.read.help('ascii', out)
    doc = out.getvalue()

    # Check a smattering of expected content
    assert "Table.read general documentation" in doc
    assert "The available built-in formats" in doc
    assert "Table.read(format='ascii') documentation" in doc
    assert "delimiter : str" in doc
    assert "ASCII reader 'ascii' details" in doc
    assert "Character-delimited table with a single header line" in doc


def test_table_write_help_hdf5():
    """
    Test dynamically created documentation help via the I/O registry for 'hdf5'.
    """
    out = StringIO()
    Table.write.help('hdf5', out)
    doc = out.getvalue()

    # Check a smattering of expected content
    assert "Table.write general documentation" in doc
    assert "The available built-in formats" in doc
    assert "Table.write(format='hdf5') documentation" in doc
    assert "Write a Table object to an HDF5 file" in doc
    assert "compression : bool or str or int" in doc


def test_table_write_help_fits():
    """
    Test dynamically created documentation help via the I/O registry for 'fits'.
    """
    out = StringIO()
    Table.write.help('fits', out)
    doc = out.getvalue()

    # Check a smattering of expected content
    assert "Table.write general documentation" in doc
    assert "The available built-in formats" in doc
    assert "Table.write(format='fits') documentation" in doc
    assert "Write a Table object to a FITS file" in doc
