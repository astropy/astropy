"""
Tests for ECSV format with column names starting with '#'.

This tests the fix for issue #18710 where column names starting with '#'
would cause ECSV read failures.
"""

import io

import pytest

from astropy.table import Table
from astropy.io import ascii


def test_column_name_starting_with_hash():
    """Test that a column name starting with '#' can be written and read back."""
    # Create a table with a column name starting with '#'
    t = Table({"#column": [1, 2, 3], "normal": [4, 5, 6]})

    # Write to ECSV
    output = io.StringIO()
    t.write(output, format="ascii.ecsv")

    # Read it back
    ecsv_content = output.getvalue()
    t_read = Table.read(ecsv_content, format="ascii.ecsv")

    # Verify the table was read correctly
    assert t_read.colnames == ["#column", "normal"]
    assert list(t_read["#column"]) == [1, 2, 3]
    assert list(t_read["normal"]) == [4, 5, 6]


def test_column_name_starting_with_hash_and_space():
    """Test that a column name starting with '# ' (hash + space) works."""
    # Create a table with a column name starting with '# '
    t = Table({"# column": [1, 2, 3], "a": [4, 5, 6]})

    # Write to ECSV
    output = io.StringIO()
    t.write(output, format="ascii.ecsv")

    # Read it back
    ecsv_content = output.getvalue()
    t_read = Table.read(ecsv_content, format="ascii.ecsv")

    # Verify the table was read correctly
    assert t_read.colnames == ["# column", "a"]
    assert list(t_read["# column"]) == [1, 2, 3]
    assert list(t_read["a"]) == [4, 5, 6]


def test_multiple_columns_starting_with_hash():
    """Test multiple column names starting with '#'."""
    t = Table({"#col1": [1, 2], "#col2": [3, 4], "normal": [5, 6]})

    # Write to ECSV
    output = io.StringIO()
    t.write(output, format="ascii.ecsv")

    # Read it back
    ecsv_content = output.getvalue()
    t_read = Table.read(ecsv_content, format="ascii.ecsv")

    # Verify the table was read correctly
    assert t_read.colnames == ["#col1", "#col2", "normal"]
    assert list(t_read["#col1"]) == [1, 2]
    assert list(t_read["#col2"]) == [3, 4]
    assert list(t_read["normal"]) == [5, 6]


def test_column_name_hash_with_comma_delimiter():
    """Test column name starting with '#' using comma delimiter."""
    t = Table({"#column": [1, 2, 3], "b": [4, 5, 6]})

    # Write to ECSV with comma delimiter
    output = io.StringIO()
    t.write(output, format="ascii.ecsv", delimiter=",")

    # Read it back
    ecsv_content = output.getvalue()
    t_read = Table.read(ecsv_content, format="ascii.ecsv")

    # Verify the table was read correctly
    assert t_read.colnames == ["#column", "b"]
    assert list(t_read["#column"]) == [1, 2, 3]
    assert list(t_read["b"]) == [4, 5, 6]


def test_ecsv_content_format():
    """Verify that column names starting with '#' are quoted in ECSV output."""
    t = Table({"#col": [1, 2], "a": [3, 4]})

    # Write to ECSV
    output = io.StringIO()
    t.write(output, format="ascii.ecsv")
    ecsv_content = output.getvalue()

    # The column header line should have the #col quoted
    lines = ecsv_content.strip().split("\n")
    # Find the column header line (last non-comment line before data)
    header_line = None
    for i, line in enumerate(lines):
        if not line.startswith("#") and i > 0:
            # This should be the column header line
            header_line = line
            break

    assert header_line is not None
    # The column name starting with # should be quoted
    assert '"#col"' in header_line or "'#col'" in header_line


def test_original_issue_18710():
    """Test the exact scenario from issue #18710."""
    # Create a table with a column whose name starts '#'
    bad_mojo = Table({"#sdfafds": [1, 2, 3], "a": [1, 2, 3]})

    # Write it to ECSV; there should be no error at write
    output = io.StringIO()
    bad_mojo.write(output, format="ascii.ecsv")

    # Read it back in; should not fail with InconsistentTableError
    ecsv_content = output.getvalue()
    result = Table.read(ecsv_content, format="ascii.ecsv")

    # Verify the table was read correctly
    assert result.colnames == ["#sdfafds", "a"]
    assert list(result["#sdfafds"]) == [1, 2, 3]
    assert list(result["a"]) == [1, 2, 3]


def test_column_name_with_quotes_and_hash():
    """Test column names that have both quotes and hash."""
    t = Table({'#"quoted"': [1, 2], "normal": [3, 4]})

    # Write to ECSV
    output = io.StringIO()
    t.write(output, format="ascii.ecsv")

    # Read it back
    ecsv_content = output.getvalue()
    t_read = Table.read(ecsv_content, format="ascii.ecsv")

    # Verify the table was read correctly
    assert t_read.colnames == ['#"quoted"', "normal"]
    assert list(t_read['#"quoted"']) == [1, 2]
    assert list(t_read["normal"]) == [3, 4]
