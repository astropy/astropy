# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Test that all expected methods are present, before I/O tests import.

This file is weirdly named so that it's the first test of I/O.
"""

from astropy.cosmology.io import convert_registry, readwrite_registry


def test_expected_readwrite_io():
    """Test that ONLY the expected I/O is registered."""

    got = {k for k, _ in readwrite_registry._readers.keys()}
    expected = {"ascii.ecsv", "ascii.html"}

    assert got == expected


def test_expected_convert_io():
    """Test that ONLY the expected I/O is registered."""

    got = {k for k, _ in convert_registry._readers.keys()}
    expected = {
        "astropy.cosmology",
        "mapping",
        "astropy.model",
        "astropy.row",
        "astropy.table",
        "yaml",
    }

    assert got == expected
