# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""I/O subpackage for cosmology.

This subpackage contains classes and functions for reading, writing
and converting cosmology objects to and from various formats.

User access to the I/O functionality is provided through the methods
on the `~astropy.cosmology.Cosmology` class (and its subclasses) and its instances:

- |Cosmology.read| for reading from a file,
- |Cosmology.write| for writing to a file,
- |Cosmology.from_format| to construct a Cosmology from an object
- |Cosmology.to_format| to convert a Cosmology to an object
"""

__all__ = [
    "CosmologyFromFormat",
    "CosmologyRead",
    "CosmologyToFormat",
    "CosmologyWrite",
    "convert_registry",
    "readwrite_registry",
]

# Importing the I/O subpackage registers the I/O methods.
from ._src.io.connect import (
    CosmologyFromFormat,
    CosmologyRead,
    CosmologyToFormat,
    CosmologyWrite,
    convert_registry,
    readwrite_registry,
)
