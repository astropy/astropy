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

from __future__ import annotations

from . import _connect
from ._connect import *

__all__: list[str] = []
__all__ += _connect.__all__
