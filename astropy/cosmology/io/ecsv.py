# Licensed under a 3-clause BSD style license - see LICENSE.rst

import inspect
import os

from astropy.io import registry as io_registry
from astropy.cosmology import Cosmology
from astropy.table import QTable

from .core import from_table, to_table


__all__ = ["read_ecsv", "write_ecsv", "ecsv_identify"]


def read_ecsv(*args, index=None, **kwargs):
    """Read a cosmology from ECSV file.

    Parameters
    ----------
    *args
        Positional arguments passed through to `~astropy.table.QTable` reader.
    index : int, str, or None (optional, keyword-only)
        The row index in the `~astropy.table.QTable`. Required for multi-row
        tables.
    **kwargs
        Keyword arguments passed through to `~astropy.table.QTable` reader.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    See Also
    --------
    astropy.cosmology.io.from_table
    """
    table = QTable.read(*args, format="ascii.ecsv", **kwargs)
    cosmo = from_table(table, index=index)
    return cosmo


def write_ecsv(cosmology, *args, **kwargs):
    """Write a cosmology to a file in an ECSV format.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance
    *args
        Positional arguments passed through to `~astropy.table.QTable` writer.
    **kwargs
        Keyword arguments passed through to `~astropy.table.QTable` writer.
    """
    table = to_table(cosmology)

    # force ascii.ecsv format
    sig = inspect.signature(table.write)
    ba = sig.bind(*args, **kwargs)
    ba.arguments["format"] = "ascii.ecsv"

    table.write(*ba.args, **ba.kwargs)


def ecsv_identify(origin, filepath, fileobj, *args, **kwargs):
    """Identify if object uses the ECSV format."""
    return filepath is not None and filepath.endswith(".ecsv")


io_registry.register_reader("ascii.ecsv", Cosmology, read_ecsv)
io_registry.register_writer("ascii.ecsv", Cosmology, write_ecsv)
io_registry.register_identifier("ascii.ecsv", Cosmology, ecsv_identify)
