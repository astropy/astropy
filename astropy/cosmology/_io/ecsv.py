# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Reader and Writer for `~astropy.cosmology.Cosmology` in the ECSV format.

This module is NOT considered public API.

"""

import inspect
import os

from astropy.io import registry as io_registry
from astropy.cosmology import Cosmology
from astropy.table import QTable

from .table import from_table, to_table


__all__ = ["read_ecsv", "write_ecsv", "ecsv_identify"]


def read_ecsv(*args, index=None, move_to_meta=False, **kwargs):
    """Read a cosmology from :class:`~astropy.io.ascii.Ecsv` file.

    Parameters
    ----------
    *args
        Positional arguments passed through to `~astropy.table.QTable` reader.
    index : int, str, or None (optional, keyword-only)
        The row index in the `~astropy.table.QTable`. Required for multi-row
        tables.
    move_to_meta : bool (optional, keyword-only)
        Whether to move keyword arguments that are not in the Cosmology class'
        signature to the Cosmology's metadata. This will only be applied if the
        Cosmology does NOT have a keyword-only argument (e.g. ``**kwargs``).
        Arguments moved to the metadata will be merged with existing metadata,
        preferring specified metadata in the case of a merge conflict
        (e.g. for ``Cosmology(meta={'key':10}, key=42)``, the ``Cosmology.meta``
        will be ``{'key': 10}``).
    **kwargs
        Keyword arguments passed through to `~astropy.table.QTable` reader.
        If 'format' is a kwarg, it must be 'ascii.ecsv'.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    See Also
    --------
    astropy.cosmology.io.from_table
    """
    # check 'format' correctness
    format = kwargs.pop("format", "ascii.ecsv")
    if format != "ascii.ecsv":  # check that if specified, it's ECSV
        raise ValueError(f"'format', if specified, must be 'ascii.ecsv' not {format}")

    table = QTable.read(*args, format="ascii.ecsv", **kwargs)
    cosmo = from_table(table, index=index, move_to_meta=move_to_meta)
    return cosmo


def write_ecsv(cosmology, *args, **kwargs):
    """Write a cosmology to a file in an :class:`~astropy.io.ascii.Ecsv` format.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance
    *args
        Positional arguments passed through to `~astropy.table.QTable` writer.
    **kwargs
        Keyword arguments passed through to `~astropy.table.QTable` writer.
        If 'format' is a kwarg, it must be 'ascii.ecsv'.
    """
    table = to_table(cosmology)

    # force ascii.ecsv format
    sig = inspect.signature(table.write)
    ba = sig.bind(*args, **kwargs)
    format = ba.arguments.setdefault("format", "ascii.ecsv")
    if format != "ascii.ecsv":  # check that if specified, it's ECSV
        raise ValueError("'format', if an argument, must be 'ascii.ecsv'")

    table.write(*ba.args, **ba.kwargs)


def ecsv_identify(origin, filepath, fileobj, *args, **kwargs):
    """Identify if object uses the :class:`~astropy.io.ascii.Ecsv` format.

    Returns
    -------
    bool
        True if ``pyyaml`` is installed and the 'filepath' suffix is '.ecsv',
        False otherwise.
    """
    return filepath is not None and filepath.endswith(".ecsv")


# ===================================================================
# Register

io_registry.register_reader("ascii.ecsv", Cosmology, read_ecsv)
io_registry.register_writer("ascii.ecsv", Cosmology, write_ecsv)
io_registry.register_identifier("ascii.ecsv", Cosmology, ecsv_identify)
