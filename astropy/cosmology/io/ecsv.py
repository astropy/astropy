# Licensed under a 3-clause BSD style license - see LICENSE.rst

import inspect
import os

from astropy.io import registry as io_registry
from astropy.cosmology import Cosmology
from astropy.table import QTable

from .common import from_table, to_table


__all__ = ["read_ecsv", "write_ecsv", "ecsv_identify"]


def read_ecsv(*args, index=None, **kwargs):
    """Read from ECSV file.

    Parameters
    ----------
    *args
    key : str or None, optional
        If the JSON is for many cosmologies, ``key`` is needed to select
        the specific cosmology.
        See ``from_mapping``.

    """
    table = QTable.read(*args, format="ascii.ecsv", **kwargs)
    cosmo = from_table(table, index=index)
    return cosmo


def write_ecsv(cosmology, *args, **kwargs):
    """Write to ECSV file.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance
    file : path-like or file-like
    overwrite : bool, optional
        Whether to overwrite the file. Applies only if "file" is path-like.

    """
    table = to_table(cosmology)

    sig = inspect.signature(table.write)
    ba = sig.bind(*args, **kwargs)
    ba.arguments["format"] = "ascii.ecsv"  # force ascii.ecsv

    table.write(*ba.args, **ba.kwargs)


def ecsv_identify(origin, filepath, fileobj, *args, **kwargs):
    return filepath is not None and filepath.endswith(".ecsv")


io_registry.register_reader("ascii.ecsv", Cosmology, read_ecsv)
io_registry.register_writer("ascii.ecsv", Cosmology, write_ecsv)
io_registry.register_identifier("ascii.ecsv", Cosmology, ecsv_identify)
