# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see Astropy LICENSE.rst

"""
Register Read/Write methods for "myformat" (JSON) with Astropy Cosmology.

With this format registered, we can start with a Cosmology from
``mypackage``, write it to a file, and read it with Astropy to create an
astropy Cosmology instance.

    >>> from mypackage.cosmology import myplanck
    >>> from mypackage.io import file_writer
    >>> file_writer('<file name>', myplanck)

    >>> from astropy.cosmology import Cosmology
    >>> cosmo = Cosmology.read('<file name>', format="myformat")
    >>> cosmo

We can also do the reverse: start with an astropy Cosmology, save it and
read it with ``mypackage``.

    >>> from astropy.cosmology import Planck18
    >>> Planck18.write('<file name>', format="myformat")

    >>> from mypackage.io import file_reader
    >>> cosmo2 = file_reader('<file name>')
    >>> cosmo2

"""

# STDLIB
import json
import os

# THIRD PARTY
import astropy.units as u
from astropy.cosmology import Cosmology
from astropy.cosmology.connect import readwrite_registry

# LOCAL
from .core import file_reader, file_writer

__doctest_skip__ = ['*']


def read_myformat(filename, **kwargs):
    """Read files in format 'myformat'.

    Parameters
    ----------
    filename : str
    **kwargs
        Keyword arguments into `astropy.cosmology.Cosmology.from_format`
        with ``format="mypackage"``.

    Returns
    -------
    `~mypackage.cosmology.MyCosmology` instance
    """
    mycosmo = file_reader(filename)  # ← read file  ↓ build Cosmology
    return Cosmology.from_format(mycosmo, format="mypackage", **kwargs)


def write_myformat(cosmology, file, *, overwrite=False):
    """Write files in format 'myformat'.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` instance
    file : str, bytes, or `~os.PathLike`
    overwrite : bool (optional, keyword-only)
        Whether to overwrite an existing file. Default is False.
    """
    cosmo = cosmology.to_format("mypackage")  # ← convert Cosmology ↓ write file
    file_writer(file, cosmo, overwrite=overwrite)


def myformat_identify(origin, filepath, fileobj, *args, **kwargs):
    """Identify if object uses ``myformat`` (JSON)."""
    return filepath is not None and filepath.endswith(".myformat")


# -------------------------------------------------------------------
# Register read/write/identify methods with Astropy Unified I/O

readwrite_registry.register_reader("myformat", Cosmology, read_myformat, force=True)
readwrite_registry.register_writer("myformat", Cosmology, write_myformat, force=True)
readwrite_registry.register_identifier("myformat", Cosmology, myformat_identify, force=True)
