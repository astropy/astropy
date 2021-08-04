# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see Astropy LICENSE.rst

"""
Register Read/Write methods for "myformat" with Astropy Cosmology.

With this format registered, we can start with a Cosmology from
``mypackage``, write it to a file, and read it with Astropy to create an
astropy Cosmology instance.

    >>> from mypackage.cosmology import somecosmologyobject
    >>> from mypackage.io import file_writer
    >>> file_writer('<file name>', somecosmologyobject)

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

# THIRD PARTY
from astropy.cosmology import Cosmology
from astropy.io import registry as io_registry


def read_myformat(filename, **kwargs):
    if isinstance(filename, (str, bytes, os.PathLike)):
        with open(filename, "r") as file:
            data = file.read()
    else:  # file-like : this also handles errors in dumping
       data = filename.read()

    mapping = {}
    ...  # process `data`, adding to `mapping`

    return Cosmology.from_format(mapping, **kwargs)


def write_myformat(cosmology, file, *, overwrite=False, **kwargs):
    data = cosmology.to_format("mapping")  # start by turning into dict
    data["cosmology"] = data["cosmology"].__name__  # change class field to str

    output = ...  # whatever it's supposed to be
    ...  # process `data` into correct "myformat"

    if isinstance(file, (str, bytes, os.PathLike)):
        # check that file exists and whether to overwrite.
        if os.path.exists(file) and not overwrite:
            raise IOError(f"{file} exists. Set 'overwrite' to write over.")
        with open(file, "w") as write_file:
            write_file.write(output)
    else:  # file-like
        file.write(output)


def myformat_identify(origin, filepath, fileobj, *args, **kwargs):
    """Identify if object uses ``myformat``."""
    return filepath is not None and filepath.endswith(".myformat")


# -------------------------------------------------------------------
# Register read/write/identify methods with Astropy Unified I/O

io_registry.register_reader("myformat", Cosmology, read_myformat)
io_registry.register_writer("myformat", Cosmology, write_myformat)
io_registry.register_identifier("myformat", Cosmology, myformat_identify)
