# Licensed under a 3-clause BSD style license - see LICENSE.rst

import json
import os

from astropy.io import registry as io_registry
from astropy.cosmology import Cosmology

from .common import from_mapping, to_mapping


__all__ = ["read_json", "write_json", "json_identify"]


def read_json(filename, key=None, **kwargs):
    """Read from JSON file.

    Parameters
    ----------
    filename : path-like or file-like
        The JSON file name or actual file.
    key : str or None, optional
        If the JSON is for many cosmologies, ``key`` is needed to select
        the specific cosmology.
        See ``from_mapping``.
    **kwargs
        Not used.

    """
    # read file, from path-like or file-like
    if isinstance(filename, (str, bytes, os.PathLike)):  # pathlike
        with open(filename, "r") as file:
            data = file.read()
    else:  # file-like : this also handles errors in dumping
        data = filename.read()

    mapping = json.loads(data)  # parse json mappable to dict

    return from_mapping(mapping, key=key)


def write_json(cosmology, file, overwrite=False, **kwargs):
    """Write to JSON file.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance
    file : path-like or file-like
    overwrite : bool, optional
        Whether to overwrite the file. Applies only if "file" is path-like.

    """
    data = to_mapping(cosmology)  # start by turning into dict
    data["cosmology"] = data["cosmology"].__name__  # change class field to str

    if isinstance(file, (str, bytes, os.PathLike)):  # pathlike
        # check that file exists and whether to overwrite.
        if os.path.exists(file) and not overwrite:
            raise IOError(f"{file} exists. Set 'overwrite' to write over.")
        with open(file, "w") as write_file:
            json.dump(data, write_file)
    else:  # file-like or error (this handles errors in dumping)
        json.dump(data, file)


def json_identify(origin, filepath, fileobj, *args, **kwargs):
    return filepath is not None and filepath.endswith(".json")


io_registry.register_reader("json", Cosmology, read_json)
io_registry.register_writer("json", Cosmology, write_json)
io_registry.register_identifier("json", Cosmology, json_identify)
