# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Reader and Writer for `~astropy.cosmology.Cosmology` in the JSON format.

This module is NOT considered public API.

"""

import json
import os

from astropy.io import registry as io_registry
from astropy.cosmology import Cosmology

from .core import from_mapping, to_mapping


__all__ = ["read_json", "write_json", "json_identify"]


_all_pathlike = (str, bytes, os.PathLike)  # all recognized path-like types


def read_json(filename, key=None, move_to_meta=False, **kwargs):
    """Read a cosmology from a JSON file.

    Parameters
    ----------
    filename : path-like or file-like
        The JSON file name or actual file.
    key : str or None, optional
        If the JSON is a nested set of cosmologies, ``key`` is needed to select
        the specific cosmology.
    move_to_meta : bool (optional, keyword-only)
        Whether to move keyword arguments that are not in the Cosmology class'
        signature to the Cosmology's metadata. This will only be applied if the
        Cosmology does NOT have a keyword-only argument (e.g. ``**kwargs``).
        Arguments moved to the metadata will be merged with existing metadata,
        preferring specified metadata in the case of a merge conflict
        (e.g. for ``Cosmology(meta={'key':10}, key=42)``, the ``Cosmology.meta``
        will be ``{'key': 10}``).
    **kwargs
        Not used.
        If 'format' is a kwarg, it must be 'json'.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance
    """
    # check 'format' correctness
    format = kwargs.pop("format", "json")
    if format != "json":  # check that if specified, it's JSON
        raise ValueError(f"'format', if specified, must be 'json' not {format}")

    # read file, from path-like or file-like
    if isinstance(filename, _all_pathlike):  # pathlike
        with open(filename, "r") as file:
            data = file.read()
    else:  # file-like : this also handles errors in dumping
        data = filename.read()

    mapping = json.loads(data)  # parse json mappable to dict
    if key is not None:  # select from dict. Enables nested storage of cosmos
        mapping = mapping[key]

    return from_mapping(mapping, move_to_meta=move_to_meta)


def write_json(cosmology, file, *, overwrite=False, **kwargs):
    """Write a cosmology to a file in a JSON format.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance
    file : path-like or file-like
        The file location or the file itself.
    overwrite : bool (optional, keyword-only)
        Whether to overwrite the file. Applies only if "file" is path-like.
    **kwargs
        Not used.
    """
    data = to_mapping(cosmology)  # start by turning into dict
    data["cosmology"] = data["cosmology"].__name__  # change class field to str

    if isinstance(file, _all_pathlike):  # pathlike
        # check that file exists and whether to overwrite.
        if os.path.exists(file) and not overwrite:
            raise IOError(f"{file} exists. Set 'overwrite' to write over.")
        with open(file, "w") as write_file:
            json.dump(data, write_file)
    else:  # file-like or error (this handles errors in dumping)
        json.dump(data, file)


def json_identify(origin, filepath, fileobj, *args, **kwargs):
    """Identify if object uses the JSON format.

    Returns
    -------
    bool
        True if the 'filepath' suffix is '.json', False otherwise.
    """
    return filepath is not None and filepath.endswith(".json")


io_registry.register_reader("json", Cosmology, read_json)
io_registry.register_writer("json", Cosmology, write_json)
io_registry.register_identifier("json", Cosmology, json_identify)
