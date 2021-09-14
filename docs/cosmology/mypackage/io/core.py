# -*- coding: utf-8 -*-

# STDLIB
import json
import os

# THIRD PARTY
import numpy as np

# LOCAL
from mypackage.cosmology import MyCosmology


def file_reader(filename):
    """Read files in format 'myformat'.

    Parameters
    ----------
    filename : str, bytes, or `~os.PathLike`

    Returns
    -------
    `~mypackage.cosmology.MyCosmology` instance
    """
    # read
    if isinstance(filename, (str, bytes, os.PathLike)):
        with open(filename, "r") as file:
            data = file.read()
    else:  # file-like : this also handles errors in dumping
        data = filename.read()

    mapping = json.loads(data)  # parse json mappable to dict

    # deserialize list to ndarray
    for k, v in mapping.items():
        if isinstance(v, list):
            mapping[k] = np.array(v)

    return MyCosmology(**mapping)


def file_writer(file, cosmo, overwrite=False):
    """Write files in format 'myformat'.

    Parameters
    ----------
    file : str, bytes, `~os.PathLike`, or file-like
    cosmo : `~mypackage.cosmology.MyCosmology` instance
    overwrite : bool (optional, keyword-only)
        Whether to overwrite an existing file. Default is False.
    """
    output = dict(cosmo)

    # serialize ndarray
    for k, v in output.items():
        if isinstance(v, np.ndarray):
            output[k] = v.tolist()

    # write
    if isinstance(file, (str, bytes, os.PathLike)):
        # check that file exists and whether to overwrite.
        if os.path.exists(file) and not overwrite:
            raise IOError(f"{file} exists. Set 'overwrite' to write over.")
        with open(file, "w") as write_file:
            json.dump(output, write_file)
    else:  # file-like
        json.dump(output, file)
