# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Configure the tests for :mod:`astropy.cosmology`."""

##############################################################################
# IMPORTS

# STDLIB
import inspect
import json
import os

# THIRD-PARTY
import pytest

# LOCAL
import astropy.cosmology.units as cu
import astropy.units as u
from astropy.cosmology import core
from astropy.cosmology.core import Cosmology
from astropy.tests.helper import pickle_protocol


###############################################################################
# FUNCTIONS

def get_redshift_methods(cosmology, allow_private=True, allow_z2=True):
    """Get redshift methods from a cosmology.

    Parameters
    ----------
    cosmology : |Cosmology| class or instance

    Returns
    -------
    set[str]
    """
    methods = set()
    for n in dir(cosmology):
        try:  # get method, some will error on ABCs
            m = getattr(cosmology, n)
        except NotImplementedError:
            continue

        # Add anything callable, optionally excluding private methods.
        if callable(m) and (not n.startswith('_') or allow_private):
            methods.add(n)

    # Sieve out incompatible methods.
    # The index to check for redshift depends on whether cosmology is a class
    # or instance and does/doesn't include 'self'.
    iz1 = 1 if inspect.isclass(cosmology) else 0
    for n in tuple(methods):
        try:
            sig = inspect.signature(getattr(cosmology, n))
        except ValueError:  # Remove non-introspectable methods.
            methods.discard(n)
            continue
        else:
            params = list(sig.parameters.keys())

        # Remove non redshift methods:
        if len(params) <= iz1:  # Check there are enough arguments.
            methods.discard(n)
        elif len(params) >= iz1 + 1 and not params[iz1].startswith("z"):  # First non-self arg is z.
            methods.discard(n)
        # If methods with 2 z args are not allowed, the following arg is checked.
        elif not allow_z2 and (len(params) >= iz1 + 2) and params[iz1 + 1].startswith("z"):
            methods.discard(n)

    return methods


def read_json(filename, **kwargs):
    """Read JSON.

    Parameters
    ----------
    filename : str
    **kwargs
        Keyword arguments into :meth:`~astropy.cosmology.Cosmology.from_format`

    Returns
    -------
    `~astropy.cosmology.Cosmology` instance

    """
    # read
    if isinstance(filename, (str, bytes, os.PathLike)):
        with open(filename, "r") as file:
            data = file.read()
    else:  # file-like : this also handles errors in dumping
        data = filename.read()

    mapping = json.loads(data)  # parse json mappable to dict

    # deserialize Quantity
    with u.add_enabled_units(cu.redshift):
        for k, v in mapping.items():
            if isinstance(v, dict) and "value" in v and "unit" in v:
                mapping[k] = u.Quantity(v["value"], v["unit"])
        for k, v in mapping.get("meta", {}).items():  # also the metadata
            if isinstance(v, dict) and "value" in v and "unit" in v:
                mapping["meta"][k] = u.Quantity(v["value"], v["unit"])

    return Cosmology.from_format(mapping, **kwargs)


def write_json(cosmology, file, *, overwrite=False):
    """Write Cosmology to JSON.

    Parameters
    ----------
    cosmology : `astropy.cosmology.Cosmology` subclass instance
    file : path-like or file-like
    overwrite : bool (optional, keyword-only)
    """
    data = cosmology.to_format("mapping")  # start by turning into dict
    data["cosmology"] = data["cosmology"].__qualname__

    # serialize Quantity
    for k, v in data.items():
        if isinstance(v, u.Quantity):
            data[k] = {"value": v.value.tolist(), "unit": str(v.unit)}
    for k, v in data.get("meta", {}).items():  # also serialize the metadata
        if isinstance(v, u.Quantity):
            data["meta"][k] = {"value": v.value.tolist(), "unit": str(v.unit)}

    # check that file exists and whether to overwrite.
    if os.path.exists(file) and not overwrite:
        raise IOError(f"{file} exists. Set 'overwrite' to write over.")
    with open(file, "w") as write_file:
        json.dump(data, write_file)


def json_identify(origin, filepath, fileobj, *args, **kwargs):
    return filepath is not None and filepath.endswith(".json")


###############################################################################
# FIXTURES

@pytest.fixture
def clean_registry():
    # TODO! with monkeypatch instead for thread safety.
    ORIGINAL_COSMOLOGY_CLASSES = core._COSMOLOGY_CLASSES
    core._COSMOLOGY_CLASSES = {}  # set as empty dict

    yield core._COSMOLOGY_CLASSES

    core._COSMOLOGY_CLASSES = ORIGINAL_COSMOLOGY_CLASSES
