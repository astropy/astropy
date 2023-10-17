# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module provides the tools used to internally run the cosmology test suite
from the installed astropy.  It makes use of the `pytest`_ testing framework.
"""

import inspect

import pytest

from astropy.cosmology import core

__all__ = ["get_redshift_methods", "clean_registry"]

###############################################################################
# FUNCTIONS


def get_redshift_methods(cosmology, include_private=True, include_z2=True):
    """Get redshift methods from a cosmology.

    Parameters
    ----------
    cosmology : |Cosmology| class or instance
    include_private : bool
        Whether to include private methods, i.e. starts with an underscore.
    include_z2 : bool
        Whether to include methods that are functions of 2 (or more) redshifts,
        not the more common 1 redshift argument.

    Returns
    -------
    set[str]
        The names of the redshift methods on `cosmology`, satisfying
        `include_private` and `include_z2`.
    """
    # Get all the method names, optionally sieving out private methods
    methods = set()
    for n in dir(cosmology):
        try:  # get method, some will error on ABCs
            m = getattr(cosmology, n)
        except NotImplementedError:
            continue

        # Add anything callable, optionally excluding private methods.
        if callable(m) and (not n.startswith("_") or include_private):
            methods.add(n)

    # Sieve out incompatible methods.
    # The index to check for redshift depends on whether cosmology is a class
    # or instance and does/doesn't include 'self'.
    iz1 = int(isinstance(cosmology, type))
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
        elif len(params) >= iz1 + 1 and not params[iz1].startswith(
            "z"
        ):  # First non-self arg is z.
            methods.discard(n)
        # If methods with 2 z args are not allowed, the following arg is checked.
        elif (
            not include_z2
            and (len(params) >= iz1 + 2)
            and params[iz1 + 1].startswith("z")
        ):
            methods.discard(n)

    return methods


###############################################################################
# FIXTURES


@pytest.fixture
def clean_registry():
    """`pytest.fixture` for clearing and restoring ``_COSMOLOGY_CLASSES``."""
    # TODO! with monkeypatch instead for thread safety.
    ORIGINAL_COSMOLOGY_CLASSES = core._COSMOLOGY_CLASSES
    core._COSMOLOGY_CLASSES = {}  # set as empty dict

    yield core._COSMOLOGY_CLASSES

    core._COSMOLOGY_CLASSES = ORIGINAL_COSMOLOGY_CLASSES
