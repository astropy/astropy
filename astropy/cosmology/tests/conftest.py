# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Configure the tests for :mod:`astropy.cosmology`."""

##############################################################################
# IMPORTS

# STDLIB
import inspect

# THIRD-PARTY
import pytest

# LOCAL
from astropy.cosmology import core
from astropy.tests.helper import pickle_protocol  # noqa: F403


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


###############################################################################
# FIXTURES

@pytest.fixture
def clean_registry():
    # TODO! with monkeypatch instead for thread safety.
    ORIGINAL_COSMOLOGY_CLASSES = core._COSMOLOGY_CLASSES
    core._COSMOLOGY_CLASSES = {}  # set as empty dict

    yield core._COSMOLOGY_CLASSES

    core._COSMOLOGY_CLASSES = ORIGINAL_COSMOLOGY_CLASSES
