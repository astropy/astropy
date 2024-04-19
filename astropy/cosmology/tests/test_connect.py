"""Test the (deprecated) ``astropy.cosmology.connect`` module."""

import re
import warnings

import pytest

from astropy.cosmology import connect
from astropy.utils.exceptions import AstropyDeprecationWarning


@pytest.mark.parametrize(
    "name",
    [
        "CosmologyRead",
        "CosmologyWrite",
        "readwrite_registry",
        "CosmologyFromFormat",
        "CosmologyToFormat",
        "convert_registry",
    ],
)
def test_deprecated_imports(name):
    """Test that deprecated imports raise a warning."""
    match = f"{name} is deprecated (since v7.0)"
    with (
        warnings.catch_warnings(),
        pytest.warns(AstropyDeprecationWarning, match=re.escape(match)),
    ):
        warnings.simplefilter("always")
        getattr(connect, name)


def test_fail_import():
    """Test that an import that doesn't exist raises an error."""
    with pytest.raises(
        AttributeError,
        match="module 'astropy.cosmology.connect' has no attribute 'foo'",
    ):
        connect.foo
